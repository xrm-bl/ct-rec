#!/usr/bin/env python3
"""
Noise2Inverse (X:1 strategy) -- 2D / 2.5D / 3D U-Net, self-contained PyTorch version.

Input : K directories (K>=2) of 32-bit float TIFF slices, each reconstructed from a
        disjoint (interleaved) subset of projections. Same slice count/order in each.

  # 2D
  python n2i_v2.py train --mode 2d  --splits rec0 rec1 rec2 rec3 --out n2i_2d.pt
  # 2.5D (2*s+1 neighbouring slices in, centre slice out)
  python n2i_v2.py train --mode 25d --s 2 --splits rec0 rec1 rec2 rec3 --out n2i_25d.pt
  # 3D (cache strongly recommended: random 3D crops from a raw memmap)
  python n2i_v2.py train --mode 3d  --splits rec0 rec1 rec2 rec3 --cache /fast/n2i.npy --out n2i_3d.pt --amp
  # inference (mode is taken from the checkpoint)
  python n2i_v2.py infer --splits rec0 rec1 rec2 rec3 --model n2i_3d.pt --outdir denoised --fp16

Artifact countermeasures built in:
  * upsampling = bilinear/trilinear resize + 3x3 conv (no checkerboard from ConvTranspose)
  * replicate padding in convolutions
  * BatchNorm (fixed statistics at inference -> output does not depend on tile contents)
  * tiled inference with overlap and cosine cross-fade; 3D slabs read with z-margins
Large volumes:
  * train --cache stores only the slices used for training. Without --start/--stop the
    central slices are selected automatically: at most --max-train-slices (500), and few
    enough that the cache fits in --mem-frac (0.7) of the RAM and in the free disk space.
  * a cache is verified against its source files (<cache>.json) and rebuilt if stale;
    infer uses a cache only if it covers exactly the whole stack, otherwise reads TIFFs.
  * CUDA out-of-memory: training splits the batch into micro-batches, inference halves the tile.
Requires: torch, numpy, tifffile
"""
import argparse, functools, glob, itertools, json, math, os, shutil, time
import numpy as np
import tifffile
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader, Dataset

DEFAULTS = {
    "2d":  dict(patch=256, batch=16, ppr=4, depth=4, tile=[1024, 1024], ov=[64, 64]),
    "25d": dict(patch=256, batch=16, ppr=4, depth=4, tile=[1024, 1024], ov=[64, 64]),
    "3d":  dict(patch=128, batch=4, ppr=2, depth=3, tile=[64, 256, 256], ov=[16, 32, 32]),
}


def ceil_to(v, m):
    return -(-v // m) * m


# ---------------------------------------------------------------- data
def scan_splits(dirs):
    lists = [sorted(glob.glob(os.path.join(d, "*.tif*"))) for d in dirs]
    n = len(lists[0])
    if n == 0:
        raise SystemExit(f"no tif files in {dirs[0]}")
    if any(len(l) != n for l in lists):
        raise SystemExit("slice counts differ between split directories")
    H, W = tifffile.imread(lists[0][0]).shape
    return lists, n, H, W


def file_sig(lists):
    """Identifies the exact source files of a cache: (path, size, mtime) of every slice."""
    return [[os.path.abspath(f), os.path.getsize(f), int(os.path.getmtime(f))] for l in lists for f in l]


def total_ram():
    try:
        if os.name == "nt":
            import ctypes

            class MemStat(ctypes.Structure):
                _fields_ = [("dwLength", ctypes.c_ulong), ("dwMemoryLoad", ctypes.c_ulong)] + \
                           [(n, ctypes.c_ulonglong) for n in ("tot", "av", "totpf", "avpf", "totv", "avv", "avx")]
            st = MemStat()
            st.dwLength = ctypes.sizeof(MemStat)
            if ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(st)):
                return int(st.tot)
        else:
            return os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
    except Exception:
        pass
    return None


class Volume:
    """K split volumes restricted to slices [z0, z1) of the stack (default: all).
    get() returns (K, len(zs), h, w) float32 with indices relative to z0, clamped."""

    def __init__(s, dirs, cache=None, lru=0, zr=None, build=False):
        lists, ntot, s.H, s.W = scan_splits(dirs)
        s.z0, z1 = zr if zr else (0, ntot)
        s.lists = [l[s.z0:z1] for l in lists]
        s.K, s.N, s.Ntot = len(s.lists), z1 - s.z0, ntot
        s.cache_path, s.lru, s.mm = None, lru, None
        if cache:
            shape = (s.K, s.N, s.H, s.W)
            sig = file_sig(s.lists)
            ok = s._cache_ok(cache, shape, sig)
            if not ok and build:
                s._build_cache(cache, shape, sig)
                ok = True
            if ok:
                s.cache_path = cache
                s.mm = np.load(cache, mmap_mode="r")
            else:
                print(f"cache {cache} does not cover exactly these slices -- reading TIFF files instead",
                      flush=True)
        s._read = functools.lru_cache(maxsize=lru)(s._read_raw) if lru > 0 else s._read_raw

    @staticmethod
    def _cache_ok(cache, shape, sig):
        if not os.path.exists(cache):
            return False
        try:
            if tuple(np.load(cache, mmap_mode="r").shape) != shape:
                return False
        except Exception:
            return False
        meta = cache + ".json"
        if not os.path.exists(meta):
            print(f"warning: {cache} has no metadata (older version) and cannot be verified; "
                  f"delete it if the data has changed", flush=True)
            return True
        try:
            with open(meta, encoding="utf-8") as f:
                return json.load(f).get("sig") == sig
        except Exception:
            return False

    def _build_cache(s, cache, shape, sig):
        need = int(np.prod(shape)) * 4
        d = os.path.dirname(os.path.abspath(cache))
        free = shutil.disk_usage(d).free + (os.path.getsize(cache) if os.path.exists(cache) else 0)
        if need > free:
            raise SystemExit(f"not enough disk space for the cache: need {need / 1e9:.1f} GB, "
                             f"free {free / 1e9:.1f} GB in {d}")
        for f in (cache, cache + ".json"):
            if os.path.exists(f):
                os.remove(f)
        print(f"building cache {cache} ({need / 1e9:.1f} GB, slices {s.z0}-{s.z0 + s.N - 1})", flush=True)
        mm = np.lib.format.open_memmap(cache + ".tmp", mode="w+", dtype=np.float32, shape=shape)
        for k in range(s.K):
            for z in range(s.N):
                mm[k, z] = tifffile.imread(s.lists[k][z])
                if (z + 1) % 200 == 0:
                    print(f"  split {k}: {z + 1}/{s.N}", flush=True)
        print(f"  writing cache to disk ({need / 1e9:.1f} GB) -- this can take a while, please wait ...",
              flush=True)
        t0 = time.time()
        mm.flush()
        del mm
        os.replace(cache + ".tmp", cache)
        print(f"  cache written ({time.time() - t0:.0f} s)", flush=True)
        with open(cache + ".json", "w", encoding="utf-8") as f:
            json.dump({"z0": s.z0, "shape": list(shape), "sig": sig}, f)

    # Windows spawns DataLoader workers and pickles the dataset: never pickle the
    # memmap (it would be copied in full); reopen it in each worker instead.
    def __getstate__(s):
        d = s.__dict__.copy()
        d["mm"] = None
        d.pop("_read", None)
        return d

    def __setstate__(s, d):
        s.__dict__.update(d)
        if s.cache_path:
            s.mm = np.load(s.cache_path, mmap_mode="r")
        s._read = functools.lru_cache(maxsize=s.lru)(s._read_raw) if s.lru > 0 else s._read_raw

    def _read_raw(s, k, z):
        return tifffile.imread(s.lists[k][z])

    def get(s, zs, y0=0, y1=None, x0=0, x1=None):
        zs = np.clip(np.asarray(list(zs)), 0, s.N - 1)
        ys, xs = slice(y0, y1), slice(x0, x1)
        if s.mm is not None:
            return np.asarray(s.mm[:, zs, ys, xs], dtype=np.float32)
        return np.stack([np.stack([s._read(k, int(z))[ys, xs] for z in zs])
                         for k in range(s.K)]).astype(np.float32)


def train_range(a, K, N, H, W):
    """Slice range [lo, hi) used for training (and stored in the cache)."""
    if a.start is not None or a.stop is not None:
        lo, hi = slice(a.start, a.stop).indices(N)[:2]
        return lo, hi, "--start/--stop"
    if not a.cache:
        return 0, N, "all slices"
    per = K * H * W * 4                       # cache bytes per slice
    n, why = N, "all slices"
    if a.train_slices:
        n, why = a.train_slices, "--train-slices"
    else:
        if a.max_train_slices and a.max_train_slices < n:
            n, why = a.max_train_slices, f"central {a.max_train_slices} slices (--max-train-slices)"
        ram = total_ram()
        if ram:
            n_ram = int(a.mem_frac * ram // per)
            if n_ram < n:
                n, why = n_ram, f"cache limited to {a.mem_frac:.0%} of RAM ({ram / 1e9:.0f} GB)"
        d = os.path.dirname(os.path.abspath(a.cache))
        free = shutil.disk_usage(d).free + (os.path.getsize(a.cache) if os.path.exists(a.cache) else 0)
        n_disk = int(0.9 * free // per)
        if n_disk < n:
            n, why = n_disk, f"cache limited by free disk space ({free / 1e9:.0f} GB)"
    n = max(0, min(n, N))
    lo = (N - n) // 2
    return lo, lo + n, why


def footprint(mode, s, pz):
    return {"2d": 1, "25d": 2 * s + 1, "3d": pz}[mode]


def make_pair(im, k, mode, s):
    """im: (K, Lz, h, w) normalised. X:1 -> input = mean of other splits, target = split k."""
    tgt = im[k]
    inp = (im.sum(0) - tgt) / (im.shape[0] - 1)
    if mode == "2d":
        return inp, tgt                    # (1,h,w), (1,h,w)
    if mode == "25d":
        return inp, tgt[s:s + 1]           # (2s+1,h,w), (1,h,w)
    return inp[None], tgt[None]            # (1,Lz,h,w), (1,Lz,h,w)


class Patches(Dataset):
    def __init__(s, vol, starts, a, mean, std, length):
        s.vol, s.starts, s.a, s.mean, s.std, s.length = vol, starts, a, mean, std, length

    def __len__(s):
        return s.length

    def __getitem__(s, _):
        rng = np.random.default_rng()
        a, v = s.a, s.vol
        p, L = a.patch, footprint(a.mode, a.s, a.pz)
        z0 = s.starts[rng.integers(len(s.starts))]
        zs = range(z0, z0 + L)
        pos = [(rng.integers(v.H - p + 1), rng.integers(v.W - p + 1)) for _ in range(a.ppr)]
        if v.mm is None:   # read whole slices once, crop several patches
            blk = v.get(zs)
            crops = [blk[..., y:y + p, x:x + p] for y, x in pos]
        else:
            crops = [v.get(zs, y, y + p, x, x + p) for y, x in pos]
        ins, tgs = [], []
        for c in crops:
            inp, tgt = make_pair((c - s.mean) / s.std, rng.integers(v.K), a.mode, a.s)
            r = int(rng.integers(4))
            inp, tgt = np.rot90(inp, r, axes=(-2, -1)), np.rot90(tgt, r, axes=(-2, -1))
            if rng.integers(2):
                inp, tgt = inp[..., ::-1], tgt[..., ::-1]
            if rng.integers(2):            # flip along z
                if a.mode == "3d":
                    inp, tgt = inp[:, ::-1], tgt[:, ::-1]
                elif a.mode == "25d":
                    inp = inp[::-1]
            ins.append(np.ascontiguousarray(inp))
            tgs.append(np.ascontiguousarray(tgt))
        return torch.from_numpy(np.stack(ins)), torch.from_numpy(np.stack(tgs))


# ---------------------------------------------------------------- network
def block(nd, ci, co):
    Conv, BN = (nn.Conv2d, nn.BatchNorm2d) if nd == 2 else (nn.Conv3d, nn.BatchNorm3d)
    return nn.Sequential(
        Conv(ci, co, 3, padding=1, padding_mode="replicate"), BN(co), nn.LeakyReLU(0.1, True),
        Conv(co, co, 3, padding=1, padding_mode="replicate"), BN(co), nn.LeakyReLU(0.1, True))


class UNet(nn.Module):
    def __init__(s, nd=2, cin=1, ch=32, depth=4, up="resize"):
        super().__init__()
        Conv = nn.Conv2d if nd == 2 else nn.Conv3d
        s.pool = nn.MaxPool2d(2) if nd == 2 else nn.MaxPool3d(2)
        s.center = cin // 2
        s.downs, s.ups, s.decs = nn.ModuleList(), nn.ModuleList(), nn.ModuleList()
        c, chs = cin, []
        for i in range(depth):
            co = ch * 2 ** i
            s.downs.append(block(nd, c, co))
            chs.append(co)
            c = co
        s.mid = block(nd, c, 2 * c)
        c *= 2
        imode = "bilinear" if nd == 2 else "trilinear"
        for co in reversed(chs):
            if up == "resize":   # resize-conv: no checkerboard pattern
                s.ups.append(nn.Sequential(
                    nn.Upsample(scale_factor=2, mode=imode, align_corners=False),
                    Conv(c, co, 3, padding=1, padding_mode="replicate")))
            else:                # original transposed conv (kept for comparison)
                T = nn.ConvTranspose2d if nd == 2 else nn.ConvTranspose3d
                s.ups.append(T(c, co, 2, stride=2))
            s.decs.append(block(nd, 2 * co, co))
            c = co
        s.out = Conv(c, 1, 1)

    def forward(s, x):
        x0, skips = x[:, s.center:s.center + 1], []
        for d in s.downs:
            x = d(x)
            skips.append(x)
            x = s.pool(x)
        x = s.mid(x)
        for up, dec, sk in zip(s.ups, s.decs, reversed(skips)):
            x = dec(torch.cat([up(x), sk], 1))
        return s.out(x) + x0   # residual on the centre slice


# ---------------------------------------------------------------- tiled prediction
def starts1d(L, t, ov):
    if L <= t:
        return [0]
    st = list(range(0, L - t, t - ov))
    st.append(L - t)
    return st


def window1d(t, ov):
    w = torch.ones(t)
    if ov > 0:
        r = torch.linspace(0, 1, ov + 2)[1:-1]          # strictly inside (0,1)
        ramp = torch.sin(r * math.pi / 2) ** 2
        w[:ov] = ramp
        w[-ov:] = torch.minimum(w[-ov:], ramp.flip(0))
    return w


@torch.no_grad()
def predict(model, x, tile, ov, m, dev, fp16):
    """x: CPU tensor (C, *spatial), normalised. Returns CPU tensor (*spatial)."""
    sp = x.shape[1:]
    Pd = [ceil_to(v, m) for v in sp]
    td = [min(ceil_to(t, m), p) for t, p in zip(tile, Pd)]
    od = [min(o, t // 2) for o, t in zip(ov, td)]
    pad = []
    for v, p in zip(reversed(sp), reversed(Pd)):
        pad += [0, p - v]
    xp = F.pad(x[None], pad, mode="replicate")[0] if any(pad) else x
    out, wsum = torch.zeros(Pd), torch.zeros(Pd)
    W = None
    for w in (window1d(t, o) for t, o in zip(td, od)):
        W = w if W is None else W[..., None] * w
    for st in itertools.product(*[starts1d(p, t, o) for p, t, o in zip(Pd, td, od)]):
        sl = tuple(slice(s0, s0 + t) for s0, t in zip(st, td))
        xt = xp[(slice(None),) + sl][None].to(dev, non_blocking=True)
        with torch.autocast(dev.type, dtype=torch.float16, enabled=fp16 and dev.type == "cuda"):
            y = model(xt)
        out[sl] += y[0, 0].float().cpu() * W
        wsum[sl] += W
    out /= wsum
    return out[tuple(slice(0, v) for v in sp)]


# ---------------------------------------------------------------- train
def train(a):
    for k, v in DEFAULTS[a.mode].items():
        if getattr(a, k, None) is None:
            setattr(a, k, v)
    lists, Ntot, H, W = scan_splits(a.splits)
    K = len(lists)
    z0, z1, why = train_range(a, K, Ntot, H, W)
    print(f"training slices {z0}-{z1 - 1} of {Ntot} ({why})", flush=True)
    nd, cin, m = (3 if a.mode == "3d" else 2), (2 * a.s + 1 if a.mode == "25d" else 1), 2 ** a.depth
    a.patch = min(a.patch, H, W) // m * m
    if a.mode == "3d" and a.pz % m:
        raise SystemExit(f"--pz must be a multiple of {m}")
    L = footprint(a.mode, a.s, a.pz)
    vb = max(L, a.nval)                       # contiguous held-out z block
    if z1 - z0 < vb + L + 1:
        raise SystemExit(f"training range too small ({z1 - z0} slices, need >= {vb + L + 1}); "
                         f"free memory/disk, raise --mem-frac or set --train-slices")

    vol = Volume(a.splits, a.cache, zr=(z0, z1), build=True)
    N = vol.N
    lo, hi = 0, N                              # indices relative to the training range
    v0 = (lo + hi - vb) // 2
    v1 = v0 + vb
    starts = [z for z in range(lo, hi - L + 1) if z + L <= v0 or z >= v1]

    rng = np.random.default_rng(0)
    samp = rng.choice(np.arange(lo, hi), min(20, hi - lo), replace=False)
    v = np.concatenate([vol.get([z]).mean(0).ravel()[::17] for z in samp])
    mean, std = float(v.mean()), float(v.std())
    print(f"mode={a.mode} K={K} N={N} {H}x{W} train-starts={len(starts)} "
          f"val=[{v0},{v1}) mean={mean:.4g} std={std:.4g}", flush=True)
    if z0:
        print(f"  (val is indexed within the training range / cache: slices {z0 + v0}-{z0 + v1 - 1} of the stack)",
              flush=True)

    vs = min(2 * a.patch, H, W) // m * m
    vy, vx = (H - vs) // 2, (W - vs) // 2
    vz = list(range(v0, v1 - L + 1))
    vz = vz[::max(1, len(vz) // 4)][:4]
    vin, vtg = [], []
    for z0 in vz:
        im = (vol.get(range(z0, z0 + L), vy, vy + vs, vx, vx + vs) - mean) / std
        for k in range(K):
            i_, t_ = make_pair(im, k, a.mode, a.s)
            vin.append(i_)
            vtg.append(t_)
    vin, vtg = torch.from_numpy(np.stack(vin)), torch.from_numpy(np.stack(vtg))

    dev = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = UNet(nd, cin, a.ch, a.depth, a.up).to(dev)
    opt = torch.optim.Adam(model.parameters(), lr=a.lr)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, a.iters)
    bs = max(1, a.batch // a.ppr)
    dl = DataLoader(Patches(vol, starts, a, mean, std, a.iters * bs + bs), batch_size=bs,
                    num_workers=a.workers, pin_memory=True, persistent_workers=a.workers > 0)
    amp = dict(device_type=dev.type, dtype=torch.bfloat16, enabled=a.amp and dev.type == "cuda")

    def save(path):
        torch.save({"model": model.state_dict(), "mean": mean, "std": std, "mode": a.mode,
                    "s": a.s, "nd": nd, "cin": cin, "ch": a.ch, "depth": a.depth, "up": a.up}, path)

    best, run, it, micro = math.inf, 0.0, 0, 1
    model.train()
    for inp, tgt in dl:
        inp = inp.flatten(0, 1).to(dev, non_blocking=True)
        tgt = tgt.flatten(0, 1).to(dev, non_blocking=True)
        while True:   # on CUDA OOM: split the batch into micro-batches (gradient accumulation)
            try:
                opt.zero_grad(set_to_none=True)
                loss_v = 0.0
                for ci, ct in zip(inp.chunk(micro), tgt.chunk(micro)):
                    with torch.autocast(**amp):
                        pred = model(ci)
                    loss = F.mse_loss(pred.float(), ct) * (ci.shape[0] / inp.shape[0])   # L2 (N2I)
                    loss.backward()
                    loss_v += loss.item()
                break
            except torch.cuda.OutOfMemoryError:
                opt.zero_grad(set_to_none=True)
                pred = loss = None
                torch.cuda.empty_cache()
                if micro >= inp.shape[0]:
                    raise SystemExit("CUDA out of memory even with one patch per step; reduce --patch")
                micro *= 2
                print(f"CUDA out of memory: splitting each batch into {micro} micro-batches", flush=True)
        opt.step()
        sched.step()
        run += loss_v
        it += 1
        if it % a.log == 0 or it == a.iters:
            model.eval()
            with torch.no_grad(), torch.autocast(**amp):
                vl = float(np.mean([F.mse_loss(model(vin[j:j + 1].to(dev)).float(),
                                               vtg[j:j + 1].to(dev)).item() for j in range(len(vin))]))
            model.train()
            tag = ""
            if vl < best:
                best = vl
                save(a.out)
                tag = " *saved"
            print(f"it {it:6d}  train {run / a.log:.5f}  val {vl:.5f}{tag}", flush=True)
            run = 0.0
        if it >= a.iters:
            break
    save(os.path.splitext(a.out)[0] + "_last.pt")


# ---------------------------------------------------------------- inference
def infer(a):
    ck = torch.load(a.model, map_location="cpu")
    dev = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = UNet(ck["nd"], ck["cin"], ck["ch"], ck["depth"], ck["up"])
    model.load_state_dict(ck["model"])
    model.to(dev).eval()
    mode, s, mean, std, m = ck["mode"], ck["s"], ck["mean"], ck["std"], 2 ** ck["depth"]
    tile = list(a.tile or DEFAULTS[mode]["tile"])
    ov = a.ov or DEFAULTS[mode]["ov"]
    nsp = 3 if mode == "3d" else 2
    if len(tile) != nsp or len(ov) != nsp:
        raise SystemExit(f"--tile/--ov need {nsp} values for mode {mode}")
    lru = 4 * (2 * s + 1) if mode == "25d" else 0
    vol = Volume(a.splits, a.cache, lru=lru * len(a.splits), build=False)
    lo, hi = slice(a.start, a.stop).indices(vol.N)[:2]
    os.makedirs(a.outdir, exist_ok=True)

    def run(blk):   # blk: (K, Z, H, W) normalised; average over the K strategies
        acc = None
        for k in range(vol.K):
            inp, _ = make_pair(blk, k, mode, s)
            x = torch.from_numpy(np.ascontiguousarray(inp))
            while True:   # on CUDA OOM: halve the tile and retry
                try:
                    y = predict(model, x, tile, ov, m, dev, a.fp16)
                    break
                except torch.cuda.OutOfMemoryError:
                    torch.cuda.empty_cache()
                    if max(tile[-2:]) <= 128:
                        raise SystemExit("CUDA out of memory even with small tiles")
                    tile[-2:] = [max(128, t // 2) for t in tile[-2:]]
                    print(f"CUDA out of memory: tile reduced to {tile}", flush=True)
            acc = y if acc is None else acc + y
        return (acc / vol.K).numpy() * std + mean

    def write(z, img):
        tifffile.imwrite(os.path.join(a.outdir, os.path.basename(vol.lists[0][z])),
                         img.astype(np.float32))

    if mode in ("2d", "25d"):
        for z in range(lo, hi):
            zs = [z] if mode == "2d" else range(z - s, z + s + 1)
            write(z, run((vol.get(zs) - mean) / std))
            if z % 50 == 0:
                print(f"slice {z}/{hi}", flush=True)
    else:
        S, M = a.slab, ov[0]   # slab thickness, z-margin on each side
        for z in range(lo, hi, S):
            z1 = min(z + S, hi)
            out = run((vol.get(range(z - M, z1 + M)) - mean) / std)
            for j, zz in enumerate(range(z, z1)):
                write(zz, out[M + j])
            print(f"slices {z}-{z1 - 1}/{hi}", flush=True)


# ---------------------------------------------------------------- cli
if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest="cmd", required=True)
    for name in ("train", "infer"):
        q = sp.add_parser(name)
        q.add_argument("--splits", nargs="+", required=True)
        q.add_argument("--cache", default=None, help="raw float32 memmap (.npy) of the splits (train: training range only)")
        q.add_argument("--start", type=int, default=None)
        q.add_argument("--stop", type=int, default=None)
    t = sp.choices["train"]
    t.add_argument("--mode", choices=["2d", "25d", "3d"], default="2d")
    t.add_argument("--s", type=int, default=2, help="2.5D half width (2s+1 slices)")
    t.add_argument("--pz", type=int, default=32, help="3D patch depth")
    t.add_argument("--patch", type=int, default=None)
    t.add_argument("--batch", type=int, default=None)
    t.add_argument("--ppr", type=int, default=None, help="patches per read")
    t.add_argument("--depth", type=int, default=None)
    t.add_argument("--ch", type=int, default=32)
    t.add_argument("--up", choices=["resize", "transpose"], default="resize")
    t.add_argument("--out", default="n2i.pt")
    t.add_argument("--iters", type=int, default=2000)
    t.add_argument("--lr", type=float, default=5e-4)
    t.add_argument("--nval", type=int, default=8)
    t.add_argument("--workers", type=int, default=8)
    t.add_argument("--log", type=int, default=500)
    t.add_argument("--amp", action="store_true", help="bf16 autocast for training")
    t.add_argument("--mem-frac", type=float, default=0.7,
                   help="with --cache and no --start/--stop: limit the cache to this fraction of RAM")
    t.add_argument("--train-slices", type=int, default=None,
                   help="with --cache and no --start/--stop: number of central slices to train on")
    t.add_argument("--max-train-slices", type=int, default=500,
                   help="with --cache and no --start/--stop: upper limit of the automatically "
                        "selected central slices (0 = no limit)")
    f = sp.choices["infer"]
    f.add_argument("--model", required=True)
    f.add_argument("--outdir", required=True)
    f.add_argument("--tile", type=int, nargs="+", default=None, help="(y x) or (z y x)")
    f.add_argument("--ov", type=int, nargs="+", default=None, help="tile overlap, same order")
    f.add_argument("--slab", type=int, default=64, help="3D: output slices per slab")
    f.add_argument("--fp16", action="store_true")
    a = ap.parse_args()
    train(a) if a.cmd == "train" else infer(a)