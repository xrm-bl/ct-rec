#!/usr/bin/env python3
"""
n2i_v4.py -- Noise2Inverse (X:1 strategy) CT denoising, 2D / 2.5D / 3D U-Net, PyTorch.

Input : K directories (K>=2) of 32-bit float TIFF slices, each reconstructed from a
        disjoint (interleaved) subset of projections. Same slice count/order in each.

  python n2i_v4.py train
      (defaults: splits = rec/0, rec/1 (, rec/2, rec/3 ...) detected automatically,
       --mode 25d, --cache image_cache.npy, models in training/)
  python n2i_v4.py infer
      (defaults: model training/best.pt, output denoised/)

Changes from n2i_v3.py (models trained with v2/v3 can be used for inference):
  * the training crop starts by default 1/4 of the image height below the top edge
    (--crop-y0 0.25; x stays centred, --crop-x0), instead of being centred: this keeps
    the rotation centre (ring artifacts) out of large images and covers a wider range of
    noise levels. --crop-y0 -1 restores the centred crop of v3.

  * time breakdown: training prints the share of data wait / GPU training / validation /
    saving / inference check at every --log and a table at the end; inference prints the
    share of waiting for reading / GPU / waiting for writing and the s/slice
  * tiles are laid out evenly with the fewest tiles that keep the overlap (the old layout
    computed up to 2.25x the image area with large tiles)
  * bench command: times the network for each memory layout (--layout cl|nchw) and replicate-
    padding implementation (--pad-impl torch|cat, identical values) and profiles GPU kernels
  * inference prepares its inputs on the GPU (normalisation, strategy inputs, tiling; --pre)
  * --pad-impl cat has a hand-written backward (exact adjoint of the padding) for training
  * defaults chosen with bench on a Quadro RTX 8000: training --layout nchw --pad-impl torch
    (channels_last made the training-mode BatchNorm ~4x slower), inference --layout cl
    --pad-impl cat
  * for inference, replicate padding is built with torch.cat (--pad-impl cat): PyTorch's own
    replication_pad kernel plus the layout conversions it caused took ~2/3 of the GPU time
    with channels_last; same values, about 2.6x faster network on a Quadro RTX 8000
  * faster inference with identical arithmetic: BatchNorm is folded into the preceding
    convolution (--no-fold-bn to disable; the difference on a test input is printed), and
    the K splitting strategies go through the network as one batch

Changes from n2i_v2.py to n2i_v3.py:
  * mixed precision by default: training uses bf16 (Ampere or newer) or fp16 with loss
    scaling (older GPUs) for the convolutions, with fp32 weights and loss (--precision);
    inference runs in fp16 by default (--fp32 to disable)
  * training data limited to a full-resolution crop of at most 1024 x 1024 pixels
    (--crop, --crop-center; pixels are cropped, never binned) and at most 512 slices
  * channels_last memory format and cudnn.benchmark for faster convolutions
  * inference overlaps reading (prefetch thread), GPU computation and writing (writer
    threads); tile blending is accumulated on the GPU
  * train writes the models to a directory (--out, default training/): 000100.pt, 000200.pt ...
    at every --log (default 100) and best.pt (lowest val); at the same interval one untrained
    slice is denoised into infer_check/%06d.tif (average.tif = mean of the splits);
    infer --crop processes
    only a central (or --crop-center) region, e.g. to compare models quickly

Artifact countermeasures: resize-conv upsampling (no checkerboard), replicate padding,
  BatchNorm (fixed statistics at inference), overlapping tiles with cosine cross-fade.
Cache: train --cache stores only the training slices/crop, verified against its source
  files (<cache>.json) and rebuilt if stale; limited to --mem-frac (0.7) of RAM and the free
  disk space. Inference reads the TIFF files sequentially (no cache needed).
CUDA out-of-memory: training splits the batch into micro-batches, inference halves the tile.
Output: each slice inherits the TIFF tags of the corresponding slice of the first split; a
  ct-rec ImageDescription ("pixel_size\tcenter\tM\ttheta0\tmin\tmax") gets new min/max.
Logging: every run is appended to n2i.log in the current directory (--log-file, --no-log).
Requires: torch, numpy, tifffile
"""
import argparse, collections, copy, faulthandler, functools, glob, itertools, json, math, os, queue, shutil, \
    subprocess, sys, threading, time
from concurrent.futures import ThreadPoolExecutor
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


def resolve_splits(a):
    """--splits as given, or auto-detected sub-directories 0, 1, 2, ... of --rec (default: rec)."""
    if a.splits:
        return
    rec = a.rec
    if not os.path.isdir(rec):
        raise SystemExit(f"split directory '{rec}' not found; give --rec or --splits")
    nums = sorted(int(d) for d in os.listdir(rec) if d.isdigit() and os.path.isdir(os.path.join(rec, d)))
    k = 0
    while k in nums:
        k += 1
    if k < 2:
        raise SystemExit(f"need sub-directories 0, 1 (, 2, 3 ...) in '{rec}'; found {nums or 'none'}")
    extra = [n for n in nums if n >= k]
    if extra:
        raise SystemExit(f"sub-directories of '{rec}' are not consecutive from 0: {nums}; "
                         f"fix them or give --splits explicitly")
    a.splits = [os.path.join(rec, str(i)) for i in range(k)]
    print(f"splits : {' '.join(a.splits)}  (K={k}, auto-detected in {rec})", flush=True)


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


def avail_ram():
    try:
        if os.name == "nt":
            import ctypes

            class MemStat(ctypes.Structure):
                _fields_ = [("dwLength", ctypes.c_ulong), ("dwMemoryLoad", ctypes.c_ulong)] + \
                           [(n, ctypes.c_ulonglong) for n in ("tot", "av", "totpf", "avpf", "totv", "avv", "avx")]
            st = MemStat()
            st.dwLength = ctypes.sizeof(MemStat)
            if ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(st)):
                return int(st.av)
        else:
            with open("/proc/meminfo") as f:
                for line in f:
                    if line.startswith("MemAvailable:"):
                        return int(line.split()[1]) * 1024
    except Exception:
        pass
    return None


def status():
    """Short resource status for the progress lines: GPU memory peak, free RAM, elapsed time."""
    parts = []
    if torch.cuda.is_available():
        parts.append(f"gpu {torch.cuda.max_memory_allocated() / 1e9:.1f}GB")
    ram = avail_ram()
    if ram is not None:
        parts.append(f"ram_free {ram / 1e9:.0f}GB")
    parts.append(f"{(time.time() - T_START) / 60:.1f}min")
    return "  ".join(parts)


class Volume:
    """K split volumes restricted to slices [z0, z1) of the stack (default: all).
    get() returns (K, len(zs), h, w) float32 with indices relative to z0, clamped."""

    def __init__(s, dirs, cache=None, lru=0, zr=None, build=False, crop=None):
        lists, ntot, s.Hfull, s.Wfull = scan_splits(dirs)
        s.z0, z1 = zr if zr else (0, ntot)
        s.crop = list(crop) if crop else [0, s.Hfull, 0, s.Wfull]     # y0, y1, x0, x1
        s.H, s.W = s.crop[1] - s.crop[0], s.crop[3] - s.crop[2]
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

    def _cache_ok(s, cache, shape, sig):
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
                d = json.load(f)
            return d.get("sig") == sig and d.get("crop", [0, s.Hfull, 0, s.Wfull]) == s.crop
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
        print(f"building cache {cache} ({need / 1e9:.1f} GB, slices {s.z0}-{s.z0 + s.N - 1}, "
              f"crop y {s.crop[0]}-{s.crop[1] - 1} x {s.crop[2]}-{s.crop[3] - 1})", flush=True)
        mm = np.lib.format.open_memmap(cache + ".tmp", mode="w+", dtype=np.float32, shape=shape)
        for k in range(s.K):
            for z in range(s.N):
                mm[k, z] = s._read_raw(k, z)
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
            json.dump({"z0": s.z0, "crop": s.crop, "shape": list(shape), "sig": sig}, f)

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
        im = tifffile.imread(s.lists[k][z])
        if s.crop != [0, s.Hfull, 0, s.Wfull]:
            im = np.ascontiguousarray(im[s.crop[0]:s.crop[1], s.crop[2]:s.crop[3]])
        return im

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


def crop_box(a, H, W):
    """Full-resolution crop [y0, y1, x0, x1] of at most --crop pixels.
    Position: --crop-center 'y,x' (pixels) if given; otherwise the crop starts at
    --crop-y0 / --crop-x0 times the image height / width from the top / left edge
    (a negative value centres the crop along that axis)."""
    c = a.crop or 0
    if c <= 0 or (c >= H and c >= W):
        return [0, H, 0, W]
    ch, cw = min(c, H), min(c, W)
    if a.crop_center:
        try:
            cy, cx = (int(float(v)) for v in a.crop_center.split(","))
        except ValueError:
            raise SystemExit("--crop-center must be 'y,x' in pixels")
        y0, x0 = cy - ch // 2, cx - cw // 2
    else:
        fy, fx = getattr(a, "crop_y0", None), getattr(a, "crop_x0", None)
        y0 = int(round(fy * H)) if fy is not None and fy >= 0 else (H - ch) // 2
        x0 = int(round(fx * W)) if fx is not None and fx >= 0 else (W - cw) // 2
    y0 = min(max(y0, 0), H - ch)
    x0 = min(max(x0, 0), W - cw)
    return [y0, y0 + ch, x0, x0 + cw]


def pick_precision(p, dev):
    """Training autocast dtype: None (fp32), torch.bfloat16 or torch.float16."""
    if dev.type != "cuda" or p == "fp32":
        return None
    if p == "auto":
        p = "bf16" if torch.cuda.get_device_capability(0)[0] >= 8 else "fp16"
    return torch.bfloat16 if p == "bf16" else torch.float16


def make_scaler(enabled):
    try:
        return torch.amp.GradScaler("cuda", enabled=enabled)
    except (AttributeError, TypeError):
        return torch.cuda.amp.GradScaler(enabled=enabled)


def footprint(mode, s, pz):
    return {"2d": 1, "25d": 2 * s + 1, "3d": pz}[mode]


# ---------------------------------------------------------------- output (ct-rec compatible tags)
# Tags written by the TIFF writer itself; everything else is copied from the source slice.
_STRUCT_TAGS = {254, 255, 256, 257, 258, 259, 262, 270, 273, 274, 277, 278, 279, 282, 283, 284, 296, 305,
                317, 322, 323, 324, 325, 338, 339, 340, 341}


def ctrec_description(desc, vmin, vmax):
    """ct-rec ImageDescription: "pixel_size\tcenter\tM\ttheta0\tmin\tmax" (printf %f).
    Replace min/max with the values of the denoised image; keep the other fields as they are.
    Returns (description, recognised)."""
    if not desc or desc.lstrip().startswith("{"):          # none, or tifffile's own JSON
        return None, False
    f = desc.split("\t")
    if len(f) == 6:
        try:
            float(f[4]), float(f[5])
            f[4], f[5] = f"{vmin:f}", f"{vmax:f}"
            return "\t".join(f), True
        except ValueError:
            pass
    return desc, False


def write_like(src, dst, img):
    """Write img (float32) to dst, inheriting the tags of the source slice src."""
    img = np.ascontiguousarray(img, dtype=np.float32)
    kw = dict(photometric="minisblack", metadata=None, rowsperstrip=img.shape[0])
    recognised = False
    try:
        with tifffile.TiffFile(src) as t:
            pg = t.pages[0]
            tags = pg.tags
            desc, recognised = ctrec_description(pg.description, float(img.min()), float(img.max()))
            if desc is not None:
                kw["description"] = desc
            if pg.rowsperstrip:
                kw["rowsperstrip"] = min(int(pg.rowsperstrip), img.shape[0])
            if "XResolution" in tags and "YResolution" in tags:
                kw["resolution"] = (tags["XResolution"].value, tags["YResolution"].value)
                if "ResolutionUnit" in tags:
                    kw["resolutionunit"] = int(tags["ResolutionUnit"].value)
            if "Software" in tags:
                kw["software"] = tags["Software"].value
            else:
                kw["software"] = False
            extra = []
            for tg in tags.values():
                if tg.code in _STRUCT_TAGS or tg.code < 254:
                    continue
                extra.append((tg.code, int(tg.dtype), tg.count, tg.value, True))
            if extra:
                kw["extratags"] = extra
    except Exception as e:
        print(f"warning: could not read tags of {src} ({e}); writing plain TIFF", flush=True)
        kw = dict(photometric="minisblack", metadata=None, rowsperstrip=img.shape[0])
    try:
        tifffile.imwrite(dst, img, **kw)
    except Exception as e:          # an exotic tag could not be copied: retry without the extra tags
        kw.pop("extratags", None)
        print(f"warning: extra tags of {src} not copied ({e})", flush=True)
        tifffile.imwrite(dst, img, **kw)
    return recognised


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
def tile_layout(P, tmax, ov, m):
    """Fewest tiles of equal size (multiple of m, <= tmax) covering P with an overlap >= ov,
    spread evenly. Returns (tile size, start positions)."""
    tmax = max(m, (int(tmax) // m) * m)
    if P <= tmax:
        return P, [0]
    ov = min(ov, tmax // 2)
    n = math.ceil((P - ov) / (tmax - ov))
    while True:
        t = ceil_to(math.ceil((P + (n - 1) * ov) / n), m)
        if t <= tmax:
            break
        n += 1
    t = min(t, P)
    return t, [round(i * (P - t) / (n - 1)) for i in range(n)]


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


class Timer:
    """Accumulates wall-clock time per stage; report() gives the share of each stage."""

    def __init__(s):
        s.t = collections.OrderedDict()

    def add(s, key, dt):
        s.t[key] = s.t.get(key, 0.0) + dt

    def report(s, total=None):
        total = total or sum(s.t.values()) or 1.0
        return "  ".join(f"{k} {v / total * 100:.0f}%" for k, v in s.t.items())

    def table(s, title, total, extra=None):
        lines = [f"--- time breakdown: {title} (total {total / 60:.1f} min) ---"]
        for k, v in s.t.items():
            lines.append(f"  {k:<14s} {v / 60:8.2f} min  {v / total * 100:5.1f} %")
        for k, v in (extra or {}).items():
            lines.append(f"  {k:<14s} {v / 60:8.2f} min  (in background threads)")
        return "\n".join(lines)


def fold_bn(model):
    """Merge every eval-mode BatchNorm into the preceding convolution (inference only).
    y = gamma * (conv(x) + b - mean) / sqrt(var + eps) + beta  ->  conv'(x) + b'
    is the same arithmetic, rearranged; the BatchNorm layers become Identity."""
    for mod in model.modules():
        if not isinstance(mod, nn.Sequential):
            continue
        layers = list(mod.children())
        for i in range(len(layers) - 1):
            c, b = layers[i], layers[i + 1]
            if isinstance(c, (nn.Conv2d, nn.Conv3d)) and isinstance(b, (nn.BatchNorm2d, nn.BatchNorm3d)):
                with torch.no_grad():
                    scale = (b.weight.double() / torch.sqrt(b.running_var.double() + b.eps))
                    shape = (-1,) + (1,) * (c.weight.dim() - 1)
                    w = c.weight.double() * scale.view(shape)
                    bias = c.bias.double() if c.bias is not None else torch.zeros_like(scale)
                    bias = (bias - b.running_mean.double()) * scale + b.bias.double()
                    c.weight.copy_(w.to(c.weight.dtype))
                    if c.bias is None:
                        c.bias = nn.Parameter(bias.to(c.weight.dtype))
                    else:
                        c.bias.copy_(bias.to(c.bias.dtype))
                mod[i + 1] = nn.Identity()
    return model


def _rep_pad_cat(x, nd):
    """Replicate padding by 1 pixel on every spatial side, built from slices and torch.cat
    (gives exactly the same values as padding_mode='replicate')."""
    for d in range(x.dim() - nd, x.dim()):
        x = torch.cat([x.narrow(d, 0, 1), x, x.narrow(d, x.size(d) - 1, 1)], dim=d)
    return x


class _RepPad(torch.autograd.Function):
    """Replicate padding by 1 pixel (same values as padding_mode='replicate') with a cheap
    hand-written backward: the gradient of the inner part plus the edge gradients added
    onto the outermost inner rows/columns (the exact adjoint of the padding)."""

    @staticmethod
    def forward(ctx, x, nd):
        ctx.nd = nd
        return _rep_pad_cat(x, nd)

    @staticmethod
    def backward(ctx, g):
        nd = ctx.nd
        for d in reversed(range(g.dim() - nd, g.dim())):   # undo the paddings in reverse order
            n = g.size(d)
            inner = g.narrow(d, 1, n - 2).clone()
            inner.narrow(d, 0, 1).add_(g.narrow(d, 0, 1))
            inner.narrow(d, n - 3, 1).add_(g.narrow(d, n - 1, 1))
            g = inner
        return g, None


def make_inputs_t(xt, mode):
    """Torch version of make_pair for all K strategies at once (same arithmetic):
    xt (K, Z, H, W) normalised -> (K, C, ...) network inputs."""
    k = xt.shape[0]
    tot = xt.sum(0)
    x = torch.stack([(tot - xt[i]) / (k - 1) for i in range(k)])
    return x.unsqueeze(1) if mode == "3d" else x


def set_padding(model, impl):
    """impl 'torch': PyTorch's own replicate padding (padding_mode='replicate').
    impl 'cat': the same values built with torch.cat, with a hand-written backward (_RepPad).
    impl 'cat-autograd': torch.cat with PyTorch's automatic backward (for bench only).
    The parameters (and state_dict keys) are unchanged, so saved models stay compatible."""
    if impl == "torch":
        return model
    for mod in model.modules():
        if isinstance(mod, (nn.Conv2d, nn.Conv3d)) and mod.padding_mode == "replicate":
            nd = len(mod.padding)
            if any(p != 1 for p in mod.padding):
                continue
            mod.padding_mode = "zeros"
            mod.padding = (0,) * nd
            if impl == "cat":
                mod.register_forward_pre_hook(lambda _m, inp, nd=nd: (_RepPad.apply(inp[0], nd),))
            else:
                mod.register_forward_pre_hook(lambda _m, inp, nd=nd: (_rep_pad_cat(inp[0], nd),))
    return model


def mem_format(layout, nd):
    if layout == "nchw":
        return torch.contiguous_format
    return torch.channels_last if nd == 2 else torch.channels_last_3d


@torch.no_grad()
def predict(model, x, tile, ov, m, dev, fp16, mf=None):
    """x: CPU tensor (B, C, *spatial), normalised; B = the K splitting strategies.
    All B inputs go through the network in one batch; returns their mean (*spatial), CPU."""
    sp = x.shape[2:]
    Pd = [ceil_to(v, m) for v in sp]
    lay = [tile_layout(p, t, o, m) for p, t, o in zip(Pd, tile, ov)]
    td = [l[0] for l in lay]
    od = [min(o, t // 2) for o, t in zip(ov, td)]
    pad = []
    for v, p in zip(reversed(sp), reversed(Pd)):
        pad += [0, p - v]
    xp = F.pad(x, pad, mode="replicate") if any(pad) else x
    # accumulate on the GPU unless the (padded) block is large (3D slabs)
    acc = dev if dev.type == "cuda" and 8 * math.prod(Pd) < 2e9 else torch.device("cpu")
    out, wsum = torch.zeros(Pd, device=acc), torch.zeros(Pd, device=acc)
    W = None
    for w in (window1d(t, o) for t, o in zip(td, od)):
        W = w if W is None else W[..., None] * w
    W = W.to(acc)
    for st in itertools.product(*[l[1] for l in lay]):
        sl = tuple(slice(s0, s0 + t) for s0, t in zip(st, td))
        xt = xp[(slice(None), slice(None)) + sl].to(dev, non_blocking=True)
        if mf is not None:
            xt = xt.contiguous(memory_format=mf)
        with torch.autocast(dev.type, dtype=torch.float16, enabled=fp16 and dev.type == "cuda"):
            y = model(xt)
        out[sl] += y[:, 0].float().mean(0).to(acc) * W
        wsum[sl] += W
    out /= wsum
    return out[tuple(slice(0, v) for v in sp)].cpu()


# ---------------------------------------------------------------- train
def train(a):
    for k, v in DEFAULTS[a.mode].items():
        if getattr(a, k, None) is None:
            setattr(a, k, v)
    if a.no_cache:
        a.cache = None
    os.makedirs(a.out, exist_ok=True)
    print(f"model output: {os.path.join(a.out, '%06d.pt')} every {a.log} iterations (--log), "
          f"{os.path.join(a.out, 'best.pt')} (lowest val)", flush=True)
    lists, Ntot, Hf, Wf = scan_splits(a.splits)
    K = len(lists)
    box = crop_box(a, Hf, Wf)
    H, W = box[1] - box[0], box[3] - box[2]
    z0, z1, why = train_range(a, K, Ntot, H, W)
    print(f"training slices {z0}-{z1 - 1} of {Ntot} ({why})", flush=True)
    print(f"training crop  y {box[0]}-{box[1] - 1}, x {box[2]}-{box[3] - 1} of {Hf}x{Wf} "
          f"({'full frame' if [H, W] == [Hf, Wf] else 'full resolution, cropped'})", flush=True)
    nd, cin, m = (3 if a.mode == "3d" else 2), (2 * a.s + 1 if a.mode == "25d" else 1), 2 ** a.depth
    a.patch = min(a.patch, H, W) // m * m
    if a.mode == "3d" and a.pz % m:
        raise SystemExit(f"--pz must be a multiple of {m}")
    L = footprint(a.mode, a.s, a.pz)
    vb = max(L, a.nval)                       # contiguous held-out z block
    if z1 - z0 < vb + L + 1:
        raise SystemExit(f"training range too small ({z1 - z0} slices, need >= {vb + L + 1}); "
                         f"free memory/disk, raise --mem-frac or set --train-slices")

    vol = Volume(a.splits, a.cache, zr=(z0, z1), build=True, crop=box)
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
    for zz in vz:
        im = (vol.get(range(zz, zz + L), vy, vy + vs, vx, vx + vs) - mean) / std
        for k in range(K):
            i_, t_ = make_pair(im, k, a.mode, a.s)
            vin.append(i_)
            vtg.append(t_)
    vin, vtg = torch.from_numpy(np.stack(vin)), torch.from_numpy(np.stack(vtg))

    dev = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    torch.backends.cudnn.benchmark = True
    mf = mem_format(a.layout, nd)
    model = set_padding(UNet(nd, cin, a.ch, a.depth, a.up), a.pad_impl).to(dev).to(memory_format=mf)
    print(f"layout: {a.layout}, padding: {a.pad_impl}", flush=True)
    dt = pick_precision(a.precision, dev)
    scaler = make_scaler(dt == torch.float16)
    print(f"precision: {'fp32' if dt is None else str(dt).replace('torch.', '') + ' mixed (fp32 weights/loss)'}"
          f"{' + loss scaling' if dt == torch.float16 else ''}", flush=True)
    opt = torch.optim.Adam(model.parameters(), lr=a.lr)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, a.iters)
    bs = max(1, a.batch // a.ppr)
    dl = DataLoader(Patches(vol, starts, a, mean, std, a.iters * bs + bs), batch_size=bs,
                    num_workers=a.workers, pin_memory=True, persistent_workers=a.workers > 0)
    amp = dict(device_type=dev.type, dtype=dt or torch.float32, enabled=dt is not None)

    def save(path):
        torch.save({"model": model.state_dict(), "mean": mean, "std": std, "mode": a.mode,
                    "s": a.s, "nd": nd, "cin": cin, "ch": a.ch, "depth": a.depth, "up": a.up,
                    "version": 4, "crop": box, "slices": [z0, z1]}, path)

    # ---- inference check: one untrained slice (centre of the validation block) at every --log
    check = None
    if not a.no_check:
        zc = a.check_slice if a.check_slice is not None else z0 + (v0 + v1) // 2
        zc = min(max(int(zc), 0), Ntot - 1)
        if a.mode == "2d":
            zs = [zc]
        elif a.mode == "25d":
            zs = list(range(zc - a.s, zc + a.s + 1))
        else:
            zs = list(range(zc - 16, zc + 17))
        ci = zs.index(zc)
        zs = [min(max(z, 0), Ntot - 1) for z in zs]
        cy0, cy1, cx0, cx1 = box if a.mode == "3d" else (0, Hf, 0, Wf)   # 3D: training crop only
        blk = np.stack([np.stack([tifffile.imread(lists[k][z])[cy0:cy1, cx0:cx1] for z in zs])
                        for k in range(K)]).astype(np.float32)
        os.makedirs(a.check_dir, exist_ok=True)
        write_like(lists[0][zc], os.path.join(a.check_dir, "average.tif"), blk[:, ci].mean(0))
        check = ((blk - mean) / std).astype(np.float32)
        del blk
        ctile = list(DEFAULTS[a.mode]["tile"])
        cov = DEFAULTS[a.mode]["ov"]
        print(f"inference check: slice {zc} -> {os.path.join(a.check_dir, '%06d.tif')} every {a.log} iterations "
              f"(reference: average.tif = mean of the {K} splits)", flush=True)

    def run_check(it):
        model.eval()
        try:
            x = make_inputs_t(torch.from_numpy(check).to(dev), a.mode)
            while True:
                try:
                    y = predict(model, x, ctile, cov, m, dev, True, mf)
                    break
                except torch.cuda.OutOfMemoryError:
                    torch.cuda.empty_cache()
                    if max(ctile[-2:]) <= 128:
                        raise
                    ctile[-2:] = [max(128, t // 2) for t in ctile[-2:]]
            out = y.numpy() * std + mean
            write_like(lists[0][zc], os.path.join(a.check_dir, f"{it:06d}.tif"), out if a.mode != "3d" else out[ci])
            return f"  [check {it:06d}.tif]"
        except Exception as e:
            torch.cuda.empty_cache() if dev.type == "cuda" else None
            return f"  [check skipped: {type(e).__name__}]"
        finally:
            model.train()

    best, run, it, micro = math.inf, 0.0, 0, 1
    tm, tm_int = Timer(), Timer()          # whole run / current --log interval
    pc = time.perf_counter
    t_train0 = t_mark = pc()

    def lap(key, t0):
        dt = pc() - t0
        tm.add(key, dt)
        tm_int.add(key, dt)
        return pc()

    model.train()
    for inp, tgt in dl:
        t_mark = lap("data wait", t_mark)
        inp = inp.flatten(0, 1).to(dev, non_blocking=True).contiguous(memory_format=mf)
        tgt = tgt.flatten(0, 1).to(dev, non_blocking=True)
        while True:   # on CUDA OOM: split the batch into micro-batches (gradient accumulation)
            try:
                opt.zero_grad(set_to_none=True)
                loss_v = 0.0
                for ci, ct in zip(inp.chunk(micro), tgt.chunk(micro)):
                    with torch.autocast(**amp):
                        pred = model(ci)
                    loss = F.mse_loss(pred.float(), ct) * (ci.shape[0] / inp.shape[0])   # L2 (N2I), fp32
                    scaler.scale(loss).backward()
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
        scaler.step(opt)
        scaler.update()
        sched.step()
        run += loss_v
        it += 1
        t_mark = lap("gpu train", t_mark)          # loss.item() has synchronised the GPU
        if it % a.log == 0 or it == a.iters:
            model.eval()
            with torch.no_grad(), torch.autocast(**amp):
                vl = float(np.mean([F.mse_loss(model(vin[j:j + 1].to(dev).contiguous(memory_format=mf)).float(),
                                               vtg[j:j + 1].to(dev)).item() for j in range(len(vin))]))
            model.train()
            t_mark = lap("validation", t_mark)
            tag = ""
            if vl < best:
                best = vl
                save(os.path.join(a.out, "best.pt"))
                tag = " *best"
            save(os.path.join(a.out, f"{it:06d}.pt"))      # every --log (and the last iteration)
            tag += f"  -> {it:06d}.pt"
            t_mark = lap("save", t_mark)
            if check is not None:
                tag += run_check(it)
                t_mark = lap("infer check", t_mark)
            print(f"it {it:6d}  train {run / a.log:.5f}  val {vl:.5f}  [{status()}]{tag}", flush=True)
            print(f"         time: {tm_int.report()}", flush=True)
            tm_int = Timer()
            run = 0.0
            t_mark = pc()
        if it >= a.iters:
            break
    print(tm.table("training", pc() - t_train0), flush=True)


# ---------------------------------------------------------------- inference
def resolve_model(a):
    """--model as given, or <train dir>/best.pt, or the newest <train dir>/%06d.pt."""
    if a.model:
        if not os.path.isfile(a.model):
            raise SystemExit(f"model '{a.model}' not found")
        return
    best = os.path.join(a.train_dir, "best.pt")
    if os.path.isfile(best):
        a.model = best
    else:
        cands = sorted(glob.glob(os.path.join(a.train_dir, "[0-9]*.pt")))
        if not cands:
            raise SystemExit(f"no model found in '{a.train_dir}' (train first, or give --model)")
        a.model = cands[-1]
    print(f"model  : {a.model}", flush=True)


def infer(a):
    resolve_model(a)
    print(f"output : {a.outdir}", flush=True)
    ck = torch.load(a.model, map_location="cpu")
    dev = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = UNet(ck["nd"], ck["cin"], ck["ch"], ck["depth"], ck["up"])
    model.load_state_dict(ck["model"])
    model.eval()
    if not a.no_fold_bn:
        ref = copy.deepcopy(model)
        fold_bn(model)
        with torch.no_grad():             # numerical check of the folding (fp32, CPU)
            shp = (2, ck["cin"]) + ((16, 32, 32) if ck["nd"] == 3 else (64, 64))
            tx = torch.randn(shp, generator=torch.Generator().manual_seed(0))
            d = (ref(tx) - model(tx)).abs().max().item()
        print(f"BatchNorm folded into the convolutions (max |difference| on a test input: {d:.1e})", flush=True)
        del ref
    set_padding(model, a.pad_impl)
    torch.backends.cudnn.benchmark = True
    mf = mem_format(a.layout, ck["nd"])
    model.to(dev).to(memory_format=mf).eval()
    print(f"layout: {a.layout}, padding: {a.pad_impl}", flush=True)
    fp16 = not a.fp32
    print(f"precision: {'fp16' if fp16 and dev.type == 'cuda' else 'fp32'} inference", flush=True)
    mode, s, mean, std, m = ck["mode"], ck["s"], ck["mean"], ck["std"], 2 ** ck["depth"]
    tile = list(a.tile or DEFAULTS[mode]["tile"])
    ov = a.ov or DEFAULTS[mode]["ov"]
    nsp = 3 if mode == "3d" else 2
    if len(tile) != nsp or len(ov) != nsp:
        raise SystemExit(f"--tile/--ov need {nsp} values for mode {mode}")
    lru = 4 * (2 * s + 1) if mode == "25d" else 0
    _, _, Hf, Wf = scan_splits(a.splits)
    box = crop_box(a, Hf, Wf)
    vol = Volume(a.splits, a.cache if box == [0, Hf, 0, Wf] else None, lru=lru * len(a.splits),
                 build=False, crop=box)
    if box != [0, Hf, 0, Wf]:
        print(f"inference crop y {box[0]}-{box[1] - 1}, x {box[2]}-{box[3] - 1} of {Hf}x{Wf} "
              f"(output images are cropped)", flush=True)
    lo, hi = slice(a.start, a.stop).indices(vol.N)[:2]
    os.makedirs(a.outdir, exist_ok=True)

    def run(blk):   # blk: raw (K, Z, H, W) block; the K strategies in one batch, averaged
        xt = blk if torch.is_tensor(blk) else torch.from_numpy(blk)
        if gpu_pre:
            xt = xt.to(dev, non_blocking=True)
        x = make_inputs_t((xt - mean) / std, mode)
        while True:   # on CUDA OOM: halve the tile and retry
            try:
                y = predict(model, x, tile, ov, m, dev, fp16, mf)
                break
            except torch.cuda.OutOfMemoryError:
                torch.cuda.empty_cache()
                if max(tile[-2:]) <= 128:
                    raise SystemExit("CUDA out of memory even with small tiles")
                tile[-2:] = [max(128, t // 2) for t in tile[-2:]]
                print(f"CUDA out of memory: tile reduced to {tile}", flush=True)
        return y.numpy() * std + mean

    warned = []

    def write(z, img):
        ok = write_like(vol.lists[0][z], os.path.join(a.outdir, os.path.basename(vol.lists[0][z])), img)
        if not ok and not warned:
            warned.append(True)
            print("note: source ImageDescription is not in the ct-rec format "
                  "(pixel_size, center, M, theta0, min, max); it was copied unchanged", flush=True)

    # jobs: (output slices, slices to read, offset of the first output slice in the block)
    if mode in ("2d", "25d"):
        jobs = [([z], [z] if mode == "2d" else list(range(z - s, z + s + 1)), 0) for z in range(lo, hi)]
    else:
        S, M = a.slab, ov[0]   # slab thickness, z-margin on each side
        jobs = [(list(range(z, min(z + S, hi))), list(range(z - M, min(z + S, hi) + M)), M)
                for z in range(lo, hi, S)]

    # input preprocessing (normalisation, the K strategy inputs, tiling) on the GPU when the
    # block fits comfortably; otherwise (large 3D slabs) on the CPU. Same arithmetic.
    blk_bytes = vol.K * len(jobs[0][1]) * vol.H * vol.W * 4 if jobs else 0
    if a.pre == "auto":
        gpu_pre = dev.type == "cuda" and 3 * blk_bytes < 0.3 * torch.cuda.get_device_properties(0).total_memory
    else:
        gpu_pre = a.pre == "gpu" and dev.type == "cuda"
    print(f"preprocessing on the {'GPU' if gpu_pre else 'CPU'}", flush=True)

    # reading (prefetch thread) -> GPU (this thread) -> writing (writer threads), overlapped
    q = queue.Queue(maxsize=max(1, a.prefetch if mode != "3d" else 1))
    stop = threading.Event()

    tm = Timer()
    bg = {"read (thread)": 0.0, "write (threads)": 0.0}
    bg_lock = threading.Lock()
    pc = time.perf_counter

    def reader():
        try:
            for job in jobs:
                if stop.is_set():
                    return
                t0 = pc()
                blk = vol.get(job[1])
                if gpu_pre:
                    blk = torch.from_numpy(blk).pin_memory()     # fast transfer to the GPU
                with bg_lock:
                    bg["read (thread)"] += pc() - t0
                q.put((job, blk))
            q.put(None)
        except BaseException as e:     # hand the error to the main thread
            q.put(e)

    def timed_write(z, img):
        t0 = pc()
        write(z, img)
        with bg_lock:
            bg["write (threads)"] += pc() - t0

    th = threading.Thread(target=reader, daemon=True)
    t_inf0 = pc()
    th.start()
    pool = ThreadPoolExecutor(max_workers=max(1, a.writers))
    pending = collections.deque()
    try:
        while True:
            t0 = pc()
            item = q.get()
            tm.add("wait read", pc() - t0)
            if item is None:
                break
            if isinstance(item, BaseException):
                raise item
            (zout, _, off), blk = item
            t0 = pc()
            out = run(blk)
            tm.add("gpu infer", pc() - t0)
            for j, z in enumerate(zout):
                img = out if mode != "3d" else out[off + j]
                pending.append(pool.submit(timed_write, z, img))
            t0 = pc()
            while len(pending) > 2 * max(1, a.writers):
                pending.popleft().result()
            tm.add("wait write", pc() - t0)
            if mode == "3d":
                print(f"slices {zout[0]}-{zout[-1]}/{hi}  [{status()}]  time: {tm.report()}", flush=True)
            elif zout[0] % 50 == 0:
                print(f"slice {zout[0]}/{hi}  [{status()}]  time: {tm.report()}", flush=True)
        t0 = pc()
        while pending:
            pending.popleft().result()
        tm.add("wait write", pc() - t0)
        total = pc() - t_inf0
        n_sl = sum(len(j[0]) for j in jobs)
        print(tm.table("inference (main thread)", total, bg), flush=True)
        print(f"  {n_sl} slices, {total / max(1, n_sl):.2f} s/slice", flush=True)
    finally:
        stop.set()
        pool.shutdown(wait=True)


# ---------------------------------------------------------------- bench
def _profile_table(prof, title, top):
    ka = prof.key_averages()
    for key in ("self_device_time_total", "self_cuda_time_total", "device_time_total", "cuda_time_total"):
        try:
            table = ka.table(sort_by=key, row_limit=top)
        except Exception:
            continue
        print(f"--- profile ({title}), top {top} by GPU time ---", flush=True)
        print(table, flush=True)
        return


def bench_train(a, cfg, state, dev):
    """Time one training step (forward + backward + Adam, mixed precision as in train) for every
    layout / padding combination, check that the gradients agree, and profile the default."""
    nd, cin = cfg[0], cfg[1]
    p = a.patch
    shape = (a.batch, cin) + ((32, p, p) if nd == 3 else (p, p))
    g = torch.Generator().manual_seed(0)
    x0 = torch.randn(shape, generator=g)
    t0_ = torch.randn((a.batch, 1) + shape[2:], generator=g)
    dt_amp = pick_precision("auto", dev)
    print(f"training step: input {tuple(shape)}, precision "
          f"{'fp32' if dt_amp is None else str(dt_amp).replace('torch.', '')} mixed, "
          f"{a.iters} timed steps after 3 warm-up steps", flush=True)

    def fresh(pad, layout):
        mdl = UNet(*cfg)
        if state is not None:
            mdl.load_state_dict(state)
        mdl = set_padding(mdl, pad).to(dev).to(memory_format=mem_format(layout, nd)).train()
        return mdl

    # gradient check in fp32 (same weights, same input; BatchNorm in training mode)
    grads = {}
    for pad in ("torch", "cat", "cat-autograd"):
        mdl = fresh(pad, "nchw")
        x = x0.to(dev)
        loss = F.mse_loss(mdl(x), t0_.to(dev))
        loss.backward()
        grads[pad] = torch.cat([q.grad.flatten() for q in mdl.parameters()]).double().cpu()
        del mdl
    gref = grads["torch"]
    for pad in ("cat", "cat-autograd"):
        rel = ((grads[pad] - gref).abs().max() / gref.abs().max()).item()
        print(f"  gradient check (fp32) padding {pad:<12s} vs torch: max relative difference {rel:.1e}",
              flush=True)

    results = []
    for layout in ("cl", "nchw"):
        for pad in ("torch", "cat", "cat-autograd"):
            mdl = fresh(pad, layout)
            mf = mem_format(layout, nd)
            opt = torch.optim.Adam(mdl.parameters(), lr=1e-4)
            scaler = make_scaler(dt_amp == torch.float16)
            amp = dict(device_type="cuda", dtype=dt_amp or torch.float32, enabled=dt_amp is not None)
            x, tg = x0.to(dev).contiguous(memory_format=mf), t0_.to(dev)

            def step():
                opt.zero_grad(set_to_none=True)
                with torch.autocast(**amp):
                    pred = mdl(x)
                loss = F.mse_loss(pred.float(), tg)
                scaler.scale(loss).backward()
                scaler.step(opt)
                scaler.update()

            for _ in range(3):
                step()
            torch.cuda.synchronize()
            t0 = time.perf_counter()
            for _ in range(a.iters):
                step()
            torch.cuda.synchronize()
            dt = (time.perf_counter() - t0) / a.iters
            results.append((layout, pad, dt))
            print(f"  layout {layout:<4s} padding {pad:<12s}: {dt * 1000:8.1f} ms/step", flush=True)
            del mdl, opt
            torch.cuda.empty_cache()
    best = min(results, key=lambda r: r[2])
    print(f"fastest: --layout {best[0]} --pad-impl {best[1]}", flush=True)

    mdl = fresh(a.pad_impl, a.layout)
    mf = mem_format(a.layout, nd)
    opt = torch.optim.Adam(mdl.parameters(), lr=1e-4)
    scaler = make_scaler(dt_amp == torch.float16)
    amp = dict(device_type="cuda", dtype=dt_amp or torch.float32, enabled=dt_amp is not None)
    x, tg = x0.to(dev).contiguous(memory_format=mf), t0_.to(dev)
    from torch.profiler import profile, ProfilerActivity
    for _ in range(2):
        opt.zero_grad(set_to_none=True)
        with torch.autocast(**amp):
            loss = F.mse_loss(mdl(x).float(), tg)
        scaler.scale(loss).backward()
        scaler.step(opt)
        scaler.update()
    torch.cuda.synchronize()
    with profile(activities=[ProfilerActivity.CPU, ProfilerActivity.CUDA]) as prof:
        for _ in range(3):
            opt.zero_grad(set_to_none=True)
            with torch.autocast(**amp):
                loss = F.mse_loss(mdl(x).float(), tg)
            scaler.scale(loss).backward()
            scaler.step(opt)
            scaler.update()
        torch.cuda.synchronize()
    _profile_table(prof, f"training, --layout {a.layout} --pad-impl {a.pad_impl}, 3 steps", a.top)


def bench(a):
    """Time the network on one tile for every combination of memory layout and padding
    implementation (all give the same result), and list the slowest GPU kernels.
    With --train: the same for one training step (forward + backward)."""
    dev = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if dev.type != "cuda":
        raise SystemExit("bench needs a CUDA GPU")
    try:
        resolve_model(a)
        ck = torch.load(a.model, map_location="cpu")
        cfg = (ck["nd"], ck["cin"], ck["ch"], ck["depth"], ck["up"])
        state = ck["model"]
    except SystemExit:
        print("no trained model found: using a randomly initialised 2.5D network", flush=True)
        cfg, state = (2, 5, 32, 4, "resize"), None
    nd, cin = cfg[0], cfg[1]
    torch.backends.cudnn.benchmark = True
    if a.batch is None:
        a.batch = 16 if a.train else 2
    if a.layout is None:
        a.layout = "nchw" if a.train else "cl"
    if a.pad_impl is None:
        a.pad_impl = "torch" if a.train else "cat"
    if a.train:
        bench_train(a, cfg, state, dev)
        return
    base = UNet(*cfg)
    if state is not None:
        base.load_state_dict(state)
    base.eval()
    fold_bn(base)
    torch.backends.cudnn.benchmark = True
    t = a.tile
    shape = (a.batch, cin) + ((max(16, t // 16), t, t) if nd == 3 else (t, t))
    x0 = torch.randn(shape, generator=torch.Generator().manual_seed(0))
    npx = a.batch * math.prod(shape[2:])
    print(f"input {tuple(shape)} ({'fp32' if a.fp32 else 'fp16'}), {a.iters} timed runs after 3 warm-up runs", flush=True)

    results, ref = [], None
    for layout in ("cl", "nchw"):
        for pad in ("torch", "cat"):
            model = set_padding(copy.deepcopy(base), pad)
            mf = mem_format(layout, nd)
            model.to(dev).to(memory_format=mf)
            x = x0.to(dev).contiguous(memory_format=mf)
            with torch.no_grad(), torch.autocast("cuda", dtype=torch.float16, enabled=not a.fp32):
                for _ in range(3):
                    y = model(x)
                torch.cuda.synchronize()
                t0 = time.perf_counter()
                for _ in range(a.iters):
                    y = model(x)
                torch.cuda.synchronize()
                dt = (time.perf_counter() - t0) / a.iters
            y = y.float().cpu()
            diff = 0.0 if ref is None else (y - ref).abs().max().item()
            ref = y if ref is None else ref
            results.append((layout, pad, dt))
            print(f"  layout {layout:<4s} padding {pad:<5s}: {dt * 1000:8.1f} ms/run  "
                  f"{npx / dt / 1e6:7.1f} Mpx/s   max |diff| vs first {diff:.1e}", flush=True)
            del model
            torch.cuda.empty_cache()
    best = min(results, key=lambda r: r[2])
    print(f"fastest: --layout {best[0]} --pad-impl {best[1]}", flush=True)

    # profile the current default combination
    model = set_padding(copy.deepcopy(base), a.pad_impl).to(dev).to(memory_format=mem_format(a.layout, nd))
    x = x0.to(dev).contiguous(memory_format=mem_format(a.layout, nd))
    from torch.profiler import profile, ProfilerActivity
    with torch.no_grad(), torch.autocast("cuda", dtype=torch.float16, enabled=not a.fp32):
        for _ in range(2):
            model(x)
        torch.cuda.synchronize()
        with profile(activities=[ProfilerActivity.CPU, ProfilerActivity.CUDA]) as prof:
            for _ in range(3):
                model(x)
            torch.cuda.synchronize()
    _profile_table(prof, f"inference, --layout {a.layout} --pad-impl {a.pad_impl}, 3 runs", a.top)


# ---------------------------------------------------------------- logging
T_START = time.time()
_LOG_FILES = []


class _Tee:
    def __init__(s, stream, f):
        s.stream, s.f = stream, f

    def write(s, x):
        s.stream.write(x)
        s.f.write(x)
        return len(x)

    def flush(s):
        s.stream.flush()
        s.f.flush()

    def isatty(s):
        return s.stream.isatty()


def setup_log(a):
    if a.no_log:
        return None
    path = a.log_file or "n2i.log"
    f = open(path, "a", encoding="utf-8", buffering=1)
    fh = open(path, "a")                      # separate raw handle for native crash dumps
    _LOG_FILES.extend([f, fh])
    faulthandler.enable(file=fh, all_threads=True)
    sys.stdout, sys.stderr = _Tee(sys.stdout, f), _Tee(sys.stderr, f)
    f.write("\n" + "=" * 78 + "\n")             # separator between runs (log file only)
    print(f"=== {os.path.basename(__file__)} {a.cmd}  {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"command : {subprocess.list2cmdline(sys.argv)}")
    print(f"cwd     : {os.getcwd()}")
    print(f"python {sys.version.split()[0]} | torch {torch.__version__} (CUDA {torch.version.cuda}) | "
          f"numpy {np.__version__} | tifffile {tifffile.__version__}")
    if torch.cuda.is_available():
        pr = torch.cuda.get_device_properties(0)
        print(f"gpu     : {pr.name}, {pr.total_memory / 1e9:.1f} GB")
    ram = total_ram()
    if ram:
        print(f"ram     : {ram / 1e9:.0f} GB")
    print(f"args    : {vars(a)}", flush=True)
    return path


# ---------------------------------------------------------------- cli
if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        prog=os.path.basename(__file__), formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Noise2Inverse denoising of CT reconstructions (2D / 2.5D / 3D U-Net).",
        epilog="""typical use (in the data directory, splits in rec\\0, rec\\1 (, rec\\2, rec\\3)):
  %(prog)s train          train; models are written to training\\ (000100.pt ..., best.pt)
  %(prog)s infer          denoise all slices with training\\best.pt into denoised\\

details of each command:
  %(prog)s train -h
  %(prog)s infer -h
  %(prog)s bench          time the network on the GPU (memory layout / padding variants)
  %(prog)s bench --train  the same for one training step""")
    sp = ap.add_subparsers(dest="cmd", metavar="{train,infer}")
    for name, hlp in (("train", "train a denoising model"), ("infer", "denoise slices with a trained model")):
        q = sp.add_parser(name, help=hlp)
        q.add_argument("--splits", nargs="+", default=None,
                       help="split directories (default: sub-directories 0, 1, ... of --rec)")
        q.add_argument("--rec", default="rec",
                       help="directory containing the split sub-directories 0, 1 (, 2, 3 ...) (default: rec)")
        q.add_argument("--cache", default=None, help="raw float32 memmap (.npy) of the splits (train: training range only)")
        q.add_argument("--start", type=int, default=None)
        q.add_argument("--stop", type=int, default=None)
        q.add_argument("--log-file", default=None, help="log file, appended to (default: n2i.log in the current directory)")
        q.add_argument("--no-log", action="store_true", help="do not write a log file")
    # fastest combinations measured with "bench" (Quadro RTX 8000): training nchw + torch
    # (BatchNorm in training mode is very slow with channels_last), inference cl + cat
    for q, lay, pad in ((sp.choices["train"], "nchw", "torch"), (sp.choices["infer"], "cl", "cat")):
        q.add_argument("--layout", choices=["cl", "nchw"], default=lay,
                       help=f"memory layout on the GPU: cl = channels_last, nchw = standard (default: {lay})")
        q.add_argument("--pad-impl", choices=["torch", "cat"], default=pad,
                       help=f"how the replicate padding is done, same values: torch or cat (default: {pad})")
    b = sp.add_parser("bench", help="time the network for each layout / padding implementation")
    b.add_argument("--model", default=None, help="model file (default: <train-dir>/best.pt ...)")
    b.add_argument("--train-dir", default="training")
    b.add_argument("--tile", type=int, default=1024, help="tile size (pixels)")
    b.add_argument("--batch", type=int, default=None,
                   help="batch (default: 2 = K strategies for inference, 16 for --train)")
    b.add_argument("--train", action="store_true", help="time a training step instead of inference")
    b.add_argument("--patch", type=int, default=256, help="--train: patch size (default 256)")
    b.add_argument("--iters", type=int, default=10, help="timed runs per combination")
    b.add_argument("--fp32", action="store_true", help="fp32 instead of fp16")
    b.add_argument("--layout", choices=["cl", "nchw"], default=None,
                   help="combination to profile (default: that of train / infer)")
    b.add_argument("--pad-impl", choices=["torch", "cat", "cat-autograd"], default=None,
                   help="combination to profile (default: that of train / infer)")
    b.add_argument("--top", type=int, default=15, help="number of kernels listed in the profile")
    b.add_argument("--log-file", default=None)
    b.add_argument("--no-log", action="store_true")
    t = sp.choices["train"]
    t.set_defaults(cache="image_cache.npy")
    t.add_argument("--no-cache", action="store_true", help="train without a cache (read TIFF files)")
    t.add_argument("--mode", choices=["2d", "25d", "3d"], default="25d")
    t.add_argument("--s", type=int, default=2, help="2.5D half width (2s+1 slices)")
    t.add_argument("--pz", type=int, default=32, help="3D patch depth")
    t.add_argument("--patch", type=int, default=None)
    t.add_argument("--batch", type=int, default=None)
    t.add_argument("--ppr", type=int, default=None, help="patches per read")
    t.add_argument("--depth", type=int, default=None)
    t.add_argument("--ch", type=int, default=32)
    t.add_argument("--up", choices=["resize", "transpose"], default="resize")
    t.add_argument("--out", default="training",
                   help="output directory for the models: %%06d.pt (iteration) and best.pt (default: training)")
    t.add_argument("--iters", type=int, default=2000)
    t.add_argument("--lr", type=float, default=5e-4)
    t.add_argument("--nval", type=int, default=8)
    t.add_argument("--workers", type=int, default=4)
    t.add_argument("--log", type=int, default=100,
                   help="interval of loss/validation output, model saving (%%06d.pt) and the inference check")
    t.add_argument("--precision", choices=["auto", "bf16", "fp16", "fp32"], default="auto",
                   help="mixed precision for training (auto: bf16 on Ampere or newer, else fp16 + loss scaling)")
    t.add_argument("--amp", action="store_true", help=argparse.SUPPRESS)   # v2 option, now the default
    t.add_argument("--crop", type=int, default=1024,
                   help="train on a full-resolution crop of at most this many pixels (0 = full frame)")
    t.add_argument("--crop-y0", type=float, default=0.25,
                   help="crop start below the top edge as a fraction of the image height "
                        "(default 0.25; negative = centred)")
    t.add_argument("--crop-x0", type=float, default=-1.0,
                   help="crop start right of the left edge as a fraction of the image width "
                        "(default: centred)")
    t.add_argument("--crop-center", default=None,
                   help="crop centre 'y,x' in pixels (overrides --crop-y0/--crop-x0)")
    t.add_argument("--check-dir", default="infer_check",
                   help="at every --log, denoise one slice and write it here as %%06d.tif (default: infer_check)")
    t.add_argument("--check-slice", type=int, default=None,
                   help="stack index of the check slice (default: centre of the validation block)")
    t.add_argument("--no-check", action="store_true", help="disable the inference check")
    t.add_argument("--mem-frac", type=float, default=0.7,
                   help="with --cache and no --start/--stop: limit the cache to this fraction of RAM")
    t.add_argument("--train-slices", type=int, default=None,
                   help="with --cache and no --start/--stop: number of central slices to train on")
    t.add_argument("--max-train-slices", type=int, default=512,
                   help="with --cache and no --start/--stop: upper limit of the automatically "
                        "selected central slices (0 = no limit)")
    f = sp.choices["infer"]
    f.add_argument("--model", default=None,
                   help="model file (default: <train-dir>/best.pt, or the newest <train-dir>/%%06d.pt)")
    f.add_argument("--train-dir", default="training", help="directory searched for the model (default: training)")
    f.add_argument("--outdir", default="denoised", help="output directory (default: denoised)")
    f.add_argument("--tile", type=int, nargs="+", default=None, help="(y x) or (z y x)")
    f.add_argument("--ov", type=int, nargs="+", default=None, help="tile overlap, same order")
    f.add_argument("--slab", type=int, default=64, help="3D: output slices per slab")
    f.add_argument("--fp16", action="store_true", help=argparse.SUPPRESS)   # v2 option, now the default
    f.add_argument("--fp32", action="store_true", help="fp32 inference (default: fp16)")
    f.add_argument("--pre", choices=["auto", "gpu", "cpu"], default="auto",
                   help="where the inputs are prepared (normalisation, strategy inputs, tiling): "
                        "auto (GPU unless the block is very large), gpu or cpu")
    f.add_argument("--no-fold-bn", action="store_true",
                   help="do not fold BatchNorm into the convolutions (for comparison)")
    f.add_argument("--prefetch", type=int, default=3, help="slices read ahead (2D/2.5D)")
    f.add_argument("--writers", type=int, default=2, help="threads writing output files")
    f.add_argument("--crop", type=int, default=0,
                   help="process only a crop of at most this many pixels (0 = full frame, default)")
    f.add_argument("--crop-center", default=None, help="crop centre 'y,x' in pixels (default: image centre)")
    a = ap.parse_args()
    if a.cmd is None:          # no command: show the help
        ap.print_help()
        sys.exit(0)
    log = setup_log(a)
    try:
        if a.cmd == "bench":
            bench(a)
        else:
            resolve_splits(a)
            train(a) if a.cmd == "train" else infer(a)
        print(f"=== finished normally ({(time.time() - T_START) / 60:.1f} min)", flush=True)
    except KeyboardInterrupt:
        print("=== interrupted by user (Ctrl+C)", flush=True)
        sys.exit(130)
    except SystemExit as e:
        if e.code not in (None, 0):
            print(f"=== stopped: {e.code}", flush=True)
        raise
    except BaseException:
        import traceback
        traceback.print_exc()
        print("=== terminated by an error (traceback above)", flush=True)
        sys.exit(1)
