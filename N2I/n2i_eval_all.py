#!/usr/bin/env python3
"""
n2i_eval_all.py -- one script for the whole evaluation of projection filtering and Noise2Inverse:
computation (filtering, split reconstruction, training, denoising) and analysis (tables, graphs).

Run it in the data directory (where raw/ with a.his, conv.bat, output.log is) with the Python of
the n2i installation; on a machine where N2I runs nothing else is needed:
  ~/data/ct-rec/ct-rec/N2I/n2i_portable/python/bin/python3 n2i_eval_all.py
  The pixel size (--dr) and the rotation centre (--rc) are asked for when not given (in a run
  without a keyboard, e.g. nohup or a batch job, they must be given as options); the in-plane
  angle --ra0 defaults to 0.0. They are kept in <out>/params.json and must match on a restart.
  GPUs: --gpus 0,1

For every combination of
    median kernel   (--medians, default 1 3 5)          } his_mgf on a.his
    gaussian sigma  (--sigmas,  default 0 0.5 1.0 2.0)  }   -> 12 filtered his files
    reconstruction  (--progs,   default his_srec_g_c_n2i his_srec_g_r_n2i: Chesler / Ramachandran filter)
    split count K   (--splits,  default 2 4)              -> 48 cases
it
    1. filters a.his with his_mgf (<out>/his/m<M>_g<S>/raw with conv.bat, output.log)
    2. reconstructs the K splits    (<out>/cases/<case>/rec/0 ... K-1)
    3. trains n2i_v4.py and denoises the central slice
    4. writes <out>/summary/:
         summary_center_slice.tif, summary.csv   ImageJ stack of the central slices (--no-stack: skip)
         analysis_report.html                    tables, graphs and images of the analysis
         analysis_metrics.csv, transfer_curves.csv
The analysis (in a --roi-size region of the central slice) measures the independent noise before
and after denoising and the signal transfer (resolution, f50) from two independent halves of the
projections, so components common to all splits do not bias it; a normal reconstruction of the
central slice (raw/rec<slice>.tif, made with 'his_rec_g_c <layer> <Rc> <Dr> <RA0>' in raw/ when missing,
or --ref) is compared with every case.

The run can be stopped and restarted: finished steps are skipped (marker files .done_*).
--analyze-only redoes only the analysis of the finished cases.
"""
import argparse, csv, glob, json, os, queue, re, shlex, shutil, subprocess, sys, threading, time
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import tifffile

LOCK = threading.Lock()
LOGF = None


def log(msg):
    line = f"{time.strftime('%Y-%m-%d %H:%M:%S')}  {msg}"
    with LOCK:
        print(line, flush=True)
        if LOGF:
            LOGF.write(line + "\n")
            LOGF.flush()


def run(cmd, cwd, logfile, env=None, dry=False):
    """Run a command, appending its output to logfile; True on success."""
    shown = " ".join(str(c) for c in cmd)
    if dry:
        log(f"[dry-run] (cd {cwd}) {shown}")
        return True
    os.makedirs(cwd, exist_ok=True)
    with open(logfile, "a", encoding="utf-8", errors="replace") as f:
        f.write(f"\n===== {time.strftime('%Y-%m-%d %H:%M:%S')}  {shown}\n")
        f.flush()
        p = subprocess.run([str(c) for c in cmd], cwd=cwd, stdout=f, stderr=subprocess.STDOUT, env=env)
    return p.returncode == 0


def done(path, step):
    return os.path.exists(os.path.join(path, f".done_{step}"))


def mark(path, step):
    with open(os.path.join(path, f".done_{step}"), "w") as f:
        f.write(time.strftime("%Y-%m-%d %H:%M:%S\n"))


N2I_PLACES = ["~/data/ct-rec/ct-rec/N2I/n2i_portable", "~/n2i", "C:\\n2i"]


def resolve_n2i(arg):
    """Find n2i_v4.py and the Python to run it with. arg may be n2i_v4.py, the launcher (n2i,
    n2i.bat) or the installation directory. Returns (n2i_v4.py, python) or (None, None)."""
    cands = []
    if arg:
        p = os.path.abspath(os.path.expanduser(arg))
        cands.append(os.path.dirname(p) if os.path.isfile(p) else p)
    cands.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.executable)))))
    if not arg:
        cands += [os.path.expanduser(d) for d in N2I_PLACES]
    for d in cands:
        script = os.path.join(d, "n2i_v4.py")
        if os.path.isfile(script):
            for py in (os.path.join(d, "python", "bin", "python3"),        # portable bundle
                       os.path.join(d, "venv", "bin", "python"),           # setup_n2i.sh
                       os.path.join(d, "venv", "Scripts", "python.exe")):  # setup_n2i.bat
                if os.path.isfile(py):
                    return script, py
            return script, sys.executable
    return None, None


def prog_label(prog):
    m = re.search(r"_([A-Za-z0-9]+)_n2i$", os.path.splitext(os.path.basename(prog))[0])
    return m.group(1) if m else re.sub(r"[^A-Za-z0-9]+", "", os.path.basename(prog))


def tif_list(d):
    return sorted(glob.glob(os.path.join(d, "*.tif*")))


# ------------------------------------------------------------------ steps
def step_filter(a, m, s, hisdir):
    """his_mgf a.his -> <hisdir>/raw/a.his, plus copies of the companion files."""
    raw = os.path.join(hisdir, "raw")
    if done(hisdir, "filter"):
        return True
    os.makedirs(raw, exist_ok=True)
    for f in a.copy_files:
        src = os.path.join(a.raw, f)
        if not os.path.exists(src):
            log(f"ERROR {os.path.basename(hisdir)}: {src} not found")
            return False
        if not a.dry_run:
            shutil.copy2(src, os.path.join(raw, f))
    out_his = os.path.join(raw, a.his)
    if os.path.exists(out_his) and not a.dry_run:
        os.remove(out_his)        # left by an interrupted run; his_mgf does not overwrite its output
    log(f"filter   {os.path.basename(hisdir)}: his_mgf median {m} sigma {s}")
    ok = run([a.mgf, os.path.abspath(os.path.join(a.raw, a.his)), os.path.abspath(os.path.join(raw, a.his)),
              m, s], hisdir, os.path.join(hisdir, "his_mgf.log"), dry=a.dry_run)
    if ok and not a.dry_run:
        if not os.path.exists(os.path.join(raw, a.his)):
            log(f"ERROR {os.path.basename(hisdir)}: his_mgf did not write {a.his}")
            return False
        mark(hisdir, "filter")
    if not ok:
        log(f"ERROR {os.path.basename(hisdir)}: his_mgf failed (see his_mgf.log)")
    return ok


def case_pipeline(a, case, gpu):
    """reconstruction -> training -> inference of the central slice, for one case."""
    name, hisdir, prog, k = case["name"], case["hisdir"], case["prog"], case["k"]
    cdir = case["dir"]
    os.makedirs(cdir, exist_ok=True)
    env = dict(os.environ)
    if gpu is not None:
        env["CUDA_VISIBLE_DEVICES"] = str(gpu)
    tag = f"{name} [gpu {gpu}]" if gpu is not None else name

    # 1. split reconstruction
    if not done(cdir, "recon"):
        rec = os.path.join(cdir, "rec")
        if os.path.isdir(rec) and not a.dry_run:
            shutil.rmtree(rec)          # unfinished earlier attempt
        try:
            raw_rel = os.path.relpath(os.path.join(hisdir, "raw"), cdir)
        except ValueError:                      # different drives on Windows
            raw_rel = os.path.join(hisdir, "raw")
        log(f"recon    {tag}: {prog} {raw_rel} {a.dr} {a.rc} {a.ra0} rec {k}")
        t0 = time.time()
        if not run([prog, raw_rel, a.dr, a.rc, a.ra0, "rec", k], cdir, os.path.join(cdir, "recon.log"),
                   env, a.dry_run):
            log(f"ERROR {name}: reconstruction failed (see recon.log)")
            return False
        if not a.dry_run:
            counts = [len(tif_list(os.path.join(rec, str(i)))) for i in range(k)]
            if min(counts) == 0 or len(set(counts)) != 1:
                log(f"ERROR {name}: unexpected reconstruction output, slices per split {counts}")
                return False
            mark(cdir, "recon")
            log(f"recon    {tag}: done ({counts[0]} slices x {k}, {(time.time() - t0) / 60:.1f} min)")

    # 2. training
    if not done(cdir, "train"):
        log(f"train    {tag}")
        t0 = time.time()
        cmd = [a.n2i_py, a.n2i, "train", "--rec", "rec", "--out", "training", "--check-dir", "infer_check"]
        cmd += shlex.split(a.train_args, posix=os.name != "nt")
        if not run(cmd, cdir, os.path.join(cdir, "run.log"), env, a.dry_run):
            log(f"ERROR {name}: training failed (see run.log / n2i.log in {cdir})")
            return False
        if not a.dry_run:
            if not os.path.exists(os.path.join(cdir, "training", "best.pt")):
                log(f"ERROR {name}: training produced no model")
                return False
            mark(cdir, "train")
            if not a.keep_cache:
                for f in ("image_cache.npy", "image_cache.npy.json"):
                    p = os.path.join(cdir, f)
                    if os.path.exists(p):
                        os.remove(p)
            log(f"train    {tag}: done ({(time.time() - t0) / 60:.1f} min)")

    # 3. inference of the central slice(s)
    if not done(cdir, "infer"):
        files = tif_list(os.path.join(cdir, "rec", "0")) if not a.dry_run else []
        n = len(files) if files else 1
        c = n // 2
        lo, hi = max(0, c - a.infer_half), min(n, c + a.infer_half + 1)
        log(f"infer    {tag}: slices {lo}-{hi - 1} of {n}")
        cmd = [a.n2i_py, a.n2i, "infer", "--rec", "rec", "--train-dir", "training",
               "--outdir", "denoised", "--start", lo, "--stop", hi] + shlex.split(a.infer_args, posix=os.name != "nt")
        if not run(cmd, cdir, os.path.join(cdir, "run.log"), env, a.dry_run):
            log(f"ERROR {name}: inference failed (see run.log / n2i.log in {cdir})")
            return False
        if not a.dry_run:
            out = os.path.join(cdir, "denoised", os.path.basename(files[c]))
            if not os.path.exists(out):
                log(f"ERROR {name}: {out} not written")
                return False
            mark(cdir, "infer")

    if a.delete_rec and not a.dry_run:
        prune_rec(a, cdir, k)
    log(f"finished {tag}")
    return True


def prune_rec(a, cdir, k):
    """--delete-rec: delete the split reconstructions except the central slice(s), which the
    summary stack and the analysis need (the kept slices are the denoised ones)."""
    files = tif_list(os.path.join(cdir, "rec", "0"))
    n = len(files)
    if n <= 2 * a.infer_half + 1:
        return                                   # already pruned
    c = n // 2
    keep = {os.path.basename(files[z]) for z in range(max(0, c - a.infer_half), min(n, c + a.infer_half + 1))}
    for i in range(k):
        for f in tif_list(os.path.join(cdir, "rec", str(i))):
            if os.path.basename(f) not in keep:
                os.remove(f)


# ------------------------------------------------------------------ summary
def summary(a, cases):
    """ImageJ stack of the central slice of every case: denoised first, then the plain average
    of the splits. Written slice by slice, so it needs little memory."""
    items, shape = [], None
    for c in cases:
        cdir = c["dir"]
        rec0 = tif_list(os.path.join(cdir, "rec", "0"))
        d_path, parts = None, None
        if done(cdir, "infer") and rec0:
            fname = os.path.basename(rec0[len(rec0) // 2])
            p = os.path.join(cdir, "denoised", fname)
            if os.path.exists(p):
                d_path = p
                if shape is None:
                    shape = tifffile.imread(p).shape
            pp = [os.path.join(cdir, "rec", str(i), fname) for i in range(c["k"])]
            if all(os.path.exists(q) for q in pp):
                parts = pp
        items.append((c, d_path, parts))
    if shape is None:
        log("summary: no finished case yet")
        return

    def load(part, d_path, parts):
        try:
            if part == "denoised" and d_path:
                img = tifffile.imread(d_path).astype(np.float32)
            elif part == "average" and parts:
                img = np.mean([tifffile.imread(q).astype(np.float32) for q in parts], axis=0).astype(np.float32)
            else:
                return None
            return img if img.shape == shape else None
        except Exception:
            return None

    plan = [(part, c, d, pp) for part in ("denoised", "average") for c, d, pp in items]
    labels, table = [], []

    def frames():
        for n, (part, c, d, pp) in enumerate(plan, 1):
            img = load(part, d, pp)
            ok = img is not None
            lab = f"{part}_{c['name']}"
            labels.append(lab if ok else lab + "_MISSING")
            table.append([n, part, c["m"], c["s"], c["plabel"], c["prog"], c["k"], "ok" if ok else "missing", c["name"]])
            yield img if ok else np.zeros(shape, np.float32)

    sdir = os.path.join(a.out, "summary")
    os.makedirs(sdir, exist_ok=True)
    tif = os.path.join(sdir, "summary_center_slice.tif")
    # the labels are known only after the frames were read: read once for the labels, then write
    for _ in frames():
        pass
    final_labels = list(labels)
    labels.clear()
    table.clear()
    tifffile.imwrite(tif, frames(), shape=(len(plan),) + tuple(shape), dtype=np.float32,
                     imagej=True, metadata={"Labels": final_labels})
    with open(os.path.join(sdir, "summary.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["slice", "type", "median", "sigma", "recon", "program", "K", "status", "case"])
        w.writerows(table)
    nok = sum(1 for t in table if t[7] == "ok")
    log(f"summary: {tif} ({len(plan)} slices, {nok} present) and summary.csv")


# ------------------------------------------------------------------ analysis
# Quantities that are not affected by components common to all splits (sample structure,
# reconstruction-specific patterns):
#   independent noise : from the difference of two independent halves of the projections
#   signal transfer   : cross-spectrum of the two halves (noise-free signal power), relative to a
#                       reference case (unfiltered, sharpest reconstruction)
import base64, html, math, struct, zlib

MED_COLORS = ["#1f77b4", "#e8590c", "#2b8a3e", "#9c36b5", "#c92a2a"]


def _halves(sp):
    k = len(sp)
    if k == 2:
        return sp[0], sp[1]
    return np.mean(sp[0::2], axis=0), np.mean(sp[1::2], axis=0)


def _gauss(img, sigma):
    if sigma <= 0:
        return img
    h, w = img.shape
    fy = np.fft.fftfreq(h)[:, None]
    fx = np.fft.rfftfreq(w)[None, :]
    return np.fft.irfft2(np.fft.rfft2(img) * np.exp(-2 * (np.pi * sigma) ** 2 * (fy ** 2 + fx ** 2)), s=img.shape)


def _ssim(x, y, sigma=1.5):
    lo, hi = np.percentile(y, [0.1, 99.9])
    L = max(hi - lo, 1e-12)
    c1, c2 = (0.01 * L) ** 2, (0.03 * L) ** 2
    mx, my = _gauss(x, sigma), _gauss(y, sigma)
    sxx = _gauss(x * x, sigma) - mx * mx
    syy = _gauss(y * y, sigma) - my * my
    sxy = _gauss(x * y, sigma) - mx * my
    m = ((2 * mx * my + c1) * (2 * sxy + c2)) / ((mx * mx + my * my + c1) * (sxx + syy + c2))
    return float(m.mean())


class Spectra:
    """Radially averaged spectra on a square ROI (Hann window)."""

    def __init__(s, n, nbins=25):
        s.n = n
        s.w = np.outer(np.hanning(n), np.hanning(n))
        fr = np.hypot(np.fft.fftfreq(n)[:, None], np.fft.fftfreq(n)[None, :])
        s.edges = np.linspace(0, 0.5, nbins + 1)
        s.fc = (s.edges[:-1] + s.edges[1:]) / 2
        s.idx = np.digitize(fr.ravel(), s.edges) - 1
        s.ok = (s.idx >= 0) & (s.idx < nbins)
        s.cnt = np.bincount(s.idx[s.ok], minlength=nbins).astype(float)
        s.nb = nbins

    def F(s, x):
        return np.fft.fft2((x - x.mean()) * s.w)

    def radial(s, P):
        return np.bincount(s.idx[s.ok], P.ravel()[s.ok], minlength=s.nb) / np.maximum(s.cnt, 1)


def _f_at(fc, H, level=0.5):
    """Frequency where H first drops below level (linear interpolation); 0.5 if it never does."""
    for j in range(1, len(fc)):
        if fc[j] > 0.03 and H[j] < level:
            if H[j - 1] <= level or H[j - 1] == H[j]:
                return float(fc[j - 1])
            return float(fc[j - 1] + (level - H[j - 1]) * (fc[j] - fc[j - 1]) / (H[j] - H[j - 1]))
    return 0.5


def analyze_case(a, c, roi, S, refimg):
    """Spectra and noise figures of one case, accumulated over the denoised central slices."""
    name, k, cdir = c["name"], c["k"], c["dir"]
    if k % 2:
        log(f"analysis: {name}: K = {k} is odd, the two halves of the projections are not equivalent: skipped")
        return None
    files = tif_list(os.path.join(cdir, "rec", "0"))
    if not files:
        log(f"analysis: {name}: split reconstructions not found: skipped")
        return None
    n = len(files)
    zc = n // 2
    zs = [z for z in range(zc - a.infer_half, zc + a.infer_half + 1)
          if 0 <= z < n and os.path.exists(os.path.join(cdir, "denoised", os.path.basename(files[z])))]
    if not zs:
        return None
    y0, y1, x0, x1 = roi
    acc = {key: np.zeros(S.nb) for key in ("PS", "Pden", "Pn", "Pres", "Paa", "Pbb")}
    var_n = var_r = 0.0
    thumbs = None
    ssim_v, bias_v = [], []
    for z in zs:
        fname = os.path.basename(files[z])
        sp = [tifffile.imread(os.path.join(cdir, "rec", str(i), fname)).astype(np.float64)[y0:y1, x0:x1]
              for i in range(k)]
        den = tifffile.imread(os.path.join(cdir, "denoised", fname)).astype(np.float64)[y0:y1, x0:x1]
        ha, hb = _halves(sp)
        avg = np.mean(sp, axis=0)
        noise = (ha - hb) / 2                      # independent noise, same variance as the average's
        res = avg - den
        Fa, Fb = S.F(ha), S.F(hb)
        acc["PS"] += S.radial((Fa * np.conj(Fb)).real)
        acc["Paa"] += S.radial(np.abs(Fa) ** 2)
        acc["Pbb"] += S.radial(np.abs(Fb) ** 2)
        acc["Pden"] += S.radial(np.abs(S.F(den)) ** 2)
        acc["Pn"] += S.radial(np.abs(S.F(noise)) ** 2)
        acc["Pres"] += S.radial(np.abs(S.F(res)) ** 2)
        var_n += noise.var()
        var_r += res.var()
        if refimg is not None and z == zc:
            r = refimg[y0:y1, x0:x1]
            ssim_v.append(_ssim(_gauss(den, 1.0), _gauss(r, 1.0)))
            lp = np.abs(_gauss(den, 4.0) - _gauss(r, 4.0)).mean()
            bias_v.append(lp / max(np.abs(r).mean(), 1e-30))
        if z == zc:
            thumbs = (avg.astype(np.float32), den.astype(np.float32), res.astype(np.float32))
    m = len(zs)
    for key in acc:
        acc[key] /= m
    var_n /= m
    var_r /= m
    coh = acc["PS"] / np.sqrt(np.maximum(acc["Paa"] * acc["Pbb"], 1e-30))
    hb_mask = S.fc >= 0.4
    bands = ((0, 0.08), (0.08, 0.2), (0.2, 0.35), (0.35, 0.5))
    bratio = []
    for lo_, hi_ in bands:
        sel = (S.fc >= lo_) & (S.fc < hi_)
        bratio.append(float(acc["Pres"][sel].sum() / max(acc["Pn"][sel].sum(), 1e-30)))
    out = dict(case=name, median=c["m"], sigma=c["s"], recon=c["plabel"], K=k, slices=len(zs),
               noise_avg=math.sqrt(var_n), removed=math.sqrt(var_r),
               noise_left=math.sqrt(max(var_n - var_r, 0.0)),
               rr_low=bratio[0], rr_mid=bratio[1], rr_high=bratio[2], rr_top=bratio[3],
               common_hf=float(coh[hb_mask].mean()),
               spectra=acc, thumbs=thumbs)
    if ssim_v:
        out["ssim_ref"] = float(np.mean(ssim_v))
        out["bias_ref"] = float(np.mean(bias_v))
    return out


def _png(img, lo, hi):
    """8-bit greyscale PNG (no external library) as a data URI."""
    a8 = np.clip((img - lo) / max(hi - lo, 1e-30) * 255, 0, 255).astype(np.uint8)
    h, w = a8.shape
    raw = b"".join(b"\x00" + a8[i].tobytes() for i in range(h))

    def chunk(t, d):
        return struct.pack(">I", len(d)) + t + d + struct.pack(">I", zlib.crc32(t + d) & 0xFFFFFFFF)
    png = (b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", struct.pack(">IIBBBBB", w, h, 8, 0, 0, 0, 0))
           + chunk(b"IDAT", zlib.compress(raw, 9)) + chunk(b"IEND", b""))
    return "data:image/png;base64," + base64.b64encode(png).decode()


def _svg_scatter(rows, title, W=520, H=360):
    """noise (log) vs f50: hollow = before, filled = after denoising; line joins a case."""
    pts = [(r, False, r["f50_avg"], r["noise_avg"]) for r in rows] + \
          [(r, True, r["f50_den"], r["noise_left"]) for r in rows]
    pts = [p for p in pts if p[3] > 0 and np.isfinite(p[2])]
    if not pts:
        return ""
    ml, mr, mt, mb = 56, 14, 40, 40
    ys = [math.log10(p[3]) for p in pts]
    ylo, yhi = math.floor(min(ys) * 4) / 4, math.ceil(max(ys) * 4) / 4
    if yhi - ylo < 0.5:
        yhi = ylo + 0.5
    X = lambda v: ml + v / 0.52 * (W - ml - mr)
    Y = lambda v: H - mb - (math.log10(v) - ylo) / (yhi - ylo) * (H - mt - mb)
    meds = sorted({r["median"] for r in rows}, key=float)
    out = [f'<svg viewBox="0 0 {W} {H}" width="100%" style="max-width:{W}px" role="img" aria-label="{html.escape(title)}">',
           f'<text x="{ml}" y="16" class="t">{html.escape(title)}</text>']
    for e in range(int(math.floor(ylo)), int(math.ceil(yhi)) + 1):
        for mant in (1, 2, 5):
            v = mant * 10 ** e
            if ylo <= math.log10(v) <= yhi:
                out.append(f'<line x1="{ml}" x2="{W - mr}" y1="{Y(v):.1f}" y2="{Y(v):.1f}" class="g"/>'
                           f'<text x="{ml - 6}" y="{Y(v) + 4:.1f}" text-anchor="end" class="a">{v:g}</text>')
    for v in (0, 0.1, 0.2, 0.3, 0.4, 0.5):
        out.append(f'<line x1="{X(v):.1f}" x2="{X(v):.1f}" y1="{mt}" y2="{H - mb}" class="g"/>'
                   f'<text x="{X(v):.1f}" y="{H - mb + 16}" text-anchor="middle" class="a">{v:.1f}</text>')
    out.append(f'<text x="{W - mr}" y="{H - 6}" text-anchor="end" class="a">f50 (cycles/pixel, larger = sharper)</text>')
    out.append(f'<text x="{ml}" y="{mt - 8}" class="a">noise std (log scale)</text>')
    for r in rows:
        col = MED_COLORS[meds.index(r["median"]) % len(MED_COLORS)]
        if r["noise_avg"] > 0 and r["noise_left"] > 0:
            out.append(f'<line x1="{X(r["f50_avg"]):.1f}" y1="{Y(r["noise_avg"]):.1f}" x2="{X(r["f50_den"]):.1f}" '
                       f'y2="{Y(r["noise_left"]):.1f}" stroke="{col}" stroke-opacity="0.35"/>')
    for r, den, xv, yv in pts:
        col = MED_COLORS[meds.index(r["median"]) % len(MED_COLORS)]
        cx, cy = X(xv), Y(yv)
        fill = col if den else "white"
        tip = html.escape(f'{r["case"]} {"after" if den else "before"}: f50 {xv:.3f}, noise {yv:.3g}')
        if r["recon"] == sorted({q["recon"] for q in rows})[0]:
            out.append(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="4.5" fill="{fill}" stroke="{col}" stroke-width="1.5"><title>{tip}</title></circle>')
        else:
            out.append(f'<path d="M{cx:.1f},{cy - 5.5:.1f} L{cx + 5:.1f},{cy + 3.5:.1f} L{cx - 5:.1f},{cy + 3.5:.1f} Z" '
                       f'fill="{fill}" stroke="{col}" stroke-width="1.5"><title>{tip}</title></path>')
    out.append("</svg>")
    return "".join(out)


def _svg_curves(fc, curves, title, W=520, H=300):
    ml, mr, mt, mb = 44, 14, 26, 40
    ymax = min(3.0, max(1.1, max(float(np.nanmax(c[1])) for c in curves) * 1.05))
    X = lambda v: ml + v / 0.5 * (W - ml - mr)
    Y = lambda v: H - mb - v / ymax * (H - mt - mb)
    out = [f'<svg viewBox="0 0 {W} {H}" width="100%" style="max-width:{W}px" role="img" aria-label="{html.escape(title)}">',
           f'<text x="{ml}" y="16" class="t">{html.escape(title)}</text>']
    for v in np.arange(0, ymax + 1e-9, 0.25):
        out.append(f'<line x1="{ml}" x2="{W - mr}" y1="{Y(v):.1f}" y2="{Y(v):.1f}" class="g"/>'
                   f'<text x="{ml - 6}" y="{Y(v) + 4:.1f}" text-anchor="end" class="a">{v:.2f}</text>')
    for v in (0, 0.1, 0.2, 0.3, 0.4, 0.5):
        out.append(f'<text x="{X(v):.1f}" y="{H - mb + 16}" text-anchor="middle" class="a">{v:.1f}</text>')
    out.append(f'<text x="{W - mr}" y="{H - 6}" text-anchor="end" class="a">spatial frequency (cycles/pixel)</text>')
    out.append(f'<line x1="{ml}" x2="{W - mr}" y1="{Y(0.5):.1f}" y2="{Y(0.5):.1f}" stroke="#888" stroke-dasharray="4 3"/>')
    ly = mt + 6
    for label, H_, col, dash in curves:
        d = " ".join(f'{"M" if i == 0 else "L"}{X(f):.1f},{Y(min(h, ymax)):.1f}' for i, (f, h) in enumerate(zip(fc, H_)))
        out.append(f'<path d="{d}" fill="none" stroke="{col}" stroke-width="1.8"{" stroke-dasharray=" + chr(34) + "5 4" + chr(34) if dash else ""}/>')
        out.append(f'<line x1="{W - mr - 190}" x2="{W - mr - 170}" y1="{ly}" y2="{ly}" stroke="{col}" stroke-width="2"'
                   f'{" stroke-dasharray=" + chr(34) + "5 4" + chr(34) if dash else ""}/>'
                   f'<text x="{W - mr - 164}" y="{ly + 4}" class="a">{html.escape(label)}</text>')
        ly += 15
    out.append("</svg>")
    return "".join(out)


def stored_params(a):
    """Rc, Dr and RA0: the options when given, else the values recorded in <out>/params.json
    (RA0 falls back to 0.0). Returns (rc, dr, ra0); rc / dr are None when unknown."""
    try:
        with open(os.path.join(a.out, "params.json"), encoding="utf-8") as f:
            st = json.load(f)
    except (OSError, ValueError):
        st = {}
    rc = a.rc if a.rc is not None else st.get("rc")
    dr = a.dr if a.dr is not None else st.get("dr")
    ra0 = a.ra0 if a.ra0 is not None else st.get("ra0", "0.0")
    return rc, dr, ra0


def make_reference(a, fname, zc, ref_path):
    """Normal reconstruction of the central slice with the same parameters as the split
    reconstructions: <ref-prog> <layer> <Rc> <Dr> <RA0> in the raw directory. The parameters are
    recorded in <out>/summary/reference.json; an existing image that was not made by this script
    with the same parameters is kept as *.before_n2i_eval.tif and made again.
    Returns a short note for the report."""
    nums = re.findall(r"\d+", os.path.splitext(fname)[0])
    slice_no = int(nums[-1]) if nums else zc
    rc, dr, ra0 = stored_params(a)
    exists = os.path.isfile(ref_path)
    if rc is None or dr is None:
        log("analysis: the rotation centre / pixel size are not known (give --rc and --dr): "
            + ("the existing reference image is used as it is" if exists else "reference comparison skipped"))
        return "existing file, parameters not verified" if exists else ""
    sdir = os.path.join(a.out, "summary")
    os.makedirs(sdir, exist_ok=True)
    recfile = os.path.join(sdir, "reference.json")
    want = {"file": os.path.abspath(ref_path), "prog": a.ref_prog, "layer": slice_no,
            "rc": float(rc), "dr": float(dr), "ra0": float(ra0)}
    if exists:
        try:
            with open(recfile, encoding="utf-8") as f:
                old = json.load(f)
            if all(old.get(k) == v for k, v in want.items()) and old.get("mtime") == int(os.path.getmtime(ref_path)):
                return f"{a.ref_prog} {slice_no} {rc} {dr} {ra0}"
        except (OSError, ValueError):
            pass
    if shutil.which(a.ref_prog) is None:
        log(f"analysis: {a.ref_prog} not found in PATH: "
            + ("the existing reference image is used as it is" if exists else "reference comparison skipped"))
        return "existing file, parameters not verified" if exists else ""
    if exists:
        bak = os.path.splitext(ref_path)[0] + ".before_n2i_eval.tif"
        if not os.path.exists(bak):
            os.replace(ref_path, bak)
            log(f"analysis: existing {ref_path} (parameters unknown) kept as {os.path.basename(bak)}")
        else:
            os.remove(ref_path)
    log(f"analysis: making the reference image: (cd {a.raw}) {a.ref_prog} {slice_no} {rc} {dr} {ra0}")
    ok = run([a.ref_prog, slice_no, rc, dr, ra0], a.raw, os.path.abspath(os.path.join(sdir, "reference_rec.log")))
    if not ok or not os.path.isfile(ref_path):
        log(f"analysis: {a.ref_prog} did not produce {ref_path} (see summary/reference_rec.log): "
            f"reference comparison skipped")
        return ""
    want["mtime"] = int(os.path.getmtime(ref_path))
    with open(recfile, "w", encoding="utf-8") as f:
        json.dump(want, f)
    return f"{a.ref_prog} {slice_no} {rc} {dr} {ra0}"


def analyze(a, cases):
    """Detailed analysis of all finished cases: CSV files and an HTML report."""
    fin = [c for c in cases if done(c["dir"], "infer")]
    if not fin:
        log("analysis: no finished case")
        return
    first = tif_list(os.path.join(fin[0]["dir"], "rec", "0"))
    h, w = tifffile.imread(first[0]).shape
    size = min(a.roi_size, h, w)
    cy, cx = (h // 2, w // 2) if not a.roi_center else (int(float(v)) for v in a.roi_center.split(","))
    y0 = min(max(cy - size // 2, 0), h - size)
    x0 = min(max(cx - size // 2, 0), w - size)
    roi = (y0, y0 + size, x0, x0 + size)
    S = Spectra(size)

    # reference image: normal reconstruction of the central slice (made with --ref-prog if missing)
    zc = len(first) // 2
    ref_path = a.ref or os.path.join(a.raw, os.path.basename(first[zc]))
    ref_note = "given with --ref" if a.ref else "existing file, parameters not verified"
    if not a.ref and not a.no_ref_rec:
        ref_note = make_reference(a, os.path.basename(first[zc]), zc, ref_path)
    refimg = None
    if os.path.isfile(ref_path):
        refimg = tifffile.imread(ref_path).astype(np.float64)
        if refimg.shape != (h, w):
            log(f"analysis: reference {ref_path} has shape {refimg.shape}, not {(h, w)}: not used")
            refimg = None
    log(f"analysis: {len(fin)} cases, ROI y {y0}-{y0 + size - 1} x {x0}-{x0 + size - 1}, "
        f"reference image {'%s' % ref_path if refimg is not None else 'none'}")

    res = []
    for c in fin:
        try:
            r = analyze_case(a, c, roi, S, refimg)
            if r:
                res.append(r)
        except Exception as e:
            log(f"analysis: {c['name']} failed ({type(e).__name__}: {e})")
    if not res:
        log("analysis: nothing to analyse")
        return

    # transfer reference: the unfiltered case with the most high-frequency signal (smallest K)
    medmin = min(res, key=lambda r: float(r["median"]))["median"]
    sigmin = min(res, key=lambda r: float(r["sigma"]))["sigma"]
    cand0 = [r for r in res if r["median"] == medmin and r["sigma"] == sigmin]
    cand = cand0 or res
    kmin = min(r["K"] for r in cand)
    cand = [r for r in cand if r["K"] == kmin]
    band = (S.fc >= 0.15) & (S.fc < 0.35)
    tref = max(cand, key=lambda r: r["spectra"]["PS"][band].sum())
    # where the reference has (almost) no signal the ratio is meaningless: floor it, cap the ratio
    PSref = np.maximum(tref["spectra"]["PS"], 1e-4 * max(float(np.max(tref["spectra"]["PS"])), 1e-30))
    curves = []
    for r in res:
        Ha = np.sqrt(np.clip(r["spectra"]["PS"] / PSref, 0, 25))
        Hd = np.sqrt(np.clip(r["spectra"]["Pden"] / PSref, 0, 25))
        r.update(f50_avg=_f_at(S.fc, Ha), f50_den=_f_at(S.fc, Hd),
                 H02_avg=float(np.interp(0.2, S.fc, Ha)), H02_den=float(np.interp(0.2, S.fc, Hd)),
                 H03_avg=float(np.interp(0.3, S.fc, Ha)), H03_den=float(np.interp(0.3, S.fc, Hd)))
        r["Ha"], r["Hd"] = Ha, Hd
        curves += [(r["case"], "average", f, v) for f, v in zip(S.fc, Ha)]
        curves += [(r["case"], "denoised", f, v) for f, v in zip(S.fc, Hd)]

    # Pareto front: no other case is both sharper (f50) and less noisy after denoising
    for r in res:
        r["pareto"] = not any((q["f50_den"] >= r["f50_den"] and q["noise_left"] <= r["noise_left"]) and
                              (q["f50_den"] > r["f50_den"] or q["noise_left"] < r["noise_left"]) for q in res)

    key = lambda r: (float(r["median"]), float(r["sigma"]), r["recon"], r["K"])
    res.sort(key=key)
    sdir = os.path.join(a.out, "summary")
    os.makedirs(sdir, exist_ok=True)
    cols = ["case", "median", "sigma", "recon", "K", "slices", "noise_avg", "removed", "noise_left",
            "f50_avg", "f50_den", "H02_avg", "H02_den", "H03_avg", "H03_den",
            "rr_low", "rr_mid", "rr_high", "rr_top", "common_hf", "ssim_ref", "bias_ref", "pareto"]
    with open(os.path.join(sdir, "analysis_metrics.csv"), "w", newline="", encoding="utf-8") as f:
        wr = csv.writer(f)
        wr.writerow(cols)
        for r in res:
            wr.writerow([(f"{r[c]:.5g}" if isinstance(r.get(c), float) else r.get(c, "")) for c in cols])
    with open(os.path.join(sdir, "transfer_curves.csv"), "w", newline="", encoding="utf-8") as f:
        wr = csv.writer(f)
        wr.writerow(["case", "type", "freq_cycles_per_pixel", "transfer"])
        wr.writerows((cn, t, f"{fv:.4f}", f"{v:.4f}") for cn, t, fv, v in curves)
    report(a, res, S, roi, tref, refimg is not None, ref_path, os.path.join(sdir, "analysis_report.html"),
           ref_note, tref in cand0)
    log(f"analysis: {os.path.join(sdir, 'analysis_report.html')}, analysis_metrics.csv, transfer_curves.csv")


def report(a, res, S, roi, tref, has_ref, ref_path, path, ref_note="", tref_unfiltered=True):
    esc = html.escape
    recons = sorted({r["recon"] for r in res})
    ks = sorted({r["K"] for r in res})
    meds = sorted({r["median"] for r in res}, key=float)
    medmin = meds[0]
    sigmin = min({r["sigma"] for r in res}, key=float)
    unf = [r for r in res if r["median"] == medmin and r["sigma"] == sigmin]
    pareto = sorted([r for r in res if r["pareto"]], key=lambda r: -r["f50_den"])

    p = ['<!doctype html><html lang="ja"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">',
         '<title>N2I evaluation report</title><style>',
         'body{font-family:system-ui,-apple-system,"Segoe UI",sans-serif;margin:24px;max-width:1150px;color:#222}',
         'h1{font-size:22px} h2{font-size:17px;margin-top:28px;border-bottom:1px solid #ddd;padding-bottom:4px}',
         'table{border-collapse:collapse;font-size:12px;font-variant-numeric:tabular-nums}',
         'td,th{border-bottom:1px solid #e3e3e3;padding:3px 7px;text-align:right;white-space:nowrap}',
         'td:first-child,th:first-child{text-align:left} th{color:#555;font-weight:600}',
         '.wrap{overflow-x:auto} .grid{display:flex;flex-wrap:wrap;gap:18px} .note{color:#555;font-size:13px}',
         'svg .g{stroke:#e5e5e5} svg .a{fill:#555;font-size:11px} svg .t{fill:#222;font-size:13px;font-weight:600}',
         '.thumbs{display:flex;flex-wrap:wrap;gap:8px;margin-bottom:10px}',
         '.thumbs div{width:220px;font-size:11px;color:#555} .thumbs img{width:220px;height:auto;display:block}',
         'tr.pf td{background:#fff7e0}</style></head><body>',
         '<h1>N2I evaluation report</h1>',
         f'<p class="note">{esc(os.path.abspath(a.out))} — {time.strftime("%Y-%m-%d %H:%M")} — '
         f'ROI y {roi[0]}–{roi[1] - 1}, x {roi[2]}–{roi[3] - 1} ({roi[1] - roi[0]}×{roi[1] - roi[0]} px) of the central slice'
         f'{" ±%d" % a.infer_half if a.infer_half else ""}. Transfer reference: <b>{esc(tref["case"])}</b> '
         f'({"unfiltered" if tref_unfiltered else "NOT unfiltered: no unfiltered case was available"}, signal).'
         f' Reference image: {esc(ref_path) + " (" + esc(ref_note) + ")" if has_ref else "none"}.</p>']

    p.append("<h2>要点</h2><ul>")
    for r in sorted(unf, key=lambda r: (r["recon"], r["K"])):
        p.append(f'<li>{esc(r["case"])}（フィルタなし）: ノイズ {r["noise_avg"]:.3g} → 約 {r["noise_left"]:.3g}、'
                 f'f50 {r["f50_avg"]:.3f} → {r["f50_den"]:.3f}</li>')
    if pareto:
        p.append("<li>ノイズ除去後に「より鮮明で、かつノイズも少ない」ケースが他にない条件（パレート最適、鮮明な順"
                 + (f"に上位 8 件、全 {len(pareto)} 件は表の黄色の行" if len(pareto) > 8 else "") + "）: "
                 + "、".join(f'{esc(r["case"])}（f50 {r["f50_den"]:.3f}、ノイズ {r["noise_left"]:.3g}）' for r in pareto[:8])
                 + "。この中から、目的に合わせて分解能とノイズの兼ね合いで選びます。</li>")
    p.append("</ul>")

    p.append("<h2>残ったノイズと分解能</h2><p class=\"note\">白抜き＝ノイズ除去前、塗りつぶし＝ノイズ除去後（同じケースを線で結んでいます）。"
             "色はメディアン、丸と三角は再構成の違い。右下ほど良い条件です。点にマウスを重ねると値が出ます。"
             f'色: {"、".join(f"<span style=color:{MED_COLORS[i % len(MED_COLORS)]}>■ median {esc(m)}</span>" for i, m in enumerate(meds))}；'
             f'丸 = {esc(recons[0])}{"、三角 = " + "・".join(esc(x) for x in recons[1:]) if len(recons) > 1 else ""}</p><div class="grid">')
    for k in ks:
        p.append(_svg_scatter([r for r in res if r["K"] == k], f"K = {k}"))
    p.append("</div>")

    p.append("<h2>信号の伝達（フィルタなしのケース）</h2><p class=\"note\">基準ケースの信号に対する、周波数ごとの伝わり方。"
             "実線＝ノイズ除去前（2 つの独立な半分のデータの相互スペクトル）、破線＝ノイズ除去後（残ったノイズの分だけわずかに大きく出ます）。"
             "f50 はこの値が 0.5 に落ちる周波数です。</p><div class=\"grid\">")
    cols = ["#1f77b4", "#e8590c", "#2b8a3e", "#9c36b5"]
    cv = []
    for i, r in enumerate(sorted(unf, key=lambda r: (r["recon"], r["K"]))):
        cv.append((f'{r["case"]} before', r["Ha"], cols[i % 4], False))
        cv.append((f'{r["case"]} after', r["Hd"], cols[i % 4], True))
    if cv:
        p.append(_svg_curves(S.fc, cv, "unfiltered cases"))
    p.append("</div>")

    # thumbnails: unfiltered cases and the strongest filter, smallest K
    sel = [r for r in res if r["K"] == min(ks) and (r in unf)]
    strong = [r for r in res if r["K"] == min(ks) and r["median"] == meds[-1] and
              r["sigma"] == max({q["sigma"] for q in res}, key=float)]
    sel += strong[:1]
    if sel:
        p.append("<h2>画像（ROI、中央スライス）</h2><p class=\"note\">左から ノイズ除去前（分割の平均）、ノイズ除去後、取り除いた成分（平均 − 除去後）。"
                 "取り除いた成分に試料の輪郭が見えなければ、構造を削っていません。</p>")
        for r in sel:
            if not r.get("thumbs"):
                continue
            avg, den, rs = r["thumbs"]
            lo, hi = np.percentile(den, [1, 99])
            q = max(abs(np.percentile(rs, 1)), abs(np.percentile(rs, 99)), 1e-30)
            p.append(f'<div class="thumbs"><div>{esc(r["case"])}<br>before<img src="{_png(avg, lo, hi)}"></div>'
                     f'<div><br>after<img src="{_png(den, lo, hi)}"></div>'
                     f'<div><br>removed<img src="{_png(rs, -q, q)}"></div></div>')

    head = ["case", "median", "sigma", "recon", "K", "noise before", "removed", "noise after (est.)", "f50 before", "f50 after",
            "H(0.3) after", "removed/noise low", "mid", "high", "top", "common HF"] + (["SSIM ref", "bias ref"] if has_ref else [])
    p.append("<h2>全ケース</h2><p class=\"note\">黄色の行はパレート最適。removed/noise は帯域ごとの「取り除いた成分 ÷ 独立なノイズ」のパワー比"
             "（1 ＝ ちょうどノイズ分、1 より大きい ＝ ノイズ以外も除去、小さい ＝ ノイズが残った）。"
             "common HF は分割間で共通の高周波成分の強さ（0.4–0.5 cycles/pixel のコヒーレンス。大きいと FRC は分解能の指標になりません）。"
             "SSIM ref・bias ref は通常の再構成（ノイズを含む）との比較で、正解との一致度ではありません。</p>"
             '<div class="wrap"><table><tr>' + "".join(f"<th>{esc(x)}</th>" for x in head) + "</tr>")
    for r in res:
        vals = [r["case"], r["median"], r["sigma"], r["recon"], r["K"], f'{r["noise_avg"]:.3g}', f'{r["removed"]:.3g}',
                f'{r["noise_left"]:.3g}', f'{r["f50_avg"]:.3f}', f'{r["f50_den"]:.3f}', f'{r["H03_den"]:.2f}',
                f'{r["rr_low"]:.2f}', f'{r["rr_mid"]:.2f}', f'{r["rr_high"]:.2f}', f'{r["rr_top"]:.2f}', f'{r["common_hf"]:.2f}']
        if has_ref:
            vals += [f'{r.get("ssim_ref", float("nan")):.3f}', f'{r.get("bias_ref", float("nan")):.3g}']
        p.append(f'<tr class="{"pf" if r["pareto"] else ""}">' + "".join(f"<td>{esc(str(v))}</td>" for v in vals) + "</tr>")
    p.append("</table></div>")
    p.append("<h2>方法と注意</h2><ul class=\"note\">"
             "<li>独立なノイズ：投影の半分ずつ（K=2 は分割 0 と 1、K=4 は 0+2 と 1+3）の差 ÷ 2。分割の単純平均に含まれるノイズと同じ分散になります。</li>"
             "<li>ノイズ除去後のノイズは「除去前のノイズの分散 − 取り除いた成分の分散」からの推定です（取り除いた成分がノイズだけの場合に正確）。</li>"
             "<li>信号の伝達：2 つの半分の相互スペクトル（独立なノイズを含まない信号のパワー）の、基準ケースに対する比の平方根。メディアンフィルタは非線形なので目安です。</li>"
             "<li>f50 が 0.5 のケースは、測定範囲内で伝達が 0.5 を下回らなかったことを示します。</li>"
             "<li>基準画像（SSIM ref、bias ref）は通常の再構成プログラム（" + esc(a.ref_prog) + "）で作った 1 種類です。"
             "再構成の方法が違う系統のケースは、それだけで基準との差が大きく出ます。同じ系統の中で比べてください。</li>"
             "<li>除去後のノイズが 0 と表示されるケースは、取り除いた成分がノイズの量を上回った（構造も削った）ことを示し、散布図には表示されません。</li>"
             "<li>ROI 1 か所での評価です。別の位置（--roi-center）でも確かめると確実です。</li></ul></body></html>")
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(p))


# ------------------------------------------------------------------ parameters
def _number(text, what):
    try:
        float(text)
    except (TypeError, ValueError):
        raise SystemExit(f"{what}: '{text}' is not a number")
    return text.strip()


def reconstruction_params(a):
    """Pixel size and rotation centre must be given (option or keyboard); they are recorded in
    <out>/params.json and must be the same when a run is restarted."""
    for key, label in (("dr", "pixel size Dr (um)"), ("rc", "rotation centre Rc (pixels)")):
        v = getattr(a, key)
        if v is None:
            if not sys.stdin.isatty():
                raise SystemExit(f"{label} is required: give --{key} (no keyboard input in this run)")
            while True:
                try:
                    v = input(f"{label}: ").strip()
                except EOFError:
                    raise SystemExit(f"{label} is required: give --{key}")
                try:
                    float(v)
                    break
                except ValueError:
                    print("  please enter a number")
        setattr(a, key, _number(v, label))
    if a.ra0 is None:
        a.ra0 = "0.0"
    a.ra0 = _number(a.ra0, "in-plane rotation angle")
    cur = {"dr": a.dr, "rc": a.rc, "ra0": a.ra0}
    pfile = os.path.join(a.out, "params.json")
    if os.path.exists(pfile):
        with open(pfile, encoding="utf-8") as f:
            old = json.load(f)
        diff = {k: (old.get(k), v) for k, v in cur.items() if old.get(k) is not None and float(old[k]) != float(v)}
        if diff:
            raise SystemExit("these values differ from the earlier run in this directory "
                             + ", ".join(f"{k}: before {o}, now {n}" for k, (o, n) in diff.items())
                             + f" -- use the same values, or another --out (or delete {a.out})")
    elif not a.dry_run:
        with open(pfile, "w", encoding="utf-8") as f:
            json.dump(cur, f)
    log(f"reconstruction parameters: Dr {a.dr} um, Rc {a.rc}, RA0 {a.ra0} deg")


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser(
        description="Evaluation of projection filtering and Noise2Inverse in one run: filter a.his (median x "
                    "gaussian), reconstruct K splits with each program, train and denoise with n2i_v4.py, then "
                    "analyse noise and resolution and write a report (HTML, CSV) and an ImageJ stack.",
        formatter_class=argparse.HelpFormatter)
    ap.add_argument("--raw", default="raw", help="directory with a.his, conv.bat, output.log (default: raw)")
    ap.add_argument("--his", default="a.his", help="his file name in --raw (default: a.his)")
    ap.add_argument("--copy-files", nargs="+", default=["conv.bat", "output.log"],
                    help="companion files copied next to every filtered his (default: conv.bat output.log)")
    ap.add_argument("--out", default="n2i_eval", help="output directory (default: n2i_eval)")
    ap.add_argument("--dr", default=None, help="pixel size (um); asked for when not given")
    ap.add_argument("--rc", default=None, help="rotation centre (pixels); asked for when not given")
    ap.add_argument("--ra0", default=None, help="in-plane rotation angle in deg (default: 0.0)")
    ap.add_argument("--medians", nargs="+", default=["1", "3", "5"],
                    help="median kernel sizes for his_mgf; 0 or 1 = no median filter (default: 1 3 5)")
    ap.add_argument("--sigmas", nargs="+", default=["0", "0.5", "1.0", "2.0"],
                    help="gaussian sigmas for his_mgf; 0 = no gaussian filter (default: 0 0.5 1.0 2.0)")
    ap.add_argument("--progs", nargs="+", default=["his_srec_g_c_n2i", "his_srec_g_r_n2i"],
                    help="split-reconstruction programs, run as '<prog> <raw> <Dr> <Rc> <RA0> rec <K>' "
                         "(default: his_srec_g_c_n2i his_srec_g_r_n2i)")
    ap.add_argument("--splits", nargs="+", type=int, default=[2, 4], help="split counts K (default: 2 4)")
    ap.add_argument("--mgf", default="his_mgf", help="his_mgf command (default: his_mgf)")
    ap.add_argument("--n2i", default=None,
                    help="n2i_v4.py, its launcher or installation directory (default: searched in the "
                         "installation of this Python, then " + ", ".join(N2I_PLACES) + ")")
    ap.add_argument("--gpus", default="",
                    help="GPUs used in parallel, e.g. 0,1 (one case per GPU at a time). Numbers count among "
                         "the GPUs already visible (CUDA_VISIBLE_DEVICES, e.g. set by Slurm). "
                         "Default: one case at a time on the GPU(s) as given by the environment")
    ap.add_argument("--train-args", default="", help="extra options for n2i train, e.g. \"--iters 1000\"")
    ap.add_argument("--infer-args", default="", help="extra options for n2i infer")
    ap.add_argument("--infer-half", type=int, default=0,
                    help="denoise the central slice +- this many slices (default: only the central one)")
    ap.add_argument("--keep-cache", action="store_true", help="keep image_cache.npy after training")
    ap.add_argument("--delete-rec", action="store_true",
                    help="when a case is finished, delete its split reconstructions except the central "
                         "slice(s) needed for the summary and the analysis")
    ap.add_argument("--summary-only", action="store_true", help="only rebuild the summary stack and the analysis")
    ap.add_argument("--analyze-only", action="store_true", help="only run the analysis of finished cases")
    ap.add_argument("--no-analyze", action="store_true", help="skip the analysis")
    ap.add_argument("--no-stack", action="store_true", help="skip the summary stack (large file)")
    ap.add_argument("--ref", default=None,
                    help="reference image (normal reconstruction of the central slice; default: raw/rec<slice>.tif, "
                         "made with --ref-prog when missing)")
    ap.add_argument("--ref-prog", default="his_rec_g_c",
                    help="program that makes the reference: run as '<prog> <layer> <Rc> <Dr> <RA0>' in --raw "
                         "(default: his_rec_g_c)")
    ap.add_argument("--no-ref-rec", action="store_true", help="do not make a missing reference image")
    ap.add_argument("--roi-size", type=int, default=512, help="analysis region size in pixels (default: 512)")
    ap.add_argument("--roi-center", default=None, help="analysis region centre 'y,x' (default: image centre)")
    ap.add_argument("--dry-run", action="store_true", help="print the commands without running them")
    a = ap.parse_args()

    global LOGF
    os.makedirs(a.out, exist_ok=True)
    LOGF = open(os.path.join(a.out, "eval.log"), "a", encoding="utf-8")
    log(f"=== n2i_eval_all.py  {' '.join(sys.argv[1:])}")

    a.n2i, a.n2i_py = resolve_n2i(a.n2i)
    if a.n2i:
        log(f"n2i: {a.n2i}  (python: {a.n2i_py})")
    # --gpus counts among the GPUs that are already visible (a job scheduler such as Slurm sets
    # CUDA_VISIBLE_DEVICES to the GPUs of the job); without --gpus the environment is left alone
    sel = [g.strip() for g in a.gpus.split(",") if g.strip()]
    vis = [v.strip() for v in os.environ.get("CUDA_VISIBLE_DEVICES", "").split(",") if v.strip()]
    if sel and vis:
        try:
            sel = [vis[int(g)] for g in sel]
        except (ValueError, IndexError):
            raise SystemExit(f"--gpus {a.gpus} does not match the visible GPUs {vis} (CUDA_VISIBLE_DEVICES)")
    gpus = sel or [None]

    cases, his_jobs = [], []
    for m in a.medians:
        for s in a.sigmas:
            hisdir = os.path.abspath(os.path.join(a.out, "his", f"m{m}_g{s}"))
            his_jobs.append((m, s, hisdir))
            for prog in a.progs:
                for k in a.splits:
                    pl = prog_label(prog)
                    name = f"m{m}_g{s}_{pl}_K{k}"
                    cases.append(dict(name=name, m=m, s=s, prog=prog, plabel=pl, k=k, hisdir=hisdir,
                                      dir=os.path.abspath(os.path.join(a.out, "cases", name))))
    log(f"{len(his_jobs)} filtered his files, {len(cases)} cases, "
        f"GPUs {gpus if gpus != [None] else 'as given by the environment'}"
        f"{' (CUDA_VISIBLE_DEVICES=' + ','.join(vis) + ')' if vis else ''}")

    unfinished = False
    if a.analyze_only:
        a.summary_only = True
        a.no_stack = True
    if not a.summary_only:
        missing = [c for c in [a.mgf] + a.progs if shutil.which(c) is None]
        if missing and not a.dry_run:
            raise SystemExit(f"not found in PATH: {' '.join(missing)}")
        if not a.n2i:
            raise SystemExit("n2i_v4.py not found; give its location with --n2i "
                             "(e.g. --n2i ~/data/ct-rec/ct-rec/N2I/n2i_portable/n2i)")
        if not os.path.isfile(os.path.join(a.raw, a.his)) and not a.dry_run:
            raise SystemExit(f"{os.path.join(a.raw, a.his)} not found (use --raw / --his)")
        odd = [k for k in a.splits if k % 2 or k < 2]
        if odd:
            log(f"note: K = {odd} is not an even number >= 2: such cases are not analysed "
                f"(his_srec_*_n2i itself accepts K = 2 or 4)")
        reconstruction_params(a)          # asked for only after the checks above have passed
        if not a.dry_run:
            hs = os.path.getsize(os.path.join(a.raw, a.his))
            free = shutil.disk_usage(os.path.abspath(a.out)).free
            log(f"disk: {len(his_jobs)} filtered his files need about {len(his_jobs) * hs / 1e9:.0f} GB; "
                f"free {free / 1e9:.0f} GB (reconstructions come on top of this)")

        # 1. filtering (sequential: disk bound)
        ok_his = {}
        for m, s, hisdir in his_jobs:
            ok_his[hisdir] = step_filter(a, m, s, hisdir)

        # 2-3. reconstruction, training, inference: one case per GPU at a time
        todo = [c for c in cases if ok_his.get(c["hisdir"])]
        free_gpus = queue.Queue()
        for g in gpus:
            free_gpus.put(g)
        failed = []

        def work(c):
            g = free_gpus.get()
            try:
                if not case_pipeline(a, c, g):
                    failed.append(c["name"])
            except Exception as e:
                log(f"ERROR {c['name']}: {type(e).__name__}: {e}")
                failed.append(c["name"])
            finally:
                free_gpus.put(g)

        with ThreadPoolExecutor(max_workers=len(gpus)) as pool:
            list(pool.map(work, todo))
        skipped = [c["name"] for c in cases if not ok_his.get(c["hisdir"])]
        if failed or skipped:
            log(f"not finished: {len(failed)} failed {failed}, {len(skipped)} skipped (filter failed) {skipped}")
            unfinished = True

    # 4. summary stack and analysis
    if not a.dry_run:
        if not a.no_stack:
            summary(a, cases)
        if not a.no_analyze:
            analyze(a, cases)
    log("=== done" + (" (some cases did not finish: see above; run again to retry them)" if unfinished else ""))
    if unfinished:
        sys.exit(1)


if __name__ == "__main__":
    main()
