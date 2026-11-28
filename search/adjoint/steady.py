"""steady.py -- the adjoint's base flow: the production state continued until it is steady (owner 2026-11-28:
"Continue to steady"), on the same mesh.

The production score (1000 iterations from cold, linearUpwind, mean of 200-1000) never quite settles, and
DAFoam's Krylov adjoint stalls on an unsteady primal (adjoint repo ROUTES.md 2026-11-20). So the adjoint is
taken on a steady continuation of the production state:

  rung 1: the production scheme (linearUpwind) continued `extra` more iterations;
  rung 2 (only if rung 1 is not steady): first-order upwind for div(phi,U) (the "steady model" that passed
          on the coarse mesh), continued from rung 1's state.

STEADY (fixed 2026-11-21 in the adjoint repo, fine_campaign/steady_base.py, kept): over the last 1,000
iterations the 100-iteration block means of the half-car drag span <= 0.10 mN (coarse mesh: 0.02 mN). A rung
that passes is the base; if neither does, the last rung's state is still handed on, flagged `steady: false`
(the adjoint may then stall; the job reports it and the gradient is not used).

    python steady.py <production state dir> <work dir> [np] [extra] [span_mN]
"""
from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

BASHRC = os.environ.get("OF_BASHRC", "/usr/lib/openfoam/openfoam2412/etc/bashrc")
SCHEMES = (("linearUpwind", None), ("upwind", "bounded Gauss upwind"))


def _times(d: Path):
    return sorted((int(p.name) for p in d.iterdir() if re.fullmatch(r"\d+", p.name)), key=int)


def drag_mN(run: Path):
    """Half-car drag per iteration (mN) over every forces file of a run, in time order (restarts append)."""
    rows = []
    for f in sorted((run / "postProcessing" / "forces").glob("*/force*.dat"), key=lambda p: float(p.parent.name)):
        a = np.loadtxt([l.replace("(", " ").replace(")", " ") for l in open(f) if l.strip() and not l.startswith("#")])
        rows.append(a[:, [0, 1]])
    a = np.vstack(rows)
    _, i = np.unique(a[:, 0], return_index=True)       # a restart re-writes its first time: keep one
    return a[i, 1] * 1e3


def flatness(D, n=1000, block=100):
    last = D[-n:]
    b = last[: len(last) // block * block].reshape(-1, block).mean(1)
    return float(b.max() - b.min()), b


def continue_run(src: Path, dst: Path, n_proc: int, extra: int, div_u: str | None):
    """Copy src (a reconstructed OpenFOAM case at its latest time) to dst and run `extra` more simpleFoam
    iterations, with div(phi,U) replaced by div_u if given."""
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst, ignore=shutil.ignore_patterns("processor*", "logs"))
    t0 = _times(dst)[-1]
    if div_u:
        fs = dst / "system" / "fvSchemes"
        s = fs.read_text()
        s2 = re.sub(r"div\(phi,U\)\s+[^;]+;", f"div(phi,U)      {div_u};", s, count=1)
        if s2 == s and div_u not in s:
            raise SystemExit("div(phi,U) not found in fvSchemes")
        fs.write_text(s2)
    cd = dst / "system" / "controlDict"
    c = cd.read_text()
    for k, v in (("startFrom", "latestTime"), ("endTime", str(t0 + extra)), ("writeInterval", str(extra)),
                 ("purgeWrite", "0")):
        c = re.sub(rf"^{k}\s+[^;]+;", f"{k}       {v};", c, flags=re.M)
    cd.write_text(c)
    dp = dst / "system" / "decomposeParDict"
    if dp.exists():
        dp.write_text(re.sub(r"numberOfSubdomains\s+\d+;", f"numberOfSubdomains {n_proc};", dp.read_text()))
    sh = f"""source {BASHRC} || true; set -e
cd {dst}; decomposePar -force -latestTime > log.decompose.{t0} 2>&1
mpirun -np {n_proc} simpleFoam -parallel > log.simpleFoam.{t0} 2>&1
reconstructPar -latestTime > log.reconstruct.{t0} 2>&1; rm -rf processor*"""
    subprocess.run(["bash", "-c", sh], check=True)
    return t0


def run(prod: Path, work: Path, n_proc=4, extra=2000, span_mN=0.10) -> dict:
    work.mkdir(parents=True, exist_ok=True)
    src, rungs = Path(prod), []
    for name, div_u in SCHEMES:
        dst = work / f"steady_{name}"
        t0 = continue_run(src, dst, n_proc, extra, div_u)
        D = drag_mN(dst)
        span, blocks = flatness(D)
        rungs.append({"scheme": name, "from_time": t0, "iterations": extra, "span_mN": span,
                      "block_means_mN": blocks.round(3).tolist(), "D_half_mN_last1000": float(D[-1000:].mean()),
                      "state": str(dst.resolve())})
        print(json.dumps(rungs[-1]), flush=True)
        if span <= span_mN:
            break
        src = dst
    last = rungs[-1]
    out = {"steady": last["span_mN"] <= span_mN, "span_mark_mN": span_mN, "scheme": last["scheme"],
           "state": last["state"], "D_half_mN": last["D_half_mN_last1000"], "rungs": rungs}
    (work / "steady.json").write_text(json.dumps(out, indent=1))
    print("STEADY_PASS" if out["steady"] else "STEADY_FAIL", last["scheme"], round(last["span_mN"], 4), flush=True)
    return out


if __name__ == "__main__":
    a = sys.argv
    run(Path(a[1]), Path(a[2]), int(a[3]) if len(a) > 3 else 4, int(a[4]) if len(a) > 4 else 2000,
        float(a[5]) if len(a) > 5 else 0.10)
