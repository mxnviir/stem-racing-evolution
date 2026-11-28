"""fine_adjoint.py -- one car's drag gradient with respect to the free-form body's numbers: the whole chain, the
same on an AWS Batch node (aws/, the 2.1 M-cell wall-resolved mesh, owner 2026-11-28) and on a GitHub runner
(the coarse end-to-end test).

    python fine_adjoint.py <job.json> <work dir> --adjoint-root <stem-racing-adjoint checkout> [--np N]

job.json (written by requests.py from a stage-2 design; every key but the car is optional):
  {"id": "...",
   "car": {"genes": [...layout...], "freeform": {...}, "support": {...}, "fwing": {...}, "rwing": {...}}
          or {"cfd_car_dir": "<folder with car_half.stl, wheelF/R.stl, assembly.json>", "layout_car_json": "..."},
   "parametrisation": {"freeform_seed": {...}, "v": [...CMA vector...], "coarse": [5,2,3,4], "amp_mm": [6,4,6]},
   "mesh": {"res": "resolved", "wall_resolved": true, "frame_k": 0},
   "steady": {"extra": 2000, "span_mN": 0.10},
   "dafoam": {"primal_iters": 300, "env": {"DAF_ADJ": "Krylov", "DAF_PC_FILL": "1", ...}},
   "link": {"tol_mm": 1.0},
   "stop_after": "steady"   (optional: production + steady only -- the forward checks of the end-to-end test)}

Steps (each one's output is kept in <work>, and a step whose output exists is skipped on a re-run):
  1. car: built from its genes (layout_car.build, the referee) and written as the CFD car folder
     (openfoam_eval.write_cfd_car), unless a CFD car folder is given;
  2. production: the route's own run (openfoam_eval.run_cfd: fresh snappy mesh, k-omega SST, linearUpwind,
     1000 iterations from cold), its run folder kept -- so the adjoint's mesh is the scored mesh;
  3. steady: steady.py (linearUpwind continued; upwind if that does not settle);
  4. DAFoam: the adjoint repo's dafoam/make_case.py on the steady state, then its runScript.py (primal
     `primal_iters` from the steady flow, then the discrete adjoint): dD/dX for every point of the car patch;
  5. link: link.py -> gradient.json (fine: per control point; coarse: per CMA number; half and whole car).
Writes <work>/result.json: the drag at each stage, the steadiness, the adjoint's health and the gradients.

DAFoam runs natively when this script is inside the DAFoam image (aws/Dockerfile: DAFOAM_RUNNER=native) and in
DAFoam's Docker image otherwise (DAFOAM_RUNNER=docker, as on GitHub runners).
"""
from __future__ import annotations

import sys as _sys, pathlib as _pl  # the repo root, then every code folder (repo_paths.py)
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2]))
import repo_paths  # noqa: E402,F401
import argparse
import json
import os
import re
import shutil
import subprocess
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
DAFOAM_DEFAULT_ENV = {                       # the coarse settings that converged on runners (adjoint repo, coarse_dafoam.yml)
    "DAF_ADJ": "Krylov", "DAF_PC_FILL": "0", "DAF_PC_LEVEL": "1", "DAF_GMRES_ITERS": "2000",
    "DAF_GMRES_RESTART": "200", "DAF_GMRES_TOL": "1e-6"}


class _Stop(Exception):
    """The job asked to stop after a step (job["stop_after"])."""


def log(msg):
    print(f"[fine_adjoint {time.strftime('%H:%M:%S')}] {msg}", flush=True)


def step_car(job, work: Path) -> tuple[Path, dict]:
    import layout_car as LC
    import openfoam_eval as OE
    car = job["car"]
    if car.get("cfd_car_dir"):                    # relative paths: to the job file's folder (batch.py package)
        rel = lambda q: q if Path(q).is_absolute() else str(Path(job["_dir"]) / q)
        lc = json.loads(Path(rel(car["layout_car_json"])).read_text()) if car.get("layout_car_json") else {}
        return Path(rel(car["cfd_car_dir"])), lc
    cfd = work / "cfd_car"
    if (cfd / "car_half.stl").exists():
        return cfd, json.loads((work / "car" / "layout_car.json").read_text())
    rec = LC.build(np.asarray(car["genes"], float), (work / "car").resolve(), freeform=car.get("freeform"),
                   support=car.get("support"), fwing=car.get("fwing"), rwing=car.get("rwing"))
    if not rec.get("legal"):
        raise SystemExit(f"car not legal: {rec.get('failed')} {rec.get('why', '')[:300]}")
    OE.write_cfd_car(work / "car", cfd)
    return cfd, json.loads((work / "car" / "layout_car.json").read_text())


def step_production(job, cfd: Path, work: Path, adjoint_root: Path, n_proc: int) -> dict:
    import openfoam_eval as OE
    pj = work / "production.json"
    if pj.exists():
        return json.loads(pj.read_text())
    m = job.get("mesh", {})
    out = work / "production"
    out.mkdir(parents=True, exist_ok=True)
    rec = OE.run_cfd(cfd, out, adjoint_root, int(m.get("frame_k", 0)), n_proc, res=m.get("res", "resolved"),
                     wall_resolved=bool(m.get("wall_resolved", True)))
    rec["state"] = str((out / "state_prod").resolve())
    if rec.get("error"):
        raise SystemExit(f"production run failed: {rec['error']}")
    pj.write_text(json.dumps(rec, indent=1, default=float))
    return rec


def step_steady(job, prod: dict, work: Path, n_proc: int) -> dict:
    import steady
    sj = work / "steady" / "steady.json"
    if sj.exists():
        return json.loads(sj.read_text())
    s = job.get("steady", {})
    return steady.run(Path(prod["state"]), work / "steady", n_proc, int(s.get("extra", 2000)), float(s.get("span_mN", 0.10)))


def dafoam_health(logtext: str) -> dict:
    """What DAFoam's log says about its adjoint (the Krylov solve can stall: 'Adjoint solution failed')."""
    h = {"failed": bool(re.search(r"[Aa]djoint solution failed|NaN|Floating point exception|Killed", logtext))}
    res = re.findall(r"(?:Main iteration|KSP Residual norm)\s*\d*\s*([0-9.eE+-]+)", logtext)
    if res:
        try:
            h["first_residual"], h["last_residual"] = float(res[0]), float(res[-1])
        except ValueError:
            pass
    m = re.search(r"DAFOAM primal D_half = ([0-9.eE+-]+)", logtext)
    h["primal_D_half_N"] = float(m.group(1)) if m else None
    return h


def step_dafoam(job, st: dict, work: Path, adjoint_root: Path, n_proc: int) -> dict:
    dj = work / "dafoam.json"
    if dj.exists():
        return json.loads(dj.read_text())
    dd = adjoint_root / "research" / "sidepod_opt" / "dafoam"
    case = work / "dafoam_case"
    d = job.get("dafoam", {})
    subprocess.run([sys.executable, str(dd / "make_case.py"), "--src", st["state"], "--out", str(case),
                    "--iters", str(int(d.get("primal_iters", 300)))], check=True)
    env = dict(os.environ, **DAFOAM_DEFAULT_ENV, **{k: str(v) for k, v in (d.get("env") or {}).items()})
    t0 = time.time()
    if os.environ.get("DAFOAM_RUNNER", "docker") == "native":
        shutil.copy(dd / "runScript.py", case / "runScript.py")
        sh = f""". {os.environ.get('DAFOAM_LOAD', '/home/dafoamuser/dafoam/loadDAFoam.sh')} >/dev/null; set -e; cd {case}
renumberMesh -overwrite > log.renumberMesh 2>&1
foamFormatConvert -constant > log.formatConvert 2>&1
rm -rf processor*
mpirun -np {n_proc} python runScript.py > log.dafoam 2>&1"""
        r = subprocess.run(["bash", "-c", sh], env=env)
    else:
        with open(case / "log.dafoam", "w") as f:
            r = subprocess.run(["bash", str(dd / "run_docker.sh"), str(case), str(n_proc)], env=env, stdout=f, stderr=subprocess.STDOUT)
    txt = (case / "log.dafoam").read_text(errors="replace") if (case / "log.dafoam").exists() else ""
    rec = {"returncode": r.returncode, "dafoam_s": round(time.time() - t0), "env": {k: env[k] for k in DAFOAM_DEFAULT_ENV},
           **dafoam_health(txt)}
    sj = case / "dafoam_sens.json"
    rec["ok"] = r.returncode == 0 and (case / "dafoam_sens.npz").exists() and not rec["failed"]
    if sj.exists():
        rec.update(json.loads(sj.read_text()))
    rec["sens"] = str((case / "dafoam_sens.npz").resolve())
    dj.write_text(json.dumps(rec, indent=1, default=float))
    return rec


def step_link(job, lc: dict, daf: dict, work: Path) -> dict:
    import freeform as fr
    import link
    gj = work / "gradient.json"
    if gj.exists():
        return json.loads(gj.read_text())
    d = lc["decoded"]
    ff = fr.Freeform.from_dict(d["freeform"])
    frame = link.Frame.from_layout(d["W"], d["x_front"], d["rear"], d["freeform"].get("nose_foam", False), d.get("skin_mm", 0.0))
    p = job.get("parametrisation") or {}
    seed = fr.Freeform.from_dict(p["freeform_seed"]) if p.get("freeform_seed") else None
    z = np.load(daf["sens"])
    g = link.gradient(z["dD_dXs"], z["Xs"], ff, frame, seed=seed, v=p.get("v"), amp=p.get("amp_mm"),
                      coarse=p.get("coarse"), tol_mm=float((job.get("link") or {}).get("tol_mm", link.TOL_MM)))
    gj.write_text(json.dumps(g, indent=1))
    return g


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("job")
    ap.add_argument("work")
    ap.add_argument("--adjoint-root", required=True)
    ap.add_argument("--np", type=int, default=int(os.environ.get("NP", os.cpu_count() or 4)))
    a = ap.parse_args()
    job = json.loads(Path(a.job).read_text())
    job["_dir"] = str(Path(a.job).resolve().parent)
    work = Path(a.work).resolve(); work.mkdir(parents=True, exist_ok=True)
    adj = Path(a.adjoint_root).resolve()
    res = {"id": job.get("id"), "np": a.np, "started": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    t = time.time()
    try:
        log("1 car"); cfd, lc = step_car(job, work); res["car_s"] = round(time.time() - t)
        log("2 production"); prod = step_production(job, cfd, work, adj, a.np)
        res["production"] = {k: prod.get(k) for k in ("D_half_mN", "D_half_se_mN", "L_half_mN", "cells", "y_plus_max", "cfd_s")}
        log("3 steady"); st = step_steady(job, prod, work, a.np)
        res["steady"] = {k: st[k] for k in ("steady", "scheme", "D_half_mN", "span_mark_mN")}
        res["steady"]["span_mN"] = st["rungs"][-1]["span_mN"]
        if job.get("stop_after") == "steady":         # a forward check of the steady model (end-to-end test)
            res["usable"] = bool(st["steady"])
            raise _Stop()
        log("4 DAFoam"); daf = step_dafoam(job, st, work, adj, a.np)
        res["dafoam"] = {k: daf.get(k) for k in ("ok", "failed", "D_half_N", "primal_D_half_N", "adjoint_s", "primal_s",
                                                  "dafoam_s", "n_surface_points", "first_residual", "last_residual", "env")}
        if daf.get("D_half_N") and st.get("D_half_mN"):
            res["dafoam"]["primal_vs_openfoam_pct"] = 100 * (1e3 * daf["D_half_N"] - st["D_half_mN"]) / st["D_half_mN"]
        if not daf.get("ok"):
            raise SystemExit("DAFoam adjoint failed (log: dafoam_case/log.dafoam)")
        log("5 link"); g = step_link(job, lc, daf, work)
        res["gradient"] = g
        res["usable"] = bool(st["steady"] and daf.get("ok"))
    except _Stop:
        pass
    except SystemExit as e:
        res["error"] = str(e)
        res["usable"] = False
    res["wall_s"] = round(time.time() - t)
    (work / "result.json").write_text(json.dumps(res, indent=1, default=float))
    log(f"done: usable={res['usable']} {res.get('error', '')}")


if __name__ == "__main__":
    main()
