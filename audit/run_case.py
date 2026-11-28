"""One drag-audit case on a GitHub runner, scored exactly as seed evolution scores (seedevo/campaign_case.py):
fresh wall-resolved snappyHexMesh, meshing frame k (domain shifted k x 0.37 mm), production scheme (linearUpwind),
cold start, 1000 iterations, score = mean of iterations 200-1000 (four 200-iteration blocks, SE from the blocks).

    python audit/run_case.py <adjoint repo checkout> <case json> <out dir>

case: {"name", "design": team_rm | hybrid_final | seedevo:<id>, "k", "bg_div" (12 = the scoring mesh), "turb"}
- team_rm:       the team pod surface the hybrid loop morphed from (gate/original_rm, rebuilt from gate_results/adjoint_*.npy)
- hybrid_final:  hybrid/final_design/car_half.stl (team_rm + the six accepted chained moves)
- seedevo:<id>:  a design of seedevo/campaign_poc, rebuilt from its stored field by seedevo.campaign build
bg_div sets the background cell (car length / bg_div); every refinement level, wake box and layer stack scales with it,
except the first-layer height (fixed y+ target). It needs the one-line env hook the workflow adds to openfoam_case.py.
"""
from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np


def forces(path):
    rows = [[float(x) for x in re.sub(r"[()]", " ", l).split()] for l in open(path) if l.strip() and not l.startswith("#")]
    return np.array(rows)


def score(Fx_mN):
    b = np.asarray(Fx_mN)[-800:].reshape(4, 200).mean(1)
    return float(b.mean()), float(b.std(ddof=1) / 2), [float(x) for x in b]


def main(adj, case, out):
    adj, out = Path(adj).resolve(), Path(out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    sp = adj / "research" / "sidepod_opt"
    c = json.loads(case) if case.strip().startswith("{") else json.loads(Path(case).read_text())
    k, bg, turb = int(c["k"]), float(c.get("bg_div", 12)), c.get("turb", "kOmegaSST")
    rec = {"case": c, "converged": False, "adjoint_repo_sha": os.environ.get("ADJ_SHA", ""),
           "protocol": "fresh wall-resolved mesh, production linearUpwind, cold, 1000 its, score its 200-1000"}
    py = sys.executable
    design = c["design"]
    if design.startswith("seedevo:"):
        i = design.split(":")[1]
        car = out / "design"
        r = subprocess.run([py, "-W", "ignore", "-m", "seedevo.campaign", "build", "seedevo/campaign_poc", i, str(car)],
                           cwd=sp, capture_output=True, text=True)
        (out / "build.log").write_text(r.stdout[-5000:] + r.stderr[-5000:])
        if r.returncode != 0:
            rec["error"] = "build: " + (r.stdout + r.stderr)[-400:]
            (out / "result.json").write_text(json.dumps(rec, indent=1)); return 1
        folder = car / "car"
    elif design.startswith("folder:"):                       # an already-built car folder (seedhybrid confirm stage)
        folder = Path(design.split(":", 1)[1])
    elif design == "hybrid_region":                          # the hybrid's in-region moves only (audit/region_only.py)
        folder = sp / "gate" / "hybrid_region"
        shutil.copytree(sp / "gate" / "original_rm", folder, ignore=shutil.ignore_patterns("_body_half.stl", "car_half.stl"))
        r = subprocess.run([py, "-W", "ignore", str(Path(__file__).resolve().parent / "region_only.py"), str(adj),
                            str(folder / "car_half.stl")], capture_output=True, text=True)
        (out / "region_only.log").write_text(r.stdout[-3000:] + r.stderr[-3000:])
        if r.returncode != 0:
            rec["error"] = "region_only: " + (r.stdout + r.stderr)[-400:]
            (out / "result.json").write_text(json.dumps(rec, indent=1)); return 1
    elif design.startswith("legal:"):                        # a morphed design put back through the legal builder
        name = design.split(":", 1)[1]
        if name == "hybrid_region":
            src = out / "region_only.stl"
            r = subprocess.run([py, "-W", "ignore", str(Path(__file__).resolve().parent / "region_only.py"), str(adj), str(src)],
                               capture_output=True, text=True)
            (out / "region_only.log").write_text(r.stdout[-3000:] + r.stderr[-3000:])
        else:
            src = {"hybrid_final": sp / "hybrid" / "final_design" / "car_half.stl",
                   "team_rm": sp / "gate" / "original_rm" / "car_half.stl"}[name]
        r = subprocess.run([py, "-W", "ignore", str(Path(__file__).resolve().parents[1] / "seedhybrid" / "legalise.py"),
                            str(adj), str(src), str(out / "legal"), design.split(":", 1)[1]], capture_output=True, text=True)
        (out / "legalise.log").write_text(r.stdout[-5000:] + r.stderr[-5000:])
        rec["legalise"] = json.loads((out / "legal" / "legalise.json").read_text()) if (out / "legal" / "legalise.json").exists() else None
        if r.returncode != 0:
            rec["error"] = "legalise: " + (r.stdout + r.stderr)[-400:]
            (out / "result.json").write_text(json.dumps(rec, indent=1, default=float)); return 1
        folder = out / "legal" / "car"
    else:
        folder = sp / "gate" / design
        if design == "hybrid_final":
            shutil.copytree(sp / "gate" / "original_rm", folder, ignore=shutil.ignore_patterns("_body_half.stl", "car_half.stl"))
            shutil.copy(sp / "hybrid" / "final_design" / "car_half.stl", folder / "car_half.stl")
        elif design != "team_rm":
            raise SystemExit(f"unknown design {design}")
        else:
            folder = sp / "gate" / "original_rm"
    env = dict(os.environ, GATE_DOMAIN_DX_MM=f"{k * 0.37:.4f}", PRODCHECK_WALL_RESOLVED="1", PRODCHECK_RES="resolved",
               PRODCHECK_RESULTS=str(out), PRODCHECK_KEEP=str(out / "state"), PRODCHECK_ITERS="1000",
               PRODCHECK_TIMEOUT=str(c.get("timeout_s", 19000)), P2_BG_DIV=f"{bg:g}", PRODCHECK_TURB=turb)
    t = time.time()
    r = subprocess.run([py, "-W", "ignore", "-u", "prodcheck.py", str(folder), c["name"], "4", "cold"], cwd=sp, env=env,
                       capture_output=True, text=True)
    rec["wall_s"] = round(time.time() - t, 1)
    (out / "prodcheck.log").write_text(r.stdout[-20000:] + r.stderr[-20000:])
    st = out / "state"
    for sub in ("logs", "postProcessing"):
        if (st / sub).exists():
            shutil.copytree(st / sub, out / sub, dirs_exist_ok=True)
    pj = out / f"prodcheck_{c['name']}.json"
    if pj.exists():
        j = json.loads(pj.read_text())
        rec.update(prodcheck=j, cells=j.get("n_cells"), y_plus_min=j.get("y_plus_min"), y_plus_max=j.get("y_plus_max"),
                   solver_converged=j.get("converged"))
    fs = sorted((st / "postProcessing" / "forces").glob("*/force.dat"))
    if r.returncode != 0 or not fs:
        rec["error"] = "solve: " + (r.stdout + r.stderr)[-400:]
        (out / "result.json").write_text(json.dumps(rec, indent=1)); return 1
    F = forces(fs[-1])
    rec["iterations_run"] = int(len(F))
    if len(F) < 1000:
        rec["error"] = f"only {len(F)} iterations"
        (out / "result.json").write_text(json.dumps(rec, indent=1)); return 1
    m, se, blocks = score(F[:, 1] * 1e3)                       # total Fx (half car, N -> mN), as campaign_case
    rec.update(drag_mN=m, drag_block_se_mN=se, drag_blocks_mN=blocks,
               lift_mN=float(F[-800:, 3].mean() * 1e3) if F.shape[1] > 3 else None, converged=True)
    for name in ("forces_car",):                               # the body alone (wheels = total - body)
        fc = sorted((st / "postProcessing" / name).glob("*/force.dat"))
        if fc:
            Fc = forces(fc[-1])
            rec["body_drag_mN"] = float(Fc[-800:, 1].mean() * 1e3)
            rec["wheels_drag_mN"] = m - rec["body_drag_mN"]
    shutil.rmtree(st, ignore_errors=True)
    (out / "result.json").write_text(json.dumps(rec, indent=1))
    print(json.dumps({k: rec.get(k) for k in ("case", "drag_mN", "drag_block_se_mN", "cells", "y_plus_max", "wall_s")}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(*sys.argv[1:4]))
