"""The layout search's CFD on the team's own OpenFOAM pipeline (owner 2026-11-28: "OpenFOAM on Actions"), not SimScale.

One case on one machine:

    python openfoam_eval.py case <x.json> <out dir> <adjoint repo> [k]

  1. builds the car from its genes (layout_car.build: Parts 1-5, the referee, the CFD solid);
  2. if it is legal, writes it as the pipeline's car folder (<out>/cfd_car: car_half.stl, an ASCII solid named
     `car`, the right half of the merged one-piece car; wheelF/wheelR.stl and assembly.json with their rotation);
  3. runs research/sidepod_opt/prodcheck.py of the adjoint repo (github.com/mxnviir/stem-racing-adjoint): a fresh
     wall-resolved snappyHexMesh (~2 M cells for the team car), k-omega SST, the production scheme (linearUpwind),
     1000 iterations from a cold start, meshing frame k (the domain shifted k x 0.37 mm: a re-run sees another mesh);
  4. scores drag and lift as the mean of iterations 200-1000 (four 200-iteration blocks; the first 200 of a cold
     start are transient), whole car (x 2), and the race time with this car's own layout and mass (race_time);
  5. writes <out>/result.json.

The pipeline's noise on this mesh (the adjoint repo's test T1, four frames of the team car): 0.57 mN half car.
"""
from __future__ import annotations

import sys as _sys, pathlib as _pl  # the repo root, then every code folder (repo_paths.py)
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2]))
import repo_paths  # noqa: E402,F401
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROUTE = "openfoam_v2412_resolved_prod1000_layoutframe"    # the setup every result of this route shares
FRAME_MM = 0.37                               # meshing frame k shifts the domain by k x this
# ONE fixed meshing frame for every layout car (the pipeline's own frame is the team car's, 257 mm): longest legal car
# is Ref A + the 223 mm block + the cartridge behind it (~276 mm from the nose tip); widest 42.5 mm (T3.4); the rear
# wing up to 65 mm (T9.4.3). Metres, car frame.
DOMAIN = [[-0.005, 0.0, 0.0], [0.285, 0.046, 0.070]]
ITERS = 1000
SKIP = 200                                    # cold-start transient left out


def forces(path):
    """(Fx, Fz) half car, N, per iteration, from an OpenFOAM forces function object's force.dat (total columns)."""
    rows = [[float(v) for v in re.sub(r"[()]", " ", l).split()] for l in open(path) if l.strip() and not l.startswith("#")]
    a = np.array(rows)
    return a[:, 1], a[:, 3]


def score(F):
    """Mean of iterations SKIP..ITERS and its block standard error (four blocks)."""
    x = np.asarray(F)[-(ITERS - SKIP):]
    b = x.reshape(4, -1).mean(1)
    return float(b.mean()), float(b.std(ddof=1) / 2.0)


def write_cfd_car(car_dir: Path, dest: Path) -> dict:
    """The built car as the pipeline's car folder: the merged solid's right half, the wheels and their rotation."""
    import trimesh
    dest.mkdir(parents=True, exist_ok=True)
    m = trimesh.load(str(car_dir / "body_full.stl"), force="mesh")            # metres, car frame, both sides
    half = m.slice_plane([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], cap=True)
    txt = trimesh.exchange.stl.export_stl_ascii(half)
    txt = re.sub(r"^solid.*$", "solid car", txt, count=1, flags=re.M)
    txt = re.sub(r"^endsolid.*$", "endsolid car", txt, flags=re.M)
    (dest / "car_half.stl").write_text(txt)
    asm = json.loads((car_dir / "parts" / "assembly.json").read_text())
    ext = []
    for s in asm["extra_surfaces"]:
        if s["name"] in ("wheelF", "wheelR"):
            (dest / f"{s['name']}.stl").write_bytes((car_dir / f"{s['name']}.stl").read_bytes())
            ext.append({"name": s["name"], "stl": f"{s['name']}.stl", "rotating": s["rotating"]})
    (dest / "assembly.json").write_text(json.dumps({"extra_surfaces": ext, "label": car_dir.name}, indent=1))
    return {"faces_half": int(len(half.faces)), "watertight": bool(half.is_watertight)}


def run_cfd(cfd_car: Path, out: Path, adjoint_root: Path, k: int, n_proc: int = 4, res: str = "resolved",
            wall_resolved: bool = True) -> dict:
    """res / wall_resolved: the route's mesh (resolved, wall-resolved) unless a test asks for another (the coarse
    end-to-end adjoint test, search/adjoint)."""
    sp = Path(adjoint_root) / "research" / "sidepod_opt"
    keep = out / "state_prod"
    env = dict(os.environ, GATE_DOMAIN=json.dumps(DOMAIN), GATE_DOMAIN_DX_MM=f"{k * FRAME_MM:.4f}", PRODCHECK_WALL_RESOLVED="1" if wall_resolved else "0", PRODCHECK_RES=res,
               PRODCHECK_RESULTS=str(out), PRODCHECK_KEEP=str(keep), PRODCHECK_ITERS=str(ITERS),
               PRODCHECK_TIMEOUT=os.environ.get("PRODCHECK_TIMEOUT", "16000"))
    env.pop("CFD_DIV_U", None)
    t0 = time.time()
    r = subprocess.run([sys.executable, "-W", "ignore", "-u", "prodcheck.py", str(cfd_car), f"{out.name}_k{k}", str(n_proc), "cold"],
                       cwd=sp, env=env, capture_output=True, text=True)
    (out / "prodcheck.log").write_text(r.stdout[-30000:] + r.stderr[-30000:])
    rec = {"cfd_s": round(time.time() - t0)}
    f = sorted(keep.glob("postProcessing/forces/*/force.dat"))
    if r.returncode != 0 or not f:
        rec["error"] = "cfd: " + (r.stdout + r.stderr)[-400:]
        return rec
    Fx, Fz = forces(f[-1])
    if len(Fx) < ITERS:
        rec["error"] = f"cfd: only {len(Fx)} iterations"
        return rec
    d, dse = score(Fx); l, lse = score(Fz)
    rec.update(D20_N=2 * d, L20_N=2 * l, D_half_mN=1e3 * d, D_half_se_mN=1e3 * dse, L_half_mN=1e3 * l, L_half_se_mN=1e3 * lse,
               iterations=int(len(Fx)))
    pj = next(out.glob("prodcheck_*.json"), None)
    if pj:
        j = json.loads(pj.read_text())
        rec.update(cells=j.get("n_cells"), y_plus_max=j.get("y_plus_max"), converged=j.get("converged"))
    return rec


def case(x_json, out, adjoint_root, k=0, cid=None):
    import layout_car as LC
    import simscale_layout as SL
    out = Path(out).resolve(); out.mkdir(parents=True, exist_ok=True)
    spec = json.loads(Path(x_json).read_text()) if len(str(x_json)) < 4096 and Path(x_json).exists() else json.loads(x_json)
    # a layout car (the genes), or a stage-2 car (simscale_shape.ShapeCampaign.job: genes + free-form body, wings,
    # supports)
    if isinstance(spec, dict):
        x = np.array(spec["x"], float)
        kw = {k: spec.get(k) for k in ("freeform", "support", "fwing", "rwing")}
    else:
        x, kw = np.array(spec, float), {}
    car = out / "car"
    rec = LC.build(x, car, **kw)
    res = {"id": cid or out.name, "k": int(k), "route": ROUTE, "genes": [float(v) for v in x], "legal": bool(rec.get("legal")),
           "v": spec.get("v") if isinstance(spec, dict) else None,
           "failed": rec.get("failed"), "why": (rec.get("why") or "")[:500], "build_s": rec.get("build_s"),
           "mass": rec.get("mass"), "layout": rec.get("layout")}
    if res["legal"]:
        res["cfd_car"] = write_cfd_car(car, out / "cfd_car")
        res.update(run_cfd(out / "cfd_car", out, Path(adjoint_root), int(k)))
        if res.get("D20_N") is not None:
            res["T"] = float(SL.race_time(rec, res["D20_N"], res["L20_N"]))
    (out / "result.json").write_text(json.dumps(res, indent=1, default=float))
    print(json.dumps({k_: res.get(k_) for k_ in ("id", "k", "legal", "failed", "D_half_mN", "L_half_mN", "T", "error")}, default=float))
    return res


if __name__ == "__main__":
    if sys.argv[1] == "case":
        case(sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]) if len(sys.argv) > 5 else 0,
             sys.argv[6] if len(sys.argv) > 6 else None)
