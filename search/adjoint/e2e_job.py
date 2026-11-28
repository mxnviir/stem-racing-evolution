"""Job files for the end-to-end test (E2E_TEST.md) and for a stage-2 car in general.

    python e2e_job.py base <out.json> [--res coarse] [--wall 0]               # the team car, v = 0, full chain
    python e2e_job.py fd <gradient.json> <out dir> [--res coarse] [--wall 0]  # top-3 |dD/dv| numbers at +-delta
                                                                              # (stop after steady); prints the list
"""
from __future__ import annotations

import sys as _sys, pathlib as _pl  # the repo root, then every code folder (repo_paths.py)
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2]))
import repo_paths  # noqa: E402,F401
import argparse
import json
from pathlib import Path

import numpy as np

import freeform as fr
import layout_genes as LG
import simscale_shape as SS

COARSE_STEADY = {"extra": 2000, "span_mN": 0.02}
COARSE_DAFOAM = {"primal_iters": 300, "env": {"DAF_ADJ": "Krylov", "DAF_PC_FILL": "0", "DAF_PC_LEVEL": "1",
                                               "DAF_GMRES_RESTART": "200", "DAF_GMRES_ITERS": "2000"}}
TARGET_MN, DELTA = 3.0, (0.1, 0.6)


def car_job(jid, genes, seed: dict, v, res="coarse", wall=False, stop_after=None):
    """A fine_adjoint job for stage-2 design v on layout `genes` (as ShapeCampaign.job builds the car)."""
    n = fr.n_genes(*SS.COARSE)
    v = np.zeros(n + len(LG.SHAPE_GENES)) if v is None else np.asarray(v, float)
    ff = fr.from_unit(fr.Freeform.from_dict(seed), v[:n], SS.AMP_MM, *SS.COARSE)
    sh = LG.shape_from_unit(v[n:], stage1=np.asarray(genes, float))
    j = {"id": jid,
         "car": {"genes": list(map(float, genes)), "freeform": dict(ff.as_dict(), nose_foam=sh["nose_foam"]),
                 "support": sh["support"], "fwing": sh["fwing"], "rwing": sh["rwing"]},
         "parametrisation": {"freeform_seed": seed, "v": v.tolist(), "coarse": list(SS.COARSE), "amp_mm": list(SS.AMP_MM)},
         "mesh": {"res": res, "wall_resolved": bool(wall), "frame_k": 0},
         "steady": COARSE_STEADY if res == "coarse" else {"extra": 2000, "span_mN": 0.10},
         "dafoam": COARSE_DAFOAM, "link": {"tol_mm": 1.0}}
    if stop_after:
        j["stop_after"] = stop_after
    return j


def team():
    x = LG.x_default()
    seed = LG.decode(x)["freeform"]
    seed = {k: v for k, v in seed.items() if k in fr.Freeform.__dataclass_fields__}
    return x, seed


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd"); ap.add_argument("args", nargs="*")
    ap.add_argument("--res", default="coarse"); ap.add_argument("--wall", type=int, default=0)
    a = ap.parse_args()
    x, seed = team()
    if a.cmd == "base":
        Path(a.args[0]).write_text(json.dumps(car_job("e2e-base", x, seed, None, a.res, a.wall), indent=1))
    elif a.cmd == "fd":
        g = json.loads(Path(a.args[0]).read_text())
        gh = np.asarray(g["coarse_half_mN_per_unit"], float)
        top = [int(k) for k in np.argsort(-np.abs(gh))[:3]]
        out = Path(a.args[1]); out.mkdir(parents=True, exist_ok=True)
        n = fr.n_genes(*SS.COARSE)
        cases = []
        for k in top:
            d = float(np.clip(TARGET_MN / max(abs(gh[k]), 1e-9), *DELTA))
            for sgn, tag in ((1, "p"), (-1, "m")):
                v = np.zeros(n + len(LG.SHAPE_GENES)); v[k] = sgn * d
                jid = f"e2e-fd-g{k}{tag}"
                (out / f"{jid}.json").write_text(json.dumps(car_job(jid, x, seed, v, a.res, a.wall, "steady"), indent=1))
                cases.append({"id": jid, "gene": k, "sign": sgn, "delta": d, "pred_half_mN_per_unit": float(gh[k])})
        (out / "fd_cases.json").write_text(json.dumps(cases, indent=1))
        print(json.dumps([c["id"] for c in cases]))


if __name__ == "__main__":
    main()
