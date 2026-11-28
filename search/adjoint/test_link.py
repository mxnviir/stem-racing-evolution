"""The adjoint link's exact test, no CFD (link.py docstring).

With DAFoam's dD/dX_i replaced by f(x_i) A_i n_i (A_i n_i: the body surface's area vector at vertex i), the link
must return d/du of  integral_body f dV . That is measured here by rebuilding the body (param_body.build, the
very build the searches use, fixed skin) at u_k +- delta and integrating f over each rebuilt solid. f = 1 (the
volume), f = x and f = z (first moments: they test WHERE on the body the link puts each number's movement).

    python test_link.py <built car dir> [--genes 7,19,33] [--delta 0.25] [--out result.json]

Marks (fixed before the first run, 2026-11-28): for every tested number and every f, the predicted derivative
within 10 % of the rebuilt one (or within 2 % of the largest of the three |values| for that f, for a number
that barely moves the body), and the same sign. 10 % allows for the 1 mm build grid, the milling and the
inactive (rule-cut) points; a wrong frame, unit, sign, normal or chain would be off by far more.

Amendment (2026-11-28 14:35, after gene 7's result only; genes 19 and 33 unseen): gene 7 came out 1.30 / 1.36 / 1.46
x the rebuilt values. Diagnosis: a third of its predicted change sat on the spine (|y| < 13, z 24-43 mm), where the
fitted field lies on the chamber wall and the halo pocket and the build's hard masks hold the surface in both
directions (per-point displacements). The link now excludes points a mask holds (link.Masks, link.build_masks), with
the hold distance set on gene 7 alone. Genes 19 and 33 are the held-out check of that fix, under the same marks.
"""
from __future__ import annotations

import sys as _sys, pathlib as _pl  # the repo root, then every code folder (repo_paths.py)
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2]))
import repo_paths  # noqa: E402,F401
import argparse
import json
import time
from pathlib import Path

import numpy as np

import freeform as fr
import link

COARSE = (5, 2, 3, 4)
AMP_MM = (6.0, 4.0, 6.0)


CACHE = None                                                    # folder of built half solids (--cache)


def build_body(d, ff, key=None):
    """Part 1's body for a decoded layout and a Freeform, at the searches' fixed skin, as run_car.build_body does;
    returns the closed half solid (mm) and its full volume, first moments (mm^3, mm^4)."""
    import trimesh
    if CACHE and key and (CACHE / f"{key}.stl").exists():
        return trimesh.load(str(CACHE / f"{key}.stl"), force="mesh")
    import param_body as pb
    import unified_phi as up
    from coarse import use_spacing
    import layout_genes as LG
    use_spacing(1.0)
    bp = pb.BodyParams.from_dict(LG.body_params(d))
    g = pb.build(d["W"], d["x_front"], d["d_halo"], bp, skin_offset_mm=d["skin_mm"], rear_len_mm=d["rear"],
                 cargo_placement=d["cargo"], freeform=ff)
    half = up.extract_half_surface(g)
    half.vertices[half.vertices[:, 1] < 0, 1] = 0.0
    half.apply_scale(1e3)                                      # m -> mm
    if CACHE and key:
        CACHE.mkdir(parents=True, exist_ok=True)
        half.export(str(CACHE / f"{key}.stl"))
    return half


def moments(m):
    """Volume and first moments (x, z) of a closed half solid (mm), doubled to the whole body."""
    V = m.volume
    c = m.center_mass
    return np.array([2 * V, 2 * V * c[0], 2 * V * c[2]])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("car")
    ap.add_argument("--genes", default="7,19,33")
    ap.add_argument("--delta", type=float, default=0.25)
    ap.add_argument("--out", default=None)
    ap.add_argument("--cache", default=None, help="keep the built bodies here (a re-run reuses them)")
    a = ap.parse_args()
    global CACHE
    CACHE = Path(a.cache) if a.cache else None
    import layout_genes as LG
    r = json.loads((Path(a.car) / "layout_car.json").read_text())
    d = r["decoded"]
    seed = fr.Freeform.from_dict(d["freeform"])                   # the built car's body is the seed (u = 0)
    nf = d["freeform"].get("nose_foam", False)
    frame = link.Frame.from_layout(d["W"], d["x_front"], d["rear"], nf, d["skin_mm"])
    genes = [int(g) for g in a.genes.split(",")]
    n = fr.n_genes(*COARSE)
    t0 = time.time()
    base = build_body(d, fr.from_unit(seed, np.zeros(n), AMP_MM, *COARSE), "base")
    print(f"base body {time.time() - t0:.0f} s, {len(base.vertices)} vertices, watertight {base.is_watertight}", flush=True)
    fs = {"V": None, "x": lambda P: P[:, 0], "z": lambda P: P[:, 2]}
    # the half solid's cap on y = 0 is inside the body: its points are inactive (phi < level), as they must be
    masks = link.build_masks(d)                                   # the build's hard masks (amendment)
    pred = {}
    for name, f in fs.items():
        _, co, info = link.predicted_moment(base, seed, frame, f=f, seed=seed, v=np.zeros(n), amp=AMP_MM, coarse=COARSE,
                                            masks=masks)
        pred[name] = 2 * co                                     # half -> whole body
    print("active vertices", info, flush=True)
    rows = []
    for k in genes:
        e = np.zeros(n); e[k] = a.delta
        t1 = time.time()
        mp = moments(build_body(d, fr.from_unit(seed, e, AMP_MM, *COARSE), f"g{k}_p{a.delta}"))
        mm = moments(build_body(d, fr.from_unit(seed, -e, AMP_MM, *COARSE), f"g{k}_m{a.delta}"))
        meas = (mp - mm) / (2 * a.delta)
        rows.append({"gene": k, "measured": meas.tolist(), "predicted": [float(pred[nm][k]) for nm in fs],
                     "build_s": round(time.time() - t1)})
        print(json.dumps(rows[-1]), flush=True)
    verdict = []
    for j, nm in enumerate(fs):
        big = max(abs(r_["measured"][j]) for r_ in rows)
        for r_ in rows:
            m_, p_ = r_["measured"][j], r_["predicted"][j]
            ok = (np.sign(m_) == np.sign(p_) or abs(m_) < 0.02 * big) and \
                 (abs(p_ - m_) <= 0.10 * abs(m_) or abs(p_ - m_) <= 0.02 * big)
            verdict.append({"gene": r_["gene"], "f": nm, "measured": m_, "predicted": p_,
                            "ratio": p_ / m_ if m_ else None, "pass": bool(ok)})
    res = {"car": str(a.car), "delta": a.delta, "coarse": COARSE, "amp_mm": AMP_MM, "active": info, "rows": rows,
           "verdict": verdict, "PASS": all(v["pass"] for v in verdict)}
    for v in verdict:
        print(f"gene {v['gene']:2d} {v['f']}: measured {v['measured']:12.3f} predicted {v['predicted']:12.3f} "
              f"ratio {v['ratio'] if v['ratio'] is None else round(v['ratio'], 3)} {'ok' if v['pass'] else 'FAIL'}")
    print("PASS" if res["PASS"] else "FAIL")
    if a.out:
        Path(a.out).write_text(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    main()
