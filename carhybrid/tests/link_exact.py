"""The gradient map's exact test, no CFD (the other agent's test_link.py, run on carhybrid/grad.py's map).

With DAFoam's dD/dX replaced by f(x) A n (the body surface's area vector), grad.py's chain must return
d/du of the integral of f over the body. That integral is measured by rebuilding the body (param_body.build, the build
the searches use) at u_k +- delta. f = 1 (the volume), x and z (first moments: they test WHERE the map moves each number).

    python carhybrid/tests/link_exact.py <evolution checkout> <layout_car.json> [--genes 7,19,33] [--delta 0.25]
                                         [--out result.json] [--cache dir]

Marks (the other agent's, fixed before their first run): for every number and f, predicted within 10 % of the rebuilt
value (or within 2 % of the largest |value| for that f) and the same sign. Reported with the held surface excluded
(grad.py's default) and without it (--no-hold's map), so the fix's effect is measured, not assumed.
"""
from __future__ import annotations

import argparse, json, sys, time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
import grad as G  # noqa: E402

CACHE = None


def build_body(d, ff, key):
    import trimesh
    if CACHE and (CACHE / f"{key}.stl").exists():
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
    half.apply_scale(1e3)
    if CACHE:
        CACHE.mkdir(parents=True, exist_ok=True)
        half.export(str(CACHE / f"{key}.stl"))
    return half


def moments(m):
    V, c = m.volume, m.center_mass
    return np.array([2 * V, 2 * V * c[0], 2 * V * c[2]])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("evo"); ap.add_argument("rec")
    ap.add_argument("--genes", default="7,19,33")
    ap.add_argument("--delta", type=float, default=0.25)
    ap.add_argument("--out"); ap.add_argument("--cache")
    a = ap.parse_args()
    global CACHE
    CACHE = Path(a.cache) if a.cache else None
    fr, SL, SS = G.evo_import(a.evo)
    d = json.loads(Path(a.rec).read_text())["decoded"]
    seed = fr.Freeform.from_dict(d["freeform"])
    grid, amp = tuple(SS.COARSE), SS.AMP_MM
    n = fr.n_genes(*grid)
    axles = (d["x_front"], d["x_front"] + d["W"], 14.12)
    if "sidepod_rules.py S1" in (Path(fr.__file__).parent / "param_body.py").read_text():
        axles = None
    hk = {"halo": fr.halo_span(d["x_front"] - 16.0, d["d_halo"])} if hasattr(fr, "halo_span") else {}
    F = G.FieldMap(fr, seed, np.zeros(n), amp, grid, axles, hk)
    t0 = time.time()
    base = build_body(d, fr.from_unit(seed, np.zeros(n), amp, *grid), "base")
    print(f"base body {time.time() - t0:.0f} s, {len(base.faces)} faces", flush=True)
    c, A = base.triangles_center, base.area_faces
    keep = np.abs(c[:, 1]) > 1e-3
    c, A = c[keep], A[keep]
    fs = {"V": np.ones(len(c)), "x": c[:, 0], "z": c[:, 2]}
    held = G.HeldMasks(d)
    pred = {}
    for tag, h in (("hold", held), ("no_hold", None)):
        F.held = h
        mask, Q, _, gn = F.on_body(c)
        W = np.stack([A[mask] * f[mask] for f in fs.values()], 1)
        pred[tag] = 2.0 * F.speeds_dot(Q, gn, W)             # (n, 3), whole body
        print(f"{tag}: {int(mask.sum())} of {len(c)} faces move with u", flush=True)
    rows, verdict = [], []
    for k in [int(g) for g in a.genes.split(",")]:
        e = np.zeros(n); e[k] = a.delta
        t1 = time.time()
        mp = moments(build_body(d, fr.from_unit(seed, e, amp, *grid), f"g{k}_p{a.delta}"))
        mm = moments(build_body(d, fr.from_unit(seed, -e, amp, *grid), f"g{k}_m{a.delta}"))
        rows.append({"gene": k, "measured": ((mp - mm) / (2 * a.delta)).tolist(),
                     **{f"predicted_{t}": pred[t][k].tolist() for t in pred}, "build_s": round(time.time() - t1)})
        print(json.dumps(rows[-1]), flush=True)
    for t in pred:
        for j, nm in enumerate(fs):
            big = max(abs(r["measured"][j]) for r in rows)
            for r in rows:
                m_, p_ = r["measured"][j], r[f"predicted_{t}"][j]
                ok = (np.sign(m_) == np.sign(p_) or abs(m_) < 0.02 * big) and \
                     (abs(p_ - m_) <= 0.10 * abs(m_) or abs(p_ - m_) <= 0.02 * big)
                verdict.append({"map": t, "gene": r["gene"], "f": nm, "measured": m_, "predicted": p_,
                                "ratio": p_ / m_ if m_ else None, "pass": bool(ok)})
    for v in verdict:
        print(f"{v['map']:8s} gene {v['gene']:2d} {v['f']}: measured {v['measured']:12.1f} predicted "
              f"{v['predicted']:12.1f} ratio {v['ratio'] if v['ratio'] is None else round(v['ratio'], 3)} "
              f"{'ok' if v['pass'] else 'FAIL'}")
    res = {"rows": rows, "verdict": verdict,
           "PASS_hold": all(v["pass"] for v in verdict if v["map"] == "hold"),
           "PASS_no_hold": all(v["pass"] for v in verdict if v["map"] == "no_hold")}
    print("hold:", "PASS" if res["PASS_hold"] else "FAIL", " no_hold:", "PASS" if res["PASS_no_hold"] else "FAIL")
    if a.out:
        Path(a.out).write_text(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    main()
