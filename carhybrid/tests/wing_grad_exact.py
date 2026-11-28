"""wing_grad.py's exact test, no CFD: with DAFoam's dD/dX replaced by the area vectors of points spread over the wing and
its struts (dD/dX_i = A_i n_i), the chain must return d(volume)/ds_j of the rebuilt parts, which the same builds give
directly. Tests the point selection, the normals, the signed distance's sign and the factor 2 / units, number by number.

    python carhybrid/tests/wing_grad_exact.py <evolution checkout> <car folder (make_car.py's out, W8-patched builder)>
                                             [--names rw_chord,rw_rise,...] [--n 20000] [--out result.json]

Mark: for each number without a branch change, |predicted - rebuilt| <= 10 % of |rebuilt| or <= 2 % of the largest
|rebuilt| over the numbers, and the same sign (as test_link's marks).
"""
from __future__ import annotations

import json, pickle, sys, tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))


def main():
    a = sys.argv[1:]
    opt = lambda k: a[a.index(k) + 1] if k in a else None
    evo, car = a[0], Path(a[1])
    import grad as G
    G.evo_import(evo)
    import wing_grad as WG
    import rear_wing_engine as rwe
    import trimesh
    ctx = pickle.loads((car / "evo" / "parts" / "rwing_ctx.pkl").read_bytes())
    base = WG.Built(rwe, ctx["design"], ctx)
    n = int(opt("--n") or 20000)
    P, A, N = [], [], []
    for q in base.parts:                                   # area-weighted points on each part, own normals
        k = max(200, int(n * q.area / sum(p.area for p in base.parts)))
        pts, tri = trimesh.sample.sample_surface_even(q, k, seed=1)
        P.append(pts); N.append(q.face_normals[tri]); A.append(np.full(len(pts), q.area / len(pts)))
    P, N, A = np.concatenate(P), np.concatenate(N), np.concatenate(A)
    # drop points inside another part (the strut's root inside the wing): not on the union's surface
    inside = np.zeros(len(P), bool)
    for q in base.parts:
        inside |= trimesh.proximity.signed_distance(q, P) > 0.05
    P, N, A = P[~inside], N[~inside], A[~inside]
    with tempfile.TemporaryDirectory() as td:
        npz = Path(td) / "sens.npz"
        np.savez(npz, Xs=P * 1e-3, dD_dXs=N * A[:, None], D_N=0.18)          # 'drag' = volume (mm^3 per mm)
        names = opt("--names").split(",") if opt("--names") else None
        out = WG.main(evo, car, str(npz), str(Path(td) / "g.json"), names=names, quiet=True)
    rho = None
    rows, big = [], 0.0
    import component_contract as cc, assembly as AS
    rho = cc.DENSITY_G_CM3[AS.PART_MATERIAL] * 1e-3
    for nm, dD, dm in zip(out["names"], out["dD_ds_N"], out["dm_ds_g"]):
        if dD is None:
            rows.append({"name": nm, "branch_change": True}); continue
        pred = dD / 2.0 * 1e3 * 2.0                          # whole-car 'N' -> mm^3 per unit: 2 x sum A V x 1e-3
        meas = dm / rho                                      # 2 x d(volume)/ds, mm^3 per unit
        rows.append({"name": nm, "predicted_mm3": pred, "rebuilt_mm3": meas}); big = max(big, abs(meas))
    ok_all = True
    for r in rows:
        if "predicted_mm3" not in r:
            print(f"{r['name']:20s} branch change (no derivative)"); continue
        p_, m_ = r["predicted_mm3"], r["rebuilt_mm3"]
        ok = (np.sign(p_) == np.sign(m_) or abs(m_) < 0.02 * big) and (abs(p_ - m_) <= 0.1 * abs(m_) or abs(p_ - m_) <= 0.02 * big)
        r["pass"] = bool(ok); ok_all &= ok
        print(f"{r['name']:20s} predicted {p_:10.2f} rebuilt {m_:10.2f} mm3/unit ratio "
              f"{(p_ / m_ if m_ else float('nan')):6.3f} {'ok' if ok else 'FAIL'}")
    print("PASS" if ok_all else "FAIL")
    if opt("--out"):
        Path(opt("--out")).write_text(json.dumps({"rows": rows, "PASS": bool(ok_all)}, indent=1, default=float))


if __name__ == "__main__":
    main()
