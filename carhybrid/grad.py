"""The DAFoam drag gradient -> the race-time gradient in stage 2's free-form numbers (carhybrid/DESIGN.md, decision 4).

    python carhybrid/grad.py <evolution checkout> <car record json> <dafoam_sens.npz> <u json or "zero"> <out json>
         [--body <car folder with body_half.stl>]

car record: layout_car.build's record (genes, decoded layout incl. "freeform" = the stage-2 seed, mass).
u: stage 2's coarse displacement numbers (45 = top 5x2 + bottom 5x4 + side 5x3, simscale_shape.COARSE / AMP_MM);
   the car the gradient was computed on is from_unit(seed, u).

For each number u_j the body+nose surface moves along its normal at speed V_j = -(dphi/du_j) / |grad phi| (mm per
unit), phi = freeform.field (explicit, so the derivative is a finite difference of a function, no building).
Drag: dD/du_j = 2 * sum_i (G_i . n_i) V_ij (half car -> whole car), over DAFoam's surface points on the body or
nose (|phi| < PHI_TOL_MM); the wings, struts, supports and wheels do not move with u.
Mass: dm/du_j = 2 * sum over the body's own half surface of V_j dA times the density (foam behind Ref A, printed nose
ahead of it); counted only while the car carries no lead (with lead, the ballast absorbs it: race mass fixed).
Race time: dT/du = dT/dD * dD/du + dT/dm * dm/du, the partials by finite differences of the evolution race model at
this car.
"""
import json, sys
from pathlib import Path

import numpy as np

PHI_TOL_MM = 1.0
EPS_U = 0.02                       # unit step for the finite difference of the field (0.08-0.12 mm of surface)
H_MM = 0.05                        # spatial step for grad phi
RHO_FOAM_G_MM3 = 0.163e-3
RHO_NOSE_G_MM3 = 1.04e-3           # printed ABS nose


def main(evo, rec_p, sens_p, u_s, out_p, body_dir=None):
    evo = Path(evo).resolve()
    sys.path.insert(0, str(evo)); import repo_paths; repo_paths.add()
    import freeform as fr
    import simscale_layout as SL
    import simscale_shape as SS
    rec = json.loads(Path(rec_p).read_text())
    d = rec["decoded"]
    seed = fr.Freeform.from_dict(d["freeform"])
    coarse, amp = SS.COARSE, SS.AMP_MM
    n = fr.n_genes(*coarse)
    u = np.zeros(n) if u_s == "zero" else np.asarray(json.loads(u_s) if u_s.strip().startswith("[") else json.loads(Path(u_s).read_text()), float)
    ref_a = d["x_front"] - 16.0
    axles = (d["x_front"], d["x_front"] + d["W"], 14.12)

    def phi(uu, P):
        ff = fr.from_unit(seed, uu, amp, *coarse)
        return fr.field(ff, P[:, 0], P[:, 1], P[:, 2], ff.x_tip_mm, ff.x_rear_mm, axles=axles)

    S = np.load(sens_p)
    X, G = S["Xs"] * 1e3, S["dD_dXs"]                       # mm; N/m per point (dD = sum G . dX, dX in m)
    p0 = phi(u, X)
    on = np.abs(p0) < PHI_TOL_MM
    Xo, Go = X[on], G[on]
    grad = np.stack([(phi(u, Xo + H_MM * e) - phi(u, Xo - H_MM * e)) / (2 * H_MM) for e in np.eye(3)], 1)
    gn = np.linalg.norm(grad, axis=1)
    ok = gn > 1e-6
    Xo, Go, grad, gn = Xo[ok], Go[ok], grad[ok], gn[ok]
    nrm = grad / gn[:, None]                                 # outward (phi increases outward)
    Gn = (Go * nrm).sum(1)                                   # N/m per point along the outward normal
    pb = phi(u, Xo)
    V = np.empty((n, len(Xo)))
    for j in range(n):
        e = np.zeros(n); e[j] = EPS_U
        V[j] = -(phi(u + e, Xo) - pb) / EPS_U / gn          # mm of outward surface move per unit u_j
    dD = 2.0 * (V * Gn[None, :]).sum(1) * 1e-3               # whole car, N per unit u_j
    # mass: the body's own half surface (area-weighted), both halves
    dm = np.zeros(n)
    if body_dir:
        import trimesh
        b = trimesh.load(Path(body_dir) / "body_half.stl", force="mesh")
        if b.bounds[1][0] < 1.0:                             # metres -> mm
            b.apply_scale(1e3)
        nose_p = Path(body_dir) / "parts" / "nose.stl"      # the printed nose ahead of Ref A moves with the same field
        if nose_p.exists():
            nz = trimesh.load(nose_p, force="mesh")
            if nz.bounds[1][0] < 1.0:
                nz.apply_scale(1e3)
            nz = nz.slice_plane([0, 0, 0], [0, 1, 0]) if nz.bounds[0][1] < -1.0 else nz     # one side, as body_half
            b = trimesh.util.concatenate([b, nz])
        c, A = b.triangles_center, b.area_faces
        keep = (np.abs(c[:, 1]) > 1e-3) & (np.abs(phi(u, c)) < PHI_TOL_MM)       # not the closing face at y = 0
        c, A = c[keep], A[keep]
        gc = np.stack([(phi(u, c + H_MM * e) - phi(u, c - H_MM * e)) / (2 * H_MM) for e in np.eye(3)], 1)
        g_n = np.maximum(np.linalg.norm(gc, axis=1), 1e-6)
        pc = phi(u, c)
        rho = np.where(c[:, 0] >= ref_a, RHO_FOAM_G_MM3, RHO_NOSE_G_MM3)
        for j in range(n):
            e = np.zeros(n); e[j] = EPS_U
            Vc = -(phi(u + e, c) - pc) / EPS_U / g_n
            dm[j] = 2.0 * (Vc * A * rho).sum()
    m = rec["mass"]
    lead = float(m.get("lead_g", 0.0) or 0.0)
    D20, L20 = float(rec.get("D20_N", 0.36)), float(rec.get("L20_N", 0.0))
    T0 = SL.race_time(rec, D20, L20)
    dT_dD = (SL.race_time(rec, D20 + 1e-3, L20) - SL.race_time(rec, D20 - 1e-3, L20)) / 2e-3      # s per N
    r2 = json.loads(json.dumps(rec)); r2["mass"]["race_mass_g"] = m["race_mass_g"] + 0.1
    dT_dm = (SL.race_time(r2, D20, L20) - T0) / 0.1 if lead <= 0 else 0.0                          # s per g
    dT = dT_dD * dD + dT_dm * dm
    gN = float(np.linalg.norm(dT))
    out = {"n_points_on_body": int(len(Xo)), "n_points_total": int(len(X)), "u": u.tolist(),
           "D20_N_used": D20, "T0_s": T0, "dT_dD_s_per_N": dT_dD, "dT_dm_s_per_g": dT_dm, "lead_g": lead,
           "dD_du_N": dD.tolist(), "dm_du_g": dm.tolist(), "dT_du_s": dT.tolist(), "norm_dT_du_s": gN,
           "largest_share": float(np.max(np.abs(dT)) / gN) if gN > 0 else None,
           "direction": (-dT / gN).tolist() if gN > 0 else None,
           "predicted_dT_ms_per_unit_step": -gN * 1e3}
    Path(out_p).write_text(json.dumps(out, indent=1))
    print(json.dumps({k: out[k] for k in ("n_points_on_body", "T0_s", "dT_dD_s_per_N", "dT_dm_s_per_g", "norm_dT_du_s",
                                          "largest_share")}))


if __name__ == "__main__":
    a = sys.argv[1:]
    body = a[a.index("--body") + 1] if "--body" in a else None
    main(*a[:5], body_dir=body)
