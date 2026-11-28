"""The DAFoam drag (and lift) gradient -> the race-time gradient in stage 2's free-form numbers (carhybrid/DESIGN.md,
decision 4; carhybrid/FINE_ADJOINT.md).

    python carhybrid/grad.py <evolution checkout> <car record json> <dafoam_sens.npz>[,<frame 2 npz>,...]
         <u json or "zero"> <out json>
         [--body <car folder with body_half.stl>] [--d20 <N>] [--l20 <N>] [--grid gu,gv,gw,gb] [--smooth <mm>]
         [--method lin|fd]

car record: layout_car.build's record (genes, decoded layout incl. "freeform" = the stage-2 seed, mass).
u: the free-form numbers on the grid --grid (default simscale_shape.COARSE = 5,2,3,4: 45 numbers = top 5x2 + bottom
   5x4 + side 5x3; "24,9,12,9" = the seed's own grid, 720 numbers, one per seed control point), mm per unit
   simscale_shape.AMP_MM; the car the gradient was computed on is freeform.from_unit(seed, u, AMP_MM, *grid).
Several npz files (frames of the same car, comma separated): each is mapped on its own and the mapped gradients are
averaged; the per-frame cosines with the mean (and pairwise) are reported.

For each number u_j the body+nose surface moves along its normal at speed V_j = -(dphi/du_j) / |grad phi| (mm per
unit), phi = freeform.field (explicit, so the derivative needs no building). --method lin (default): phi depends on u
only through the three interpolated surfaces (top, bottom, side heights at the point), each LINEAR in the grid, and
the grid is linear in u (freeform.from_unit, up to its clip of the side at 0); so dphi/du_j = (dphi/ds) (ds/du_j),
dphi/ds by a finite difference of the field in the surface (3 x 2 field calls), ds/du_j by interpolating the grid
response to u_j (one interpolation per number): 720 numbers on 10^5 points run in about a minute.
--method fd: the original finite difference of the whole field per number (one field call per number; kept as the
reference the tests compare lin against).
Drag: dD/du_j = 2 * sum_i (G_i . n_i) V_ij (half car -> whole car), over DAFoam's surface points on the body or
nose (|phi| < PHI_TOL_MM); the wings, struts, supports and wheels do not move with u.
Lift (if the npz has dL_dXs): dL/du_j exactly as drag, from dL_dXs.
Mass: dm/du_j = 2 * sum over the body's own half surface of V_j dA times the density (foam behind Ref A, printed nose
ahead of it); counted only while the car carries no lead (with lead, the ballast absorbs it: race mass fixed).
Race time: dT/du = dT/dD * dD/du + dT/dL * dL/du + dT/dm * dm/du, the partials by central finite differences of the
evolution race model (simscale_layout.race_time) at this car's own D20 and L20.

Lift sign convention (checked end to end, 2026-11-28):
- race model: simscale_layout.race_time(rec, D20, L20) -> objective.RaceModel -> race_objective's parameter
  lift_20_n = "lift at 20 m/s, N, + up" (whole car); the wheels carry m g - L, so L20 > 0 (upforce) lowers the
  rolling and guide friction: dT/dL20 < 0 (the team car: -2.0 ms per N, against dT/dD20 = +185 ms per N).
- score: audit/run_case.py's lift_mN = mean of column 3 of postProcessing/forces/*/force.dat = total_z (OpenFOAM
  v2412 columns: Time total_x total_y total_z ...), HALF car, in the car frame (x from the nose tip to the rear =
  the flow direction, z up from the track at z = 0: part2-simulation's domain clamps z_min to the track). So
  L20 = 2 * lift_mN / 1e3 with NO sign change (as simscale_step_eval's L20_N = Fz * scale). P0's team car: lift_mN
  +97 / +89 (frames 0 / 1): an upforce of ~0.19 N whole car.
- DAFoam: L is the force along +z (fixedDirection [0, 0, 1], scale rho: N, half car, on car + wheelF + wheelR as D),
  so dL_dXs is d(upforce)/dX in N per m, the same sign as lift_mN; whole car x 2 here, as for drag.
--l20 is therefore 2 x the score's lift_mN / 1e3 (whole car, + up).

--smooth <mm>: a Gaussian smoothing of the CFD parts of the gradient (dD/du, dL/du; not dm/du, which is exact) over
the control grid: each number's control point sits at (x, y) (top, bottom) or (x, z) (side) in mm (from_unit's
grid over the seed's x_tip..x_rear and freeform's Y_MAX_MM / Z_MAX_MM); g_s = K g / K 1 with
K_jk = exp(-|p_j - p_k|^2 / 2 s^2) within the same surface. Needed when the grid is finer than the adjoint's mesh
noise (the 720-point grid's spacing is ~11 mm along the car, 4-5 mm across); default off. The unsmoothed parts are
kept as *_raw.
"""
import dataclasses, json, sys, time
from pathlib import Path

import numpy as np

PHI_TOL_MM = 1.0
EPS_U = 0.02                       # unit step for the finite difference of the field (0.08-0.12 mm of surface)
H_MM = 0.05                        # spatial step for grad phi
DS_MM = 0.01                       # surface-height step for dphi/ds (method lin)
RHO_FOAM_G_MM3 = 0.163e-3
RHO_NOSE_G_MM3 = 1.04e-3           # printed ABS nose
SURF = ("top", "bottom", "side")


def evo_import(evo):
    """Put the pinned evolution checkout's code folders first on sys.path (its own repo_paths.py, loaded by file: this
    repo has a repo_paths module of its own, which a test session may already have imported) and return
    (freeform, simscale_layout, simscale_shape) from it."""
    import importlib, importlib.util
    evo = Path(evo).resolve()
    if str(evo) not in sys.path:
        sys.path.insert(0, str(evo))
    spec = importlib.util.spec_from_file_location("_evo_repo_paths", evo / "repo_paths.py")
    rp = importlib.util.module_from_spec(spec); spec.loader.exec_module(rp); rp.add()
    mods = [importlib.import_module(m) for m in ("freeform", "simscale_layout", "simscale_shape")]
    for m in mods:
        if not Path(m.__file__).resolve().is_relative_to(evo):
            raise ImportError(f"{m.__name__} came from {m.__file__}, not the evolution checkout {evo}")
    return mods


def parse_grid(s, SS):
    if s is None:
        return tuple(SS.COARSE)
    g = tuple(int(v) for v in str(s).replace(" ", "").split(","))
    return g if len(g) == 4 else g + (None,)


def cos(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    d = np.linalg.norm(a) * np.linalg.norm(b)
    return float(a @ b / d) if d > 0 else None


class FieldMap:
    """The free-form field at numbers u on a grid, and the normal speeds of its surface per number."""

    def __init__(self, fr, seed, u, amp, grid, axles, hk):
        self.fr, self.seed, self.u, self.amp, self.grid, self.axles, self.hk = fr, seed, u, amp, grid, axles, hk
        self.n = fr.n_genes(*grid)
        self.ff = fr.from_unit(seed, u, amp, *grid)
        self._resp = None

    def phi(self, ff, P):
        return self.fr.field(ff, P[:, 0], P[:, 1], P[:, 2], ff.x_tip_mm, ff.x_rear_mm, axles=self.axles, **self.hk)

    def phi_u(self, uu, P):
        return self.phi(self.fr.from_unit(self.seed, uu, self.amp, *self.grid), P)

    def on_body(self, P):
        """(mask of P on the body (|phi| < tol, grad phi > 0), points, outward normals, |grad phi|)."""
        p0 = self.phi(self.ff, P)
        idx = np.flatnonzero(np.abs(p0) < PHI_TOL_MM)
        Q = P[idx]
        grad = np.stack([(self.phi(self.ff, Q + H_MM * e) - self.phi(self.ff, Q - H_MM * e)) / (2 * H_MM)
                         for e in np.eye(3)], 1)
        gn = np.linalg.norm(grad, axis=1)
        ok = gn > 1e-6
        mask = np.zeros(len(P), bool); mask[idx[ok]] = True
        return mask, Q[ok], grad[ok] / gn[ok, None], gn[ok]

    def responses(self):
        """Per number j: (surface index, d grid / d u_j) of from_unit at u (central difference: exact for the linear
        map, half the slope exactly at the side's clip)."""
        if self._resp is None:
            fr, n = self.fr, self.n
            sizes = [a * b for a, b in fr.coarse_grid_shapes(*self.grid).values()]
            surf_of = np.repeat(np.arange(3), sizes)
            R = []
            for j in range(n):
                e = np.zeros(n); e[j] = EPS_U
                fp, fm = fr.from_unit(self.seed, self.u + e, self.amp, *self.grid), fr.from_unit(self.seed, self.u - e, self.amp, *self.grid)
                s = surf_of[j]
                R.append((s, (fp.grids()[s] - fm.grids()[s]) / (2 * EPS_U)))
            self._resp = R
        return self._resp

    def speeds_dot(self, Q, gn, W, method="lin"):
        """sum_i V_ij W_i (n x m): V_ij = -(dphi_i/du_j) / |grad phi_i| (mm per unit u_j), Q (N, 3) on the body."""
        W = np.asarray(W, float).reshape(len(Q), -1)
        out = np.zeros((self.n, W.shape[1]))
        if method == "fd":
            pb = self.phi(self.ff, Q)
            Wn = W / gn[:, None]
            for j in range(self.n):
                e = np.zeros(self.n); e[j] = EPS_U
                out[j] = -((self.phi_u(self.u + e, Q) - pb) / EPS_U) @ Wn
            return out
        fr, ff = self.fr, self.ff
        A = []                                               # -(dphi/ds) / |grad phi| * W, per surface
        for name in SURF:
            g = np.asarray(getattr(ff, name), float)
            fp = dataclasses.replace(ff, **{name: tuple(g + DS_MM)})
            fm = dataclasses.replace(ff, **{name: tuple(g - DS_MM)})
            dps = (self.phi(fp, Q) - self.phi(fm, Q)) / (2 * DS_MM)
            A.append(-(dps / gn)[:, None] * W)
        # freeform.surfaces' coordinates: a along the car, b across (top, bottom: |y| / Y_MAX) or up (side: z / Z_MAX)
        a = (Q[:, 0] - ff.x_tip_mm) / (ff.x_rear_mm - ff.x_tip_mm)
        b = (np.abs(Q[:, 1]) / fr.Y_MAX_MM, np.abs(Q[:, 1]) / fr.Y_MAX_MM, Q[:, 2] / fr.Z_MAX_MM)
        live = [np.any(A[s] != 0, axis=1) for s in range(3)]  # points the surface actually moves (smax weights)
        for j, (s, R) in enumerate(self.responses()):
            if not np.any(R):
                continue
            L = live[s]
            ds = fr._interp(R, a[L], b[s][L])                # linear in the grid: d surface / d u_j at the points
            out[j] = ds @ A[s][L]
        return out


def control_positions(fr, seed, grid):
    """(n, 2) mm positions of the numbers' control points ((x, y) top/bottom, (x, z) side) and their surface index."""
    pos, surf = [], []
    for s, (a, b) in enumerate(fr.coarse_grid_shapes(*grid).values()):
        span = fr.Z_MAX_MM if s == 2 else fr.Y_MAX_MM
        xs = seed.x_tip_mm + (np.arange(a) / (a - 1) if a > 1 else np.full(a, 0.5)) * (seed.x_rear_mm - seed.x_tip_mm)
        ys = (np.arange(b) / (b - 1) if b > 1 else np.full(b, 0.5)) * span
        X, Y = np.meshgrid(xs, ys, indexing="ij")
        pos.append(np.stack([X.ravel(), Y.ravel()], 1)); surf += [s] * (a * b)
    return np.concatenate(pos), np.asarray(surf)


def smoother(fr, seed, grid, s_mm):
    p, sf = control_positions(fr, seed, grid)
    d2 = ((p[:, None, :] - p[None, :, :]) ** 2).sum(-1)
    K = np.exp(-d2 / (2 * s_mm ** 2)) * (sf[:, None] == sf[None, :])
    return lambda g: (K @ np.asarray(g, float)) / K.sum(1)


def prolongation(fr, seed, coarse, fine, amp):
    """P (n_fine x n_coarse): the fine grid's numbers of the shape that coarse numbers make (exactly linear while the
    side is not clipped): the coarse space as a subspace of the fine one (same mm per unit). Fine must be the seed's
    own grid (one number per seed control point)."""
    nc = fr.n_genes(*coarse)
    g0 = np.concatenate([g.ravel() for g in fr.from_unit(seed, np.zeros(nc), amp, *coarse).grids()])
    amps = np.concatenate([np.full(a * b, amp[s]) for s, (a, b) in enumerate(fr.coarse_grid_shapes(*fine).values())])
    P = np.zeros((len(g0), nc))
    for j in range(nc):
        e = np.zeros(nc); e[j] = 1.0
        g = np.concatenate([g.ravel() for g in fr.from_unit(seed, e, amp, *coarse).grids()])
        P[:, j] = (g - g0) / amps
    return P


def load_u(u_s, n):
    if u_s == "zero":
        return np.zeros(n)
    u = np.asarray(json.loads(u_s) if u_s.strip().startswith("[") else json.loads(Path(u_s).read_text()), float)
    if len(u) != n:
        raise SystemExit(f"u has {len(u)} numbers, the grid {n}")
    return u


def main(evo, rec_p, sens_p, u_s, out_p, body_dir=None, d20=None, l20=None, grid=None, smooth=None, method="lin",
         quiet=False):
    t_start = time.time()
    fr, SL, SS = evo_import(evo)
    rec = json.loads(Path(rec_p).read_text())
    d = rec["decoded"]
    seed = fr.Freeform.from_dict(d["freeform"])
    grid, amp = parse_grid(grid, SS), SS.AMP_MM
    n = fr.n_genes(*grid)
    u = load_u(u_s, n)
    ref_a = d["x_front"] - 16.0
    axles = (d["x_front"], d["x_front"] + d["W"], 14.12)
    # the halo's visibility region (evo_patches/soft_cap.py v2): the builder's field compresses the top there too
    hk = {"halo": fr.halo_span(ref_a, d["d_halo"])} if hasattr(fr, "halo_span") else {}
    F = FieldMap(fr, seed, u, amp, grid, axles, hk)

    paths = [p for p in (sens_p if isinstance(sens_p, (list, tuple)) else str(sens_p).split(",")) if p]
    frames, dDs, dLs, D_N, L_N = [], [], [], [], []
    for p in paths:
        S = np.load(p)
        X = S["Xs"] * 1e3                                    # mm
        G = [S["dD_dXs"]] + ([S["dL_dXs"]] if "dL_dXs" in S.files else [])   # N/m per point (dF = sum G . dX, dX in m)
        mask, Q, nrm, gn = F.on_body(X)
        Wn = np.stack([(g[mask] * nrm).sum(1) for g in G], 1)                  # along the outward normal
        r = 2.0 * F.speeds_dot(Q, gn, Wn, method) * 1e-3     # whole car, N per unit u_j
        dDs.append(r[:, 0])
        if len(G) > 1:
            dLs.append(r[:, 1])
        fD = float(S["D_N"]) if "D_N" in S.files else None
        fL = float(S["L_N"]) if "L_N" in S.files else None
        D_N.append(fD); L_N.append(fL)
        frames.append({"npz": str(p), "n_points_on_body": int(mask.sum()), "n_points_total": int(len(X)),
                       "D_half_N": fD, "L_half_N": fL, "has_lift": len(G) > 1})
    have_lift = len(dLs) == len(paths)
    if 0 < len(dLs) < len(paths):
        raise SystemExit("some frames have dL_dXs and some do not")
    dD = np.mean(dDs, 0)
    dL = np.mean(dLs, 0) if have_lift else np.zeros(n)
    for f, a, b in zip(frames, dDs, dLs if have_lift else [None] * len(frames)):
        f["cos_dD_vs_mean"] = cos(a, dD)
        if b is not None:
            f["cos_dL_vs_mean"] = cos(b, dL)
    pair = [[cos(a, b) for b in dDs] for a in dDs]
    raw = {}
    if smooth:
        sm = smoother(fr, seed, grid, float(smooth))
        raw = {"dD_du_N_raw": dD.tolist(), "dL_du_N_raw": dL.tolist() if have_lift else None}
        dD, dL = sm(dD), (sm(dL) if have_lift else dL)
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
        keep = np.abs(c[:, 1]) > 1e-3                        # not the closing face at y = 0
        c, A = c[keep], A[keep]
        mask, Qc, _, gc = F.on_body(c)
        rho = np.where(Qc[:, 0] >= ref_a, RHO_FOAM_G_MM3, RHO_NOSE_G_MM3)
        dm = 2.0 * F.speeds_dot(Qc, gc, (A[mask] * rho)[:, None], method)[:, 0]
    m = rec["mass"]
    lead = float(m.get("lead_g", 0.0) or 0.0)
    # the race model's partials at this car's own drag and lift (whole car, N; lift + up): --d20 / --l20 (2 x the audit
    # protocol's half-car drag_mN / lift_mN / 1e3), else the record's, else DAFoam's own primal (2 x the frames' mean
    # D_N / L_N), else 0.36 N / 0 (P0 used 0.36 for a car measured at 0.467 N: dT/dD 0.178 vs 0.186 ms/mN, 4 % low)
    def pick(cli, key, daf, default):
        if cli is not None:
            return float(cli), "cli"
        if rec.get(key) is not None:
            return float(rec[key]), "rec"
        v = [x for x in daf if x is not None]
        if v and len(v) == len(daf):
            return 2.0 * float(np.mean(v)), "dafoam_primal"
        return default, "default"
    D20, D20_src = pick(d20, "D20_N", D_N, 0.36)
    L20, L20_src = pick(l20, "L20_N", L_N, 0.0)
    T0 = SL.race_time(rec, D20, L20)
    h = 1e-3
    dT_dD = (SL.race_time(rec, D20 + h, L20) - SL.race_time(rec, D20 - h, L20)) / (2 * h)       # s per N
    dT_dL = (SL.race_time(rec, D20, L20 + h) - SL.race_time(rec, D20, L20 - h)) / (2 * h) if have_lift else 0.0
    r2 = json.loads(json.dumps(rec)); r2["mass"]["race_mass_g"] = m["race_mass_g"] + 0.1
    dT_dm = (SL.race_time(r2, D20, L20) - T0) / 0.1 if lead <= 0 else 0.0                          # s per g
    tD, tL, tm = dT_dD * dD, dT_dL * dL, dT_dm * dm
    dT = tD + tL + tm
    gN = float(np.linalg.norm(dT))
    nrm_ = lambda v: float(np.linalg.norm(v))
    out = {"grid": [g for g in grid], "n": n, "method": method, "smooth_mm": float(smooth) if smooth else None,
           "n_points_on_body": int(np.mean([f["n_points_on_body"] for f in frames])),
           "n_points_total": int(np.mean([f["n_points_total"] for f in frames])),
           "n_frames": len(frames), "frames": frames, "cos_dD_frames": pair, "u": u.tolist(),
           "D20_N_used": D20, "D20_source": D20_src, "L20_N_used": L20, "L20_source": L20_src, "has_lift": have_lift,
           "T0_s": T0, "dT_dD_s_per_N": dT_dD, "dT_dL_s_per_N": dT_dL, "dT_dm_s_per_g": dT_dm, "lead_g": lead,
           "dD_du_N": dD.tolist(), "dL_du_N": dL.tolist() if have_lift else None, "dm_du_g": dm.tolist(),
           **raw, "dT_du_s": dT.tolist(), "norm_dT_du_s": gN,
           "norm_terms_s": {"drag": nrm_(tD), "lift": nrm_(tL), "mass": nrm_(tm)},
           "largest_share": float(np.max(np.abs(dT)) / gN) if gN > 0 else None,
           "direction": (-dT / gN).tolist() if gN > 0 else None,
           "predicted_dT_ms_per_unit_step": -gN * 1e3}
    out["run_s"] = round(time.time() - t_start, 2)
    Path(out_p).write_text(json.dumps(out, indent=1))
    if not quiet:
        print(json.dumps({k: out[k] for k in ("n", "n_frames", "n_points_on_body", "T0_s", "dT_dD_s_per_N", "dT_dL_s_per_N",
                                              "dT_dm_s_per_g", "norm_dT_du_s", "largest_share", "run_s")}))
    return out


if __name__ == "__main__":
    a = sys.argv[1:]
    opt = lambda k: a[a.index(k) + 1] if k in a else None
    main(*a[:5], body_dir=opt("--body"), d20=opt("--d20"), l20=opt("--l20"), grid=opt("--grid"),
         smooth=opt("--smooth"), method=opt("--method") or "lin")
