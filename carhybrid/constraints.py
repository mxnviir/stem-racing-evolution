"""The rules as linear constraints on the stage-2 numbers: the search never proposes a body that breaks a rule
(carhybrid/CONSTRAINTS.md). Owner, 2026-11-28: "it shouldn't invent surface that breaks a rule in the first place".

    python carhybrid/constraints.py <builder checkout> <layout: team | genes json | Lxxxx> [--grid gu,gv,gw,gb]
           [--out c.npz] [--check-random N]

Why it is exact and linear
    The body is the solid { bottom(x,|y|) <= z <= top(x,|y|), |y| <= side(x,z) } blended by a smooth maximum
    (freeform.field). The three surfaces are cubic splines of their control grids, and the grids are the seed's
    plus from_unit's cubic interpolation of the numbers u: both linear. So at any point each surface is s0 + J u,
    exactly (the one non-linearity, from_unit's clamp of the side grid at 0, becomes the constraint "side grid >= 0").
    The smooth maximum only ever removes material (smax >= max), so a rule that keeps the SHARP body out of a region
    keeps the built body out of it too.

Which rules (the builder's own masks, built as param_body.build builds them)
    - forced solid (hard_mask_solid + build_rule_masks' solid): T4.2 cargo and its marking room, the axle supports,
      T5.5's 3 mm safety zone round the chamber, the halo -> canister loft, the nose stub, the halo pocket bridge:
      in every column through them top >= their top and bottom <= their bottom, and in every row side >= their |y|;
    - CLIPPING air (forced air that bounds where a surface may go): outside the rule envelope (the block, T3.4/T3.5,
      the region bounds), T7.9's zones, T7.2's wheel keep-clear, T4.4.2's view air round the halo, T8.5.1's nose box,
      the halo region's 27.175 mm (freeform.field; T4.4.2-3). The team's sidepod cap was dropped (owner,
      2026-11-28; evo_patches/sidepod_rules.py): in each column the body stays in its allowed interval
      (top <= the air above it, bottom >= the air below it), or the column is empty (top < bottom);
    - FIXED CUTS (forced air that is always cut, whatever the numbers): the halo pocket and seat, the chamber bore,
      axle holes and wheel discs (fixed_hardware's voids), the block's underside slot: NOT constraints. They hold
      their surface; the gradient must treat it as held (carhybrid/grad.py, the other agent's link.Masks).
    Margins: half a build cell (the masks are cells) + 0.05 mm for air; + 0.8 mm for solid (the smooth maximum can
    take up to blend/4 per blend off the required material).

What it gives
    build()    -> Constraints: A u <= b (rows tagged by rule), with the box |u| <= 1.
    project()  -> the nearest legal u (Euclidean, in the numbers): ADMM, exact to ~1e-6.
    check()    -> the builder's own field on the builder's own grid at u: cells of clipping air the body fills and
                  cells of required solid it leaves empty (both must be 0 on a projected u).
Not covered (stays in the builder, measured): the milling's closing of concavities under the 6.35 mm ball, the mass
skin (0 on every sweep5 car; check() reports it as an assumption), and the parts (wings, supports: their engines).
"""
from __future__ import annotations

import json
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
AIR_MARGIN_MM = 0.05
SOLID_MARGIN_MM = 0.8


def builder(evo):
    sys.path.insert(0, str(HERE))
    from grad import evo_import
    evo_import(evo)


@dataclass
class Constraints:
    A: np.ndarray
    b: np.ndarray
    tags: list                      # per row: (rule family, surface, x, y or z)
    n: int
    grid: tuple
    info: dict = field(default_factory=dict)

    def slack(self, u):
        return self.b - self.A @ np.asarray(u, float)

    def violated(self, u, tol=1e-6):
        s = self.slack(u)
        return np.nonzero(s < -tol)[0]

    def summary(self, u):
        s = self.slack(u)
        bad = s < -1e-6
        out = {}
        for k in np.nonzero(bad)[0]:
            key = f"{self.tags[k][0]}:{self.tags[k][1]}"
            o = out.setdefault(key, {"rows": 0, "worst_mm": 0.0})
            o["rows"] += 1
            o["worst_mm"] = max(o["worst_mm"], float(-s[k]))
        return out


# ------------------------------------------------------------------ masks
def masks_of(d, spacing_mm=1.0, cap=True):
    """(xs, ys, zs, solid, clip, fixed) on the builder's grid (mm), y >= 0 half, as param_body.build makes them.
    cap=False leaves out the halo region's 27.175 mm (freeform.field's own soft compression, not a mask)."""
    import unified_phi as up
    import freeform as fr
    from coarse import use_spacing
    use_spacing(spacing_mm)
    ugk = {"with_cargo": True, "cargo_placement": d.get("cargo")}
    if d.get("rear") is not None:
        from dataclasses import replace as _replace
        from bounding_volumes import default_rule_envelope
        ugk["rule_envelope"] = _replace(default_rule_envelope(), rearpod_max_length_m=(d["rear"] - 16.0) / 1e3)
    g = up.build_unified_geometry(d["W"], d["x_front"], d["d_halo"], init_mode="full", **ugk)
    g.phi.apply_hard_constraints()
    air = np.asarray(g.phi.hard_mask_air, bool).copy()
    solid = np.asarray(g.phi.hard_mask_solid, bool).copy()
    rair, rsolid = up.build_rule_masks(g)
    air |= rair; solid |= rsolid; solid &= ~air
    r, h = g.region, up.GRID_SPACING_M * 1e3
    xs = (r.origin_m[0] * 1e3 + h * np.arange(r.shape[0]))
    ys = (r.origin_m[1] * 1e3 + h * np.arange(r.shape[1]))
    zs = (r.origin_m[2] * 1e3 + h * np.arange(r.shape[2]))
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")
    # fixed cuts: always cut whatever the numbers (not constraints)
    fixed = np.zeros_like(air)
    fh = getattr(g, "fixed_hardware", None) or getattr(g.phi, "fixed_hardware", None)
    if fh is not None and getattr(fh, "combined_void_mask", None) is not None:
        fixed |= np.asarray(fh.combined_void_mask, bool)
    fixed |= up.build_halo_pocket_forbidden_mask(r.origin_m, r.shape, g.landmarks["ref_plane_A_m"], d["d_halo"])
    try:
        from fixed_hardware import halo_seat_mask
        import halo_seat as hs
        fixed |= halo_seat_mask(r.shape, r.origin_m, g.landmarks["ref_plane_A_m"] * 1e3 + d["d_halo"], field=hs.seat_air,
                                views=False)
    except Exception:                                    # noqa: BLE001  (older builders: the seat is in the pocket mask)
        pass
    fixed |= up.model_block_slot_mask(r)
    # the halo region's 27.175 mm (freeform.field; T4.4.2-3), as clipping air. The team's sidepod cap is dropped
    # (owner, 2026-11-28; evo_patches/sidepod_rules.py S1): only rules bound the sidepods
    T = fr.CAP_TAPER_MM
    xa, xh0, xh1, hw = fr.halo_span(g.landmarks["ref_plane_A_m"] * 1e3, d["d_halo"])
    cap_air = (X >= xa - T) & (X <= xh1) & ((np.abs(Y) <= hw + T) | (X >= xh0 - T)) & (Z > fr.SIDEPOD_TOP_MM)
    if cap:
        air |= cap_air & ~solid
    clip = air & ~fixed
    # border cells: the region's edge (clipping)
    j0 = int(np.searchsorted(ys, -1e-9))
    sl = slice(j0, None)
    return xs, ys[sl], zs, solid[:, sl], clip[:, sl], fixed[:, sl], g, h


# ------------------------------------------------------------------ surfaces as linear maps of u
def surface_maps(seed, grid, amp, xs, ys, zs, frame):
    """(top0, J_top), (bot0, J_bot) on the (x, y) columns and (side0, J_side) on the (x, z) rows, and (S0, J_S) the
    side control grid (its >= 0 clamp): each = value at u = 0 and d value / d u (exact, linear)."""
    import freeform as fr
    gu, gv, gw, gb = grid
    n = fr.n_genes(gu, gv, gw, gb)
    x_tip, x_rear = frame
    Xc, Yc = np.meshgrid(xs, ys, indexing="ij")
    Xs, Zs = np.meshgrid(xs, zs, indexing="ij")

    def at(u):
        T, B, S = fr.from_unit(seed, u, amp, gu, gv, gw, gb).grids()
        uu = (Xc - x_tip) / (x_rear - x_tip); vv = np.abs(Yc) / fr.Y_MAX_MM
        side = fr._interp(S, (Xs - x_tip) / (x_rear - x_tip), Zs / fr.Z_MAX_MM)
        return fr._interp(T, uu, vv).ravel(), fr._interp(B, uu, vv).ravel(), side.ravel()

    base = at(np.zeros(n))
    J = [np.zeros((len(v), n)) for v in base]
    h = 1.0
    for k in range(n):
        e = np.zeros(n); e[k] = h
        p = at(e)
        for q in range(3):
            J[q][:, k] = (p[q] - base[q]) / h
    return base, J, n


def side_grid_map(seed, grid, amp):
    """(S0, J): the side control grid before from_unit's clamp, linear in u: its rows give 'side grid >= 0'."""
    import freeform as fr
    from scipy import ndimage
    gu, gv, gw, gb = grid
    n = fr.n_genes(gu, gv, gw, gb)
    shapes = list(fr.coarse_grid_shapes(gu, gv, gw, gb).values())
    T, B, S = seed.grids()
    a, b_ = shapes[2]
    off = shapes[0][0] * shapes[0][1] + shapes[1][0] * shapes[1][1]

    def d_of(c):
        ca = np.linspace(0, a - 1, S.shape[0]); cb = np.linspace(0, b_ - 1, S.shape[1])
        CA, CB = np.meshgrid(ca, cb, indexing="ij")
        return ndimage.map_coordinates(c, [CA.ravel(), CB.ravel()], order=min(3, max(a, b_) - 1) if min(a, b_) > 1 else 1,
                                       mode="nearest")
    J = np.zeros((S.size, n))
    for k in range(a * b_):
        c = np.zeros(a * b_); c[k] = 1.0
        J[:, off + k] = amp[2] * d_of(c.reshape(a, b_))
    return S.ravel(), J


# ------------------------------------------------------------------ the constraints
def build(evo, layout="team", grid=None, spacing_mm=1.0) -> Constraints:
    builder(evo)
    import layout_genes as LG
    import freeform as fr
    import simscale_shape as SS
    sys.path.insert(0, str(HERE))
    from make_car import genes_of, grid_of
    t0 = time.time()
    d = LG.decode(genes_of(Path(evo), layout))
    seed = fr.Freeform.from_dict(d["freeform"])
    g4 = tuple(int(v) for v in grid.split(",")) if isinstance(grid, str) else (grid or tuple(SS.COARSE))
    g4 = g4 if len(g4) == 4 else g4 + (None,)
    amp = SS.AMP_MM
    xs, ys, zs, solid, clip, fixed, geom, h = masks_of(d, spacing_mm)
    t_masks = time.time() - t0
    frame = (seed.x_tip_mm, seed.x_rear_mm)
    inx = (xs >= frame[0]) & (xs <= frame[1])
    (top0, bot0, side0), (Jt, Jb, Js), n = surface_maps(seed, g4, amp, xs, ys, zs, frame)
    nx, ny, nz = len(xs), len(ys), len(zs)
    rows_A, rows_b, tags = [], [], []

    def add(J, rhs, tag):
        rows_A.append(J); rows_b.append(rhs); tags.extend(tag)

    ci_all, cj_all = np.meshgrid(np.arange(nx), np.arange(ny), indexing="ij")
    # columns (x, y): top / bottom
    t_max = np.full((nx, ny), np.nan); t_min = np.full((nx, ny), np.nan)
    b_min = np.full((nx, ny), np.nan); b_max = np.full((nx, ny), np.nan)
    empty = np.zeros((nx, ny), bool)
    for i in np.nonzero(inx)[0]:
        for j in range(ny):
            S_ = solid[i, j]; C_ = clip[i, j]
            if S_.any():
                k = np.nonzero(S_)[0]; lo, hi = k.min(), k.max()
                t_min[i, j] = zs[hi] + 0.5 * h + SOLID_MARGIN_MM
                b_max[i, j] = zs[lo] - 0.5 * h - SOLID_MARGIN_MM
                above = np.nonzero(C_[hi:])[0]
                below = np.nonzero(C_[:lo + 1])[0]
                if len(above):
                    t_max[i, j] = zs[hi + above.min()] - 0.5 * h - AIR_MARGIN_MM
                if len(below):
                    b_min[i, j] = zs[below.max()] + 0.5 * h + AIR_MARGIN_MM
            elif C_.any():
                free = ~C_
                if not free.any():
                    empty[i, j] = True
                    continue
                # the body keeps to the longest allowed run of this column (or the column is left empty)
                runs, k = [], 0
                while k < nz:
                    if free[k]:
                        k0 = k
                        while k < nz and free[k]:
                            k += 1
                        runs.append((k0, k - 1))
                    else:
                        k += 1
                k0, k1 = max(runs, key=lambda r: r[1] - r[0])
                if k1 < nz - 1:
                    t_max[i, j] = zs[k1 + 1] - 0.5 * h - AIR_MARGIN_MM
                if k0 > 0:
                    b_min[i, j] = zs[k0 - 1] + 0.5 * h + AIR_MARGIN_MM
    # where required solid meets clipping air in the same column the surface must sit ON the boundary: the margins
    # would otherwise contradict each other (top >= solid + 0.8 and top <= air - 0.55 a cell apart)
    for lo_, hi_ in ((t_min, t_max), (b_max, b_min)):
        clash = ~np.isnan(lo_) & ~np.isnan(hi_) & ((lo_ > hi_) if lo_ is t_min else (lo_ < hi_))
        mid = 0.5 * (lo_ + hi_)
        if lo_ is t_min:
            t_min[clash] = mid[clash] - 0.05; t_max[clash] = mid[clash] + 0.05
        else:
            b_max[clash] = mid[clash] + 0.05; b_min[clash] = mid[clash] - 0.05
    idx = lambda i, j: i * ny + j
    for arr, sgn, J0, s0, fam in ((t_max, 1, Jt, top0, "air_above"), (t_min, -1, Jt, top0, "solid_top"),
                                  (b_min, -1, Jb, bot0, "air_below"), (b_max, 1, Jb, bot0, "solid_bottom")):
        I, Jj = np.nonzero(~np.isnan(arr))
        if len(I):
            r = idx(I, Jj)
            # sgn * (s0 + J u) <= sgn * bound
            add(sgn * J0[r], sgn * (arr[I, Jj] - s0[r]), [(fam, "top" if J0 is Jt else "bottom", float(xs[a]), float(ys[c]))
                                                           for a, c in zip(I, Jj)])
    I, Jj = np.nonzero(empty)
    if len(I):
        r = idx(I, Jj)
        add(Jt[r] - Jb[r], -(top0[r] - bot0[r]) - 0.1, [("air_column", "top-bottom", float(xs[a]), float(ys[c])) for a, c in zip(I, Jj)])
    # rows (x, z): side >= the required solid's |y|
    s_min = np.full((nx, nz), np.nan)
    for i in np.nonzero(inx)[0]:
        for k in range(nz):
            jj = np.nonzero(solid[i, :, k])[0]
            if len(jj):
                s_min[i, k] = ys[jj.max()] + 0.5 * h + SOLID_MARGIN_MM
    I, K = np.nonzero(~np.isnan(s_min))
    if len(I):
        r = I * nz + K
        add(-Js[r], -(s_min[I, K] - side0[r]), [("solid_side", "side", float(xs[a]), float(zs[c])) for a, c in zip(I, K)])
    # the side control grid stays >= 0 (from_unit's clamp: linear only there)
    S0, JS = side_grid_map(seed, g4, amp)
    add(-JS, S0.copy(), [("side_grid_nonneg", "side grid", float(q), 0.0) for q in range(len(S0))])
    # the box
    add(np.eye(n), np.ones(n), [("box", "u<=1", float(q), 0.0) for q in range(n)])
    add(-np.eye(n), np.ones(n), [("box", "u>=-1", float(q), 0.0) for q in range(n)])
    A = np.vstack(rows_A); b = np.concatenate(rows_b)
    keep = np.abs(A).sum(1) > 1e-12                      # rows u cannot move: kept only if violated at every u
    const_bad = (~keep) & (b < -1e-9)
    C = Constraints(A[keep], b[keep], [t for t, k in zip(tags, keep) if k], n, g4,
                    {"layout": layout, "grid": list(map(lambda v: v if v is None else int(v), g4)), "n_rows": int(keep.sum()),
                     "rows_u_cannot_fix": int(const_bad.sum()),
                     "rows_u_cannot_fix_tags": sorted({f"{t[0]}:{t[1]}" for t, k in zip(tags, const_bad) if k}),
                     "masks_s": round(t_masks, 1), "build_s": round(time.time() - t0, 1),
                     "cells": {"solid": int(solid.sum()), "clip": int(clip.sum()), "fixed": int(fixed.sum())}})
    C._ctx = (d, seed, g4, amp)
    return C


# ------------------------------------------------------------------ projection
def project(C: Constraints, v, rho=1.0, iters=4000, tol=1e-7):
    """argmin ||u - v|| s.t. A u <= b: HiGHS' QP solver (active set, exact; highspy) when installed, else OSQP, else
    ADMM (which does not converge on the 720-number set: CONSTRAINTS.md)."""
    try:
        return _project_highs(C, v)
    except ImportError:
        pass
    try:
        import osqp
        import scipy.sparse as sp
        v = np.asarray(v, float)
        m = osqp.OSQP()
        m.setup(P=sp.identity(C.n, format="csc"), q=-v, A=sp.csc_matrix(C.A), l=np.full(len(C.b), -np.inf), u=C.b,
                eps_abs=1e-7, eps_rel=1e-7, max_iter=200000, polishing=True, verbose=False)
        r = m.solve()
        if r.x is not None and r.info.status in ("solved", "solved inaccurate"):
            u = np.asarray(r.x, float)
            return u, {"solver": "osqp", "status": r.info.status, "iterations": int(r.info.iter),
                       "max_violation_mm": float(max(0.0, -C.slack(u).min())), "moved": float(np.linalg.norm(u - v))}
    except ImportError:
        pass
    return _project_admm(C, v, rho, iters, tol)


_HIGHS_OPTS = {}


def _project_highs(C: Constraints, v):
    import highspy
    import scipy.sparse as sp
    v = np.asarray(v, float)
    n, m = C.n, len(C.b)
    inf = highspy.kHighsInf
    A = sp.csc_matrix(C.A)
    model = highspy.HighsModel()
    lp = model.lp_
    lp.num_col_, lp.num_row_ = n, m
    lp.col_cost_ = -v
    lp.col_lower_, lp.col_upper_ = np.full(n, -inf), np.full(n, inf)
    lp.row_lower_, lp.row_upper_ = np.full(m, -inf), np.asarray(C.b, float)
    lp.a_matrix_.format_ = highspy.MatrixFormat.kColwise
    lp.a_matrix_.start_, lp.a_matrix_.index_, lp.a_matrix_.value_ = A.indptr, A.indices, A.data
    H = model.hessian_                                    # 1/2 u'u (lower triangle of I)
    H.dim_, H.format_ = n, highspy.HessianFormat.kTriangular
    H.start_, H.index_, H.value_ = np.arange(n + 1), np.arange(n), np.ones(n)
    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    for k_, v_ in _HIGHS_OPTS.items():
        h.setOptionValue(k_, v_)
    h.setOptionValue("time_limit", 1800.0)
    h.passModel(model)
    h.run()
    status = h.modelStatusToString(h.getModelStatus())
    u = np.asarray(h.getSolution().col_value, float)
    if len(u) != n or status != "Optimal":
        raise RuntimeError(f"HiGHS QP: {status}")
    return u, {"solver": "highs-qp", "status": status, "iterations": int(h.getInfo().qp_iteration_count),
               "max_violation_mm": float(max(0.0, -C.slack(u).min())), "moved": float(np.linalg.norm(u - v))}


def _project_admm(C: Constraints, v, rho=1.0, iters=4000, tol=1e-7):
    """argmin ||u - v|| s.t. A u <= b (ADMM on z = A u; the u-step's matrix factored once)."""
    from scipy.linalg import cho_factor, cho_solve
    A, b = C.A, C.b
    v = np.asarray(v, float)
    K = cho_factor(np.eye(C.n) + rho * A.T @ A)
    z = np.minimum(A @ v, b); w = np.zeros(len(b)); u = v.copy()
    for it in range(iters):
        u = cho_solve(K, v + rho * A.T @ (z - w))
        Au = A @ u
        z_new = np.minimum(Au + w, b)
        w += Au - z_new
        r = np.abs(Au - z_new).max(); s = rho * np.abs(A.T @ (z_new - z)).max()
        z = z_new
        if r < tol and s < tol:
            break
    # the last bit of infeasibility (ADMM's primal residual) taken out by a few exact active-set steps
    for _ in range(20):
        s = C.slack(u)
        bad = s < 0
        if not bad.any():
            break
        Ab = C.A[bad]
        try:
            u = u + Ab.T @ np.linalg.lstsq(Ab @ Ab.T, s[bad] - 1e-9, rcond=1e-10)[0]
        except np.linalg.LinAlgError:
            break
    return u, {"iterations": it + 1, "max_violation_mm": float(max(0.0, -C.slack(u).min())), "moved": float(np.linalg.norm(u - v))}


# ------------------------------------------------------------------ the builder's own check
def feasibility(C: Constraints):
    """The smallest worst violation any u can reach (mm; 0 = feasible): min t s.t. A u - t <= b, t >= 0 (HiGHS LP),
    and the rows that bind at that optimum, by rule family. Tells which rules the numbers cannot meet."""
    from scipy.optimize import linprog
    import scipy.sparse as sp
    n, m = C.n, len(C.b)
    is_box = np.array([t[0] == "box" for t in C.tags])
    A = sp.hstack([sp.csr_matrix(C.A), sp.csr_matrix(np.where(is_box, 0.0, -1.0)[:, None])])
    c = np.zeros(n + 1); c[-1] = 1.0
    r = linprog(c, A_ub=A, b_ub=C.b, bounds=[(None, None)] * n + [(0, None)], method="highs")
    if r.status != 0:
        return {"status": r.message}
    u, t = r.x[:n], r.x[-1]
    s = C.slack(u)
    bind = {}
    for k in np.nonzero((s < 1e-6 - t) | (s < -1e-6))[0]:
        key = f"{C.tags[k][0]}:{C.tags[k][1]}"
        bind.setdefault(key, []).append((round(C.tags[k][2], 1), round(C.tags[k][3], 1)))
    return {"min_worst_violation_mm": float(t), "u": u,
            "binding": {k: {"rows": len(v), "x_range": [min(q[0] for q in v), max(q[0] for q in v)]} for k, v in bind.items()}}


def split_sharp(C: Constraints, tol_mm=0.05):
    """The rules the numbers can meet and the ones they cannot (sharper than the control grid can draw): an L1 LP,
    min sum s_i s.t. A u - s <= b, s >= 0 (box rows hard). Rows left with s_i > tol are SHARP: the builder's mask stays
    in charge there and the gradient holds that surface. Returns (Constraints without them, report)."""
    from scipy.optimize import linprog
    import scipy.sparse as sp
    n, m = C.n, len(C.b)
    soft = np.array([t[0] != "box" for t in C.tags])
    ms = int(soft.sum())
    S_ = sp.lil_matrix((m, ms)); S_[np.nonzero(soft)[0], np.arange(ms)] = -1.0
    A = sp.hstack([sp.csr_matrix(C.A), S_.tocsr()])
    c = np.r_[np.zeros(n), np.ones(ms)]
    r = linprog(c, A_ub=A, b_ub=C.b, bounds=[(None, None)] * n + [(0, None)] * ms, method="highs")
    if r.status != 0:
        raise RuntimeError(r.message)
    s_all = np.zeros(m); s_all[soft] = r.x[n:]
    sharp = s_all > tol_mm
    rep = {}
    for k in np.nonzero(sharp)[0]:
        key = f"{C.tags[k][0]}:{C.tags[k][1]}"
        o = rep.setdefault(key, {"rows": 0, "worst_mm": 0.0, "x": [1e9, -1e9]})
        o["rows"] += 1; o["worst_mm"] = max(o["worst_mm"], float(s_all[k]))
        o["x"] = [min(o["x"][0], C.tags[k][2]), max(o["x"][1], C.tags[k][2])]
    keep = ~sharp
    # rows kept with a slack of at most tol_mm are met to within it: their bound is relaxed by that slack, so the
    # reachable set is feasible exactly (the L1 program's u meets it) and a QP can project onto it
    D = Constraints(C.A[keep], C.b[keep] + s_all[keep] + 1e-9, [t for t, k in zip(C.tags, keep) if k], C.n, C.grid,
                    dict(C.info, sharp_rows=int(sharp.sum()), sharp=rep, constrained_rows=int(keep.sum())))
    D._ctx = C._ctx
    D._sharp_tags = [t for t, k in zip(C.tags, sharp) if k]
    D._u_l1 = r.x[:n]
    return D, rep


def check(C: Constraints, evo, u, spacing_mm=1.0):
    """The builder's field (freeform.field as param_body calls it) on the builder's grid at u: cells of clipping air
    it fills, cells of required solid it leaves empty, before any mask acts. Both 0 = the masks have nothing to do."""
    import freeform as fr
    d, seed, g4, amp = C._ctx
    ff = fr.from_unit(seed, np.asarray(u, float), amp, *g4)
    xs, ys, zs, solid, clip, fixed, geom, h = masks_of(d, spacing_mm)
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")
    xa = geom.landmarks["ref_plane_A_m"] * 1e3
    phi = fr.field(ff, X, Y, Z, ff.x_tip_mm, ff.x_rear_mm, axles=None,     # as param_body (sidepod_rules S1)
                   halo=fr.halo_span(xa, d["d_halo"]))
    body = phi < 0
    return {"clip_air_filled_cells": int((body & clip).sum()), "required_solid_empty_cells": int((~body & solid).sum()),
            "fixed_cut_cells_filled": int((body & fixed).sum()), "body_cells": int(body.sum()),
            "note": "fixed cuts are always cut (not constraints); the mass skin is assumed 0 (sweep5: 0 on every car)"}


def main(argv):
    if len(argv) < 2 or argv[0] in ("-h", "--help"):
        print(__doc__); return 0
    opt = lambda k, dflt=None: argv[argv.index(k) + 1] if k in argv else dflt
    evo, layout = argv[0], argv[1]
    C = build(evo, layout, opt("--grid"))
    print(json.dumps(C.info, indent=1))
    u0 = np.zeros(C.n)
    print("seed (u = 0) violations:", json.dumps(C.summary(u0), indent=1))
    print("builder check at u = 0:", check(C, evo, u0), flush=True)
    t0 = time.time()
    C, rep = split_sharp(C)                              # the reachable rules (hard); the sharp ones stay with the masks
    print(f"split ({time.time() - t0:.0f} s): {C.info['constrained_rows']} rows reachable, {C.info['sharp_rows']} sharp:",
          json.dumps(rep, indent=1), flush=True)
    t0 = time.time()
    p0, inf0 = project(C, u0)
    inf0["s"] = round(time.time() - t0, 1)
    print("projected seed:", inf0, "violations", C.summary(p0))
    print("builder check at the projected seed:", check(C, evo, p0))
    rng = np.random.default_rng(1)
    for k in range(int(opt("--check-random", 2))):
        v = rng.uniform(-1, 1, C.n)
        p, inf = project(C, v)
        print(f"random {k}: violations before {len(C.violated(v))} rows; projected {inf}; builder check", check(C, evo, p))
    if opt("--out"):
        np.savez(opt("--out"), A=C.A, b=C.b, info=json.dumps(C.info), u_seed_projected=p0)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
