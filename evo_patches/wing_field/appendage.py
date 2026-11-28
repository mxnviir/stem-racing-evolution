"""appendage.py -- a wing's supports as a FIELD the search designs, and the plastic-foam joint made exact.

Written into the builder checkout as car/parts/appendage.py by evo_patches/wing_field.py. Design note:
evo_patches/wing_field/WING_FIELD.md (why, the rules, what is measured).

Supports as a field
    The supports of a wing are the solid where P(p) < 0, P = seed(p) - AMP_MM * sum_j u_j g_j(p): seed is the
    signed distance (mm, + outside) of a starting support (the loft/strut engines' own, or none), g_j are Gaussian
    bumps on a regular grid of control points over the wing's support ZONE, u_j in [-1, 1] the numbers the search
    moves. u_j > 0 grows material round its point, < 0 takes it away: struts, pods, fairings, bridges, two legs merging
    into one, a support with a hole -- any topology the numbers reach. Nothing here judges drag: what the numbers
    should be is decided by the evolution (moves, recombination) and the adjoint (DAFoam's dD/dX on the support's
    surface, mapped to dD/du_j exactly as carhybrid/grad.py maps the body), on the CFD score. Material outside the
    zone is cut; pieces not joined to both the wing and the body (or nose) are dropped; walls thinner than the print
    wall are opened away; strength is a hard limit (fe_check), never an objective.

The plastic-foam joint (joint())
    Wherever printed plastic meets milled foam:
      1. the foam the plastic takes (its pocket) is what a 6.35 mm ball end mill can cut, 3-axis from the side the
         support comes from: the morphological opening, by a ball of radius R_MM, of the pocket the plastic asks
         for plus the open air it is reached through, with the ball's centre on a vertical line open to the air
         (no undercut). Every concave corner the foam keeps is therefore at least R_MM round;
      2. the foam that must stay foam is never taken: the 3 mm wall round the cartridge chamber (T5.5) and the bore
         (T5.5 "not intersected by any object"); above that wall the zone may be air or plastic;
      3. where the plastic's outside meets the foam's outside the corner is filled with plastic as a rolling-ball
         fillet of radius R_MM (the closing of their union by the ball, kept only round the joint line);
      4. foam = body - pocket and plastic = (support + fillet - body) + pocket, from ONE pair of Boolean operations on
         the same meshes: every face the two share is the same face, so the printed part is the exact positive of
         the milled cavity (to floating point), and the cavity's concave radii are the ball's.
    R_MM is 3.175 mm (the 6.35 mm ball) plus BUILD_MARGIN_MM, so a fillet measured on the voxel grid is never under
    3.175.

Grids: millimetres, x from the nose tip, y >= 0 (right half; the car is mirrored), z up from the track.
"""
from __future__ import annotations

import numpy as np

R_BALL_MM = 3.175            # the 6.35 mm ball end mill
BUILD_MARGIN_MM = 0.15       # built this much rounder, so the measured radius is never under R_BALL_MM
R_MM = R_BALL_MM + BUILD_MARGIN_MM
H_MM = 0.25                  # joint grid
CANISTER_WALL_MM = 3.0       # T5.5: foam this thick round the chamber, always
PRINT_WALL_MM = 1.2          # thinnest printed wall kept (the opening by a ball of half this)
AMP_MM = 4.0                 # mm of field per unit u
BUMP_SIGMA_MM = 2.5          # a control point's reach
CTRL_SPACING_MM = 4.0        # control-point spacing over a zone
E_ABS_MPA, NU_ABS, ALLOW_MPA = 2000.0, 0.35, 20.0     # ABS, as wing_strut / support_engine


# ------------------------------------------------------------------ grids and fields
class Grid:
    """A regular grid of cell centres over [lo, hi] (mm), spacing h."""

    def __init__(self, lo, hi, h=H_MM):
        self.lo, self.h = np.asarray(lo, float), float(h)
        n = np.ceil((np.asarray(hi, float) - self.lo) / h).astype(int) + 1
        self.shape = tuple(int(v) for v in n)
        self.xs, self.ys, self.zs = (self.lo[i] + h * np.arange(n[i]) for i in range(3))

    def points(self, mask=None):
        X, Y, Z = np.meshgrid(self.xs, self.ys, self.zs, indexing="ij")
        P = np.stack([X, Y, Z], -1)
        return P[mask] if mask is not None else P


def occupancy(mesh_mm, g: Grid):
    """Cells of g inside a closed mesh (mm): vertical rays through every (x, y) column, inside between pairs of
    crossings (parity), so a 3-million-cell grid takes seconds."""
    X, Y = np.meshgrid(g.xs, g.ys, indexing="ij")
    O = np.c_[X.ravel(), Y.ravel(), np.full(X.size, g.zs[0] - 5.0)]
    D = np.tile([0.0, 0.0, 1.0], (len(O), 1))
    occ = np.zeros(g.shape, bool)
    loc, ri, _ = mesh_mm.ray.intersects_location(O, D, multiple_hits=True)
    if len(ri) == 0:                                     # the mesh is nowhere over this grid
        return occ
    loc = np.asarray(loc).reshape(-1, 3)
    order = np.lexsort((loc[:, 2], ri))
    loc, ri = loc[order], ri[order]
    starts = np.searchsorted(ri, np.arange(len(O)))
    ends = np.searchsorted(ri, np.arange(len(O)), side="right")
    nx, ny = g.shape[0], g.shape[1]
    for k in np.nonzero(ends > starts)[0]:
        zs = np.unique(np.round(loc[starts[k]:ends[k], 2], 6))
        i, j = divmod(k, ny)
        for a, b in zip(zs[0::2], zs[1::2]):
            occ[i, j] |= (g.zs >= a) & (g.zs <= b)
    return occ


def sdf(occ, h):
    """Signed distance (mm, + outside) of a voxel set, to half a cell."""
    from scipy.ndimage import distance_transform_edt as edt
    if not occ.any():
        return np.full(occ.shape, 1e3)
    return (edt(~occ) - edt(occ)) * h + np.where(occ, 0.5 * h, -0.5 * h)


def erode(occ, r, h):
    from scipy.ndimage import distance_transform_edt as edt
    return edt(occ) * h >= r


def dilate(occ, r, h):
    from scipy.ndimage import distance_transform_edt as edt
    return edt(~occ) * h <= r


def opening(occ, r, h):
    return dilate(erode(occ, r, h), r, h)


def closing(occ, r, h):
    return erode(dilate(occ, r, h), r, h)


def components_touching(occ, *touch):
    """The connected pieces of occ that meet every one of the masks in touch (6-connected)."""
    from scipy.ndimage import label
    lab, n = label(occ)
    keep = np.zeros(n + 1, bool)
    for k in range(1, n + 1):
        piece = lab == k
        keep[k] = all((piece & t).any() for t in touch)
    return keep[lab]


# ------------------------------------------------------------------ the support field
def control_points(zone, g: Grid, spacing=CTRL_SPACING_MM):
    """Control points on a regular lattice over the zone's cells (mm), in a fixed order (x, then y, then z)."""
    P = g.points(zone)
    if not len(P):
        return np.zeros((0, 3))
    lo, hi = P.min(0), P.max(0)
    axes = [np.arange(lo[i] + 0.5 * spacing, hi[i], spacing) if hi[i] - lo[i] > spacing else np.array([0.5 * (lo[i] + hi[i])])
            for i in range(3)]
    C = np.stack(np.meshgrid(*axes, indexing="ij"), -1).reshape(-1, 3)
    # only the points inside (or within a cell of) the zone
    idx = np.clip(np.round((C - g.lo) / g.h).astype(int), 0, np.array(g.shape) - 1)
    return C[zone[idx[:, 0], idx[:, 1], idx[:, 2]]]


def field(seed_sdf, u, ctrl, g: Grid, amp=AMP_MM, sigma=BUMP_SIGMA_MM):
    """P = seed - amp * sum_j u_j exp(-|p - c_j|^2 / 2 sigma^2) on g (mm)."""
    P = np.array(seed_sdf, float, copy=True)
    u = np.asarray(u, float)
    r = int(np.ceil(3 * sigma / g.h))
    for c, uj in zip(ctrl, u):
        if abs(uj) < 1e-9:
            continue
        i0 = np.round((c - g.lo) / g.h).astype(int)
        sl = tuple(slice(max(i0[k] - r, 0), min(i0[k] + r + 1, g.shape[k])) for k in range(3))
        xs, ys, zs = g.xs[sl[0]], g.ys[sl[1]], g.zs[sl[2]]
        d2 = ((xs - c[0])[:, None, None] ** 2 + (ys - c[1])[None, :, None] ** 2 + (zs - c[2])[None, None, :] ** 2)
        P[sl] -= amp * uj * np.exp(-d2 / (2 * sigma * sigma))
    return P


def support_solid(P, zone, wing, base, g: Grid, wall=PRINT_WALL_MM):
    """The support the field draws, made buildable: inside the zone, walls of at least `wall` (the opening by a ball
    of half of it), only the pieces that join the wing to the base (the body or the nose, both as masks)."""
    s = (P < 0) & zone
    s = opening(s, 0.5 * wall, g.h) & zone
    near_w = dilate(wing, g.h * 1.5, g.h)
    near_b = dilate(base, g.h * 1.5, g.h)
    return components_touching(s | (wing & near_w & s), near_w, near_b)


# ------------------------------------------------------------------ the joint
def chamber_keep(g: Grid, rear_face_mm, depth_mm, bore_r_mm, axis_z_mm, wall=CANISTER_WALL_MM):
    """Cells that must stay foam or air as they are: the bore and the 3 mm wall round it (T5.5), over the bore's
    depth from the rear face (and the wall's own thickness ahead of its end)."""
    X, Y, Z = np.meshgrid(g.xs, g.ys, g.zs, indexing="ij")
    rr = np.hypot(Y, Z - axis_z_mm)
    return (rr <= bore_r_mm + wall) & (X >= rear_face_mm - depth_mm - wall) & (X <= rear_face_mm + 1e-6)


def millable_pocket(want, body, keep, g: Grid, r=R_MM, from_top=True):
    """The pocket a ball of radius r can cut, 3-axis from above (or below): the union of the balls whose centres
    lie where a whole ball fits in (want or air) without touching kept foam, on a vertical line open to the air.
    Returns the pocket (cells of the body), a subset of want & body & ~keep."""
    free = (want | ~body) & ~keep
    ctr = erode(free, r, g.h)
    # the shank: a centre is reachable only if every centre position above it (below, from underneath) is too
    if from_top:
        ok = np.flip(np.logical_and.accumulate(np.flip(ctr, axis=2), axis=2), axis=2)
    else:
        ok = np.logical_and.accumulate(ctr, axis=2)
    return dilate(ok, r, g.h) & body & want & ~keep


def joint(body, support, keep, g: Grid, r=R_MM, from_top=True):
    """The plastic-foam joint (module docstring, 1-4), on voxels: (foam, plastic, pocket, fillet) masks. The meshes
    that are made (mesh_joint) take these for their Booleans, so foam and plastic are exact complements there."""
    want = support & body
    pocket = millable_pocket(want, body, keep, g, r, from_top)
    outside = support & ~body & ~keep
    plastic0 = outside | pocket
    foam = body & ~pocket
    # the rolling-ball fillet where the plastic's outside meets the foam's outside, filled with plastic
    both = foam | plastic0
    closed = closing(both, r, g.h)
    near = dilate(plastic0, r + g.h, g.h) & dilate(foam, r + g.h, g.h)
    fillet = closed & ~both & near & ~keep
    return foam, plastic0 | fillet, pocket, fillet


def concave_radius_check(solid, g: Grid, r=R_BALL_MM, region=None):
    """mm3 of air the ball of radius r cannot reach in a concave corner of the solid (its complement minus the
    complement's opening by the ball), inside region: 0 means every concave corner there is at least r round."""
    air = ~solid
    left = air & ~opening(air, r - g.h, g.h)
    if region is not None:
        left &= region
    return float(left.sum()) * g.h ** 3


# ------------------------------------------------------------------ meshes
def mesh_of(occ, g: Grid, smooth=0.6):
    """A closed mesh (mm) of a voxel set: marching cubes on its smoothed signed distance (sub-cell surface)."""
    import trimesh
    from scipy.ndimage import gaussian_filter
    from skimage.measure import marching_cubes
    if not occ.any():
        return None
    d = gaussian_filter(sdf(np.pad(occ, 2), g.h), smooth)
    v, f, _, _ = marching_cubes(d, 0.0, spacing=(g.h, g.h, g.h))
    v += g.lo - 2 * g.h
    m = trimesh.Trimesh(v, f[:, ::-1], process=True)
    trimesh.repair.fix_normals(m)
    return m


def mesh_joint(body_mesh_mm, plastic, pocket, g: Grid):
    """The two parts as made, from one pair of Booleans on the same meshes: foam = body - pocket, plastic = (plastic
    voxels' mesh - body) + pocket. The pocket's faces in the foam and in the plastic are the same faces."""
    import trimesh
    pk = mesh_of(pocket, g) if pocket.any() else None
    pl = mesh_of(plastic & ~pocket, g)
    foam = body_mesh_mm if pk is None else trimesh.boolean.difference([body_mesh_mm, pk], engine="manifold")
    parts = [p for p in (trimesh.boolean.difference([pl, body_mesh_mm], engine="manifold") if pl is not None else None, pk)
             if p is not None and not p.is_empty]
    plastic_mesh = trimesh.boolean.union(parts, engine="manifold") if len(parts) > 1 else (parts[0] if parts else None)
    return foam, plastic_mesh


# ------------------------------------------------------------------ strength (hard limit, never an objective)
def _h8_stiffness(E, nu, a):
    """8-node brick (cube of side a) stiffness, 24 x 24, by 2 x 2 x 2 Gauss quadrature."""
    C = E / ((1 + nu) * (1 - 2 * nu)) * np.array(
        [[1 - nu, nu, nu, 0, 0, 0], [nu, 1 - nu, nu, 0, 0, 0], [nu, nu, 1 - nu, 0, 0, 0],
         [0, 0, 0, (1 - 2 * nu) / 2, 0, 0], [0, 0, 0, 0, (1 - 2 * nu) / 2, 0], [0, 0, 0, 0, 0, (1 - 2 * nu) / 2]])
    nodes = np.array([[x, y, z] for z in (-1, 1) for y in (-1, 1) for x in (-1, 1)], float)
    K = np.zeros((24, 24)); Bs = []
    g = 1 / np.sqrt(3)
    for xi in (-g, g):
        for eta in (-g, g):
            for ze in (-g, g):
                dN = 0.125 * np.array([nodes[:, 0] * (1 + eta * nodes[:, 1]) * (1 + ze * nodes[:, 2]),
                                       nodes[:, 1] * (1 + xi * nodes[:, 0]) * (1 + ze * nodes[:, 2]),
                                       nodes[:, 2] * (1 + xi * nodes[:, 0]) * (1 + eta * nodes[:, 1])]) * (2 / a)
                B = np.zeros((6, 24))
                B[0, 0::3], B[1, 1::3], B[2, 2::3] = dN[0], dN[1], dN[2]
                B[3, 0::3], B[3, 1::3] = dN[1], dN[0]
                B[4, 1::3], B[4, 2::3] = dN[2], dN[1]
                B[5, 0::3], B[5, 2::3] = dN[2], dN[0]
                K += B.T @ C @ B * (a / 2) ** 3
                Bs.append(B)
    return K, C, np.mean(Bs, axis=0)


def fe_check(solid, fixed, loaded, load_n, h, E=E_ABS_MPA, nu=NU_ABS):
    """Linear elastic voxel FE of a printed solid (mm, N, MPa): cells `fixed` held (glued to the foam or the nose),
    load_n (N, x y z) spread over the cells `loaded` (the wing). Returns (max von Mises MPa, its cell, max |u| mm)."""
    import scipy.sparse as sp
    from scipy.sparse.linalg import spsolve
    idx = np.argwhere(solid)
    if not len(idx):
        return np.inf, None, np.inf
    Ke, C, Bc = _h8_stiffness(E, nu, h)
    off = np.array([[x, y, z] for z in (0, 1) for y in (0, 1) for x in (0, 1)])
    nshape = np.array(solid.shape) + 1
    nid = np.ravel_multi_index((idx[:, None, :] + off[None]).reshape(-1, 3).T, nshape).reshape(-1, 8)
    used, inv = np.unique(nid, return_inverse=True)
    inv = inv.reshape(-1, 8)
    dof = (inv[:, :, None] * 3 + np.arange(3)).reshape(-1, 24)
    rows = np.repeat(dof, 24, axis=1).ravel(); cols = np.tile(dof, (1, 24)).ravel()
    K = sp.csr_matrix((np.tile(Ke.ravel(), len(idx)), (rows, cols)), shape=(3 * len(used),) * 2)
    f = np.zeros(3 * len(used))
    li = np.nonzero(loaded[tuple(idx.T)])[0]
    if not len(li):
        return np.inf, None, np.inf
    ln = np.unique(inv[li])
    for k in range(3):
        f[ln * 3 + k] += load_n[k] / len(ln)
    fi = np.nonzero(fixed[tuple(idx.T)])[0]
    if not len(fi):
        return np.inf, None, np.inf
    fixed_dof = np.unique((inv[fi][:, :, None] * 3 + np.arange(3)).ravel())
    free = np.setdiff1d(np.arange(3 * len(used)), fixed_dof)
    u = np.zeros(3 * len(used))
    Kf = K[free][:, free].tocsr()
    if len(free) < 60000:
        u[free] = spsolve(Kf.tocsc(), f[free])
    else:                                                # memory: conjugate gradients, Jacobi preconditioned
        from scipy.sparse.linalg import cg, LinearOperator
        dinv = 1.0 / Kf.diagonal()
        u[free], _ = cg(Kf, f[free], rtol=1e-8, maxiter=20000, M=LinearOperator(Kf.shape, lambda v: dinv * v))
    s = (C @ (Bc @ u[dof].T)).T                         # stress at each cell's centre, MPa
    vm = np.sqrt(0.5 * ((s[:, 0] - s[:, 1]) ** 2 + (s[:, 1] - s[:, 2]) ** 2 + (s[:, 2] - s[:, 0]) ** 2)
                 + 3 * (s[:, 3] ** 2 + s[:, 4] ** 2 + s[:, 5] ** 2))
    k = int(np.argmax(vm))
    return float(vm[k]), idx[k].tolist(), float(np.abs(u).max())
