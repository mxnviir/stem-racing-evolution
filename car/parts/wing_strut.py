"""wing_strut.py -- a wing's strut as an organic swept solid that the search designs (user 2026-11-23: "actual
struts which are designed by the evolution algorithm itself. Things that smoothly loft into the sides of the
aerofoil ... fully organically designed and lofts perfectly. It should only be throttled by the minimum
thickness and general structural stability of the ABS").

A strut LEAVES the wing through its side face -- its first section IS the wing's aerofoil there, in the
face's own plane, so the wing runs on into it with no seam -- and follows a smooth path (a cubic Bezier) to
the body or the nose. Its sections morph from the aerofoil into its own section (a superellipse, long along
the flow) and, at the end, flare into its foot. What the search sets (all continuous):

  aim_out_deg, aim_fwd_deg   where it heads: from straight down (rear) / straight in (front), leaning outward
                             or inward, forward or aft -- the first surface that way is its foot
  leave, arrive              how far it carries on out of the face, and into the foot, before it bends
                             (shares of its length): straight and stiff, or long sweeping curves
  mid_x, mid_y, mid_z        a free point its path is drawn towards half way (shares of its length): S-bends,
                             arcs up and over, loops fore or aft, swan necks -- any path that leaves the face
                             square to it and arrives square to its foot
  swell                      its section along its length: fatter (> 1) or slimmer (< 1) half way
  chord_mm, thick_mm         its own section between the wing and the foot (streamlined: long along x)
  flare                      its foot, times that section
  round                      its section's squareness (superellipse exponent)

Throttled only by: the print wall (no section thinner than cc.BUILD_WALL_MM) and the ABS's strength -- every
section along it is checked for bending and axial stress under LOAD_N at the wing (ALLOW_MPA, as the wheel
supports), and a strut that fails is thickened until it holds -- and by the rules, which the callers measure.
"""
from __future__ import annotations

import math

import numpy as np

N_RING = 192            # points round each section
STRAIGHT_MM = 1.0       # straight out of the wing's face before it curves
LEAVE_AROUND = 0.5      # go_around's repair: how far out of the face (share of its length) before it turns
FACE_GROW_MM = 0.1      # the first section 0.1 mm round the aerofoil: the strut holds the whole face
N_PATH = 40             # sections along the path
ALLOW_MPA = 20.0        # ABS, as support_engine (printed, with margin)
LOAD_N = (1.0, 0.5, 2.0)    # per strut at the wing, N: drag (x), side (y), down / bump (z): assumed, stated


def _resample(P2, n=N_RING):
    """A closed outline (k, 2) as n points evenly along its perimeter, anticlockwise, starting at its point
    furthest along +u (an aerofoil's trailing edge): any outline, however thin, keeps its shape."""
    P = np.asarray(P2, float)
    x, y = P[:, 0], P[:, 1]
    if 0.5 * (x * np.roll(y, -1) - np.roll(x, -1) * y).sum() < 0:
        P = P[::-1]
    k = int(np.argmax(P[:, 0]))
    P = np.roll(P, -k, axis=0)
    Q = np.vstack([P, P[:1]])
    seg = np.linalg.norm(np.diff(Q, axis=0), axis=1)
    t = np.r_[0.0, np.cumsum(seg)]
    want = np.linspace(0.0, t[-1], n, endpoint=False)
    return np.c_[np.interp(want, t, Q[:, 0]), np.interp(want, t, Q[:, 1])]


def _superellipse(a, b, p, n=N_RING):
    """Its outline evenly along the perimeter from +u, anticlockwise (as _resample)."""
    th = np.linspace(-np.pi, np.pi, 4 * n, endpoint=False)
    c, s = np.cos(th), np.sin(th)
    P = np.c_[a * np.sign(c) * np.abs(c) ** (2 / p), b * np.sign(s) * np.abs(s) ** (2 / p)]
    return _resample(P, n)


def _bezier(P, s):
    """A Bezier curve of any degree through its control points P (k, 3) at parameters s (de Casteljau)."""
    s = np.asarray(s, float)[:, None, None]
    Q = np.repeat(np.asarray(P, float)[None], len(s), axis=0)
    while Q.shape[1] > 1:
        Q = (1 - s) * Q[:, :-1] + s * Q[:, 1:]
    return Q[:, 0]


def _rmf(points, t0, u0):
    """Rotation-minimising frames (double reflection) along a polyline: tangents T, section axes U, V."""
    T = np.gradient(points, axis=0)
    T /= np.linalg.norm(T, axis=1, keepdims=True)
    T[0] = t0 / np.linalg.norm(t0)
    U = np.zeros_like(points)
    U[0] = u0 - (u0 @ T[0]) * T[0]
    U[0] /= np.linalg.norm(U[0])
    for i in range(len(points) - 1):
        v1 = points[i + 1] - points[i]
        c1 = v1 @ v1
        if c1 < 1e-12:
            U[i + 1] = U[i]
            continue
        rL = U[i] - (2 / c1) * (v1 @ U[i]) * v1
        tL = T[i] - (2 / c1) * (v1 @ T[i]) * v1
        v2 = T[i + 1] - tL
        c2 = v2 @ v2
        U[i + 1] = rL - (2 / c2) * (v2 @ rL) * v2 if c2 > 1e-12 else rL
        U[i + 1] /= np.linalg.norm(U[i + 1])
    V = np.cross(T, U)
    return T, U, V


def _section_props(R):
    """Area and second moments (about the centroid, along the section's u and v) of a closed polygon (n, 2)."""
    x, y = R[:, 0], R[:, 1]
    x1, y1 = np.roll(x, -1), np.roll(y, -1)
    cr = x * y1 - x1 * y
    A = 0.5 * cr.sum()
    cx = ((x + x1) * cr).sum() / (6 * A)
    cy = ((y + y1) * cr).sum() / (6 * A)
    Iyy = ((x ** 2 + x * x1 + x1 ** 2) * cr).sum() / 12 - A * cx ** 2      # about the v axis (bending in u)
    Ixx = ((y ** 2 + y * y1 + y1 ** 2) * cr).sum() / 12 - A * cy ** 2      # about the u axis (bending in v)
    return abs(A), abs(Ixx), abs(Iyy), (cx, cy)


def sweep(face_xz, face_point, n_out, foot, n_foot, leave, arrive, chord_mm, thick_mm, flare, roundness,
          into_wing=0.3, into_foot=1.5, wall_mm=1.2, mid=(0.0, 0.0, 0.0), swell=1.0, around=0.0, face_v=None):
    """The strut as rings (N_PATH + 2, N_RING, 3) and its path, before any check.
    face_xz: the face's outline, (k, 2), absolute: (x, z) in the wing's side face (n_out along +-y), or, with
    face_v (a unit world vector square to n_out and to x), (x, coordinate along face_v) -- e.g. a pylon's
    footprint under the wing: n_out = -z, face_v = +y; face_point: the face's centre (x, y, z); n_out: its
    outward normal; foot: where it lands; n_foot: the surface's outward normal there."""
    face_point, n_out, foot, n_foot = (np.asarray(v, float) for v in (face_point, n_out, foot, n_foot))
    fv = np.array([0.0, 0.0, 1.0]) if face_v is None else np.asarray(face_v, float) / np.linalg.norm(face_v)
    L = float(np.linalg.norm(foot - face_point))
    # straight out of the face for STRAIGHT_MM (its sections square to it: the aerofoil itself), then the curve
    Ps = face_point + STRAIGHT_MM * n_out
    P1 = Ps + max(leave, 0.05) * L * n_out
    P3 = foot - into_foot * n_foot
    P2 = foot + max(arrive, 0.05) * L * n_foot
    s_c = np.linspace(0.0, 1.0, N_PATH)
    # the free mid point; `around` (0-1) draws it toward the corner that takes the strut out of the face
    # first, then down to its foot's level, then across to it -- clear of the wing it leaves (a repair, below)
    Pm = 0.5 * (P1 + P2) + L * np.asarray(mid, float)
    corner = P1 + np.dot(P3 - P1, n_foot) * n_foot           # out of the face, down to its foot's level
    Pm = (1.0 - around) * Pm + around * corner
    pts = np.vstack([face_point + np.outer([-into_wing, 0.0, 0.5 * STRAIGHT_MM], n_out),
                     _bezier(np.array([Ps, P1, Pm, P2, P3]), s_c)])
    s = np.r_[0.0, 0.0, 0.0, s_c]                          # the morph's parameter: the straight part is s = 0
    T, U, V = _rmf(pts, n_out, np.array([1.0, 0.0, 0.0]))
    # the aerofoil, in its face's (u, v) = (x, z) about the face's centre; the strut's own section; the foot
    from shapely.geometry import Polygon
    grown = Polygon(np.asarray(face_xz, float) - np.array([face_point[0], face_point @ fv])).buffer(FACE_GROW_MM, quad_segs=16)
    air = _resample(np.asarray(grown.exterior.coords)[:-1])
    a, b = max(chord_mm / 2, wall_mm / 2), max(thick_mm / 2, wall_mm / 2)
    own = _superellipse(a, b, roundness)
    foot_sec = _superellipse(a * flare, b * flare, roundness)
    hold = 0.0
    vsign = 1.0 if V[0] @ fv >= 0 else -1.0                      # the face's v axis is face_v (z for a side face)
    rings = []
    for i, si in enumerate(s):
        w1 = np.clip((si - hold) / 0.35, 0.0, 1.0)            # aerofoil -> own over the next 35 %
        w1 = w1 * w1 * (3 - 2 * w1)
        w2 = np.clip((si - 0.7) / 0.3, 0.0, 1.0)               # own -> foot over the last 30 %
        w2 = w2 * w2 * (3 - 2 * w2)
        sw = 1.0 + (swell - 1.0) * math.sin(math.pi * si)      # swelling (or slimming) along its length
        sec = (1 - w1) * air + w1 * own * sw
        sec = (1 - w2) * sec + w2 * foot_sec
        # u along the flow (U), v up at the face
        rings.append(pts[i] + np.outer(sec[:, 0], U[i]) + np.outer(sec[:, 1], vsign * V[i]))
    return np.array(rings), pts, (T, U, V), s


def stress_mpa(rings, pts, frames, load_at, load=LOAD_N):
    """The highest bending + axial stress (MPa) over the sections, the strut held at its foot and loaded at
    `load_at` (the wing) by `load` (N): a cantilever, each section carrying the moment of the load about it."""
    T, U, V = frames
    F = np.asarray(load, float)
    # (the section's own axes: rings were laid out along U and +-V)
    worst = 0.0
    for i in range(1, len(pts) - 1):
        R = rings[i] - pts[i]
        P2 = np.c_[R @ U[i], R @ V[i]]
        A, Iu, Iv, (cu, cv) = _section_props(P2)
        if A <= 0 or Iu <= 0 or Iv <= 0:
            return np.inf
        M = np.cross(np.asarray(load_at) - pts[i], F)                         # N mm about the section's centre
        Mu, Mv = M @ U[i], M @ V[i]                                         # bending about u, about v
        du = np.abs(P2[:, 0] - cu).max(); dv = np.abs(P2[:, 1] - cv).max()
        sigma = abs(Mu) * dv / Iu + abs(Mv) * du / Iv + abs(F @ T[i]) / A
        worst = max(worst, sigma)
    return worst


def mesh(rings):
    """A closed solid (manifold3d) from the rings: the union of the convex hulls of each run of three neighbouring
    rings (every section is convex: the aerofoil, the superellipses and the blends between them). Whatever
    the path does -- a bend tighter than the section, a loop through itself -- this is one valid solid, the
    swept section's own volume; a ruled skin through the rings would fold through itself there instead.
    None if it is not one piece (zero-volume slivers of the union are dropped)."""
    import manifold3d as mf
    # three rings to a hull, so neighbouring hulls overlap (two that only touch on a shared ring stay apart)
    Ms = [mf.Manifold.hull_points(np.vstack(rings[i:i + 3])) for i in range(len(rings) - 2)]
    M = mf.Manifold.batch_boolean(Ms, mf.OpType.Add)
    if M.is_empty():
        return None
    pcs = sorted(M.decompose(), key=lambda q: -q.volume())
    if len(pcs) > 1 and pcs[1].volume() > 1e-3 * pcs[0].volume():
        return None
    return pcs[0]


def aim(n_out, out_deg, fwd_deg, base):
    """The direction a strut heads: `base` (straight down for the rear, straight in for the front) leaned
    toward n_out by out_deg and toward +x (forward is -x: fwd_deg leans toward -x) by fwd_deg."""
    base = np.asarray(base, float) / np.linalg.norm(base)
    n_out = np.asarray(n_out, float)
    side = n_out if abs(n_out @ base) < 0.9 else (np.array([0.0, 0.0, -1.0]) if abs(base[2]) < 0.9
                                                  else np.array([0.0, 1.0, 0.0]))   # a pylon straight down: out is +y
    d = math.cos(math.radians(out_deg)) * base + math.sin(math.radians(out_deg)) * side   # inward: lean down/up
    d = math.cos(math.radians(fwd_deg)) * d + math.sin(math.radians(fwd_deg)) * np.array([-1.0, 0.0, 0.0])
    return d / np.linalg.norm(d)


def first_hit(m, o, d):
    """(point, outward normal) of the first surface of m from o along d (mm), or (None, None)."""
    d = np.asarray(d, float) / np.linalg.norm(d)
    loc, ri, tri = m.ray.intersects_location(np.asarray(o, float)[None], d[None], multiple_hits=True)
    if not len(loc):
        return None, None
    k = int(np.argmin(np.linalg.norm(loc - o, axis=1)))
    return loc[k], m.face_normals[tri[k]]


def nearest(m, o):
    """(point, outward normal, distance) of m's nearest surface point to o."""
    import trimesh
    q, dist, tri = trimesh.proximity.closest_point(m, np.asarray(o, float)[None])
    return q[0], m.face_normals[tri[0]], float(dist[0])


def design(face_xz, face_point, n_out, target, base_dir, g, reach_mm=60.0, accept=None, wall_mm=1.2, load_at=None,
           go_around=False, face_v=None):
    """A strut (manifold3d, mm) for genes g (dict: aim_out_deg, aim_fwd_deg, leave, arrive, chord_mm, thick_mm,
    flare, round) from the face to `target` (trimesh, mm), or (None, info). It heads along its aim (leaned back
    towards base_dir when the result breaks a rule: accept(manifold) -> bool), lands on the first surface that
    way (else the target's nearest point), and is thickened until no section is over ALLOW_MPA."""
    info = {}
    o = np.asarray(face_point, float) + 0.5 * np.asarray(n_out, float)
    if target.is_volume and bool(target.contains(o[None])[0]):
        info["in_target"] = True                        # the face is (all but) in what it would join: held
        return None, info
    out, fwd = float(g["aim_out_deg"]), float(g["aim_fwd_deg"])
    chord, thick = max(float(g["chord_mm"]), wall_mm), max(float(g["thick_mm"]), wall_mm)
    for k in np.linspace(1.0, 0.0, 8):
        d = aim(n_out, out * k, fwd * k, base_dir)
        foot, nf = first_hit(target, np.asarray(face_point) + 0.5 * np.asarray(n_out), d)
        if foot is None or np.linalg.norm(foot - face_point) > reach_mm:
            foot, nf, dist = nearest(target, face_point)
            if dist > reach_mm:
                continue
        c, t = chord, thick
        M = None
        for _ in range(12):
            # a strut that breaks a rule is pulled back towards a plain one: its lean, how far it carries on out
            # of the face and into the foot, and its flare, all by the same share k (k = 1: as the search asked);
            # go_around: its mid point drawn by 1 - k toward the corner out of the face then down (clear of the wing)
            rings, pts, frames, _ = sweep(face_xz, face_point, n_out, foot, nf, (k * float(g["leave"]) + (1.0 - k) * LEAVE_AROUND) if go_around
                                          else float(g["leave"]) * (0.1 + 0.9 * k),
                                          float(g["arrive"]) * ((0.1 + 0.9 * k) if go_around else (0.4 + 0.6 * k)), c, t, 1.0 + (float(g["flare"]) - 1.0) * k,
                                          float(g["round"]), wall_mm=wall_mm,
                                          mid=k * np.array([float(g.get("mid_x", 0.0)), float(g.get("mid_y", 0.0)),
                                                            float(g.get("mid_z", 0.0))]),
                                          swell=max(1.0 + (float(g.get("swell", 1.0)) - 1.0) * k, wall_mm / max(min(c, t), 1e-6)),
                                          around=(1.0 - k) if go_around else 0.0, face_v=face_v)
            sig = stress_mpa(rings, pts, frames, load_at if load_at is not None else face_point)
            if sig <= ALLOW_MPA:
                break
            c, t = c * 1.12, t * 1.18                       # thicker (mostly), until the ABS holds
        info.update(stress_mpa=round(float(sig), 2), chord_mm=round(c, 2), thick_mm=round(t, 2))
        if sig > ALLOW_MPA:
            info["too_weak"] = True
            continue
        M = mesh(rings)
        if M is None:
            info["folds"] = True
            continue
        if accept is not None and not accept(M):
            info["leaned_back"] = round(float(k), 2)
            continue
        info.update(foot_mm=[round(float(v), 2) for v in foot], aim_used=[round(out * k, 1), round(fwd * k, 1)],
                    share=round(float(k), 2))
        return M, info
    return None, info


def whole(half_mm):
    """The whole part (mm) from its right half: the strut lands on the real surface, not on the half's cap on
    y = 0 (a strut cast inwards met that cap and stood on the centre line, gallery 7 G00)."""
    import manifold3d as mf
    import support_engine as se
    import wing_loft as wl
    m = half_mm.copy()
    mir = m.copy(); mir.vertices[:, 1] *= -1.0; mir.invert()
    try:
        U = se._to_manifold(m) + se._to_manifold(mir)
        t = wl.to_trimesh(U)
        if t.is_volume:
            return t
    except Exception:  # noqa: BLE001
        pass
    return half_mm
