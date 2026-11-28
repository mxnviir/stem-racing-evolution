"""Runtime patch of the pinned evolution builder (mxnviir/stem-racing-evolution @ e9c896b): the program decides
which parts of the nose are milled foam and which are printed ABS, and the nose meets the body flush.

    import nose_material; nose_material.apply(<evolution checkout>)      # or: python nose_material.py <checkout>

Before: the whole nose ahead of Ref A was one printed part (freeform.build_nose, weighed solid at the printed
density once the front wing's pylons are on it: 6.69 g of the team car's 13.07 g of print), the milled body ended
at Ref A, and the body's mesh rounded its own front edge (marching cubes + Taubin smoothing on the 1 mm grid), so
the printed nose stood up to 2.6 mm proud of the body just behind Ref A.
After (car/body/nose_split.py, written by this patch): the nose's material is decided per car by one plane x_cut:
printed ABS ahead of it (where the front wing runs through the nose or its supports land, plus a load-checked
collar; anything outside the official block; anything under the foam's 3 mm wall), milled foam behind it, as part
of the body's own field (same masks, milling, skin), so the foam nose IS the body and has no joint at Ref A. The
body's field runs on LAP_MM ahead of x_cut and the body's mesh is cut exactly at x_cut, so the face the print is
glued to is sharp and the print, from the same field, meets it flush. A tenon (joints.nose_part's 8 x 8 mm) runs
from the print into a pocket milled from the underside. Masses: the foam nose is in the body (0.163 g/cm3 in
unified_phi.compute_mass_com and in layout_car.mass's body mesh), the print in parts_mass_g["nose"] at the printed
density. The referee measures T8 on the union (the body ahead of Ref A was already measured; T8.7 / T8.8's
"within 5 mm of the nose" now also counts the foam nose) and checks the glued joint.
EVO_NOSE_SPLIT=0 in the environment turns the split off at run time (the builder then behaves as before).
Only fixed-skin builds (--skin-mm: every search build) are split; a mass-sized build keeps the printed nose.

Idempotent: each anchor must be found exactly once (or its patched text already present), else it stops.
See evo_patches/NOSE.md.
"""
from pathlib import Path

TAG = "evo_patches/nose_material.py"

NOSE_SPLIT_SRC = r'''"""nose_split.py -- (written by evo_patches/nose_material.py) the nose's material, decided per car: milled foam (one piece with the body, from the same
official block) wherever it can be, printed ABS only where it must be. Written into the checkout by
evo_patches/nose_material.py; the why, the rules and the numbers are in evo_patches/NOSE.md.

The nose is still ONE drawn shape (the free-form field ahead of Ref A, car/body/freeform.py): the front wing is
designed against it exactly as before (front_wing_engine.build on freeform.build_nose), and every T8 rule is
measured on the union. Only the material inside it is split, by one plane x = x_cut (car frame, mm):

  ahead of x_cut   printed ABS: it holds the front wing (where the wing runs through the nose, where its struts
                   or underside pylons land, plus COLLAR_MM of ABS behind them, the joint checked against the
                   strut generator's load case), and everything outside the official block (the block's front is
                   223 mm ahead of the rear face) or thinner than the foam's 3 mm wall.
  behind x_cut     foam: the body's own field, milled with the body (same masks, same 3-axis ball, same skin),
                   so between x_cut and Ref A there is no joint at all; the foam meets the print at x_cut, flush.

x_cut = the most rearward of those limits, never behind Ref A (x_cut = Ref A: the whole nose is printed, as
before, but the body's front face is cut exactly and the print fills to it, so the junction is flush).
The print is glued to the foam on the whole x_cut face and by a tenon (joints.nose_part's 8 x 8 mm, from the
underside up to the section's centre) in a pocket milled from the bottom.

plan()          the decision (pure numbers, JSON): x_cut, why, the tenon, the joint's stresses.
field_front()   how far ahead of x_cut the body's field runs (LAP_MM), so that the mesh's rounding of the
                field's front edge is cut off and the foam's face at x_cut is sharp.
cut_body()      the milled body's half mesh cut exactly at x_cut, the tenon's pocket taken out (and the tenon).
printed_nose()  the printed piece: the same field ahead of x_cut (OVERLAP_MM into the foam), with the tenon.
foam_weight()   per-cell weight for unified_phi.compute_mass_com: the field's cells ahead of x_cut and in the
                pocket are not foam (they are printed, weighed at the printed density).
"""
from __future__ import annotations

import math

import numpy as np

H_MM = 0.25                  # grid the nose is evaluated on (as freeform.build_nose)
LAP_MM = 5.0                 # the body's field runs this far ahead of x_cut; the mesh is cut at x_cut
OVERLAP_MM = 0.3             # the printed piece runs this far into the foam (one closed solid in the union)
BLOCK_MARGIN_MM = 1.0        # plan(): foam this far inside the block's DESIGN front (rear face - 223) ...
REF_BLOCK_MARGIN_MM = 0.5    # fit_to_block(): ... and this far inside the block as the REFEREE places it on the built
                             # body (the measured bore end + 52 - 223: 1.2-2.3 mm aft of the design front on the
                             # 1 mm grid, sweep3; the bore is drawn a cell long and its closed end gauged at 0.6 r)
COLLAR_MM = 1.2              # ABS behind the last point the wing / its supports touch (load-bearing minimum)
TENON_LEN_MM = 8.0           # joints.nose_part: 8 mm aft of the joint face ...
TENON_HALF_W_MM = 4.0        # ... 8 mm wide (the 6.35 mm ball needs >= 6.35) ...
TENON_MIN_HALF_W_MM = 3.175 + 0.15   # ... never under the ball plus the glue gap
TENON_MIN_H_MM = 3.0
FOAM_WALL_MM = 3.0           # rules_constants.MACHINED_MIN_WALL_MM / geometry_contract.MIN_MACHINED_WALL_MM
THIN_TOL_MM2 = 1.0           # a section with more than this under the 3 mm wall is thin (rounded edges: less)
Y_CAP_MM, Z_CAP_MM = 14.5, 24.5      # unified_phi.build_rule_masks' T8.5.1 cap on the body's field
Z_TRACK_MM = 1.7             # T3.7's 1.5 mm + freeform.MARGIN_MM
NOSE_OVERHANG_MM = 40.0      # T8.2
LOAD_N = (1.0, 0.5, 2.0)     # wing_strut.LOAD_N: per strut / half wing at the wing: drag, side, down
N_HALVES = 2                 # both halves of the wing load the nose (side loads taken the same way: crosswind)
FOAM_BOND_MPA = 0.5          # support_engine.FOAM_BEARING_MPA: the foam side of a glued face (assumption)
ABS_ALLOW_MPA = 20.0         # wing_strut.ALLOW_MPA
BLOCK_LENGTH_MM = 223.0      # rules_constants.BLOCK_LENGTH_MM


def _block():
    from geometry_contract import (CHAMBER_AXIS_Z_MM, MODEL_BLOCK_BORE_ABOVE_BOTTOM_MM,
                                   MODEL_BLOCK_SLOT_DEPTH_MM, MODEL_BLOCK_SLOT_WIDTH_MM)
    z0 = CHAMBER_AXIS_Z_MM - MODEL_BLOCK_BORE_ABOVE_BOTTOM_MM
    return z0, MODEL_BLOCK_SLOT_WIDTH_MM / 2.0, z0 + MODEL_BLOCK_SLOT_DEPTH_MM


def _grid(x0, x1, h=H_MM):
    xs = np.arange(x0, x1 + 0.5 * h, h)
    ys = np.arange(-Y_CAP_MM - 1.0, Y_CAP_MM + 1.0 + 0.5 * h, h)
    ys = ys - ys[np.argmin(np.abs(ys))]                    # a cell centred on y = 0
    zs = np.arange(0.0, Z_CAP_MM + 1.5 + 0.5 * h, h)
    return xs, ys, zs


def drawn_phi(ff, ref_a, X, Y, Z, skin_mm=0.0):
    """The drawn nose and body (mm, negative inside) as the body's field sees it ahead of Ref A: the free-form
    field less the skin, inside T8.5.1's cap (|y| 14.5, z 24.5), T3.7 and T8.2."""
    import freeform as fr
    phi = fr.field(ff, X, Y, Z, ff.x_tip_mm, ff.x_rear_mm) - skin_mm
    cap = np.where(X <= ref_a, np.maximum(np.abs(Y) - Y_CAP_MM, Z - Z_CAP_MM), -1e3)
    box = np.maximum.reduce([cap, Z_TRACK_MM - Z, (ref_a - NOSE_OVERHANG_MM + 0.2) - X])
    return np.maximum(phi, box)


def _foam_ok(Y, Z):
    """Inside the block's section (above its underside, out of its 6 x 6 mm slot): signed, mm, negative inside."""
    zb, sh, st = _block()
    slot = np.minimum(sh - np.abs(Y), st - Z)                 # > 0 inside the slot
    return np.maximum(zb - Z, slot)


def _sections(S, ys, zs, h=H_MM):
    """Per x: area (mm2), centroid z, I about y and z (mm4), extreme fibres (mm), from a boolean (x, y, z) solid."""
    Yg, Zg = np.meshgrid(ys, zs, indexing="ij")
    out = []
    for k in range(S.shape[0]):
        m = S[k]
        n = int(m.sum())
        if n == 0:
            out.append(None)
            continue
        y, z = Yg[m], Zg[m]
        A = n * h * h
        zc = float(z.mean())
        Iy = float(((z - zc) ** 2).sum() * h * h + A * h * h / 12.0)
        Iz = float((y ** 2).sum() * h * h + A * h * h / 12.0)
        out.append(dict(A=A, zc=zc, Iy=Iy, Iz=Iz, cz=float(np.abs(z - zc).max() + h / 2), cy=float(np.abs(y).max() + h / 2),
                        zmin=float(z.min() - h / 2), zmax=float(z.max() + h / 2), ymax=float(np.abs(y).max() + h / 2)))
    return out


def joint_stress(sec, x_joint, load_at, n_halves=N_HALVES):
    """Stresses (MPa) on a glued face (a section `sec` at x_joint) from the wing's load at load_at = (x, z):
    N_HALVES x LOAD_N (drag, side, down). Normal: axial + both bendings; shear: both forces + the torsion."""
    F = np.array(LOAD_N) * n_halves
    F[2] = -F[2]                                               # down
    r = np.array([load_at[0] - x_joint, 0.0, load_at[1] - sec["zc"]])
    M = np.cross(r, F)                                         # N mm about the face's centroid
    A, Iy, Iz = sec["A"], sec["Iy"], sec["Iz"]
    sigma = abs(F[0]) / A + abs(M[1]) * sec["cz"] / Iy + abs(M[2]) * sec["cy"] / Iz
    tau = math.hypot(F[1], F[2]) / A + abs(M[0]) * math.hypot(sec["cy"], sec["cz"]) / (Iy + Iz)
    return dict(sigma_mpa=float(sigma), tau_mpa=float(tau), util_mpa=float(sigma + tau), area_mm2=float(A),
                lever_mm=float(-r[0]), load_n=F.tolist())


def _wing_contact_x(wing_mm, nose_mm, tol=0.3):
    """The most rearward x at which the wing touches the nose (or runs inside it), or None."""
    import trimesh
    pts = np.vstack([wing_mm.vertices, trimesh.sample.sample_surface(wing_mm, 6000, seed=7)[0]])
    near = np.abs(pts[:, 1]) <= nose_mm.bounds[1, 1] + tol
    pts = pts[near & (pts[:, 0] <= nose_mm.bounds[1, 0] + tol)]
    if not len(pts):
        return None
    sd = np.concatenate([trimesh.proximity.signed_distance(nose_mm, pts[i:i + 256]) for i in range(0, len(pts), 256)])
    hit = pts[sd > -tol]
    return float(hit[:, 0].max()) if len(hit) else None


def attachment(ff, x_front_mm, fwing, nose_kw, wheel_design="team_stl"):
    """The front wing as assembly.build makes it (front_wing_engine.build on freeform.build_nose: the same
    calls, so the same wing), and how far back it holds on to the nose."""
    import component_contract as cc
    import freeform as fr
    import front_wing_engine as fwe
    import wheel as wh
    from geometry_contract import FRONT_WHEEL_INNER_Y_MM as FY
    ref_a = cc.ref_plane_A(x_front_mm)
    nm = fr.build_nose(ff, ff.x_tip_mm, ff.x_rear_mm, ref_a, blend_after_mm=nose_kw.get("blend_after_ref_a_mm", 1.0))
    if nm is None:
        return None
    f, r = wh.design(wheel_design)
    R = 0.5 * (f.R + r.R)
    res = fwe.build(fwe.from_dict(fwing), nm["nose"], x_front_mm, (x_front_mm, FY, R), None)
    nose = nm["nose"].copy(); nose.apply_scale(1e3)
    wing = res["wing"].copy(); wing.apply_scale(1e3)
    xs = {}
    if res["pylons"] is not None:
        xs["supports"] = float(res["pylons"].bounds[1, 0] * 1e3)
    xw = _wing_contact_x(wing, nose)
    if xw is not None:
        xs["wing"] = xw
    info = res["info"]
    return dict(x_mm=max(xs.values()) if xs else None, by=xs,
                load_at_mm=[float(info["x_le_root_mm"] + 0.5 * fwing["chord_mm"]), float(info["z_mm"])],
                held_by_nose=bool(info.get("held_by_nose")), root_in_nose=bool(info.get("root_in_nose")),
                struts=len(info.get("struts", [])))


def plan(ff, x_front_mm, W_mm, rear_mm, fwing, nose_kw, skin_mm=0.0, wheel_design="team_stl") -> dict:
    """Where the foam ends and the print begins (x_cut), and why. Deterministic: the same car, the same plan."""
    from scipy import ndimage
    ref_a = x_front_mm - 16.0
    rear_face = x_front_mm + W_mm + rear_mm
    block_front = rear_face - BLOCK_LENGTH_MM
    att = attachment(ff, x_front_mm, fwing, nose_kw, wheel_design)
    x0 = max(ff.x_tip_mm, ref_a - NOSE_OVERHANG_MM) - 1.0
    xs, ys, zs = _grid(x0, ref_a + TENON_LEN_MM + 2.0)
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")
    phi = drawn_phi(ff, ref_a, X, Y, Z, skin_mm)
    S_all = phi < 0
    S = S_all & (_foam_ok(Y, Z) < 0)                           # what foam the block can give
    del X, Y, Z, phi
    # (c) thin: the drawn nose where no ball of the 3 mm wall fits (its opening by a 1.5 mm ball). On the drawn
    # shape: the block's slot and underside are (b), and cut the body behind Ref A the same way
    r = 0.5 * FOAM_WALL_MM / H_MM
    core = ndimage.distance_transform_edt(np.pad(S_all, 1)) > r
    opened = ndimage.distance_transform_edt(~core) <= r
    thin = S_all & ~opened[1:-1, 1:-1, 1:-1]
    thin_area = thin.sum(axis=(1, 2)) * H_MM * H_MM
    del core, opened, thin
    sel = (xs < ref_a) & (thin_area > THIN_TOL_MM2) & S_all.any(axis=(1, 2))
    x_thin = float(xs[sel].max()) + H_MM if sel.any() else float(xs[0])
    limits = {"block": block_front + BLOCK_MARGIN_MM, "thin": x_thin,
              "overhang": ref_a - NOSE_OVERHANG_MM + 0.2}
    if att and att["x_mm"] is not None:
        limits["attachment"] = att["x_mm"] + COLLAR_MM
    secs = _sections(S, ys, zs)
    secs_all = _sections(S_all, ys, zs)
    load_at = att["load_at_mm"] if att else [ref_a - 10.0, 10.0]
    x_c = max(limits.values())
    why = max(limits, key=limits.get)
    bond = None
    if x_c < ref_a - 1e-9:
        x_first = x_c
        x_c, why_c = ref_a, "bond"                             # unless a face below carries the loads
        for k in range(int(np.searchsorted(xs, x_first - 1e-9)), len(xs)):
            if xs[k] >= ref_a - 1e-9:
                break
            sec = secs[k]
            if sec is None or secs_all[k] is None:
                continue
            b = joint_stress(sec, xs[k], load_at)              # (d) the glued face carries the wing's loads
            if b["util_mpa"] <= FOAM_BOND_MPA:
                x_c, bond, why_c = float(xs[k]), b, why
                break
        if why_c == "bond":
            limits["bond"] = float(x_c)
        elif x_c > x_first + H_MM + 1e-9:
            limits["bond"] = float(x_c); why_c = "bond"
        why = why_c
    else:
        x_c = ref_a
    k = int(np.argmin(np.abs(xs - x_c)))
    sec = secs_all[k] if secs_all[k] is not None else secs[k]
    if bond is None and sec is not None:
        bond = joint_stress(sec, x_c, load_at)
    bond = dict(bond or {}, abs_sigma_mpa=(bond or {}).get("sigma_mpa"))
    # the tenon: from the underside to the section's centre (joints.nose_part), its pocket's walls >= 3 mm
    tenon, walls = None, {}
    kk = (xs >= x_c - 1e-9) & (xs <= x_c + TENON_LEN_MM + 1e-9)
    if sec is not None and kk.any():
        zc = 0.5 * (sec["zmin"] + sec["zmax"])
        hw = TENON_HALF_W_MM
        Sk = S[kk]
        Yg, Zg = np.meshgrid(ys, zs, indexing="ij")
        top = min(float(Zg[m & (np.abs(Yg) <= hw)].max()) for m in Sk if (m & (np.abs(Yg) <= hw)).any())
        bot = max(float(Zg[m & (np.abs(Yg) <= hw)].min()) for m in Sk if (m & (np.abs(Yg) <= hw)).any())
        zc = min(zc, top - FOAM_WALL_MM)
        band = (Zg >= bot) & (Zg <= zc)
        side = min(float(np.abs(Yg[m & band]).max()) for m in Sk if (m & band).any()) if zc > bot else 0.0
        hw = min(hw, side - FOAM_WALL_MM)
        walls = {"top_mm": top - zc, "side_mm": side - hw, "height_mm": zc - bot}
        if hw >= TENON_MIN_HALF_W_MM and zc - bot >= TENON_MIN_H_MM:
            tenon = {"x0": float(x_c), "x1": float(x_c + TENON_LEN_MM), "hw": float(hw),
                     "z0": float(sec["zmin"] - 1.0), "z1": float(zc)}
    foam_cm3 = float(S[(xs >= x_c) & (xs < ref_a)].sum() * H_MM ** 3 / 1e3)
    return {"version": 1, "ref_a_mm": float(ref_a), "x_cut_mm": float(x_c), "foam_len_mm": float(ref_a - x_c),
            "x_field_front_mm": float(x_c - LAP_MM), "block_front_mm": float(block_front), "skin_mm": float(skin_mm),
            "decided_by": why, "limits_mm": {k_: float(v) for k_, v in limits.items()},
            "attachment": att, "bond": bond, "tenon": tenon, "tenon_walls_mm": walls,
            "foam_nose_drawn_cm3": foam_cm3}


def field_front(p) -> float:
    return p["x_field_front_mm"]


def _box_sdf(X, aY, Z, t):
    return np.maximum.reduce([t["x0"] - X, X - t["x1"], aY - t["hw"], t["z0"] - Z, Z - t["z1"]])


def foam_weight(geom, p, dx_mm):
    """Weight (0..1) of each cell of geom's grid as foam: 0 ahead of x_cut and in the tenon's pocket."""
    r = geom.region
    xs = r.origin_m[0] * 1e3 + np.arange(r.shape[0]) * dx_mm
    ys = r.origin_m[1] * 1e3 + np.arange(r.shape[1]) * dx_mm
    zs = r.origin_m[2] * 1e3 + np.arange(r.shape[2]) * dx_mm
    wx = np.clip((xs - p["x_cut_mm"]) / dx_mm + 0.5, 0.0, 1.0)[:, None, None]
    w = np.broadcast_to(wx, r.shape).copy()
    t = p.get("tenon")
    if t:
        X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")
        inside = np.clip(-_box_sdf(X, np.abs(Y), Z, t) / dx_mm + 0.5, 0.0, 1.0)
        w *= 1.0 - inside
    return w


def _man(mesh_m):
    import manifold3d as mf
    v = np.asarray(mesh_m.vertices, np.float64) * 1e3
    return mf.Manifold(mf.Mesh(vert_properties=v.astype(np.float32), tri_verts=np.asarray(mesh_m.faces, np.uint32)))


def _tm(M, scale=1e-3):
    import trimesh
    g = M.to_mesh()
    t = trimesh.Trimesh(np.asarray(g.vert_properties)[:, :3], np.asarray(g.tri_verts), process=True)
    t.apply_scale(scale)
    return t


def _cube(x0, x1, y0, y1, z0, z1):
    import manifold3d as mf
    return mf.Manifold.cube([x1 - x0, y1 - y0, z1 - z0]).translate([x0, y0, z0])


def _closed(M, scale=1e-3):
    """The largest piece of a manifold as a trimesh that stays one closed solid when written to an ASCII STL and read
    back (edges under 0.01 mm merged, as freeform.build_nose: Boolean cuts leave 2 um slivers an STL collapses)."""
    import io
    import trimesh
    for tol in (0.0, 0.01, 0.02, 0.05):             # 0: the Boolean's own mesh (the body untouched away from the cut)
        Q = max((M.simplify(tol) if tol > 0 else M).decompose(), key=lambda q: q.volume())
        t = _tm(Q, scale)
        back = trimesh.load(io.BytesIO(t.export(file_type="stl_ascii").encode()), file_type="stl")
        if t.is_volume and back.is_volume and len(back.split(only_watertight=False)) == 1:
            return t
        if tol == 0.0 and t.is_volume:
            # what the file's rounding leaves is the solid and a few zero-volume two-face slivers at the cut:
            # the solid alone, as read back (no re-meshing of the rest of the body)
            pcs = sorted(back.split(only_watertight=False), key=lambda q: -abs(q.volume))
            if pcs and pcs[0].is_volume and sum(abs(q.volume) for q in pcs[1:]) < 1e-3 * scale ** 3:   # 0.001 mm3
                again = trimesh.load(io.BytesIO(pcs[0].export(file_type="stl_ascii").encode()), file_type="stl")
                if again.is_volume and len(again.split(only_watertight=False)) == 1:
                    return pcs[0]
    return t


def referee_block_front(half_m):
    """The official block's front (mm) exactly as the referee places it (legality.manufactured_checks): the
    chamber measured on the body, both halves (legality.chamber_measures on as_drawn_solids' union), its rear face
    at the bore's end + 52 mm, the block 223 mm long. None if the referee's module or the chamber is not found."""
    try:
        import legality
        import rules_constants as rc
        import trimesh
    except ImportError:
        return None
    m = half_m.copy()
    m.apply_scale(1e3)
    mir = m.copy()
    mir.vertices[:, 1] *= -1.0
    mir.invert()
    try:
        full = trimesh.boolean.union([m, mir], engine="manifold")
    except ValueError:
        return None
    ch = legality.chamber_measures(full)
    if ch is None:
        return None
    return float(ch["rim_x_mm"] - ch["depth_mm"] + rc.BLOCK_BORE_DEPTH_MM - rc.BLOCK_LENGTH_MM)


def fit_to_block(half_m, p):
    """The plan made to agree with the referee's block on the built body (half_m, the milled half before the cut,
    m): the foam's face moved aft to REF_BLOCK_MARGIN_MM behind the block's front as the referee measures it (on
    the H_MM grid, never behind Ref A), the tenon with it. Returns the (new) plan; the old one is not changed."""
    f = referee_block_front(half_m)
    q = dict(p)
    if f is None:
        return q
    q["block_front_design_mm"] = p.get("block_front_design_mm", p["block_front_mm"])
    q["block_front_mm"] = f
    need = f + REF_BLOCK_MARGIN_MM
    q["limits_mm"] = dict(p["limits_mm"], block_referee=float(need))
    if p["x_cut_mm"] < need - 1e-9 and p["x_cut_mm"] < p["ref_a_mm"] - 1e-9:
        x = float(min(math.ceil(need / H_MM - 1e-9) * H_MM, p["ref_a_mm"]))
        dx = x - p["x_cut_mm"]
        q.update(x_cut_mm=x, foam_len_mm=float(p["ref_a_mm"] - x), decided_by="block_referee",
                 moved_for_referee_block_mm=float(dx))
        if p.get("tenon"):
            q["tenon"] = dict(p["tenon"], x0=p["tenon"]["x0"] + dx, x1=p["tenon"]["x1"] + dx)
    return q


def cut_body(half_m, p):
    """(the milled body's right half cut exactly at x_cut with the tenon's pocket, the tenon's right half) in
    metres; the tenon is the body's own shape inside the pocket box (None if there is no tenon)."""
    M = _man(half_m)
    M = M ^ _cube(p["x_cut_mm"], 1e3, -1.0, 1e3, -1e3, 1e3)
    ten = None
    t = p.get("tenon")
    if t:
        box = _cube(t["x0"] - 0.01, t["x1"], -1.0, t["hw"], t["z0"], t["z1"])
        T = M ^ box
        if T.volume() > 1.0:
            ten = _closed(T)
            M = M - box
    return _closed(M), ten


def printed_nose(ff, p, tenon_half_m=None):
    """The printed piece: the field ahead of x_cut + OVERLAP_MM (freeform.build_nose's mesh, the body's own
    limits), with the tenon. {"nose": right half (m), "_full": both halves (m)} or None."""
    import manifold3d as mf
    import trimesh
    from skimage import measure
    ref_a, x1 = p["ref_a_mm"], p["x_cut_mm"] + OVERLAP_MM
    x0 = max(ff.x_tip_mm, ref_a - NOSE_OVERHANG_MM) - 1.0
    xs, ys, zs = _grid(x0, x1 + H_MM)
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")
    phi = np.maximum(drawn_phi(ff, ref_a, X, Y, Z, p.get("skin_mm", 0.0)), X - x1)
    del X, Y, Z
    if (phi < 0).sum() < 10:
        return None
    phi = np.pad(phi, 1, constant_values=1.0)
    v, f, _, _ = measure.marching_cubes(phi, 0.0, spacing=(H_MM,) * 3)
    v = v - H_MM + np.array([xs[0], ys[0], zs[0]])
    full = trimesh.Trimesh(v, f[:, ::-1], process=True)
    if full.volume < 0:
        full.invert()
    M = mf.Manifold(mf.Mesh(vert_properties=np.asarray(full.vertices, np.float32), tri_verts=np.asarray(full.faces, np.uint32)))
    M = sorted(M.decompose(), key=lambda q: -q.volume())[0].simplify(0.01)
    if tenon_half_m is not None:
        T = _man(tenon_half_m)
        T = T + T.mirror([0.0, 1.0, 0.0])
        M = (M + T).simplify(0.01)
    half = M ^ _cube(-500.0, 500.0, 0.0, 500.0, -500.0, 500.0)
    return {"nose": _tm(half), "_full": _tm(M)}


def printed_mass_com(half_m, p, rho_g_cm3, shell_wall_mm=None):
    """(kg, COM m) of a printed nose (right half, m: the piece, its tenon, any pylons): what is ahead of x_cut and
    in the tenon's pocket, both halves; the OVERLAP_MM into the foam is the foam's. shell_wall_mm: weighed as a
    shell of that wall (nose.nose_mass_com's rule for a nose with no pylons), the tenon solid."""
    M = _man(half_m)
    ahead = M ^ _cube(-1e3, p["x_cut_mm"], -1e3, 1e3, -1e3, 1e3)
    parts = [ahead]
    t = p.get("tenon")
    if t:
        parts.append(M ^ _cube(p["x_cut_mm"], t["x1"] + 0.01, -1.0, t["hw"] + 0.01, t["z0"] - 0.01, t["z1"] + 0.01))
    vols = [q.volume() for q in parts]                         # mm3, right half
    if shell_wall_mm is not None and vols[0] > 0:
        a = _tm(ahead, 1.0)
        side = a.face_normals[:, 0] < 0.99
        side &= ~(np.abs(a.triangles_center[:, 1]) < 1e-6)     # nor the symmetry cap
        vols[0] = min(float(a.area_faces[side].sum()) * shell_wall_mm, vols[0])
    vol = sum(vols)
    if vol <= 0:
        return 0.0, (p["x_cut_mm"] / 1e3, 0.0, 0.0)
    cs = []
    for q in parts:
        b = _tm(q, 1.0)
        cs.append(b.center_mass if b.is_volume and len(b.faces) else np.zeros(3))
    c = sum(v * np.asarray(cc_) for v, cc_ in zip(vols, cs)) / vol
    return 2 * vol / 1e3 * rho_g_cm3 * 1e-3, (float(c[0]) / 1e3, 0.0, float(c[2]) / 1e3)


def foam_mass_com(body_half_m, p, rho_g_cm3=0.163):
    """(kg, COM m) of the foam nose: the milled body ahead of Ref A, both halves (it is in the body's mass)."""
    M = _man(body_half_m) ^ _cube(-1e3, p["ref_a_mm"], -1e3, 1e3, -1e3, 1e3)
    v = M.volume()
    if v <= 0:
        return 0.0, (p["ref_a_mm"] / 1e3, 0.0, 0.0)
    c = _tm(M, 1.0).center_mass
    return 2 * v / 1e3 * rho_g_cm3 * 1e-3, (float(c[0]) / 1e3, 0.0, float(c[2]) / 1e3)


def gates(p) -> dict:
    """Margins (+ passes) the builder measures on its own plan (the referee copies them in as P4.*)."""
    g = {}
    b = p.get("bond") or {}
    if b.get("util_mpa") is not None:
        g["nose_joint_foam_stress"] = FOAM_BOND_MPA - b["util_mpa"]
        g["nose_joint_abs_stress"] = ABS_ALLOW_MPA - b["util_mpa"]
    if p["foam_len_mm"] > 0:
        g["nose_foam_inside_block"] = p["x_cut_mm"] - p["block_front_mm"]       # (the referee's, once fitted)
        g["nose_foam_not_thin"] = p["x_cut_mm"] - p["limits_mm"]["thin"]
        if p.get("attachment") and p["attachment"].get("x_mm") is not None:
            g["nose_foam_behind_wing_attachment"] = p["x_cut_mm"] - p["attachment"]["x_mm"] - COLLAR_MM
    w = p.get("tenon_walls_mm") or {}
    if p.get("tenon") and w:
        g["nose_tenon_pocket_walls"] = min(w["top_mm"], w["side_mm"]) - FOAM_WALL_MM
    return {k: float(v) for k, v in g.items()}
'''

PATCHES = [
    # ---------------------------------------------------------------- the body: its field runs past Ref A
    ("car/body/param_body.py",
     """          cargo_placement: dict | None = None, lattice: dict | None = None, freeform=None):
    \"\"\"UnifiedGeometry whose body is `bp`, with every rule mask applied.
""",
     """          cargo_placement: dict | None = None, lattice: dict | None = None, freeform=None,
          nose_split=None):
    \"\"\"UnifiedGeometry whose body is `bp`, with every rule mask applied.

    nose_split: the nose's material plan (car/body/nose_split.py, evo_patches/nose_material.py): the body's field
      runs on to the plan's field front ahead of Ref A (the foam nose is body), labelled and weighed as foam.
"""),
    ("car/body/param_body.py",
     """    x_rear_mm = base.landmarks["rear_face_m"] * 1e3
""",
     """    x_rear_mm = base.landmarks["rear_face_m"] * 1e3
    # evo_patches/nose_material.py: the foam part of the nose is the body's own field; it ends at the plan's field
    # front (LAP_MM ahead of x_cut, where export_body cuts the mesh), foam (main body) not nose (printed) density
    x_front_cut = x_ref_a
    if nose_split is not None:
        import nose_split as _nsx
        x_front_cut = _nsx.field_front(nose_split)
        for _g in (base, mill_base):
            if _g is not None:
                _g.labels[_g.labels == up.LABEL_NOSE] = up.LABEL_MAIN_BODY
                _g.landmarks["nose_split"] = nose_split
"""),
    ("car/body/param_body.py",
     """        phi0 = np.maximum(phi0, np.maximum(x_ref_a - X, X - x_rear_mm))      # the milled body: Ref A .. rear face
""",
     """        phi0 = np.maximum(phi0, np.maximum(x_front_cut - X, X - x_rear_mm))  # the milled body: Ref A (or the foam nose) .. rear face
"""),
    ("car/body/param_body.py",
     """        field = np.maximum(field, np.maximum(x_ref_a - X, X - x_rear_mm))
        if bar_to is not None:
""",
     """        field = np.maximum(field, np.maximum(x_front_cut - X, X - x_rear_mm))   # (evo_patches/nose_material.py)
        if bar_to is not None:
"""),
    # ---------------------------------------------------------------- T8.5.1's cap on the body: ahead of Ref A only
    ("car/body/param_body.py",
     """    ugk = {"with_cargo": with_cargo, "cargo_placement": cargo_placement}
""",
     """    ugk = {"with_cargo": with_cargo, "cargo_placement": cargo_placement}
    # evo_patches/nose_material.py: with the nose's material split, no printed cone takes its root from the body's
    # section behind Ref A, so T8.5.1's cap (25 / 15 mm, 0.5 inside) stops at Ref A, where the rule stops; aft of
    # it, it made a 0.8 mm step in the body's top 2-3 mm behind Ref A
    up.T851_CAP_AFT_MM = 0.0 if nose_split is not None else 1.0
"""),
    ("car/body/unified_phi.py",
     """def apply_build_rules(geom: "UnifiedGeometry") -> None:
""",
     """T851_CAP_AFT_MM = 1.0     # evo_patches/nose_material.py: how far aft of Ref A the T8.5.1 cap runs (+ a cell; 0: to Ref A)


def apply_build_rules(geom: "UnifiedGeometry") -> None:
"""),
    ("car/body/unified_phi.py",
     """    air = (X <= x_ref_a + 1.0 + d * 1e3) & ((np.abs(Y) > 14.5) | (Z > 24.5))
""",
     """    aft = T851_CAP_AFT_MM + d * 1e3 if T851_CAP_AFT_MM > 0 else 0.0      # (evo_patches/nose_material.py)
    air = (X <= x_ref_a + aft) & ((np.abs(Y) > 14.5) | (Z > 24.5))
"""),
    ("car/body/exact_zones.py",
     """    return np.minimum(ref_A_mm + 1.0 - X, np.maximum(aY - 14.5, Z - 24.5))
""",
     """    import unified_phi as _up                    # (evo_patches/nose_material.py: T851_CAP_AFT_MM)
    return np.minimum(ref_A_mm + getattr(_up, "T851_CAP_AFT_MM", 1.0) - X, np.maximum(aY - 14.5, Z - 24.5))
"""),
    # ---------------------------------------------------------------- the milling: no relief of T8.5.1's corner
    ("car/body/machining.py",
     """    Hm = H & ~drilled & ~border & ~pocket
""",
     """    Hm = H & ~drilled & ~border & ~pocket
    if (getattr(geom, "landmarks", None) or {}).get("nose_split"):
        # evo_patches/nose_material.py: T8.5.1's ceilings over the nose (the envelope's 25 / 15 mm, the body's cap
        # 0.5 inside) end at Ref A in a square corner no ball fits. Relieved "for every body" (below), that corner
        # cut a ball-sized dip into the top of whatever body stands there: the 2-3 mm lip under the old printed nose,
        # and the same dip in a foam nose. The drawn nose stands clear of the ceilings, so they stay air (S2 &= ~H)
        # but are not relieved; a body that did reach one is cut flat to it.
        _r = geom.region
        _xs = (_r.origin_m[0] + np.arange(_r.shape[0]) * up.GRID_SPACING_M) * 1e3
        _ys = (_r.origin_m[1] + np.arange(_r.shape[1]) * up.GRID_SPACING_M) * 1e3
        _zs = (_r.origin_m[2] + np.arange(_r.shape[2]) * up.GRID_SPACING_M) * 1e3
        _xa = geom.landmarks["ref_plane_A_m"] * 1e3
        Hm &= ~((_xs <= _xa + d_mm + 1e-6)[:, None, None]
                & ((np.abs(_ys) > 14.0)[None, :, None] | (_zs > 24.0)[None, None, :]))
"""),
    # ---------------------------------------------------------------- the body's mass: not the printed part
    ("car/body/unified_phi.py",
     """    fill = fill_fraction(geom.phi.grid, dx)

    for name in ("nose", "sidepod", "rearpod", "main_body"):
""",
     """    fill = fill_fraction(geom.phi.grid, dx)
    # evo_patches/nose_material.py: the field's cells ahead of the foam nose's face and in the tenon's pocket are
    # cut off the milled body (run_car.export_body): they are printed, and weighed with the print
    _ns = (getattr(geom, "landmarks", None) or {}).get("nose_split")
    if _ns:
        import nose_split as _nsx
        fill = fill * _nsx.foam_weight(geom, _ns, dx * 1e3)

    for name in ("nose", "sidepod", "rearpod", "main_body"):
"""),
    # ---------------------------------------------------------------- run_car: plan, body, mesh cut, parts
    ("car/build/run_car.py",
     """                 freeform=_freeform(a))
    return g, {"param_body": bp.as_dict()""",
     """                 freeform=_freeform(a), nose_split=getattr(a, "_nose_split", None))   # (evo_patches/nose_material.py)
    return g, {"param_body": bp.as_dict()"""),
    ("car/build/run_car.py",
     """def part4_kwargs(a) -> dict:
""",
     """def _nose_split(a, nk):
    \"\"\"evo_patches/nose_material.py: the nose's material plan for this car (car/body/nose_split.py), once per run:
    a free-form car with the searches' front wing and a fixed skin (every search build). None: all printed.\"\"\"
    import os
    if hasattr(a, "_nose_split"):
        return a._nose_split
    plan = None
    ff = _freeform(a)
    fj = json.loads(a.fwing_json) if getattr(a, "fwing_json", None) else {}
    if (os.environ.get("EVO_NOSE_SPLIT", "1") != "0" and ff is not None and "pieces" in fj
            and getattr(a, "skin_mm", None) is not None):
        import nose_split as nsx
        rear = a.rear_mm if getattr(a, "rear_mm", None) is not None else 56.0      # Ref B + 40
        plan = nsx.plan(ff, a.x_front, a.W, rear, fj, nk, a.skin_mm)
        (Path(a.out) / "nose_split.json").write_text(json.dumps(plan, indent=1, default=float))
        log(f"  nose: printed ahead of x {plan['x_cut_mm']:.2f} mm ({plan['decided_by']}), foam "
            f"{plan['foam_len_mm']:.2f} mm to Ref A, tenon {'yes' if plan.get('tenon') else 'no'}")
    a._nose_split = plan
    return plan


def part4_kwargs(a) -> dict:
"""),
    ("car/build/run_car.py",
     """    kw = {"support": sup, "nose": ns.NoseCone(**nk), "freeform": _freeform(a),
          "supports": not getattr(a, "no_supports", False), "tethers": not getattr(a, "no_tethers", False)}
""",
     """    kw = {"support": sup, "nose": ns.NoseCone(**nk), "freeform": _freeform(a),
          "supports": not getattr(a, "no_supports", False), "tethers": not getattr(a, "no_tethers", False)}
    kw["nose_split"] = _nose_split(a, nk)            # evo_patches/nose_material.py
"""),
    ("car/build/run_car.py",
     """    half = up.extract_half_surface(geom)
    half.vertices[half.vertices[:, 1] < 0, 1] = 0.0
""",
     """    half = up.extract_half_surface(geom)
    half.vertices[half.vertices[:, 1] < 0, 1] = 0.0
    _ns = (getattr(geom, "landmarks", None) or {}).get("nose_split")
    if _ns:     # evo_patches/nose_material.py: the foam nose's face cut exactly at x_cut, the tenon's pocket
        import nose_split as _nsx
        # the foam inside the block as the referee places it on THIS body (its measured bore), the plan updated
        _ns = _nsx.fit_to_block(half, _ns)
        geom.landmarks["nose_split"] = _ns
        path.with_name("nose_split.json").write_text(json.dumps(_ns, indent=1, default=float))
        half, _ten = _nsx.cut_body(half, _ns)
        _tp = path.with_name("nose_tenon_half.stl")
        if _ten is not None:
            _ten.export(str(_tp))
        elif _tp.exists():
            _tp.unlink()
"""),
    # ---------------------------------------------------------------- assembly: the printed piece
    ("car/parts/assembly.py",
     """          tether_offsets_mm: tuple = (5.0, 5.0), freeform=None) -> dict:
""",
     """          tether_offsets_mm: tuple = (5.0, 5.0), freeform=None, nose_split=None) -> dict:
"""),
    ("car/parts/assembly.py",
     """    nm = None
    if nose is not None:
""",
     """    nm = nm_drawn = None             # evo_patches/nose_material.py: nm_drawn, the whole drawn nose (foam + print)
    if nose is not None:
"""),
    ("car/parts/assembly.py",
     """        gates.update(ns.nose_gates(nose, nm, x_front_mm))
        _export_half(nm["nose"], out / "nose.stl")
        surfaces.append({"name": "nose", "stl": str(out / "nose.stl"), "rotating": None})
        parts_mass["nose"] = ns.nose_mass_com(nose, nm["_full"])
""",
     """        gates.update(ns.nose_gates(nose, nm, x_front_mm))
        nm_drawn = nm
        if nose_split is not None and freeform is not None:
            # evo_patches/nose_material.py: the printed piece only (ahead of x_cut, and the tenon); the foam nose
            # is the body's (car/body/nose_split.py). The wing is designed against the whole drawn nose, as before.
            import nose_split as nsx
            _pj = Path(body_half_stl).with_name("nose_split.json")      # as fitted to the referee's block
            if _pj.exists():
                nose_split = json.loads(_pj.read_text())
            _tp = Path(body_half_stl).with_name("nose_tenon_half.stl")
            nm = nsx.printed_nose(freeform, nose_split, trimesh.load(str(_tp), force="mesh") if _tp.exists() else None)
            if nm is None:
                raise ValueError("illegal car: no printed nose ahead of the foam")
            gates.update(nsx.gates(nose_split))
        _export_half(nm["nose"], out / "nose.stl")
        surfaces.append({"name": "nose", "stl": str(out / "nose.stl"), "rotating": None})
        if nm is not nm_drawn:       # (evo_patches/nose_material.py) the print, as nose_mass_com weighs a nose
            parts_mass["nose"] = nsx.printed_mass_com(nm["nose"], nose_split, cc.DENSITY_G_CM3[nose.material],
                                                      shell_wall_mm=nose.wall_mm)
        else:
            parts_mass["nose"] = ns.nose_mass_com(nose, nm["_full"])
"""),
    ("car/parts/assembly.py",
     """        res = fwe.build(front, nm["nose"], x_front_mm, (x_front_mm, FY, R), body)
""",
     """        res = fwe.build(front, (nm_drawn or nm)["nose"], x_front_mm, (x_front_mm, FY, R), body)   # (nose_material)
"""),
    ("car/parts/assembly.py",
     """            parts_mass["nose"] = _full_mass_com([nose_with], cc.DENSITY_G_CM3[PART_MATERIAL])
""",
     """            parts_mass["nose"] = _full_mass_com([nose_with], cc.DENSITY_G_CM3[PART_MATERIAL])
            if nm is not nm_drawn:   # evo_patches/nose_material.py: the print only, not its 0.3 mm into the foam
                parts_mass["nose"] = nsx.printed_mass_com(nose_with, nose_split, cc.DENSITY_G_CM3[PART_MATERIAL])
"""),
    ("car/parts/assembly.py",
     """    # -- rear wing ---------------------------------------------------------
    forward = [body, halo] + sup
""",
     """    if nm is not None and nm is not nm_drawn:
        # evo_patches/nose_material.py: the nose's two materials, for the reports and the race model's breakdown
        fk, fc = nsx.foam_mass_com(body, nose_split)
        info["nose_split"] = dict(nose_split, foam_g=fk * 1e3, foam_com_mm=[v * 1e3 for v in fc],
                                  printed_g=parts_mass["nose"][0] * 1e3,
                                  printed_com_mm=[v * 1e3 for v in parts_mass["nose"][1]])

    # -- rear wing ---------------------------------------------------------
    forward = [body, halo] + sup
"""),
    # ---------------------------------------------------------------- the referee
    ("car/build/legality.py",
     """def _visibility_checks(assembly: dict, solids: dict) -> dict:
""",
     """def _foam_nose_half_width(assembly: dict, solids: dict) -> float:
    \"\"\"evo_patches/nose_material.py: the widest point of the body ahead of Ref A (a foam nose), mm; 0 if none.
    T8.7's and T8.8's "within 5 mm of the nose" is measured on the whole nose, printed and foam.\"\"\"
    body = solids.get("machined_body")
    if body is None:
        return 0.0
    v = np.asarray(body.vertices)
    a = v[v[:, 0] < assembly["x_front_mm"] - 16.0 - 0.5]
    return float(np.abs(a[:, 1]).max()) if len(a) else 0.0


NOSE_JOINT_LOAD_N = (1.0, 0.5, 2.0)    # the strut generator's case (wing_strut.LOAD_N) per half wing: drag, side, down
NOSE_JOINT_HALVES = 2
NOSE_JOINT_FOAM_MPA = 0.5              # the foam side of a glued face (support_engine.FOAM_BEARING_MPA): assumption
NOSE_TENON_MIN_MM = 4.0                # joints.cut's embed: a tenon at least this far into the foam


def _yz_section(m_mm, x):
    \"\"\"The (y, z) section of a closed mesh (mm) at x, both halves, as a shapely geometry (None if empty).\"\"\"
    import manifold3d as mf
    from shapely.affinity import scale
    from shapely.geometry import Polygon
    M = mf.Manifold(mf.Mesh(vert_properties=np.asarray(m_mm.vertices, np.float32),
                            tri_verts=np.asarray(m_mm.faces, np.uint32)))
    cs = M.warp_batch(lambda V: np.c_[V[:, 1], V[:, 2], V[:, 0]]).slice(float(x))
    g = None
    for ring in cs.to_polygons():
        if len(ring) < 3:
            continue
        p = Polygon(np.asarray(ring)).buffer(0)
        g = p if g is None else g.symmetric_difference(p)
    if g is None or g.is_empty:
        return None
    return g.union(scale(g, xfact=-1.0, origin=(0, 0))).buffer(0)


def nose_joint_checks(assembly: dict, body_mm, nose_mm, fwing_mm) -> dict:
    \"\"\"evo_patches/nose_material.py: the printed nose glued to the foam nose, measured on the files: the face
    they share (the body's section just behind its front face, inside the print's just ahead of it), the stress on
    it from the strut generator's load case at the wing (NOSE_JOINT_HALVES x NOSE_JOINT_LOAD_N), and the tenon's
    reach into the foam.\"\"\"
    from shapely import contains_xy
    r = {}
    xf = float(body_mm.bounds[0, 0])
    a, b = _yz_section(body_mm, xf + 0.2), _yz_section(nose_mm, xf - 0.2)
    face = a.intersection(b) if a is not None and b is not None else None
    if face is None or face.is_empty or fwing_mm is None:
        r["nose_joint_bond_stress"] = _entry(-1.0, "no glued face between the printed nose and the foam")
        return r
    y0, z0, y1, z1 = face.bounds
    h = 0.1
    Yg, Zg = np.meshgrid(np.arange(y0, y1, h) + h / 2, np.arange(z0, z1, h) + h / 2, indexing="ij")
    inside = contains_xy(face, Yg, Zg)
    y, z = Yg[inside], Zg[inside]
    A = inside.sum() * h * h
    zc = float(z.mean())
    Iy, Iz = float(((z - zc) ** 2).sum() * h * h), float((y ** 2).sum() * h * h)
    cz, cy = float(np.abs(z - zc).max()), float(np.abs(y).max())
    wc = fwing_mm.center_mass if fwing_mm.is_volume else fwing_mm.bounds.mean(axis=0)
    F = np.array(NOSE_JOINT_LOAD_N) * NOSE_JOINT_HALVES * np.array([1.0, 1.0, -1.0])
    M = np.cross(np.array([wc[0] - xf, 0.0, wc[2] - zc]), F)
    sigma = abs(F[0]) / A + abs(M[1]) * cz / Iy + abs(M[2]) * cy / Iz
    tau = float(np.hypot(F[1], F[2])) / A + abs(M[0]) * float(np.hypot(cy, cz)) / (Iy + Iz)
    r["nose_joint_bond_stress"] = _entry(NOSE_JOINT_FOAM_MPA - (sigma + tau),
                                         f"glued face at x {xf:.2f} mm, {A:.0f} mm2: {sigma:.3f} MPa normal + "
                                         f"{tau:.3f} shear under {NOSE_JOINT_HALVES} x {NOSE_JOINT_LOAD_N} N at the wing "
                                         f"(foam {NOSE_JOINT_FOAM_MPA} MPa)")
    if (assembly.get("info") or {}).get("nose_split", {}).get("tenon"):
        import trimesh
        aft = trimesh.intersections.slice_mesh_plane(nose_mm, [1.0, 0, 0], [xf + 0.5, 0, 0], cap=True)
        L = float(aft.bounds[1, 0] - xf) if aft is not None and len(aft.vertices) else 0.0
        r["nose_has_tenon"] = _entry(L - NOSE_TENON_MIN_MM, f"tenon {L:.2f} mm into the foam, in a pocket from the bottom")
    return r


def _visibility_checks(assembly: dict, solids: dict) -> dict:
"""),
    ("car/build/legality.py",
     """        nose_y = half("nose").bounds[1, 1] if "nose" in stl else 0.0
""",
     """        nose_y = half("nose").bounds[1, 1] if "nose" in stl else 0.0
        nose_y = max(nose_y, _foam_nose_half_width(assembly, solids))      # evo_patches/nose_material.py
"""),
    ("car/build/legality.py",
     """        nose_y = half("nose").bounds[1, 1] if "nose" in stl else 1.5
""",
     """        nose_y = half("nose").bounds[1, 1] if "nose" in stl else 1.5
        nose_y = max(nose_y, _foam_nose_half_width(assembly, solids))      # evo_patches/nose_material.py
"""),
    ("car/build/legality.py",
     """        r["T8.5.1_nose_half_width"] = _entry(REGS["T8.5.1_nose_half_width"] - max(ys), f"{max(ys):.2f} mm out")
""",
     """        r["T8.5.1_nose_half_width"] = _entry(REGS["T8.5.1_nose_half_width"] - max(ys), f"{max(ys):.2f} mm out")
    if (assembly.get("info") or {}).get("nose_split") and "nose" in parts:      # evo_patches/nose_material.py
        r.update(nose_joint_checks(assembly, mm(body), mm(parts["nose"]), mm(parts["fwing"]) if "fwing" in parts else None))
"""),
    # ---------------------------------------------------------------- the finalists' manufacturing files
    ("car/parts/joints.py",
     """    keep = trimesh.creation.box(bounds=[[x_ref_a_mm / 1e3, -0.01, -0.01], [1.0, 0.2, 0.2]])
    milled = _bool("intersection", [body_half, keep])          # aft of Ref A
""",
     """    # evo_patches/nose_material.py: a split nose's foam is milled with the body (from x_cut), its print already
    # carries the tenon and the body the pocket: the print is kept as it is, glued in (not cut again here)
    split = (assembly.get("info") or {}).get("nose_split")
    x_keep = split["x_cut_mm"] if split else x_ref_a_mm
    keep = trimesh.creation.box(bounds=[[x_keep / 1e3, -0.01, -0.01], [1.0, 0.2, 0.2]])
    milled = _bool("intersection", [body_half, keep])          # aft of Ref A (or of a foam nose's face)
"""),
    ("car/parts/joints.py",
     """    cone = parts.get("nose")                  # the CFD cone (Part 4 nose.py), if any
    np_ = nose_part(body_half, x_ref_a_mm)
""",
     """    split_nose = parts.pop("nose", None) if split else None      # (evo_patches/nose_material.py)
    cone = parts.get("nose")                  # the CFD cone (Part 4 nose.py), if any
    np_ = nose_part(body_half, x_ref_a_mm) if not split else None
"""),
    ("car/parts/joints.py",
     """    mb = machined_body(milled, joints)
""",
     """    if split_nose is not None:    # evo_patches/nose_material.py: glued on its face and by its tenon in the pocket
        parts["nose"] = split_nose
        t = split.get("tenon")
        report["parts"]["nose"] = ({"joint": "tenon glued in a pocket from the bottom, face glued at x_cut (nose_split)",
                                    "pocket_cm3_per_side": (t["x1"] - t["x0"]) * t["hw"] * (t["z1"] - t["z0"]) / 1e3}
                                   if t else {"joint": "face glued at x_cut, no tenon fits (nose_split)"})
        report["mass_delta_kg"]["nose"] = 0.0
    mb = machined_body(milled, joints)
"""),
    # ---------------------------------------------------------------- the race model's roll-up: the breakdown
    ("search/layout/layout_car.py",
     """    tot = sum(i[1] for i in items)
    return dict(car_g=car_g, foam_g=foam, lead_g=lead,""",
     """    tot = sum(i[1] for i in items)
    # evo_patches/nose_material.py: the nose's foam is in "foam" (the body's mesh), its print in "printed"
    _ns = (asm.get("info") or {}).get("nose_split") or {}
    nose_split = {k: _ns.get(k) for k in ("x_cut_mm", "foam_len_mm", "decided_by", "foam_g", "printed_g")} if _ns else None
    return dict(nose_split=nose_split, car_g=car_g, foam_g=foam, lead_g=lead,"""),
]

MODULE = "car/body/nose_split.py"


def apply(evo):
    """Patch the checkout `evo` in place. Returns a list of (file, status)."""
    evo = Path(evo)
    out = []
    mod = evo / MODULE
    if mod.exists() and mod.read_text() == NOSE_SPLIT_SRC:
        out.append((MODULE, "already written"))
    else:
        assert not mod.exists() or "evo_patches/nose_material.py" in mod.read_text(), \
            f"{MODULE} exists and is not this patch's"
        mod.write_text(NOSE_SPLIT_SRC)
        out.append((MODULE, "written"))
    texts = {}
    for f, old, new in PATCHES:
        p = evo / f
        assert p.exists(), f"nose_material: {p} not found (is {evo} an evolution checkout?)"
        s = texts.get(f, p.read_text())
        if s.count(new) == 1:
            out.append((f, "already applied"))
        else:
            n = s.count(old)
            assert n == 1, f"nose_material: anchor found {n} times in {f} (builder changed?):\n{old[:200]}"
            s = s.replace(old, new, 1)
            out.append((f, "applied"))
        texts[f] = s
    for f, s in texts.items():
        if (evo / f).read_text() != s:
            (evo / f).write_text(s)
    return out


if __name__ == "__main__":
    import sys
    for f, st in apply(sys.argv[1]):
        print(f"nose_material: {f}: {st}")
