"""front_wing_engine.py -- the searches' front wing: one piece or two, with the supports it designs itself.

  pieces      1: one wing across the centre line; 2: a left and a right wing (T8.6.1: two segments of at
              least 25 mm each where another part, the nose or a support, passes through the span).
  element     a flat NACA 0015 (user 2026-11-22: low drag, low lift), chord 15.3-25 mm (T8.6.2),
              STRAIGHT (user 2026-11-23: no sweep, no dihedral, no tilt), so every vertical cut is the
              true section (T8.6.2 and T8.6.3 are measured on such cuts).
  supports    lofts (wing_loft.py, like the halo's loft): each grows from a patch of the wing -- anywhere on
              the span the rules let a support stand, over part or all of the chord -- to the nearest point
              of the nose, above or below, its outline morphing from the patch's to its foot's; printed with
              the nose (nose.stl). A wing passing through the nose is held by the nose (T8.6.1's figures).
              Numbers: width, where on the span, how much of the chord, the foot's scale, the waist.

Legal by construction where that is cheap, measured always (gates()):
  T8.5.1  supports inside |y| <= 15, z <= 25, ahead of Ref A.
  T8.5.2  the wing ahead of Ref A; no higher than 20 mm beyond 15 mm out (it is level: its height).
  T8.6.1  each element's span from the nose or support that cuts it (or from the centre) >= 25 (the
          tip is pushed out until it is), one piece >= 50.
  T8.6.2/3  chord 15-25, thickness 2-6 (the section: 15 % of a 15.3-25 mm chord).
  T8.7    5 mm of clear air to the track, the wheels, the body and the nose, except within 5 mm of the
          wing's joints with its supports (the fillet the rule allows): a wing under the nose is lowered
          to keep it, and refused if that would take it under the track's 5 mm.
  T8.8    nothing ahead of the wing outside the nose (the nose is all there is ahead of it).
"""
from __future__ import annotations

import dataclasses
import math
from dataclasses import dataclass, field

import numpy as np

import component_contract as cc

T_FRAC = 0.15                  # NACA 0015
CLEAR_MM = 5.0                 # T8.7
SUPPORT_BAND_MM = 15.0         # T8.5.1
HIGH_OUTBOARD_MM = 20.0        # T8.5.2
SEGMENT_MIN_MM = 25.0          # T8.6.1
OVERHANG_MM = 40.0             # T8.2: the car's front no more than this ahead of Ref A
MARGIN_MM = 0.3                # kept off every limit
TRACK_CLEAR_MM = 5.0            # T8.7: 5 mm over the track
EMBED_MM = 1.0                 # a pylon runs this far into the nose and into the wing
N_SPAN = 24                    # sections along an element


@dataclass
class Support:
    """A split wing's side struts (wing_strut.py): each leaves its half through the inner end -- its first
    section the half's own aerofoil -- and sweeps, as the search shapes it, into the nose. Shape numbers only."""
    aim_out_deg: float = 0.0       # from straight in towards the nose: leaning down (+, the nose's lower part) or up
    aim_fwd_deg: float = 0.0       # leaning forward (+) or aft (-)
    leave: float = 0.35
    arrive: float = 0.35
    chord_mm: float = 10.0
    thick_mm: float = 3.0
    flare: float = 1.4
    round: float = 2.5
    mid_x: float = 0.0          # a free point its path is drawn towards half way (share of its length) ...
    mid_y: float = 0.0
    mid_z: float = 0.0
    swell: float = 1.0          # its section half way, times its own


@dataclass
class FrontWingDesign:
    """STRAIGHT (user 2026-11-23: no tilt, no angling): a level NACA 0015, square to the car."""
    pieces: int = 1                # 1 or 2
    chord_mm: float = 15.3         # T8.6.2 minimum + 0.3 for the print
    te_gap_mm: float = 6.0         # trailing edge this far ahead of Ref A (T8.7: >= 5 to the body)
    z_mm: float = 8.0              # chord line height
    tip_y_mm: float = 36.5         # tip: just outboard of the front wheel by default
    root_y_mm: float = 8.0         # two pieces: where each half starts (inside the support band)
    support: Support = field(default_factory=Support)


def from_dict(d: dict) -> FrontWingDesign:
    keep = {f.name for f in dataclasses.fields(FrontWingDesign)}
    d = dict(d)
    sp = d.pop("support", None) or {}
    return FrontWingDesign(**{k: v for k, v in d.items() if k in keep and k not in ("support", "supports")},
                           support=Support(**sp) if isinstance(sp, dict) else sp)


def params_dict(d: FrontWingDesign) -> dict:
    return dataclasses.asdict(d)


def naca00(chord, t_frac=T_FRAC, n=48):
    """Closed symmetric NACA 4-digit section (x from 0 at the LE to chord, z), anticlockwise, with the
    trailing edge closed to the printer's minimum."""
    beta = np.linspace(0.0, math.pi, n)
    xc = 0.5 * (1 - np.cos(beta))
    yt = 5 * t_frac * (0.2969 * np.sqrt(xc) - 0.1260 * xc - 0.3516 * xc ** 2 + 0.2843 * xc ** 3 - 0.1036 * xc ** 4)
    te = cc.BUILD_WALL_MM / 2 / chord
    yt = np.maximum(yt, te * xc)                  # a printable trailing edge, never thinner than the wall
    upper = np.c_[xc[::-1], yt[::-1]]
    lower = np.c_[xc[1:-1], -yt[1:-1]]
    return np.vstack([upper, lower]) * chord     # TE (top) -> LE -> TE (bottom): anticlockwise in (x, z)


def _loft(stations):
    import support_engine as se
    return se.loft_sections(stations)


def element(d: FrontWingDesign, y0, y1, x_le_root, z_root):
    """One straight wing element from y0 to y1 (y0 < y1) (manifold, mm)."""
    sec = naca00(d.chord_mm)
    return _loft([(float(y), np.c_[x_le_root + sec[:, 0], z_root + sec[:, 1]]) for y in (y0, y1)])


def _loft_z(stations):
    """Closed solid through (x, y) sections stacked along z: [(z, (n, 2) array)] (manifold, mm)."""
    import trimesh
    import support_engine as se
    rings = [np.c_[P[:, 0], P[:, 1], np.full(len(P), z)] for z, P in stations]
    n, N = len(rings[0]), len(rings)
    V = np.vstack(rings + [np.mean(rings[0], 0)[None], np.mean(rings[-1], 0)[None]])
    F = []
    for i in range(N - 1):
        for j in range(n):
            a, b = i * n + j, i * n + (j + 1) % n
            c_, d_ = (i + 1) * n + j, (i + 1) * n + (j + 1) % n
            F += [(a, b, d_), (a, d_, c_)]
    c0, c1 = N * n, N * n + 1
    for j in range(n):
        F.append((c0, (j + 1) % n, j))
        F.append((c1, (N - 1) * n + j, (N - 1) * n + (j + 1) % n))
    m = trimesh.Trimesh(V, np.array(F), process=False)
    if m.volume < 0:
        m.invert()
    return se._to_manifold(m)


def _nose_hits(nose_tm, x, y, z0):
    """(z of the nose's surface above z0, below z0) along the vertical through (x, y), or None."""
    o = np.array([[x, y, z0], [x, y, z0]])
    dirs = np.array([[0.0, 0.0, 1.0], [0.0, 0.0, -1.0]])
    loc, ri, _ = nose_tm.ray.intersects_location(o, dirs, multiple_hits=True)
    loc = np.asarray(loc).reshape(-1, 3); ri = np.asarray(ri).reshape(-1)
    up = loc[ri == 0, 2]
    dn = loc[ri == 1, 2]
    return (float(up.min()) if len(up) else None), (float(dn.max()) if len(dn) else None)


def _ray_hit(m, o, d):
    """The first surface of m from o along d (mm), or None."""
    loc, _, _ = m.ray.intersects_location(o[None], (d / np.linalg.norm(d))[None], multiple_hits=True)
    if not len(loc):
        return None
    return loc[np.argmin(np.linalg.norm(loc - o, axis=1))]


def build(d: FrontWingDesign, nose_half_m, x_front_mm, wheel, body_half_m=None):
    """The wing and its pylons, right halves (trimesh, metres), for this nose. wheel = (x_axle, y_inner,
    R) of the front wheel. Returns {"wing": trimesh, "pylons": trimesh or None, "gates": {...}, "info": {...}}.
    A wing the nose does not hold, closer than 5 mm to it, is lowered by what it lacks and built again
    (down to the track's own 5 mm)."""
    r = _build_once(d, nose_half_m, x_front_mm, wheel, body_half_m)
    if d.pieces == 1 and not r["info"].get("held_by_nose") and not r["info"].get("struts"):
        # one wing across the centre line is held only by the nose it runs through (T1.7): moved to the
        # nose's own height at its chord, then built again
        nose = nose_half_m.copy(); nose.apply_scale(1e3)
        xm = r["info"]["x_le_root_mm"] + 0.5 * d.chord_mm
        sx = nose.section(plane_origin=[xm, 0, 0], plane_normal=[1.0, 0, 0])
        if sx is not None:
            zc = 0.5 * float(sx.vertices[:, 2].min() + sx.vertices[:, 2].max())
            r2 = _build_once(dataclasses.replace(d, z_mm=zc), nose_half_m, x_front_mm, wheel, body_half_m)
            r2["info"]["moved_to_the_nose_mm"] = round(zc - r["info"]["z_mm"], 2)
            r = r2
    for _ in range(3):
        short = -r["gates"].get("T8.7_fw_nose", 0.0)
        if short <= 0 or r["info"].get("held_by_nose"):
            break
        d = dataclasses.replace(d, z_mm=r["info"]["z_mm"] - short - MARGIN_MM)
        r2 = _build_once(d, nose_half_m, x_front_mm, wheel, body_half_m)
        r2["info"]["lowered_for_clear_air_mm"] = r["info"].get("lowered_for_clear_air_mm", 0.0) + short + MARGIN_MM
        r = r2
    return r


def _build_once(d: FrontWingDesign, nose_half_m, x_front_mm, wheel, body_half_m=None):
    import support_engine as se
    import trimesh
    ref_a = x_front_mm - 16.0
    nose = nose_half_m.copy(); nose.apply_scale(1e3)
    info = dict(design=params_dict(d))
    t = T_FRAC * d.chord_mm
    # the whole wing between the car's front limit (T8.2: 40 mm ahead of Ref A) and te_gap ahead of Ref A
    gap = float(np.clip(d.te_gap_mm, CLEAR_MM, OVERHANG_MM - MARGIN_MM - d.chord_mm))
    x_le = ref_a - gap - d.chord_mm
    info.update(x_le_root_mm=x_le, front_x_mm=x_le, te_gap_used_mm=gap)
    # level, under 20 mm (T8.5.2 outboard of 15 mm) and 5 mm off the track
    z = float(np.clip(d.z_mm, CLEAR_MM + t / 2 + MARGIN_MM, HIGH_OUTBOARD_MM - t / 2 - MARGIN_MM))
    d = dataclasses.replace(d, z_mm=z)
    info.update(z_mm=z)
    # where the nose cuts the span, at the wing's height (its half-width there)
    xs = np.linspace(x_le, x_le + d.chord_mm, 7)
    cut, cut_min = 0.0, np.inf
    for x in xs:
        s = nose.section(plane_origin=[x, 0, 0], plane_normal=[1, 0, 0])
        near = s.vertices[np.abs(s.vertices[:, 2] - z) < t / 2 + 0.5] if s is not None else np.zeros((0, 3))
        w = float(near[:, 1].max()) if len(near) else 0.0
        cut, cut_min = max(cut, w), min(cut_min, w)
    info["nose_cut_y_mm"] = cut                           # widest along the chord: where the span is cut
    info["nose_cut_min_y_mm"] = cut_min                   # narrowest: a root inside it is held all along
    # the nose's half-width anywhere over the wing's length (where a pylon can reach it from above or below)
    hw_any = 0.0
    for x in xs:
        s_ = nose.section(plane_origin=[x, 0, 0], plane_normal=[1, 0, 0])
        if s_ is not None:
            hw_any = max(hw_any, float(s_.vertices[:, 1].max()))
    info["nose_half_width_mm"] = hw_any
    if d.pieces == 2:
        # each half starts where the search says, inside T8.5.1's band: in the nose (the nose holds it) or
        # off it (its inner end is lofted into the nose)
        y0 = float(np.clip(d.root_y_mm, 0.5, SUPPORT_BAND_MM - 1.0))
    else:
        y0 = -0.5                                         # across the centre line (mirrored: one piece)
    seg_from = max(cut, max(y0, 0.0))
    tip = max(d.tip_y_mm, seg_from + SEGMENT_MIN_MM + MARGIN_MM + 0.3,      # (+ the side support's 0.3 mm into it)
              (SEGMENT_MIN_MM + MARGIN_MM) if d.pieces == 1 and cut == 0.0 else 0.0)
    info["tip_y_mm"] = tip
    el = element(d, y0, tip, x_le, z)
    el = el ^ se._box(-1e3, 1e3, 0.0, 1e3, -1e3, 1e3)     # right half
    wing = se._to_trimesh(el)
    # is the wing held by the nose (does it run into it)?
    inside = trimesh.boolean.intersection([wing, nose], engine="manifold") if nose.is_volume else None
    held = inside is not None and not inside.is_empty and abs(inside.volume) > 1.0
    info["held_by_nose"] = bool(held)
    # side struts (a split wing): each half's inner end swept into the nose, printed with it
    import wing_loft as wl
    import wing_strut as ws
    pylons, placed = [], []
    sp = d.support
    root_in_nose = d.pieces == 2 and nose.is_volume and bool(nose.contains(np.array([[x_le + 0.5 * d.chord_mm, y0 + 0.2, z]]))[0])
    if d.pieces == 2 and not root_in_nose:
        sec = naca00(d.chord_mm)
        face_xz = np.c_[x_le + sec[:, 0], z + sec[:, 1]]
        face_pt = np.array([x_le + 0.5 * d.chord_mm, y0, z])
        Mnose = se._to_manifold(nose)
        half = se._box(-1e3, 1e3, 0.0, 1e3, -1e3, 1e3)

        wb = wing.bounds

        def accept(M):
            # (owner 2026-11-28: no "strut must not run back into its wing's span" check on the FRONT wing; it
            # blocked every strut on split wings set high or on narrow noses. The referee still measures T8.6.1,
            # each element's span from what cuts it, on every built car, so a strut that shortens one is refused.)
            t = wl.to_trimesh(M ^ half)
            b = t.bounds
            return (b[1, 1] <= SUPPORT_BAND_MM - MARGIN_MM and b[1, 2] <= 25.0 - MARGIN_MM
                    and b[1, 0] <= ref_a - MARGIN_MM and b[0, 0] >= ref_a - OVERHANG_MM + MARGIN_MM
                    and b[0, 2] >= TRACK_CLEAR_MM)
        # "out" for the front: leaning from straight in (-y) towards down (-z)
        M, sinfo = ws.design(face_xz, face_pt, np.array([0.0, -1.0, 0.0]), ws.whole(nose), np.array([0.0, -1.0, 0.0]),
                             dict(dataclasses.asdict(sp)), reach_mm=30.0, accept=accept, wall_mm=cc.BUILD_WALL_MM)
        if M is not None:
            pylons.append(M)
            placed.append(dict(y_mm=round(y0, 2), t_mm=0.0, **sinfo))
        elif sinfo.get("in_target"):
            # its end is all but in the nose (within 0.5 mm): run on 1.5 mm into it, held by it
            y0 = max(y0 - 1.5, 0.5)
            wing = se._to_trimesh(element(d, y0, tip, x_le, z) ^ se._box(-1e3, 1e3, 0.0, 1e3, -1e3, 1e3))
            root_in_nose = True
            info["root_run_into_the_nose"] = True
        else:
            info["supports_not_placed"] = sinfo
    info["root_in_nose"] = root_in_nose
    held = held or root_in_nose
    if not held and not placed:
        # (owner 2026-11-28: supports not only struts from the root) a pylon out of the wing's underside, straight
        # down into the nose, anywhere along the span inside T8.5.1's 15 mm band
        import wing_loft as wl
        import wing_strut as ws
        sp_ = dict(dataclasses.asdict(d.support))
        w = max(float(sp_["thick_mm"]), cc.BUILD_WALL_MM)
        lo_y = max(y0, 0.0) + 0.5 * w + 0.5
        hi_y = min(tip, SUPPORT_BAND_MM) - 0.5 * w - MARGIN_MM - 0.3
        half_ = se._box(-1e3, 1e3, 0.0, 1e3, -1e3, 1e3)

        def accept_p(M):
            b = wl.to_trimesh(M ^ half_).bounds
            return (b[1, 1] <= SUPPORT_BAND_MM - MARGIN_MM and b[1, 2] <= 25.0 - MARGIN_MM
                    and b[1, 0] <= ref_a - MARGIN_MM and b[0, 0] >= ref_a - OVERHANG_MM + MARGIN_MM
                    and b[0, 2] >= TRACK_CLEAR_MM)
        tried = []
        for y_s in (np.linspace(lo_y, hi_y, 3) if hi_y > lo_y else []):
            fxz = np.array([[x_le + 0.2 * d.chord_mm, y_s - 0.5 * w], [x_le + 0.8 * d.chord_mm, y_s - 0.5 * w],
                            [x_le + 0.8 * d.chord_mm, y_s + 0.5 * w], [x_le + 0.2 * d.chord_mm, y_s + 0.5 * w]])
            M, sinfo = ws.design(fxz, np.array([x_le + 0.5 * d.chord_mm, y_s, z]), np.array([0.0, 0.0, -1.0]), ws.whole(nose),
                                 np.array([0.0, 0.0, -1.0]), sp_, reach_mm=30.0, accept=accept_p, wall_mm=cc.BUILD_WALL_MM,
                                 face_v=np.array([0.0, 1.0, 0.0]))
            tried.append(round(float(y_s), 2))
            if M is not None:
                pylons.append(M)
                placed.append(dict(y_mm=round(float(y_s), 2), t_mm=0.0, pylon="underside", **sinfo))
                info.pop("supports_not_placed", None)
                break
        info["underside_pylons_tried_y_mm"] = tried
    info["struts"] = placed
    info["pylon_y_placed_mm"] = [p["y_mm"] for p in placed]
    pyl = None
    if pylons:
        import manifold3d as mf
        Pm = mf.Manifold.batch_boolean(pylons, mf.OpType.Add) ^ se._box(-1e3, 1e3, 0.0, 1e3, -1e3, 1e3)
        pyl = wl.to_trimesh(Pm)
    info["supported"] = bool(held or placed)
    gates = measure(d, wing, pyl, nose, ref_a, wheel, body_half_m, info)
    wing.apply_scale(1e-3)
    if pyl is not None:
        pyl.apply_scale(1e-3)
    return {"wing": wing, "pylons": pyl, "gates": gates, "info": info}


def _chunked(f, pts, n=256):
    """f over the points a few hundred at a time: trimesh's closest point pairs every point with every
    triangle near it, which on a fine nose ran to 35 million pairs (a 10 GB build killed, gallery G01)."""
    return np.concatenate([np.asarray(f(pts[i:i + n]), float) for i in range(0, len(pts), n)]) if len(pts) else np.zeros(0)


def measure(d, wing, pyl, nose, ref_a, wheel, body_half_m, info) -> dict:
    """Gates (+ = passes, mm) measured on the built wing and pylons (all mm, right halves)."""
    import trimesh
    g = {}
    b = wing.bounds
    g["fw_supported"] = 0.0 if info["supported"] else -1.0
    g["T8.5.2_fw_ahead_of_ref_A"] = ref_a - float(b[1, 0])
    g["T8.2_fw_inside_overhang"] = OVERHANG_MM - (ref_a - float(b[0, 0]))
    out = wing.vertices[wing.vertices[:, 1] > SUPPORT_BAND_MM]
    if len(out):
        g["T8.5.2_fw_height_outboard"] = HIGH_OUTBOARD_MM - float(out[:, 2].max())
    seg_from = max(info["nose_cut_y_mm"], max(float(b[0, 1]), 0.0) if d.pieces == 2 else 0.0)
    if pyl is not None:
        seg_from = max(seg_from, float(pyl.bounds[1, 1]))
    seg = float(b[1, 1]) - seg_from
    g["T8.6.1_fw_segment"] = seg - SEGMENT_MIN_MM if (seg_from > 0 or d.pieces == 2) else 2 * float(b[1, 1]) - 2 * SEGMENT_MIN_MM
    g["T8.6.2_fw_chord_min"] = d.chord_mm - 15.0
    g["T8.6.2_fw_chord_max"] = 25.0 - d.chord_mm
    g["T8.6.3_fw_thickness_min"] = T_FRAC * d.chord_mm - 2.0
    g["T8.6.3_fw_thickness_max"] = 6.0 - T_FRAC * d.chord_mm
    if pyl is not None:
        pb = pyl.bounds
        g["T8.5.1_fw_support_half_width"] = SUPPORT_BAND_MM - float(pb[1, 1])
        g["T8.5.1_fw_support_height"] = 25.0 - float(pb[1, 2])
        g["T8.5.1_fw_support_ahead_of_ref_A"] = ref_a - float(pb[1, 0])
    # T8.7: clear air from the wing's surface to everything else, the joints excepted
    pts, _ = trimesh.sample.sample_surface(wing, 4000, seed=5)
    # the joints: where the wing meets the nose (it runs into it) and its pylons. Beside each, the corner
    # the 5 mm ball on its rod cannot enter is the fillet the rule allows: a spanwise band 5 mm wide
    # outboard of the nose where it cuts the wing, and 5 mm either side of each pylon (as the referee).
    near_joint = np.zeros(len(pts), bool)
    # as the referee (legality: points within 5 mm sideways of the nose's widest are its joint's fillet)
    near_joint |= pts[:, 1] < float(nose.bounds[1, 1]) + CLEAR_MM
    if info.get("held_by_nose"):
        near_joint |= pts[:, 1] < max(info["nose_cut_y_mm"], info.get("nose_half_width_mm", 0.0)) + CLEAR_MM
    for p in info.get("struts", []):
        near_joint |= np.abs(pts[:, 1] - p["y_mm"]) < p["t_mm"] / 2 + CLEAR_MM
    free = pts[~near_joint]
    g["T8.7_fw_track"] = float(wing.bounds[0, 2]) - CLEAR_MM
    xa, y_in, R = wheel
    g["T8.7_fw_wheel"] = (xa - R) - float(b[1, 0]) - CLEAR_MM
    if len(free):
        if nose.is_volume:
            sd = _chunked(lambda q: trimesh.proximity.signed_distance(nose, q), free)
            dn = np.where(sd > 0, 0.0, np.abs(sd))
        else:
            dn = np.full(len(free), 99.0)
        g["T8.7_fw_nose"] = float(dn.min()) - CLEAR_MM
        info["fw_nose_clear_mm"] = float(dn.min())
    if body_half_m is not None:
        body = body_half_m.copy(); body.apply_scale(1e3)
        if len(free):
            db = _chunked(lambda q: trimesh.proximity.closest_point(body, q)[1], free)
            g["T8.7_fw_body"] = float(db.min()) - CLEAR_MM
    return g
