"""The rear wing as the regulations allow it, with zero angle of attack, and the strict reading of T9.5.1
(carhybrid/PLAN_WINGS_FLATNESS.md section 2; docs/rules/README.md).

    python evo_patches/apply.py <evolution checkout>      (applied with the other modules, after legality_fixes)

Owner, 2026-11-28 evening: "you can have the rear wing do the things according to the regs, like have it bend ...
it's just that the wing has to have no AoA"; "it allows a cambered section with level chord line".

What changes in the builder (car/parts/rear_wing_engine.py, search/layout/layout_params.json, layout_genes.py):
W1 The wing's shape along its span, every section's chord line level (AoA 0 everywhere, so no twist). For the
   normalised span position s = |y| / half_span:
     height   z(s)  = z_root + rise_mm * s^rise_p             (gull, droop, dihedral; T9.5.4: |rise| <= 13.5 mm)
     LE       x(s)  = x_le   + sweep_mm * s^sweep_p           (swept forward or back)
     chord    c(s)  = chord + (tip_chord - chord) * s         (T9.5.2: every section 15.3-24.7 mm)
     thick    t(s)  = t + (tip_t - t) * s                     (T9.5.3: every section 2.3-5.7 mm)
     section  NACA 4-digit, camber of either sign (+ lifts, - presses down), at camber_pos; chord line level.
   With every new number at its default (rise 0, sweep 0, tip ratios 1, camber 0) the wing is today's.
W2 Placement over the whole span, not the root: every section between Ref B and the body's back (and T9.4.2), the
   T7.9 zones at every section's leading edge, the L1 chord cut on the furthest-aft section.
W3 Supports from the wing as built: the tip struts leave through the TIP section (its own chord, thickness, camber,
   height and LE); a pylon's footprint sits under the section at its own y.
W4 T9.5.1, strict (docs/rules/README.md): a part joined to the rear wing breaks its span there; the wing needs one
   unbroken stretch of >= 50 mm. Tip struts (they run 0.3 mm into the tip) leave 2 (half_span - 0.35); a pylon at y
   leaves 2 (y - w/2) in the middle: a pylon that would leave less than 50.3 mm is refused (the next root, then L5's
   move forward). The half span's floor rises to 25.6 mm (2 x (25.6 - 0.35) = 50.5).
W5 The engine's own gates measured per section (chord, thickness, top-surface deviation T9.5.4, strict span), not on
   the bounding box (which a bent or swept wing would make meaningless).
W7 The rise is cut when the wing, placed at the lowest height T9.6 allows, would top T9.4.3's 65 mm (wing1:
   wg_L08 and wg_L12, rise 11-13 mm, broke T9.4.3/T3.5); the height is then found again. Recorded as
   info["rise_cut_for_T9.4.3_mm"].
W6 The parameters: rw_camber (+-0.06) is released (it was fixed at 0); new shape parameters rw_camber_pos, rw_rise,
   rw_rise_p, rw_sweep, rw_sweep_p, rw_tip_chord_ratio, rw_tip_t_ratio; rw_aoa stays fixed at 0.
L7 (in car/build/legality.py, the referee) T9.5.1 strict: the rear span is the longest stretch of the wing that no
   other part touches or passes through -- its own supports included -- measured on the full car (the stretch
   through the centre counts twice the half). Replaces the reading "its own pylon under it does not break it".
"""
from __future__ import annotations

import json
from pathlib import Path

RW = "car/parts/rear_wing_engine.py"
LG = "search/layout/layout_genes.py"
LP = "search/layout/layout_params.json"
REF = "car/build/legality.py"

EDITS = [
    dict(id="W1a_design_fields", file=RW,
         old='''    lift_mm: float = 0.0        # above the lowest height T9.6 and T9.7 allow
    support: Support = field(default_factory=Support)''',
         new='''    lift_mm: float = 0.0        # above the lowest height T9.6 and T9.7 allow
    # evo_patches/rear_wing_shape.py W1: the shape along the span (every section's chord line level: AoA 0)
    camber: float = 0.0         # NACA camber, fraction of chord: + lifts, - presses down
    camber_pos: float = 0.4     # where the camber peaks, share of the chord
    rise_mm: float = 0.0        # the tip's height above the root (- droops): T9.5.4
    rise_p: float = 2.0         # how the height runs out to the tip: s^rise_p
    sweep_mm: float = 0.0       # the tip's leading edge behind the root's (- forward)
    sweep_p: float = 1.0
    tip_chord_ratio: float = 1.0
    tip_t_ratio: float = 1.0
    support: Support = field(default_factory=Support)'''),
    dict(id="W1b_span_shape", file=RW,
         old='''def _wing(d, x_le, z):''',
         new='''def span_shape(d, y):
    """evo_patches W1: (dx, dz, chord, t) of the section at span position y (mm, |y| used): the leading edge's and
    the chord line's offsets from the root's, the local chord and thickness."""
    hs = max(float(d.half_span_mm), 1e-6)
    s = np.clip(np.abs(np.asarray(y, float)) / hs, 0.0, 1.0)
    dx = float(d.sweep_mm) * s ** max(float(d.sweep_p), 0.5)
    dz = float(d.rise_mm) * s ** max(float(d.rise_p), 0.5)
    c = d.chord_mm + (tip_chord(d) - d.chord_mm) * s
    t = d.t_mm + (tip_t(d) - d.t_mm) * s
    return dx, dz, c, t


def tip_chord(d):
    return float(np.clip(d.chord_mm * float(d.tip_chord_ratio), cc.REAR_CHORD_MIN + MARGIN_MM, cc.REAR_CHORD_MAX - MARGIN_MM))


def tip_t(d):
    return float(np.clip(d.t_mm * float(d.tip_t_ratio), cc.WING_THICK_MIN + MARGIN_MM, cc.WING_THICK_MAX - MARGIN_MM))


def section_at(d, y):
    """evo_patches W1: the section at span position y, (x from its LE, z about its chord line): NACA 4-digit with the
    local chord and thickness, camber of either sign, chord line level."""
    import wings as wg
    _dx, _dz, c, t = span_shape(d, y)
    c, t = float(c), float(t)
    m = float(getattr(d, "camber", 0.0))
    tf = t / c
    P = None
    for _ in range(3):
        P = wg.naca4(c, tf, abs(m), float(getattr(d, "camber_pos", 0.4)), te_mm=cc.BUILD_WALL_MM)
        th = _thickness_of(P)
        if th <= 0 or abs(th - t) < 1e-3:
            break
        tf *= t / th
    if m < 0:
        P = P * np.array([1.0, -1.0])
        P = P[::-1]                                       # keep the loop's orientation
    return P


def _thickness_of(P):
    """The section's thickness: its largest vertical extent at one x (its upper minus lower surface)."""
    xs = np.linspace(P[:, 0].min(), P[:, 0].max(), 60)[1:-1]
    n = len(P) // 2
    up, lo = P[: n + 1], P[n:]
    zu = np.interp(xs, up[:, 0], up[:, 1]) if up[0, 0] < up[-1, 0] else np.interp(xs, up[::-1, 0], up[::-1, 1])
    zl = np.interp(xs, lo[:, 0], lo[:, 1]) if lo[0, 0] < lo[-1, 0] else np.interp(xs, lo[::-1, 0], lo[::-1, 1])
    return float(np.max(np.abs(zu - zl)))


def stations(d, n=24):
    return np.linspace(-0.5, float(d.half_span_mm), max(N_SPAN, n))


def _wing(d, x_le, z):'''),
    dict(id="W1c_wing_loft", file=RW,
         old='''    sec = section(d)
    st = [(float(y), np.c_[x_le + sec[:, 0], z + sec[:, 1]]) for y in np.linspace(-0.5, d.half_span_mm, N_SPAN)]
    return se.loft_sections(st) ^ se._box(-1e3, 1e3, 0.0, 1e3, -1e3, 1e3)''',
         new='''    st = []
    for y in stations(d):                                 # evo_patches W1: each section its own
        dx, dz, _c, _t = span_shape(d, y)
        sec = section_at(d, y)
        st.append((float(y), np.c_[x_le + float(dx) + sec[:, 0], z + float(dz) + sec[:, 1]]))
    return se.loft_sections(st) ^ se._box(-1e3, 1e3, 0.0, 1e3, -1e3, 1e3)'''),
    dict(id="W2a_place_whole_span", file=RW,
         old='''    room = (x_aft_max - ref_b) - MARGIN_MM - chord
    x_le = float(np.clip(ref_b + d.x_le_mm, ref_b + MARGIN_MM, ref_b + MARGIN_MM + max(room, 0.0)))''',
         new='''    # evo_patches W1/W2: the shape numbers clipped to the rules, and the whole span placed (not the root alone)
    hs = float(max(hs, cc.REAR_SPAN_MIN / 2 + 0.6))                     # W4: 2 x (hs - 0.35) >= 50.5
    d = dataclasses.replace(d, half_span_mm=hs, chord_mm=chord, t_mm=t,
                            camber=float(np.clip(getattr(d, "camber", 0.0), -0.06, 0.06)),
                            camber_pos=float(np.clip(getattr(d, "camber_pos", 0.4), 0.2, 0.6)),
                            rise_mm=float(np.clip(getattr(d, "rise_mm", 0.0), -13.5, 13.5)),
                            sweep_mm=float(np.clip(getattr(d, "sweep_mm", 0.0), -15.0, 15.0)))

    def _ext(dd):
        ys = stations(dd)
        sh = [span_shape(dd, y) for y in ys]
        return min(float(a) for a, _, _, _ in sh), max(float(a) + float(c) for a, _, c, _ in sh)
    for _ in range(30):                                   # a sweep the room cannot take is scaled down
        lo_dx, hi_ext = _ext(d)
        room = (x_aft_max - ref_b) - MARGIN_MM - (hi_ext - lo_dx)
        if room >= 0 or abs(d.sweep_mm) < 1e-6:
            break
        d = dataclasses.replace(d, sweep_mm=0.8 * d.sweep_mm if abs(d.sweep_mm) > 0.2 else 0.0)
        info["sweep_cut_for_room"] = True
    lo_dx, hi_ext = _ext(d)
    room = (x_aft_max - ref_b) - MARGIN_MM - (hi_ext - lo_dx)
    x_le = float(np.clip(ref_b + d.x_le_mm, ref_b + MARGIN_MM - lo_dx, ref_b + MARGIN_MM - lo_dx + max(room, 0.0)))'''),
    dict(id="W2b_t79_every_section", file=RW,
         old='''        zv = ez.t79_zones(np.full(Y_.size, x_le - ZONE_MARGIN_MM), Y_.ravel(), Z_.ravel(), W_mm, x_front_mm)''',
         new='''        zv = ez.t79_zones(x_le + np.asarray(span_shape(d, Y_.ravel())[0], float) - ZONE_MARGIN_MM, Y_.ravel(), Z_.ravel(),
                          W_mm, x_front_mm)                    # evo_patches W2: every section's LE'''),
    dict(id="W2c_chord_cut_every_section", file=RW,
         old='''        _sx = float(section(dataclasses.replace(d, half_span_mm=hs, chord_mm=chord, t_mm=t))[:, 0].max())
        over = x_le + _sx - x_aft_max''',
         new='''        _dd = dataclasses.replace(d, half_span_mm=hs, chord_mm=chord, t_mm=t)
        over = x_le + _ext(_dd)[1] - x_aft_max                # evo_patches W2: the furthest-aft section'''),
    dict(id="W2d_keep_shape", file=RW,
         old='''    d = dataclasses.replace(d, half_span_mm=hs, chord_mm=chord, t_mm=t)
    info.update(half_span_mm=hs, chord_mm=chord, t_mm=t, x_le_mm=x_le)''',
         new='''    d = dataclasses.replace(d, half_span_mm=hs, chord_mm=chord, t_mm=t)
    info.update(half_span_mm=hs, chord_mm=chord, t_mm=t, x_le_mm=x_le,
                shape={k: getattr(d, k) for k in ("camber", "camber_pos", "rise_mm", "rise_p", "sweep_mm", "sweep_p")},
                tip_chord_mm=tip_chord(d), tip_t_mm=tip_t(d))'''),
    dict(id="W3a_tip_face", file=RW,
         old='''    sec = section(d)
    face_xz = np.c_[x_le + sec[:, 0], z + sec[:, 1]]
    face_pt = np.array([x_le + 0.5 * chord, hs, z])''',
         new='''    # evo_patches W3: the tip section as built (its own chord, thickness, camber, height and LE)
    _tdx, _tdz, _tc, _tt = span_shape(d, hs)
    sec = section_at(d, hs)
    face_xz = np.c_[x_le + float(_tdx) + sec[:, 0], z + float(_tdz) + sec[:, 1]]
    face_pt = np.array([x_le + float(_tdx) + 0.5 * float(_tc), hs, z + float(_tdz)])'''),
    dict(id="W3b_tip_attempt", file=RW,
         old='''            fxz, fpt, n_out, fv, around, y_j = face_xz, np.array([x_le + 0.5 * chord, hs, z]), np.array([0.0, 1.0, 0.0]), None, True, hs''',
         new='''            fxz, fpt, n_out, fv, around, y_j = face_xz, face_pt, np.array([0.0, 1.0, 0.0]), None, True, hs'''),
    dict(id="W3c_W4_pylon", file=RW,
         old='''            y_s = float(np.clip(root * hs, CENTRE_CLEAR_MM + 0.5 * w + 0.5, hs - 0.5 * w - 0.5))
            xa, xb = x_le + 0.2 * chord, x_le + 0.8 * chord
            fxz = np.array([[xa, y_s - 0.5 * w], [xb, y_s - 0.5 * w], [xb, y_s + 0.5 * w], [xa, y_s + 0.5 * w]])
            fpt, n_out, fv, around, y_j = np.array([x_le + 0.5 * chord, y_s, z]), np.array([0.0, 0.0, -1.0]), np.array([0.0, 1.0, 0.0]), False, y_s''',
         new='''            y_s = float(np.clip(root * hs, CENTRE_CLEAR_MM + 0.5 * w + 0.5, hs - 0.5 * w - 0.5))
            # evo_patches W4 (T9.5.1 strict): a pylon joined to the wing breaks its span there; the stretch between
            # the two pylons must stay >= 50 mm unbroken, or this root is refused
            if 2.0 * (y_s - 0.5 * w) < cc.REAR_SPAN_MIN + 0.3:
                return None, {"refused": "T9.5.1 strict: the pylons would leave %.2f mm unbroken" % (2.0 * (y_s - 0.5 * w))}, y_s
            _pdx, _pdz, _pc, _pt = span_shape(d, y_s)          # evo_patches W3: the section above the pylon
            xa, xb = x_le + float(_pdx) + 0.2 * float(_pc), x_le + float(_pdx) + 0.8 * float(_pc)
            fxz = np.array([[xa, y_s - 0.5 * w], [xb, y_s - 0.5 * w], [xb, y_s + 0.5 * w], [xa, y_s + 0.5 * w]])
            fpt, n_out, fv, around, y_j = np.array([x_le + float(_pdx) + 0.5 * float(_pc), y_s, z + float(_pdz)]), np.array([0.0, 0.0, -1.0]), np.array([0.0, 1.0, 0.0]), False, y_s'''),
    dict(id="W5_gates_per_section", file=RW,
         old='''    g["T9.5.1_rw_span"] = 2 * float(b[1, 1]) - cc.REAR_SPAN_MIN
    g["T9.5.2_rw_chord_min"] = float(b[1, 0] - b[0, 0]) - cc.REAR_CHORD_MIN
    g["T9.5.2_rw_chord_max"] = cc.REAR_CHORD_MAX - float(b[1, 0] - b[0, 0])
    g["T9.5.3_rw_thick_min"] = float(b[1, 2] - b[0, 2]) - cc.WING_THICK_MIN
    g["T9.5.3_rw_thick_max"] = cc.WING_THICK_MAX - float(b[1, 2] - b[0, 2])''',
         new='''    # evo_patches W5: per section (a bent or swept wing's bounding box says nothing), and T9.5.1 strict (W4)
    ys = stations(d)
    secs = [section_at(d, y) for y in ys]
    sh = [span_shape(d, y) for y in ys]
    chords = [float(P[:, 0].max() - P[:, 0].min()) for P in secs]
    thick = [_thickness_of(P) for P in secs]
    tops = [float(dz) + float(P[:, 1].max()) for P, (_dx, dz, _c, _t) in zip(secs, sh)]
    if any(p.get("root", 1.0) < 1.0 - 1e-6 for p in placed):
        unbroken = min(2.0 * (float(p["y_mm"]) - 0.5 * max(float(d.support.thick_mm), cc.BUILD_WALL_MM))
                       for p in placed if p.get("root", 1.0) < 1.0 - 1e-6)
    else:
        unbroken = 2.0 * (float(b[1, 1]) - 0.35)
    unbroken -= 0.3                                      # the referee measured 0.2 mm less on a pylon (test, 2026-11-28)
    g["T9.5.1_rw_span"] = unbroken - cc.REAR_SPAN_MIN
    g["T9.5.2_rw_chord_min"] = min(chords) - cc.REAR_CHORD_MIN
    g["T9.5.2_rw_chord_max"] = cc.REAR_CHORD_MAX - max(chords)
    g["T9.5.3_rw_thick_min"] = min(thick) - cc.WING_THICK_MIN
    g["T9.5.3_rw_thick_max"] = cc.WING_THICK_MAX - max(thick)
    g["T9.5.4_rw_height_deviation"] = 15.0 - (max(tops) - min(tops))'''),
    dict(id="W6_design_from_params", file=LG,
         old='''    return dict(chord_mm=float(v["rw_chord"]), t_mm=float(v["rw_t"]), half_span_mm=float(v["rw_half_span"]),
                x_le_mm=float(v["rw_x_le"]), lift_mm=float(v["rw_lift"]),
                support=_strut(v, "rwsup"))''',
         new='''    return dict(chord_mm=float(v["rw_chord"]), t_mm=float(v["rw_t"]), half_span_mm=float(v["rw_half_span"]),
                x_le_mm=float(v["rw_x_le"]), lift_mm=float(v["rw_lift"]),
                # evo_patches/rear_wing_shape.py W6: the shape along the span (AoA stays 0)
                camber=float(v.get("rw_camber", 0.0)), camber_pos=float(v.get("rw_camber_pos", 0.4)),
                rise_mm=float(v.get("rw_rise", 0.0)), rise_p=float(v.get("rw_rise_p", 2.0)),
                sweep_mm=float(v.get("rw_sweep", 0.0)), sweep_p=float(v.get("rw_sweep_p", 1.0)),
                tip_chord_ratio=float(v.get("rw_tip_chord_ratio", 1.0)), tip_t_ratio=float(v.get("rw_tip_t_ratio", 1.0)),
                support=_strut(v, "rwsup"))'''),
    dict(id="W7_rise_within_height", file=RW,
         old='''    base = se._to_trimesh(_wing(d, x_le, 0.0))
    pts0, _ = trimesh.sample.sample_surface_even(base, 1500, seed=1)
    keep = np.ones(len(pts0), bool)
    for a, b in bands():
        keep &= ~((pts0[:, 1] >= a) & (pts0[:, 1] <= b))
    q0 = pts0[keep] if keep.any() else pts0
    lo0 = float(base.bounds[0, 2])

    def short(zr, exact=False):
        q = q0 + np.array([0.0, 0.0, zr])
        if exact:
            gap = min(float(np.concatenate([trimesh.proximity.closest_point(m, q[i:i + 256])[1]
                                            for i in range(0, len(q), 256)]).min()) for m in [body] + others)
        else:
            gap = float(tree.query(q)[0].min()) - 0.15
        gap = min(gap, float(q[:, 2].min()))
        w_ = base.copy(); w_.apply_translation([0.0, 0.0, zr])
        return cc.CLEAR_AIR + MARGIN_MM - gap, _front_view_top(w_, fwd) + 0.1 - (lo0 + zr)
    lo, hi = 0.0, cc.REAR_Z_MAX - 1.0
    if max(short(hi)) > 0:
        z = hi
    else:
        for _ in range(16):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if max(short(mid)) > 0 else (lo, mid)
        z = hi
        for _ in range(3):
            a, b = short(z, exact=True)
            if max(a, b) <= 0:
                break
            z += max(a, b) + 0.05
    wing_tm = se._to_trimesh(_wing(d, x_le, z))''',
         new='''    import wing_strut as _ws0
    _top_max = cc.REAR_Z_MAX - MARGIN_MM - _ws0.FACE_GROW_MM - 0.05
    for _rise_try in range(8):                        # evo_patches W7: a rise the height limit (T9.4.3) cannot take
                                                      # is cut and the height found again (as W2 does for the sweep)
        base = se._to_trimesh(_wing(d, x_le, 0.0))
        pts0, _ = trimesh.sample.sample_surface_even(base, 1500, seed=1)
        keep = np.ones(len(pts0), bool)
        for a, b in bands():
            keep &= ~((pts0[:, 1] >= a) & (pts0[:, 1] <= b))
        q0 = pts0[keep] if keep.any() else pts0
        lo0 = float(base.bounds[0, 2])

        def short(zr, exact=False):
            q = q0 + np.array([0.0, 0.0, zr])
            if exact:
                gap = min(float(np.concatenate([trimesh.proximity.closest_point(m, q[i:i + 256])[1]
                                                for i in range(0, len(q), 256)]).min()) for m in [body] + others)
            else:
                gap = float(tree.query(q)[0].min()) - 0.15
            gap = min(gap, float(q[:, 2].min()))
            w_ = base.copy(); w_.apply_translation([0.0, 0.0, zr])
            return cc.CLEAR_AIR + MARGIN_MM - gap, _front_view_top(w_, fwd) + 0.1 - (lo0 + zr)
        lo, hi = 0.0, cc.REAR_Z_MAX - 1.0
        if max(short(hi)) > 0:
            z = hi
        else:
            for _ in range(16):
                mid = 0.5 * (lo + hi)
                lo, hi = (mid, hi) if max(short(mid)) > 0 else (lo, mid)
            z = hi
            for _ in range(3):
                a, b = short(z, exact=True)
                if max(a, b) <= 0:
                    break
                z += max(a, b) + 0.05
        wing_tm = se._to_trimesh(_wing(d, x_le, z))
        _over = float(wing_tm.bounds[1, 2]) - _top_max
        if _over <= 0 or abs(float(d.rise_mm)) < 1e-6:
            break
        _r = float(d.rise_mm)
        d = dataclasses.replace(d, rise_mm=float(np.sign(_r) * max(0.0, abs(_r) - _over - 0.1)))
        info["rise_cut_for_T9.4.3_mm"] = round(float(info["shape"]["rise_mm"]) - float(d.rise_mm), 3)
    info["shape"]["rise_mm"] = float(d.rise_mm)'''),
    dict(id="L7a_strict_span_fn", file=REF,
         old='''def _wing_checks(assembly: dict, body_mm, parts_mm: dict) -> dict:''',
         new='''def _rear_span_strict(elements, blockers, tol_mm=0.05):
    """evo_patches L7, T9.5.1 strict (docs/rules/README.md): the longest stretch of the rear wing that no other part
    touches or passes through, its own supports included; on the full car (the stretch through the centre counts
    its two halves). Points of each blocker within tol_mm of (or inside) the wing break the span at their y."""
    import trimesh
    W = trimesh.util.concatenate(list(elements))
    y_hi = float(min(e.bounds[1, 1] for e in elements))
    y_lo = float(max(0.0, max(e.bounds[0, 1] for e in elements)))
    hit = []
    for m in blockers:
        if m is None or not len(getattr(m, "faces", [])):
            continue
        n = int(min(80000, max(4000, m.area * 6)))
        P = m.sample(n)
        P = P[(P[:, 1] >= y_lo - 1.0) & (P[:, 1] <= y_hi + 1.0) & (P[:, 0] >= W.bounds[0, 0] - 1.0)
              & (P[:, 0] <= W.bounds[1, 0] + 1.0) & (P[:, 2] >= W.bounds[0, 2] - 1.0) & (P[:, 2] <= W.bounds[1, 2] + 1.0)]
        if not len(P):
            continue
        dist = np.concatenate([trimesh.proximity.closest_point(W, P[i:i + 512])[1] for i in range(0, len(P), 512)])
        inside = W.contains(P) if W.is_volume else np.zeros(len(P), bool)
        hit.extend(P[(dist <= tol_mm) | inside, 1].tolist())
    if not hit:
        return 2.0 * y_hi if y_lo <= 0.05 else y_hi - y_lo, []
    ys = np.sort(np.clip(np.asarray(hit), y_lo, y_hi))
    iv, a, b = [], ys[0], ys[0]
    for y in ys[1:]:
        if y - b > 0.3:
            iv.append((a, b)); a = y
        b = y
    iv.append((a, b))
    runs, prev = [], y_lo
    for a, b in iv:
        runs.append((prev, a)); prev = b
    runs.append((prev, y_hi))
    best = 0.0
    for a, b in runs:
        L = (b - a) * (2.0 if a <= y_lo + 0.05 and y_lo <= 0.05 else 1.0)
        best = max(best, L)
    return best, [(round(float(a), 2), round(float(b), 2)) for a, b in iv]


def _wing_checks(assembly: dict, body_mm, parts_mm: dict) -> dict:'''),
    dict(id="L7b_strict_span_use", file=REF,
         old='''                                    el, [body_mm] + list(others.values()), 2))
        back = float(parts_mm["rwing"].bounds[0, 0])''',
         new='''                                    el, [body_mm] + list(others.values()), 2))
            # evo_patches L7: T9.5.1 strict -- the wing's own supports break its span too
            _span, _breaks = _rear_span_strict(el, [body_mm] + list(others.values()) + list(rest))
            r["T9.5.1_rear_span"] = _entry(_span - 50.0, f"{_span:.2f} mm unbroken (strict: own supports break it; "
                                                         f"breaks at y {_breaks})")
        back = float(parts_mm["rwing"].bounds[0, 0])'''),
]

NEW_PARAMS = [
    ("rw_camber_pos", 0.4, 0.2, 0.6, "-", "where the camber peaks, share of the chord"),
    ("rw_rise", 0.0, -13.5, 13.5, "mm", "the tip's height above the root (- droops): gull, arch, dihedral; T9.5.4 <= 15"),
    ("rw_rise_p", 2.0, 1.0, 4.0, "-", "how the height runs out to the tip (s^p)"),
    ("rw_sweep", 0.0, -10.0, 10.0, "mm", "the tip's LE behind the root's (- forward)"),
    ("rw_sweep_p", 1.0, 1.0, 3.0, "-", "how the sweep runs out to the tip (s^p)"),
    ("rw_tip_chord_ratio", 1.0, 0.62, 1.6, "-", "tip chord / root chord (each section kept 15.3-24.7 mm, T9.5.2)"),
    ("rw_tip_t_ratio", 1.0, 0.5, 2.4, "-", "tip thickness / root thickness (each section kept 2.3-5.7 mm, T9.5.3)"),
]
WHY = "owner 2026-11-28: the rear wing any shape the regulations allow, zero angle of attack, camber allowed (evo_patches/rear_wing_shape.py)"


def _params(root: Path) -> str:
    p = root / LP
    P = json.loads(p.read_text())
    names = {q["name"] for q in P["params"]}
    if all(n in names for n, *_ in NEW_PARAMS) and any(q["name"] == "rw_camber" and q["mode"] == "shape" for q in P["params"]):
        return "already applied"
    for q in P["params"]:
        if q["name"] == "rw_camber":
            q.update(lo=-0.06, hi=0.06, mode="shape", why=WHY)
    i = max(k for k, q in enumerate(P["params"]) if q["name"].startswith("rw_")) + 1
    for n, v, lo, hi, unit, what in NEW_PARAMS:
        if n not in names:
            P["params"].insert(i, {"group": "rear wing", "name": n, "value": v, "lo": lo, "hi": hi, "mode": "shape",
                                   "unit": unit, "what": what, "rule": "T9.5", "why": WHY})
            i += 1
    p.write_text(json.dumps(P, indent=1) + "\n")
    return "applied"


def apply(evo) -> str:
    root = Path(evo)
    done = []
    for e in EDITS:
        p = root / e["file"]
        t = p.read_text()
        if e["new"] in t:
            done.append(f"{e['id']}: already applied")
            continue
        assert t.count(e["old"]) == 1, f"rear_wing_shape {e['id']}: anchor not found once in {e['file']}"
        p.write_text(t.replace(e["old"], e["new"], 1))
        done.append(f"{e['id']}: applied")
    done.append(f"params: {_params(root)}")
    return "rear_wing_shape: " + "; ".join(done)
