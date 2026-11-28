"""Builder fixes so that the geometry engine does not make the cars the referee then rejects (sweep2-sweep4:
half of 32 random layouts illegal). Each fix makes the PART that broke a rule respect it by construction; the
referee (legality.py) is not touched.

    python evo_patches/apply.py <evolution checkout>      (applies this module with the others)

Idempotent, as the other modules: each edit asserts that its anchor or its patched text is in the file, replaces
the anchor once and leaves a patched file alone.

L1 rear wing behind the body's back (T5.6, 7 of 32 cars in sweep2 and sweep4: L00 L02 L06 L08 L14 L16 L22).
   The rear-wing engine placed the wing anywhere from Ref B (rear axle + 16) to Ref B + 40 (T9.4.2). On a short tail
   (rear face less than ~56 mm behind the axle) that put the wing's trailing edge, and its struts, behind the body's
   back face, where T5.6 wants nothing round the protruding cartridge (the referee's ring, from the rear face to the
   cartridge's end, all round: 130-1820 mm3 of wing there). Now the wing and its supports also end before the body's
   back face (0.3 mm margin): the chord is cut to fit (never below T9.5.2's 15 mm), also after T7.9 has pushed the
   wing aft, and a support that reaches past the back face is refused like one past T9.4.2.
"""
from __future__ import annotations

from pathlib import Path

RW = "car/parts/rear_wing_engine.py"
LG = "car/build/legality.py"
RC = "car/build/run_car.py"
SE = "car/parts/support_engine.py"
FW = "car/parts/front_wing_engine.py"

EDITS = [
    dict(
        id="L1a_rear_wing_room_ends_at_body_back",
        file=RW,
        old='''    room = cc.REAR_OVERHANG_MAX - 2 * MARGIN_MM - chord
    x_le = float(np.clip(ref_b + d.x_le_mm, ref_b + MARGIN_MM, ref_b + MARGIN_MM + max(room, 0.0)))
''',
        new='''    # evo_patches L1: the wing (and its supports) also end before the body's back face: behind it, T5.6 wants
    # nothing round the protruding cartridge. A short tail therefore gets a shorter chord (never below T9.5.2).
    body_back = float(body.bounds[1, 0])
    x_aft_max = min(ref_b + cc.REAR_OVERHANG_MAX, body_back) - MARGIN_MM
    chord = float(np.clip(min(chord, x_aft_max - ref_b - MARGIN_MM - 0.05), cc.REAR_CHORD_MIN + MARGIN_MM, cc.REAR_CHORD_MAX - MARGIN_MM))
    info["x_aft_max_mm"] = x_aft_max
    room = (x_aft_max - ref_b) - MARGIN_MM - chord
    x_le = float(np.clip(ref_b + d.x_le_mm, ref_b + MARGIN_MM, ref_b + MARGIN_MM + max(room, 0.0)))
''',
    ),
    dict(
        id="L1b_rear_wing_chord_after_T79",
        file=RW,
        old='''    d = dataclasses.replace(d, half_span_mm=hs, chord_mm=chord, t_mm=t)
    info.update(half_span_mm=hs, chord_mm=chord, t_mm=t, x_le_mm=x_le)
''',
        new='''    # evo_patches L1: T7.9 may have pushed the wing aft: cut the chord so the trailing edge stays ahead of the
    # body's back face (the section's own x extent, so a cambered or pitched section is measured as built)
    for _ in range(3):
        _sx = float(section(dataclasses.replace(d, half_span_mm=hs, chord_mm=chord, t_mm=t))[:, 0].max())
        over = x_le + _sx - x_aft_max
        if over <= 0 or chord <= cc.REAR_CHORD_MIN + MARGIN_MM + 1e-6:
            break
        chord = float(max(cc.REAR_CHORD_MIN + MARGIN_MM, chord - over - 0.05))
        info["chord_cut_for_T5.6_mm"] = round(float(d.chord_mm) - chord, 3)
    d = dataclasses.replace(d, half_span_mm=hs, chord_mm=chord, t_mm=t)
    info.update(half_span_mm=hs, chord_mm=chord, t_mm=t, x_le_mm=x_le)
''',
    ),
    dict(
        id="L1c_rear_wing_support_ahead_of_body_back",
        file=RW,
        old='''            return (b[0, 0] >= ref_b + MARGIN_MM and b[1, 0] <= ref_b + cc.REAR_OVERHANG_MAX - MARGIN_MM
''',
        new='''            return (b[0, 0] >= ref_b + MARGIN_MM and b[1, 0] <= min(ref_b + cc.REAR_OVERHANG_MAX - MARGIN_MM, x_aft_max)
''',
    ),
    dict(
        id="L2a_clamp_helper",
        file=RC,
        old='''def export_body(geom, path: Path) -> dict:
''',
        new='''def _clamp_to_official_block(half, path):
    """evo_patches L2: the half body (metres) intersected with the official Model Block as the referee places it
    (legality.block_body_inside_official_block: rear face = measured bore end + 52 mm, bore on the measured axis,
    underside slot), every face 0.02 mm inside the referee's. Returns the body (unchanged if the boolean fails or
    would remove more than 5 %), and writes what it did to body_clamp.json next to path."""
    import numpy as np
    import trimesh
    import legality as _lg
    rep = {"applied": False}
    try:
        pm = half.copy(); pm.apply_scale(1e3)
        mir = pm.copy(); mir.vertices[:, 1] *= -1.0; mir.invert()
        full = trimesh.boolean.union([pm, mir], engine="manifold")
        ch = _lg.chamber_measures(full)
        if ch is None:
            rep["why"] = "no chamber found"
        else:
            rear = ch["rim_x_mm"] - ch["depth_mm"] + _lg.rc.BLOCK_BORE_DEPTH_MM
            blk = _lg.official_block_mm(rear, slack_mm=-0.02, axis_z_mm=ch["axis_z_mm"])
            out = trimesh.boolean.intersection([pm, blk], engine="manifold")
            v0 = abs(float(pm.volume))
            removed = v0 - abs(float(out.volume)) if out is not None and not out.is_empty else v0
            rep.update(block_rear_face_mm=rear, removed_mm3=round(removed, 3), body_mm3=round(v0, 1))
            if out is None or out.is_empty or not out.is_volume:
                rep["why"] = "the intersection is not a volume"
            elif removed > 0.05 * v0:
                rep["why"] = "it would remove more than 5 % of the body"
            else:
                out.apply_scale(1e-3)
                half, rep["applied"] = out, True
    except Exception as e:                                   # noqa: BLE001  (never fails a build: the referee judges)
        rep["why"] = f"error: {e!r}"[:300]
    path.with_name("body_clamp.json").write_text(json.dumps(rep, indent=1, default=float))
    return half


def export_body(geom, path: Path) -> dict:
''',
    ),
    dict(
        id="L2b_clamp_call",
        file=RC,
        old='''    if len(half.faces) > MAX_BODY_FACES:      # a 0.3 mm grid otherwise runs a 16 GB machine out of memory
''',
        new='''    half = _clamp_to_official_block(half, path)      # evo_patches L2: inside the block the referee places
    if len(half.faces) > MAX_BODY_FACES:      # a 0.3 mm grid otherwise runs a 16 GB machine out of memory
''',
    ),
    dict(
        id="L3a_support_variants_helper",
        file=SE,
        old='''def build(sd: SupportDesign, body_half_m, axles, out_dir=None, keep_out_mm=None):
''',
        new='''def _build_axle_checked(sd, body_mm, x, z, R, y_in, w, hub, tag):
    """evo_patches L3: build_axle_sized, measured (gates); if any gate of this axle is short, the design's
    variants are tried and the first meeting every gate is kept, else the one whose worst gate is least short."""
    variants = [sd]
    for kw in (dict(chord_mm=sd.chord_mm * 0.85), dict(chord_mm=sd.chord_mm * 1.15),
               dict(saddle_h_mm=sd.saddle_h_mm * 1.5), dict(saddle_h_mm=sd.saddle_h_mm * 0.75),
               dict(chord_mm=sd.chord_mm * 0.7),
               # L3 v2 (sweep5 L21: foam over the slot 0.24 mm short of the 3 mm wall at a thin tail): a saddle
               # deep enough to cut through, so no foam is left over it (the gate's "or absent")
               dict(saddle_h_mm=sd.saddle_h_mm * 2.5), dict(saddle_h_mm=sd.saddle_h_mm * 4.0)):
        variants.append(dataclasses.replace(sd, **kw))
    best = None
    for i, v in enumerate(variants):
        try:
            r = build_axle_sized(v, body_mm, x, z, R, y_in, w, hub, tag)
        except Exception:                                    # noqa: BLE001  (a variant that cannot be built)
            if i == 0 and len(variants) == 1:
                raise
            continue
        r["_struct"]["disc_r_mm"] = v.disc_r_mm
        g = gates(r, tag)
        worst = min(g.values()) if g else 0.0
        r["_struct"]["design_variant"] = i
        r["_struct"]["design_variant_params"] = {k: getattr(v, k) for k in ("chord_mm", "saddle_h_mm")}
        if worst >= 0:
            return r, g
        if best is None or worst > best[2]:
            best = (r, g, worst)
    if best is None:
        raise RuntimeError(f"support_engine: no variant of the {tag} support could be built")
    return best[0], best[1]


def build(sd: SupportDesign, body_half_m, axles, out_dir=None, keep_out_mm=None):
''',
    ),
    dict(
        id="L3b_support_variants_call",
        file=SE,
        old='''        r = build_axle_sized(sd, body_mm, x, z, R, y_in, w, hub, tag)
        r["_struct"]["disc_r_mm"] = sd.disc_r_mm
        gts.update(gates(r, tag))
''',
        new='''        r, _g = _build_axle_checked(sd, body_mm, x, z, R, y_in, w, hub, tag)     # evo_patches L3
        gts.update(_g)
''',
    ),
    dict(
        id="L4_front_wing_segment",
        file=FW,
        old='''        r2["info"]["lowered_for_clear_air_mm"] = r["info"].get("lowered_for_clear_air_mm", 0.0) + short + MARGIN_MM
        r = r2
    return r
''',
        new='''        r2["info"]["lowered_for_clear_air_mm"] = r["info"].get("lowered_for_clear_air_mm", 0.0) + short + MARGIN_MM
        r = r2
    # evo_patches L4: T8.6.1 is measured from the side strut where it leaves the wing's root (its flare), not
    # from the root: a wing short on its own segment gate is built again with the tip moved out by the shortfall
    for _ in range(2):
        short = -r["gates"].get("T8.6.1_fw_segment", 0.0)
        if short <= 0:
            break
        d = dataclasses.replace(d, z_mm=r["info"]["z_mm"], tip_y_mm=r["info"]["tip_y_mm"] + short + MARGIN_MM + 0.2)
        r2 = _build_once(d, nose_half_m, x_front_mm, wheel, body_half_m)
        r2["info"]["tip_moved_for_T8.6.1_mm"] = r["info"].get("tip_moved_for_T8.6.1_mm", 0.0) + short + MARGIN_MM + 0.2
        if min(r2["gates"].values()) < min(r["gates"].values()) and -r2["gates"].get("T8.6.1_fw_segment", 0.0) >= short:
            break                                          # no better: keep the first
        r = r2
    return r
''',
    ),
    dict(
        id="L5_rear_wing_forward_if_unsupported",
        file=RW,
        old='''def build(d: RearWingDesign, body_half_m, x_front_mm, W_mm, wheel_rear, clear_of=(), forward_parts=()):
''',
        new='''def build(d: RearWingDesign, body_half_m, x_front_mm, W_mm, wheel_rear, clear_of=(), forward_parts=()):
    """evo_patches L5: _build_at, and if none of the wing's supports could be placed, again further forward (the
    x_le gene halved, then at Ref B): the first supported wing is kept, else the first."""
    first = _build_at(d, body_half_m, x_front_mm, W_mm, wheel_rear, clear_of, forward_parts)
    if first["info"].get("supported"):
        return first
    for x_le in (0.5 * d.x_le_mm, 0.0):
        if x_le >= d.x_le_mm - 1e-6:
            continue
        r = _build_at(dataclasses.replace(d, x_le_mm=x_le), body_half_m, x_front_mm, W_mm, wheel_rear, clear_of,
                      forward_parts)
        if r["info"].get("supported"):
            r["info"]["moved_forward_for_supports_mm"] = float(d.x_le_mm - x_le)
            return r
    return first


def _build_at(d: RearWingDesign, body_half_m, x_front_mm, W_mm, wheel_rear, clear_of=(), forward_parts=()):
''',
    ),
    dict(
        id="L6_nose_joint_sliver_face",
        file=LG,
        old='''    inside = contains_xy(face, Yg, Zg)
    y, z = Yg[inside], Zg[inside]
''',
        new='''    inside = contains_xy(face, Yg, Zg)
    if not inside.any():                                  # evo_patches L6: a sliver thinner than the grid
        r["nose_joint_bond_stress"] = _entry(-1.0, f"no usable glued face (a sliver of {face.area:.3f} mm2)")
        return r
    y, z = Yg[inside], Zg[inside]
''',
    ),
]


def apply(evo) -> str:
    root = Path(evo)
    done = []
    for e in EDITS:
        p = root / e["file"]
        t = p.read_text()
        if e["new"] in t:
            done.append(f"{e['id']}: already applied")
            continue
        assert t.count(e["old"]) == 1, f"legality_fixes {e['id']}: anchor not found once in {e['file']}"
        p.write_text(t.replace(e["old"], e["new"], 1))
        done.append(f"{e['id']}: applied")
    return "legality_fixes: " + "; ".join(done)
