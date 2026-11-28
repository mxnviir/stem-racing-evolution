"""Runtime fix of the pinned evolution builder (car/body/freeform.py, field): the sidepod-top cap without hard edges or a flat deck.

Before: between the axles, nothing outboard of the 12.5 mm spine may rise above 27.175 mm (the team's sidepod top), applied
as np.where(between, np.maximum(phi, cap), phi): a hard maximum (a sharp-edged flat shelf along the sidepods: 214 mm2 of
flat cap and a crease on the team layout, carhybrid/shape_audit.py) switched on and off abruptly at the two axle planes
(a step in the body at each axle).
After: the same limit, applied with the field's own smooth maximum (blend_mm, as the top/bottom/side surfaces already
meet), and faded in over CAP_TAPER_MM OUTSIDE the axle span (fully on everywhere between the axles). Both only ever
remove material relative to the hard version (smooth max >= max; the taper adds the cap where it was not applied),
so the 27.175 mm limit still holds everywhere it held before. The re-entrant corner between the spine and the cap is
left sharp on purpose: rounding it would add material above the limit.

v2 (2026-11-28, sweep2's shape audits): the same compression where the halo's visibility rule cuts the body (T4.4.2,
halo_seat.view_air: nothing above the rim, 24 + 3.175 = 27.175 mm, ahead of the halo across its width or beside it
along its length; T4.4.3 over it). The spine (|y| < 12.5 mm) was not compressed there, so the rule's cut truncated the
seed's spine (rising to 30-33 mm ahead of the halo) into a flat deck at the rim from Ref A to the halo: 170-660 mm2 of
flat top per car in sweep2. field() takes halo = freeform.halo_span(Ref A, d_halo); param_body.build passes it, and
carhybrid/grad.py passes the same, so the gradient sees the field the builder builds. The weight is zero at Ref A and
ahead of it, so the nose (and evo_patches/nose_material.py's foam/print split there) is unchanged.
v2.1 (sweep3): the halo-region weight fades out over the 5 mm BEFORE the halo's rear face, not after it; behind the rear
face the body must stand to the halo's top (team_halo_face_joined failed on L03 and L19 with the fade behind it).

    import soft_cap; soft_cap.apply(<evolution checkout>)
"""
from pathlib import Path

TOP_OLD = """    top, bot, side = surfaces(ff, X, Y, Z, x_tip, x_rear)
    k = ff.blend_mm"""
TOP_NEW = """    top, bot, side = surfaces(ff, X, Y, Z, x_tip, x_rear)
    k = ff.blend_mm
    if axles is not None or halo is not None:
        # evo_patches/soft_cap.py: where a 27.175 mm limit applies (outboard of the spine between the axles: the team's
        # sidepod top; ahead of, beside and over the halo: T4.4.2-3, whose rim is 24 + 3.175 = the same 27.175 mm) the
        # top surface is COMPRESSED into the last CAP_SOFT_MM below it (z0 + s * tanh((top - z0) / s), s = CAP_SOFT_MM):
        # monotone, smooth, always below the limit, the identity well below it. Truncating it (the cap or the rule's cut
        # alone) left a flat deck wherever the seed's top rose above 27.175 mm; this keeps the surface's own shape.
        s_ = CAP_SOFT_MM
        z0 = SIDEPOD_TOP_MM - s_
        tc = np.where(top > z0, z0 + s_ * np.tanh((top - z0) / s_), top)
        wt = np.zeros(np.shape(top))
        if axles is not None:
            xf_, xr_ = axles[0], axles[1]
            tx = np.clip(np.minimum(X - xf_, xr_ - X) / CAP_TAPER_MM + 1.0, 0.0, 1.0)
            ty = np.clip((np.abs(Y) - SPINE_HALF_MM) / CAP_TAPER_MM + 1.0, 0.0, 1.0)
            wt = (tx * tx * (3.0 - 2.0 * tx)) * (ty * ty * (3.0 - 2.0 * ty))
        if halo is not None:
            xa_, xh0_, xh1_, hw_ = halo
            T_ = CAP_TAPER_MM

            def ss_(u):
                c = np.clip(u, 0.0, 1.0)
                return c * c * (3.0 - 2.0 * c)
            # Ref A .. the halo's rear face, fading out BEFORE it: behind the rear face the body must stand up to the
            # halo's top to back it (referee: team_halo_face_joined; sweep3 L03 / L19 lost it when the fade ran past it)
            wh = (ss_((X - xa_) / T_) * ss_((xh1_ - X) / T_) *
                  np.maximum(ss_((hw_ + T_ - np.abs(Y)) / T_), ss_((X - xh0_ + T_) / T_)))   # its width ahead; all y beside
            wt = np.maximum(wt, wh)
        top = top + wt * (tc - top)"""
OLD = """        between = (X > xf) & (X < xr)
        cap = np.minimum(Z - SIDEPOD_TOP_MM, np.abs(Y) - SPINE_HALF_MM)   # > 0 above the top AND outboard of the spine
        phi = np.where(between, np.maximum(phi, cap), phi)"""
NEW = """        cap = np.minimum(Z - SIDEPOD_TOP_MM, np.abs(Y) - SPINE_HALF_MM)   # > 0 above the top AND outboard of the spine
        # evo_patches/soft_cap.py: smooth maximum (no sharp shelf edge) and a taper just outside the axle span (no step
        # at the axles); fully applied between the axles, so the limit holds as before
        t = np.clip(np.minimum(X - xf, xr - X) / CAP_TAPER_MM + 1.0, 0.0, 1.0)
        w = t * t * (3.0 - 2.0 * t)
        phi = phi + w * (_smax(phi, cap, k) - phi)"""
SIG_OLD = "def field(ff: Freeform, X, Y, Z, x_tip, x_rear, axles=None):\n"
SIG_NEW = ('''def halo_span(ref_a_mm, d_halo_mm):
    """(evo_patches/soft_cap.py) (Ref A, halo front, halo rear, halo half-width), mm: where T4.4.2-3 keep the body
    below the halo's rim (halo_seat.view_air, seat_air); for field(..., halo=)."""
    import halo_pocket as hp
    import halo_seat as hs
    x0 = hp.compute_halo_pocket_box_m(ref_a_mm / 1e3, d_halo_mm)["x_min_m"] * 1e3
    x0, x1, hw, _top = hs.box_mm(x0)
    return (float(ref_a_mm), float(x0), float(x1), float(hw))


def field(ff: Freeform, X, Y, Z, x_tip, x_rear, axles=None, halo=None):
''')
CALL_OLD = "                        axles=(x_front_mm, x_front_mm + W_mm, None))\n"
CALL_NEW = ("                        axles=(x_front_mm, x_front_mm + W_mm, None),\n"
            "                        halo=fr.halo_span(x_ref_a, d_halo_mm))     # evo_patches/soft_cap.py v2\n")
CONST_OLD = "SPINE_HALF_MM = 12.5            # the spine: the halo pocket's half-width (Appendix ix: 25 mm wide)\n"
CONST_NEW = CONST_OLD + ("CAP_TAPER_MM = 5.0              # evo_patches/soft_cap.py: the cap fades in over this length outside the axles\n"
                         "CAP_SOFT_MM = 6.0               # evo_patches/soft_cap.py: the top surface is compressed into this band below the cap\n")


def apply(evo):
    p = Path(evo) / "car" / "body" / "freeform.py"
    q = Path(evo) / "car" / "body" / "param_body.py"
    t, u = p.read_text(), q.read_text()
    if "evo_patches/soft_cap.py" in t:
        assert SIG_NEW in t and CALL_NEW in u, "soft_cap: an older version is applied (start from a clean checkout)"
        return "soft_cap: already applied"
    assert (t.count(OLD) == 1 and t.count(CONST_OLD) == 1 and t.count(TOP_OLD) == 1 and t.count(SIG_OLD) == 1
            and u.count(CALL_OLD) == 1), "soft_cap: anchor not found (builder changed?)"
    p.write_text(t.replace(OLD, NEW, 1).replace(CONST_OLD, CONST_NEW, 1).replace(TOP_OLD, TOP_NEW, 1)
                 .replace(SIG_OLD, SIG_NEW, 1))
    q.write_text(u.replace(CALL_OLD, CALL_NEW, 1))
    return "soft_cap: applied (v2: the sidepods and the halo's visibility region)"


if __name__ == "__main__":
    import sys
    print(apply(sys.argv[1]))
