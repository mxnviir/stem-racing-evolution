"""Runtime fix of the pinned evolution builder (car/body/freeform.py, field): the sidepod-top cap without hard edges.

Before: between the axles, nothing outboard of the 12.5 mm spine may rise above 27.175 mm (the team's sidepod top), applied
as np.where(between, np.maximum(phi, cap), phi): a hard maximum (a sharp-edged flat shelf along the sidepods: 214 mm2 of
flat cap and a crease on the team layout, carhybrid/shape_audit.py) switched on and off abruptly at the two axle planes
(a step in the body at each axle).
After: the same limit, applied with the field's own smooth maximum (blend_mm, as the top/bottom/side surfaces already
meet), and faded in over CAP_TAPER_MM OUTSIDE the axle span (fully on everywhere between the axles). Both only ever
remove material relative to the hard version (smooth max >= max; the taper adds the cap where it was not applied),
so the 27.175 mm limit still holds everywhere it held before. The re-entrant corner between the spine and the cap is
left sharp on purpose: rounding it would add material above the limit.

    import soft_cap; soft_cap.apply(<evolution checkout>)
"""
from pathlib import Path

OLD = """        between = (X > xf) & (X < xr)
        cap = np.minimum(Z - SIDEPOD_TOP_MM, np.abs(Y) - SPINE_HALF_MM)   # > 0 above the top AND outboard of the spine
        phi = np.where(between, np.maximum(phi, cap), phi)"""
NEW = """        cap = np.minimum(Z - SIDEPOD_TOP_MM, np.abs(Y) - SPINE_HALF_MM)   # > 0 above the top AND outboard of the spine
        # evo_patches/soft_cap.py: smooth maximum (no sharp shelf edge) and a taper just outside the axle span (no step
        # at the axles); fully applied between the axles, so the limit holds as before
        t = np.clip(np.minimum(X - xf, xr - X) / CAP_TAPER_MM + 1.0, 0.0, 1.0)
        w = t * t * (3.0 - 2.0 * t)
        phi = phi + w * (_smax(phi, cap, k) - phi)"""
CONST_OLD = "SPINE_HALF_MM = 12.5            # the spine: the halo pocket's half-width (Appendix ix: 25 mm wide)\n"
CONST_NEW = CONST_OLD + "CAP_TAPER_MM = 5.0              # evo_patches/soft_cap.py: the cap fades in over this length outside the axles\n"


def apply(evo):
    p = Path(evo) / "car" / "body" / "freeform.py"
    t = p.read_text()
    if "evo_patches/soft_cap.py" in t:
        return "soft_cap: already applied"
    assert t.count(OLD) == 1 and t.count(CONST_OLD) == 1, "soft_cap: anchor not found (builder changed?)"
    p.write_text(t.replace(OLD, NEW, 1).replace(CONST_OLD, CONST_NEW, 1))
    return "soft_cap: applied"


if __name__ == "__main__":
    import sys
    print(apply(sys.argv[1]))
