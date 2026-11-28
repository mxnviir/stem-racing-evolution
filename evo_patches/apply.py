"""Runtime fixes to the pinned evolution builder's mass bookkeeping (mxnviir/stem-racing-evolution @ e9c896b).

    python evo_patches/apply.py <evolution checkout> [--check]

Idempotent: every patch asserts that either its exact anchor text or its already-patched text is in the file,
replaces the anchor once, and leaves an already-patched file alone. --check only reports (exit 1 if a patch is
neither applicable nor applied). See evo_patches/README.md for the why and the before/after numbers.

The race model (car/race, search/layout/layout_car.mass, simscale_layout.race_time) is NOT touched: the fixes
are in the builder (run_car.py) and the parts engine (assembly.py), so that the referee's T3.6 weighs the car
the race model weighs, and the printed parts are weighed at the material the team prints.
"""
from __future__ import annotations

import sys
from pathlib import Path

TAG = "evo_patches"

PATCHES = [
    # ------------------------------------------------------------------------------------------------ P1
    dict(
        id="P1_no_joints_mass_includes_supports",
        file="car/build/run_car.py",
        why="--no-joints (every search build): the referee's T3.6 'mass as manufactured' summed the body, the "
            "printed nose/wings/guides, the bare wheels and the halo, and left out the wheel supports and hubcaps "
            "(8.50 g on the team car), which the race model (layout_car.mass, via fixed_hardware_kwargs' wheel "
            "masses) and mass_state both carry. T3.6 read 45.29 g ('needs lead') for a 53.79 g car.",
        old='''        body_g = sum(c.mass_kg for c in up.compute_mass_com(geom)) * 1e3
        asm["joints"] = {"manufactured_g": body_g + sum(asm["parts_mass_g"].values()), "mass_delta_kg": {},
                         "manufactured": {}, "skipped": True}
        return asm''',
        new='''        body_g = sum(c.mass_kg for c in up.compute_mass_com(geom)) * 1e3
        # evo_patches P1: the wheel supports (part + hubcap - the foam its slot removes, per side, as Part 4
        # puts them in the wheels' masses the race model uses) are made and weighed too. Without them T3.6
        # read 45.29 g for the team car the race model weighs at 53.79 g (supports 8.50 g).
        sup_g = 2.0 * sum(asm.get("info", {}).get("support_mass_g_each", []))
        asm["joints"] = {"manufactured_g": body_g + sum(asm["parts_mass_g"].values()) + sup_g, "mass_delta_kg": {},
                         "manufactured": {}, "skipped": True, "supports_g": sup_g}
        return asm''',
    ),
    # ------------------------------------------------------------------------------------------------ P2
    dict(
        id="P2_fixed_skin_is_not_a_failed_sizing",
        file="car/build/run_car.py",
        why="A fixed-skin build (--skin-mm, every search build) is never sized to a mass, but build_car still "
            "compared it with the 48.2 g target and logged 'MASS NOT CONVERGED' with mass_converged False. "
            "Nothing resizes the body (param_body.build returns make(skin) at once); the line only misled.",
        old='''    info["car_mass_error_g"] = err * 1e3
    info["mass_converged"] = abs(err) <= MASS_TOL_KG
    if not info["mass_converged"]:
        log(f"  MASS NOT CONVERGED: car {(T + err)*1e3:.3f} g, {err*1e3:+.3f} g from the target")
    return geom, info, asm''',
        new='''    info["car_mass_error_g"] = err * 1e3
    # evo_patches P2: a fixed skin is not sized (param_body.build: make(skin_offset_mm), one pass): its
    # body is what the design draws, so the optimiser's mass changes reach the race model; no target applies.
    info["sized"] = getattr(a, "skin_mm", None) is None
    info["mass_converged"] = (abs(err) <= MASS_TOL_KG) if info["sized"] else None
    if info["sized"] and not info["mass_converged"]:
        log(f"  MASS NOT CONVERGED: car {(T + err)*1e3:.3f} g, {err*1e3:+.3f} g from the target")
    elif not info["sized"]:
        log(f"  fixed skin, not sized: car {(T + err)*1e3:.3f} g as built")
    return geom, info, asm''',
    ),
    # ------------------------------------------------------------------------------------------------ P3
    dict(
        id="P3_printed_parts_are_ABS",
        file="car/parts/assembly.py",
        why="The team prints ABS (team, 2026-11-14: joints.DENSITY_G_CM3, nose.NoseCone.material, "
            "beam_support/support_engine ABS_G_CM3, carhybrid/grad.py's nose density). assembly.py still weighed "
            "the nose (+ front-wing pylons), wings, wing struts and tether guides as PLA 1.24 g/cm3, so a search "
            "car's printed parts were 19 % heavier than the same parts in the joints path (the finalist's "
            "manufacturing files, 1.04).",
        old='''PART_MATERIAL = "PLA"           # wings, mounts, tether guides: printed, solid''',
        new='''# evo_patches P3: ABS, as the team prints (team, 2026-11-14) and as joints.DENSITY_G_CM3 weighs the same parts.
PART_MATERIAL = "ABS"           # nose, wings, mounts, tether guides: printed, solid''',
    ),
]


def apply(root: Path, check_only: bool = False) -> list[tuple[str, str]]:
    out = []
    for p in PATCHES:
        f = root / p["file"]
        if not f.exists():
            raise SystemExit(f"{p['id']}: {f} not found (is {root} an evolution checkout?)")
        s = f.read_text()
        if p["new"] in s:
            out.append((p["id"], "already applied"))
            continue
        n = s.count(p["old"])
        if n != 1:
            if check_only:
                out.append((p["id"], f"NOT APPLICABLE (anchor found {n} times)"))
                continue
            raise SystemExit(f"{p['id']}: anchor text found {n} times in {p['file']} (expected exactly 1): "
                             "the checkout is not the pinned e9c896b builder")
        if check_only:
            out.append((p["id"], "applicable"))
            continue
        f.write_text(s.replace(p["old"], p["new"], 1))
        out.append((p["id"], "applied"))
    return out


def main(argv):
    if not argv or argv[0] in ("-h", "--help"):
        print(__doc__)
        return 2
    root = Path(argv[0]).resolve()
    check = "--check" in argv[1:]
    res = apply(root, check_only=check)
    for pid, st in res:
        print(f"{TAG}: {pid}: {st}")
    if not check:                       # the separate patch modules (each idempotent, each with its own anchors)
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        for mod in ("soft_cap", "nose_material"):
            print(f"{TAG}: {__import__(mod).apply(root)}")
    return 1 if any(st.startswith("NOT") for _p, st in res) else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
