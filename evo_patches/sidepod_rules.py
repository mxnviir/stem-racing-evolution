"""Sidepods limited by the rules only (owner, 2026-11-28 evening: "drop it [the team's 27.175 mm sidepod cap]. Just
make the max width of the sidepod the one in the rules").

    python evo_patches/apply.py <evolution checkout>      (applied after soft_cap, which it builds on)

S1 The team's sidepod-top cap is gone. param_body no longer passes `axles` to freeform.field, so neither the soft
   compression into the 6 mm under 27.175 mm between the axles nor the hard cap outboard of the spine applies. The
   halo's own limit stays: T4.4.2-3 keep the body below the halo's rim ahead of it, beside it and over it (soft_cap's
   `halo` region and the builder's halo view air). That is a rule.
S2 (UNDER TEST, not applied yet) The body's width is the rules' width everywhere: the official Model Block (T1.14, T3.1.2; 65 mm wide, so 32.5 mm
   each side), as the envelope's sidepod corridor already was. The envelope's 28 mm body half-width was a team choice
   ("NOT itself a regulation number", bounding_volumes.default_rule_envelope: sized to hold the T4.2 cargo's wide
   end); it is raised to the block's 32.5 mm. Everything that IS a rule still cuts the body: the T7.9 wheel zones, the
   T7.2 wheel keep-clear, T8.5.1's nose box, the halo, the chamber, the block clamp (legality_fixes L2).
"""
from __future__ import annotations

from pathlib import Path

PB = "car/body/param_body.py"
BV = "car/body/bounding_volumes.py"

EDITS = [
    dict(id="S1_no_team_sidepod_cap", file=PB,
         old='''        phi0 = fr.field(freeform, X, Y, Z, freeform.x_tip_mm, freeform.x_rear_mm,
                        axles=(x_front_mm, x_front_mm + W_mm, None),
                        halo=fr.halo_span(x_ref_a, d_halo_mm))     # evo_patches/soft_cap.py v2''',
         new='''        phi0 = fr.field(freeform, X, Y, Z, freeform.x_tip_mm, freeform.x_rear_mm,
                        axles=None,                                  # evo_patches/sidepod_rules.py S1: no team cap
                        halo=fr.halo_span(x_ref_a, d_halo_mm))     # evo_patches/soft_cap.py v2: the halo's rule'''),
    # S2 (the main body's 28 mm half-width raised to the block's 32.5 mm) is held back until a test build shows the
    # widened body stays clear of the wheels; until then the sidepod corridor, already 32.5 mm, is the rules' width
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
        assert t.count(e["old"]) == 1, f"sidepod_rules {e['id']}: anchor not found once in {e['file']}"
        p.write_text(t.replace(e["old"], e["new"], 1))
        done.append(f"{e['id']}: applied")
    return "sidepod_rules: " + "; ".join(done)
