# sweep3: the 32 readiness cars on the full patch set (2026-11-28)

Run 46132129860: sweep2's 32 layouts rebuilt with every `evo_patches` fix (mass bookkeeping P1-P3, soft cap v2, the
multi-material nose). Renders of all 32: `GALLERY.png` (side and front three-quarter views; full renders per car in
`sw_*/render.png`).

## Legality: 11 / 32 (sweep2: 17 / 32), a regression

Six cars legal in sweep2 are not now (L01, L03, L11, L23, S03, S05); none gained. Causes, all traced:

| Failure | Cars | Cause | Fix |
|---|---|---|---|
| block_body_inside_official_block (new part) | L01, L11, L23, S03 (+ L13, L19 already illegal) | the nose split puts a ~1 mm full-section slab of foam ahead of the referee's block front (e.g. L01 179 mm3 at x 32.5-33.3). The plan takes the block's rear face as designed (x_front + W + rear); the referee takes the measured bore end + 52 mm, ~1.8 mm further back | nose_material.py (its author agent, in progress) |
| team_halo_face_joined | L03, L19 | soft cap v2's halo-region compression faded out over 5 mm BEHIND the halo's rear face, lowering the body that must back the halo | soft_cap v2.1 (commit 5951321): fades out before the rear face; the field behind it is now identical to the uncompressed one |
| P4.machining_rear_slot_cut_from_below | S05 | not yet traced | with the nose fix |
| T4.4.4_halo_fits_pocket | L19 | 5 mm3 of foam 0.28 mm above the pocket floor; not yet traced | with the nose fix |

The other failures are the same as sweep2's (T5.6 cartridge view 7, front-wing span 2, claw clearance 1, rear loads 1) and
the ~330 mm3 at the bore that every car shows (the 1 mm grid's bore, waived as before).

## The nose: what the program decided

x_cut decided by: wall thickness ("thin") on 19 cars, the official block on 9, the wing's attachment on 4. On most layouts
the nose ahead of Ref A is too slender for a 3 mm foam wall, so it stays printed (foam 0-0.7 g). Car mass changed by
-0.20 g on average (median +0.07, best -2.36 g): the foam nose pays on short, blunt noses, not on the slender ones most
layouts draw.

## Shape

- The flat deck is gone where it was the field's own plateau. The sidepod-cap flat is 0 mm2 on 18 cars, at most 211
  (L17); sweep2 had 126-660 mm2 on every car.
- What remains flat at the halo is the rim ring the body must meet flush (labelled `halo rim` in the audit, 100-506 mm2).
- "Unexplained" flats in these shape.json files (480-1260 mm2) were measured with the old audit. `carhybrid/shape_audit.py`
  now names the features they are: the nose's glued face and tenon pocket, the official block's 6 x 6 mm underside slot,
  the halo's rear seat and loft deck, the ballast container, the chamber end (measured ~50 mm from the rear face, not
  52), wheel-zone faces. On a rebuilt team car the unexplained area drops from 869 to 51 mm2. The next sweep reports
  with it.
- Section rectangularity (area / bounding box) median 0.74-0.76, max 0.94 (L03, one station): no box sections elsewhere.
