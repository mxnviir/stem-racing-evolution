# sweep5: the 32 readiness cars with evo_patches/legality_fixes.py (2026-11-28)

Run 46153774023: sweep4's 32 layouts rebuilt with legality fixes L1-L4. Renders: `GALLERY.png`.

## Legality: 28 / 32 (sweep4: 16 / 32); no car legal in sweep4 became illegal

| sweep4 failure | cars | sweep5 |
|---|---|---|
| T5.6 cartridge seen all round (rear wing behind the body's back) | L00 L02 L06 L08 L14 L16 L22 | all pass T5.6 (L1); L06 and L16 then had no rear-wing support that could land ahead of the back (rw_supported): L5 |
| block_body_inside_official_block | L07 S02 S07 | pass (L2) |
| front / rear support flush with the underside | L07 L12 | pass (L3) |
| foam over the rear support's slot | L13 L21 | L13 passes; L21 short by 0.24 mm: L3 v2 (deeper saddles) |
| T8.6.1 front span | L19 L20 | pass (L4) |
| rear slot cut from below | S05 | passes (L3) |
| team_halo_face_joined | S02 | still fails (also in sweep2, before any patch): open |

The four left are re-run as sweep5r with L5 and L3 v2.
