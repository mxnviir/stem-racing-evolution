# Drag audit b2: the original hybrid's final design, legalised (2026-11-28)

`legal:hybrid_final` = `hybrid/final_design/car_half.stl` through the legal builder (`seedhybrid/legalise.py`): accepted;
the rules changed 1,382 mm3 inside the region and reset 456 mm3 outside it; pod 67.8 cm3; car 48.30 g. Audit protocol,
frames 0-3. Raw: `audit/results/b2/`.

| frame | legal hybrid | CAD team pod (b1) | legal team pod (seed evolution's rebuild) | unlegalised hybrid (b1) |
|---|---|---|---|---|
| k0 | 174.95 | 180.97 | 174.18 | 171.42 |
| k1 | 173.79 | 182.86 | 175.26 | 171.23 |
| k2 | 174.73 | 181.03 | 174.63 | 171.82 |
| k3 | 173.82 | 182.92 | 174.94 | 171.76 |
| mean | **174.32** | 181.95 | 174.75 | 171.56 |

| comparison (paired by frame) | mean | SE | in SE |
|---|---|---|---|
| legal hybrid - CAD team pod | **-7.62** | 0.85 | 9.0 |
| legal hybrid - legal team pod | **-0.43** | 0.53 | 0.8 |
| legal hybrid - unlegalised hybrid | +2.77 | 0.31 | 8.9 |

Body and wheels (legal hybrid minus unlegalised hybrid): body -9.94 +- 0.31, wheels +12.70 +- 0.18 mN. The wheel-drag
reduction that made up two thirds of the hybrid's gain (b1: -6.85 mN) is gone after legalisation.

## Verdict on mark T5, and why the mark was wrong

- **As written** (legal hybrid - CAD team pod <= -7.27 mN): **met**, -7.62 mN.
- **The mark compared with the wrong reference.** The legal builder by itself moves the team pod from 181.95 (CAD
  surface) to 174.75 mN (its legal rebuild): -7.2 mN with no optimisation at all. Against the team pod built the same way,
  the legal hybrid is -0.43 +- 0.53 mN: **no measurable gain**. The fair reading: **the original hybrid's shape gain does not
  survive legalisation.** The mark should have been written against the legally rebuilt reference; it is reported as
  written and as corrected.

## What this changes

- Small morph-scale moves (peak 2 mm, RMS 0.42 mm) are not preserved by a 1 mm-grid, 6.35 mm-ball legal builder, or the
  gain lived in features the rules forbid. Batch b3 separates the two.
- For the seeded sidepod campaign (`seedhybrid/`): its confirm stage (PROTOCOL.md E4) is now expected to show a large
  legalisation loss; the in-loop region taper alone does not make moves millable.
- The whole-car pipeline (`carhybrid/`) is not exposed in the same way: every candidate is built by the legal builder
  before it is scored (legal by construction).

## Batch b3 (launched 2026-11-28, marks fixed before results)

All on frames 0-3, the audit protocol, each through the identical path:
- `hybrid_region`: the hybrid with every move outside the legal region removed (89.5 % of the move magnitude is inside the
  region), not legalised.
- `legal:hybrid_region`: the same, legalised.
- `legal:team_rm`: the team pod surface the hybrid started from, legalised by the same `legalise.py` path (the fair
  reference for both legal designs).

Readings (paired by frame):
- **R1 (did the out-of-region moves carry the gain?)** hybrid_region - CAD team pod vs the full hybrid's -10.39: if the
  region-only design keeps >= 70 % of it, the out-of-region moves were not the main source.
- **R2 (does legalisation erase in-region moves?)** legal:hybrid_region - legal:team_rm: a gain counts if <= -1.4 mN and
  >= 3 SE.
- **R3 (is the legal-builder path the same as seed evolution's rebuild?)** legal:team_rm vs seed evolution's legal team pod
  (174.75): reported.
