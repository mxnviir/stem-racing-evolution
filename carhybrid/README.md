# Whole-car hybrid (carhybrid)

The adjoint-led hybrid loop on whole cars built the stem-racing-evolution way, scored by race time.
**Design, decisions and the checks each phase must pass: [DESIGN.md](DESIGN.md).** What was tried before and why it
was dropped: [../docs/METHODS_TRIED.md](../docs/METHODS_TRIED.md).

## The car it builds (pinned: `claude/loving-keller-q72mvv` at `e9c896b`, read only)

- **Body and nose, no cone.** Three smooth free-form surfaces (top, bottom, side outline) fitted once to the team car's
  body (`search/layout/organic_seed.json`), stretched to each car's wheelbase plus tail, edges rounded by a smooth
  blend. The nose is grown from the body's own Ref A section and ends 2 mm ahead of the front wing
  (`car/body/freeform.py` `ref_a_section`, `seed_nose_from_section`; `search/layout/layout_genes.organic_freeform`).
- **Front wing**: straight, parallel NACA 0015, 15.3 mm chord, shape fixed; moved: one piece or two, height 6.5-18 mm,
  gap ahead of Ref A (within T8.2's 40 mm), root and tip position. Held by the nose it passes through, a strut into the
  nose (split halves; the span-return check removed for the front wing), or an underside pylon into the nose inside
  T8.5.1's 15 mm band.
- **Rear wing**: the rear-wing designer; moved: position behind Ref B (T9.4.2), height above T9.6 / T9.7's lowest,
  half span (<= 42.5 mm width), chord; straight, one piece. Supports: tip struts or an underside pylon at `rwsup_root`
  of the half span, with the fallback tip -> 80 % -> 60 % -> 40 %; T9.6 measured from where each support joins.
- **Strut generator**: curved path through a free midpoint, aerofoil-to-strut section, swell, flared foot, one valid
  solid for any path, thickened until ABS holds 20 MPa; start face in any plane (underside pylons).
- **Stage 1 numbers (12)**: wheelbase, rear axle to tail, halo to tail; front wing pieces, height, gap, root, tip; rear
  wing chord, half span, position, height. Fixed: cargo position and orientation, nose length, profile and tip height.
- **Rule fixes**: the T9.6 joint finds a pylon top buried in the wing (team car 5.51 mm, not 4.51); block-check specks
  under 1 mm3 waived where the free-form floor grazes the block.

## What this folder adds

| File | Job |
|---|---|
| `make_car.py` | build one car (layout + stage-2 free-form numbers) with the pinned builder; export its CFD solid; adapt it |
| `adapt.py` | evolution CFD export (one closed solid + wheels) -> OpenFOAM half-car folder (Part 2 format) |
| `grad.py` | DAFoam's surface drag gradient -> race-time gradient in the 45 free-form numbers (chain rule through the field function; mass from the same normal speeds; race-model partials at this car) |
| `request.json` | what the workflow runs next (pushing a change to it starts a run) |
| `../.github/workflows/carhybrid.yml` | per case: build -> adapt -> score (audit protocol) or DAFoam -> gradient map |
| `results/<batch>/` | every case's record, logs and gradient (committed by the workflow) |

## Status

| Phase | State |
|---|---|
| P0 plumbing | Local: team layout builds legal (9 min), adapter makes one closed half car (x 16.5-270.2 mm), gradient mapping runs on it with a synthetic sensitivity (1,930 of 4,000 surface points on the body; race-model slopes 0.178 ms/mN, 7.52 ms/g, car 53.8 g with no lead). On runners: run 46117504777 (score frames 0-1, DAFoam) |
| P1-P4 | after P0's checks (DESIGN.md) |
