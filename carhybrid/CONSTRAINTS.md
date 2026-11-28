# Rules as constraints on the stage-2 numbers (carhybrid/constraints.py)

Owner, 2026-11-28: the search must not invent surface that breaks a rule. Applying the rules after the numbers makes
the adjoint and evolution push into regions the rules then clip: dead moves, wrong predictions, flats. Guided CMA is
on hold until this is in.

## Method

The body is the solid between a top and a bottom surface, inside a side outline, blended by a smooth maximum
(`freeform.field`).
- All three surfaces are linear in the numbers u: cubic splines of control grids, which are the seed plus
  `from_unit`'s interpolation of u.
- The one exception is the side grid's clamp at 0, which becomes a constraint (side grid ≥ 0).
- The smooth maximum only removes material, so a rule that keeps the sharp body out of a region holds for the built
  body too.

`build()` takes the builder's own masks, built as `param_body.build` builds them (`build_unified_geometry`, plus
`build_rule_masks`, plus the team's sidepod and halo 27.175 mm cap from `freeform.field`), and classifies them.

| Class | What it covers | Becomes |
|---|---|---|
| Required solid | T4.2 cargo, marking room, axle supports, the T5.5 3 mm zone round the chamber, the halo–canister loft, nose stub, pocket bridge | top ≥, bottom ≤ in each column; side ≥ in each row |
| Clipping air | Rule envelope / block, T7.9 zones, T7.2 wheel keep-clear, T4.4.2 view air, T8.5.1 nose box, 27.175 mm caps | the body keeps to the column's allowed interval, or the column is empty (top < bottom) |
| Fixed cuts | Halo pocket and seat, bore, axle holes, wheel discs, the block's underside slot | not constraints: always cut; the gradient must hold that surface |

Margins: half a build cell + 0.05 mm for air; + 0.8 mm for solid (the blend). Where solid and air meet in a column,
the surface is set on the boundary.

Functions:
- `project()`: the nearest legal u, by ADMM.
- `feasibility()` and `split_sharp()`: linear programs (HiGHS), described under Findings.
- `check()`: the builder's own field on its own grid at u. It counts the clipping-air cells filled, the
  required-solid cells left empty, and the fixed-cut cells filled.

## Findings (team car, builder e9c896b + evo_patches)

1. **The seed itself breaks the rules (u = 0).**
   - The builder check finds 263 cells of clipping air filled and 5,919 cells of required solid empty.
   - The builder's masks have been quietly fixing that on every car.
   - This is why the rendered cars show flats and plateaus.
2. **No set of numbers meets every rule.** The smallest achievable worst violation (min-max LP) is 12.8 mm on the
   45-number grid and 11.3 mm on the 720-number grid. Some rule features are sharper than a height-map surface on
   these grids can draw. For example, an empty wheel-zone column lies millimetres from a column that must hold the
   halo loft 40 mm tall, while the fine grid's control points across the car are 4.5 mm apart.
3. **What the constraint set will be.** Split by an L1 linear program (`split_sharp`):
   - rows the numbers can meet become hard constraints, so the search never enters them;
   - rows they cannot meet are SHARP features: the builder's mask stays in charge there, and the gradient holds that
     surface (the held-surface fix the other agent found for `grad.py`).
   - **Result (720 numbers, team car):** 1,823 of 26,481 rows (7 %) are sharp; 24,658 can be met together (the
     L1 program finds a u meeting all of them). By family, the sharp rows are: air above the top 597, top over
     required solid 610, empty air columns 454 (the wheel zones, x 72-257), required side 95 (x 94-134, beside the
     cargo), air below 36, bottom under solid 31.
   - **Open: the projection.** `project()`'s ADMM does not converge on this set: it is 203 mm off after 4,000
     iterations, ill-conditioned with 720 numbers and 24 k rows, although the set is feasible. Next: a proper QP
     solver (OSQP, or HiGHS' QP), started from the L1 program's feasible u. Then `check()` of projected numbers and
     real builds.

## Not covered (stays in the builder, measured)

- The milling's closing of concavities under the 6.35 mm ball.
- The mass skin offset: 0 on every sweep5 car, an assumption that `check()` notes.
- The parts (wings, supports): their engines build legal by construction.

## Still to do

1. Finish the split. Then:
   - project the seed;
   - confirm with `check()` that projected numbers leave the masks nothing to do outside the sharp features;
   - build real cars from projected numbers.
2. In the builder, with projected numbers: switch off the 27.175 mm soft compression. The constraint keeps the top
   below the cap, so the compression would only distort it.
3. Wire `project()` into the adjoint step (projected gradient) and the samples of any CMA, and the held-surface masks
   into `grad.py`.
