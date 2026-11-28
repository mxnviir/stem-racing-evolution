# The method of this repository (2026-11-28)

**Goal (owner):** the best car, measured as race time on the audit protocol. The goal is not to beat another
method, and there is no CMA here. Each part below is marked **built**, **running** or **not yet**.

## 0. The car as numbers

A car is described by three groups of numbers. The builder (pinned `e9c896b` plus `evo_patches/`) turns them into a
manufacturable, legal car, or refuses.

| Group | Numbers | What they set |
|---|---|---|
| **Layout** | 12 genes in the unit cube (`search/layout/layout_params.json`, `mode: layout`) | wheelbase, axle positions, halo position, nose, wing positions and so on |
| **Body** | u: 45 numbers (coarse) or 720 (fine) | the free-form body and sidepods: top, bottom and side height maps, `u` = bumps on the stage-1 car's own organic seed |
| **Wing and struts** | rear wing (chord, thickness, span, x, lift, camber, camber position, rise, sweep, tip ratios; AoA fixed at 0) and the 13 strut numbers | `evo_patches/rear_wing_shape.py`, `wing_strut.py` |

**Rules are part of the description, not a filter afterwards:**
- The builder's masks carve rule air and force rule solid: T4.2 cargo, T5.5 3 mm chamber wall, T7.9 and T7.2 wheel
  zones, T8.5.1 nose box, T4.4.2–3 halo visibility, the official block (T1.14/T3.1.2: 32.5 mm each side).
- The wing engine places each part where the rules allow (W2–W7, L1–L7).
- The referee (`car/build/legality.py`, 185 checks, T9.5.1 read strictly) judges every car.
- Team-only limits are removed when the owner drops them; the 27.175 mm sidepod cap was dropped 2026-11-28
  (`sidepod_rules.py`).
- **Not yet:** the body numbers carry the rules as constraints, so the search never asks for a body the masks must
  cut (`constraints.py`, 93 % of rule rows reachable; the projection is not finished).

## 1. Every score is measured the same way (built)

**The audit protocol:**
- OpenFOAM v2412, half car, k-ω SST;
- a fresh wall-resolved 2.1 M-cell mesh per car, frame k of 4;
- 1,000 iterations; the mean of iterations 200–1,000 with block standard errors.

**Race time:** `simscale_layout.race_time(car, drag, lift)`, which uses the car's own mass and layout.

**Comparisons are paired by frame:**
- Δ = T(candidate) − T(base) on the same frames;
- SE from the block errors and the frame scatter;
- a step is accepted only if Δ ≤ −3 SE and Δ ≤ −0.1 ms.

## 2. Outer loop: where the parts go (Bayesian layout search; running)

`carhybrid/bo.py` uses the builder's own GP + expected-improvement search, fed with our OpenFOAM race times:
- start: 16 Sobol cars;
- then 4 batches of 8, chosen by expected improvement;
- then a confirmation run: the top 3 plus two references (L0039 and the team car) on all 4 frames, with the paired
  mark.

Cars that fail to build are kept as failures, and the search stays away from them.

**State:**
- bo_1 is done (15 of 16 cars timed).
- bo_2 is running. It uses the builder from before rear_wing_shape and before the cap drop, so its data describes an
  older car.

## 3. Inner loop: the shape at a fixed layout (trust-region hybrid; built, not yet run at full scale)

`carhybrid/loop.py`, one iteration at numbers u_i:

1. **Gradient (adjoint).**
   - DAFoam on the car gives dDrag/dX and dLift/dX on the surface.
   - `grad.py` chains them through the builder's surface map to dT/du: drag, lift and mass through the race model.
   - It agrees with the other agent's independent map (cosine 1.0000).
2. **Candidate steps.** u_i + a·r·d, with a = 0.5, 1 and 2, where d is the box-projected steepest descent and r the
   trust radius.
   - Each candidate is built (legal or refused) and scored on frames 0 and 1 against u_i itself.
3. **Judgement by measurement, not by the gradient.**
   - Take the best candidate that passes the paired mark.
   - The radius grows ×1.5 if the longest step won, and halves if nothing won.
   - Stop at r < 0.02 or after 20 iterations.

So the adjoint proposes and the CFD decides. A gradient that is wrong at this shape costs one batch of three
builds, not the car.

**Evolution where the adjoint cannot reach:**
- discrete choices (the support topology, one or two vanes, which joint type);
- moves the surface sensitivity does not see (the support field's growth move, `wing_field/appendage.py`).

**Not yet in `loop.py`:** today it has an adjoint arm and a random-direction control arm. The plan is to add the
evolution moves as extra candidates in the same step, judged by the same paired mark. The growth move itself is built
and prototyped (`evo_patches/wing_field/WING_FIELD.md`).

**Evidence so far:**
- **P1, built and measured:** the gradient predicts on the production mesh.
  - The race-time step gave −11.3 ms.
  - The drag step beat its opposite by 29.6 ± 2.8 mN (whole car) and beat both random steps.
- **The earlier sidepod hybrid:** −3.05 ± 0.80 mN half car against the team pod (audit b3).

## 4. Fidelity (built, not yet run at full scale)

- **Coarse DAFoam (GitHub runners):** gradients converge with GMRES at 6,000 iterations.
- **Fine DAFoam on the 2.1 M scoring mesh (AWS, spot, about 315 GB):**
  - `carhybrid/aws/` is built and tested dry;
  - it waits for the owner's AWS setup and the independent audit.

## 5. Open before compute tests (in order)

1. The rules as constraints on the body numbers: finish the projection, re-seed below the limits, switch off the
   halo-region squash.
2. `grad.py`: treat surface held by a mask as held, not as moving.
3. The rear wing and strut numbers in the search:
   - in the Bayesian search: as layout-level numbers;
   - in the inner loop: as numbers with a DAFoam gradient (dX/dθ by finite differences of the generator, no CFD).
4. Wing generator robustness: 5 of 16 random wing1 cars found no support. Each such car costs a search slot.
5. Re-judge the earlier cars under the strict T9.5.1, and decide what to do with bo_1/bo_2's older-builder data.
6. AWS: `aws_smoke` → `aws_scaling` → `f0`.
