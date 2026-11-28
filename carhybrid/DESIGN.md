# Whole-car hybrid: design

Written 2026-11-28, before any run of this pipeline. The adjoint-led hybrid loop (validated on the trial car's
sidepod: `audit/RESULTS_b1.md` T1, and being tested from three seeds in `seedhybrid/`) applied to the whole car the
stem-racing-evolution way: the car builder, its rules and parts, and race time as the objective.

## What it builds on

| Piece | Source (pinned) | What it brings |
|---|---|---|
| Car builder, referee, parts, race model | `stem-racing-evolution` branch `claude/loving-keller-q72mvv` at `e9c896b` (another agent's branch; read only, never committed to) | The free-form body and organic nose (`car/body/freeform.py`, `search/layout/layout_genes.organic_freeform`, `organic_seed.json`); the front wing as a straight parallel NACA 0015 at 15.3 mm chord, one piece or two, supported through the nose, by a strut into the nose, or by an underside pylon (`front_wing_engine.py`); the rear-wing designer with tip struts or underside pylons, `rwsup_root` and the fallback 1.0 -> 0.8 -> 0.6 -> 0.4 (`rear_wing_engine.py`); the strut generator with a start face in any plane (`wing_strut.py`); 12 stage-1 numbers (`layout_params.json`); the T9.6 joint fix and the sub-1 mm3 block-speck waiver (`legality.py`) |
| OpenFOAM half-car CFD, DAFoam, morphing | `stem-racing-adjoint` at `9ab9092` (read only) | The audit-validated scoring protocol (fresh wall-resolved ~2.1 M cells, production scheme, frames), the coarse DAFoam gradient, the mesh morph |
| Hybrid generation rule, analysis | this branch: `seedhybrid/` | The acceptance rule, recombination of measured moves, per-generation credit, the confirm logic |

## Decisions, and why

1. **CFD: OpenFOAM, not SimScale.** The hybrid needs an adjoint; SimScale has none. This environment has no SimScale
   keys. The audit protocol (noise, mesh ladder, transition-model check) was validated on OpenFOAM. Consequence:
   stage-1 results timed on SimScale are not mixed with OpenFOAM numbers; stage 1 is re-timed on OpenFOAM (phase 3).
2. **Objective: race time** from the evolution race model (`simscale_layout.race_time`: this car's layout, mass with
   lead to 48.2 g, COM, drag and lift). CFD gives whole-car drag and lift (half car x 2) at 20 m/s.
3. **What the adjoint moves: stage 2's free-form body numbers, not mesh vertices.** Stage 2 moves the body and nose
   by a coarse grid of control displacements (`freeform.from_unit`: top 5 x 2, bottom 5 x 4, side 5 x 3 = 45
   numbers, 6 mm per unit). Every candidate is a vector of those numbers, so **every candidate goes through the
   builder, the milling and the referee: legal by construction**. This removes the sidepod loop's largest gap (moves
   that the milling rules later undo, `seedhybrid/PROTOCOL.md` E4).
4. **The gradient in those numbers by the chain rule.** DAFoam gives drag's sensitivity to moving each surface point,
   G_i (coarse mesh, steady model, as validated on the team pod). The free-form field phi(x; u) is an explicit
   function (smooth max of the top, bottom and side surfaces), so the normal speed of the surface for number u_j is
   V_ij = -(d phi / d u_j) / |grad phi| at each body surface point, by a finite difference of the field function
   (no building). Then dD/du_j = sum_i (G_i . n_i) V_ij A_i over points on the body and nose (|phi| < 1 mm); wings,
   struts and wheels get V = 0. Race time: dT/du = (dT/dD) dD/du + (dT/dm) dm/du, with dm/du from the same normal
   speeds (foam behind Ref A at 0.163 g/cm3, the printed nose at ABS density) and dT/dD, dT/dm by finite
   differences of the race model at this car (dT/dm is 0 while the car carries lead).
5. **Measurement: auto morph.** A child is scored by morphing the parent's mesh onto the child's built surface
   (OpenFOAM's morph path) when the change is small, falling back to a fresh mesh when it fails; the parent is
   re-run on the same path every generation (as in the sidepod loop). Claims are re-scored on fresh meshes
   (the audit protocol) in the confirm stage.
6. **Stages as stem-racing-evolution.** Stage 1 (layout, 12 numbers, Bayesian) then stage 2 (shape) - with stage 2's
   CMA-ES replaced by the hybrid generation rule, and CMA-ES kept as the comparison arm (phase 4).

## Phases, each with a check fixed now

| Phase | What | Check (must pass before the next phase) |
|---|---|---|
| P0 plumbing | Build one car (team layout) with the pinned builder; adapt its CFD export (`body_full.stl`, wheels) to an OpenFOAM half-car folder; score it on the audit protocol (frame 0); compute race time; coarse DAFoam gradient on it; map the gradient to the 45 free-form numbers | **P0a** the car builds legal and the adapted half car meshes (checkMesh OK); **P0b** drag within 10 % of the SimScale reading of the same layout if one exists (a gross-error check only: different solver and mesh); **P0c** DAFoam converges; **P0d** the mapped gradient is finite and not dominated by one number (no single |dT/du_j| above 50 % of the norm) |
| P1 gradient check | The pre-registered trust test, in the numbers: build and score the car at u = -h g/|g|, +h g/|g| and two random unit directions (h = 0.25, i.e. 1.5 mm at the 6 mm/unit scale), fresh meshes, frames 0 and 1 | **P1** along-gradient race time lower than against-gradient by >= 3 SE, and lower than both random directions; predicted vs measured slope reported (the sidepod adjoint overstated size about 2x) |
| P2 stage-2 hybrid | Hybrid generations from the stage-1 car: children = parent re-run, DAFoam step, two recombined directions, three random directions (in the numbers), plus a fresh DAFoam gradient where the parent moved; accept on >= 3 SE and >= the floor in race time | `seedhybrid/PROTOCOL.md` marks E1, E3, E5 applied to this pipeline, with drag replaced by race time |
| P3 stage 1 on OpenFOAM | The evolution layout search (GP + expected improvement over the 12 numbers, 16 Sobol cars then batches of 8, about 48 cars) timed on the audit protocol | Its best layout is re-run on frames 0-3; if it beats the layout P2 started from by >= 3 SE, P2 is re-run from it |
| P4 comparison | Stage 2 by CMA-ES (evolution's method) at the same budget from the same stage-1 car | The hybrid's legal final beats CMA-ES's by >= 3 SE in race time (frames 0-3), or the claim "hybrid is better" is not made |

P0 and P1 decide whether the adjoint carries over from the sidepod to the whole car at all; everything after them
depends on it. If P1 fails, P2 runs as random + recombination only (the adjoint dropped), and that is reported.

## Starting layout for P0-P2

The team layout (`layout_genes.x_default()`) for P0 (it builds in about 8 minutes and is the reference every search
compares against). For P2, until P3 has re-timed stage 1 on OpenFOAM: **L0039**, the fastest SimScale-timed layout
that still builds legal on the pinned builder (`docs/checks/STATUS_2026-11-27.md` on that branch: 0.75721 s
adjusted). Starting P2 from a SimScale-ranked layout is a stated shortcut; P3 decides whether it was a good one.

## Known limits (stated before results)

- The free-form numbers move only the body and nose. The wings' few stage-2 numbers (rear-wing thickness,
  `rwsup_root`) and the supports' plates are moved by random and recombined directions only (no adjoint for them).
- The gradient ignores what the masks and the milling do to the surface; every step is built and scored, so this
  costs efficiency, not correctness.
- DAFoam has converged on only some shapes (sidepod: team pod yes, 5 of 6 others no). P0c tests it on a whole car.
