# The fine-mesh adjoint: the whole-car hybrid built around DAFoam on the scoring mesh

Owner, 2026-11-28: a 512 GB AWS machine will run DAFoam on the production mesh; make the whole pipeline work with it.

## Why the pipeline changes

Until now the adjoint ran on a coarse mesh (~105 k cells, first-order upwind, wall functions) while every car is judged
on the production mesh (2.1 M cells, wall-resolved, linearUpwind, fresh mesh per frame). The adjoint repo measured
what that mismatch costs: a coarse adjoint design worth -10 mN came out at -1.5 to -3 mN one mesh finer, and coarse
step sizes are overstated about 2x. With the adjoint on **the same mesh, schemes and boundary conditions as the
score**, that mismatch is gone, and three things become worth doing that were not:

1. **More design freedom.** One adjoint solve gives the derivative for every number at once, so the search can use
   the seed's own control grid (720 points: top 24 x 9, bottom 24 x 9, side 24 x 12) instead of 45 coarse numbers.
2. **Lift as well as drag.** The race model takes D20 and L20; a second function in the same DAFoam run gives dL/dX.
   Mass is exact and needs no CFD (grad.py integrates the field).
3. **Gradient steps instead of random children.** A trusted gradient makes a line search along it the main move;
   random directions remain only as the control arm the claim needs.

What it does not change (and the design keeps): frame-to-frame scoring noise (1.75 mN half car on the whole car,
2.3 mN on the fine sidepod mesh), so every accept is paired by frame and judged in SE; legality is not differentiable,
so every step is built by the evolution builder and refereed; a car with no lead is ~90 % mass in race time (P0),
so the drag/lift part matters most once a car carries lead.

## Memory and mesh size (adjoint repo, fine_campaign/README.md)

DAFoam's Krylov adjoint needs ~150 kB per cell (Jacobian colouring and preconditioner). 512 GB therefore fits about
3 M cells; the 2.1 M production mesh (~315 GB) fits with room; 5 M does not (~750 GB). The production mesh is the
target: the adjoint must be on the mesh the score uses. F0 below measures the law before anything relies on it.

## Where things run

- `bigmem`: the AWS machine, registered as a self-hosted GitHub Actions runner (labels `self-hosted`, `bigmem`). Runs
  only `fine_adjoint` jobs: build -> production mesh and solve (the score, frame k) -> DAFoam from that converged state
  (same mesh, schemes, BCs; `dafoam/make_case.py`) -> gradient mapped to the free-form numbers. One job at a time.
- GitHub-hosted runners (free, up to ~20 at once): every build, every score, the line-search candidates.
- Setup of the machine: `carhybrid/aws/README.md`.

## Interfaces (fixed so the workflow and the loop agree)

Request case (carhybrid/request.json):
`{"name", "type": "fine_adjoint", "layout", "u" (json list or "zero"), "k" (frame), "bg_div" (12.889 = production for
the whole car), "grid" ("5,2,3,4" default | "24,9,12,9" = the seed's own grid), "np" (MPI ranks, default all cores),
"runner" ("bigmem" default | "hosted" for a plumbing test on a coarse bg_div)}`

Outputs in `out/<name>/` (committed by collect to `carhybrid/results/<batch>/<name>/`):
- `result.json`: the audit-protocol score of this frame (drag_mN, lift_mN, block SE, cells, y+), as for `score` cases;
- `steady.json`: the primal's steadiness (100-iteration block means over the last iterations, span; mark A1);
- `dafoam_sens.npz`: `Xs` (m), `dD_dXs`, `dL_dXs` (N per m, half car), `D_N`, `L_N`; `dafoam_sens.json` (summary,
  GMRES convergence for each function);
- `resources.json`: cells, ranks, peak memory (GB, whole machine), primal / adjoint wall times;
- `grad.json`: carhybrid/grad.py's map (dD/du, dL/du, dm/du, dT/du, the race-model partials at this car's own D20 and L20).

## Gates (fixed before any fine result)

| Gate | Question | Passes if |
|---|---|---|
| F0 | the memory law | peak memory at three mesh sizes (bg_div ladder, team car) fits a line in cells; extrapolated 2.1 M peak < 450 GB |
| F1 (BC) | DAFoam's primal is the pipeline's flow | DAFoam's primal D within 2 % of the pipeline's drag on the same state |
| F2 (A1) | the flow is steady enough to linearise | 100-iteration block means span <= 0.10 mN over the last 500 iterations (half car) |
| F3 | the adjoint solved | GMRES reaches its tolerance for D (and L) |
| P1F | the fine gradient predicts | carhybrid/P1.md's marks, with this gradient, frames 0-3 (four, not two: the fine gradient's steps are smaller) |

A gate missed stops the line it gates and is reported; nothing downstream runs on a failed gate.

## The loop (P2F; carhybrid/loop.py)

Per iteration i, at the current numbers u_i:
1. `grad` stage: one `fine_adjoint` job on bigmem at u_i (frame k_i = i mod 4) -> g_i (race time per unit u).
2. `line` stage (hosted runners): candidates u_i + a * d_i for a in {0.5, 1, 2} x the trust radius, d_i the projected
   descent direction (unit box), plus u_i itself, each built and scored on two frames. A candidate the builder or the
   referee rejects counts as failed.
3. Accept the best legal candidate if its paired race-time change is <= -3 SE and <= -0.1 ms; the trust radius grows x1.5
   when the largest step wins and halves when nothing is accepted. Stop at radius < 0.02 or 20 iterations.
The no-adjoint control arm (same budget, random unit directions in place of d_i) runs alongside for the claim.
