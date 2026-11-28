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

- **AWS EC2, driven by the AWS CLI** (`carhybrid/aws/`, README there): `launch.sh <batch>` ships this repository at
  HEAD plus the pinned builder and adjoint code to S3, starts one instance (spot by default; m7a.32xlarge in Ohio: 128 cores, 512 GiB; aws/README.md section 5),
  which installs exactly what the GitHub workflow installs, runs the batch's `fine_adjoint` cases one at a time with
  `carhybrid/fine_adjoint.sh` (build -> production mesh and solve, the score, frame k -> DAFoam from that converged
  state: same mesh, schemes, BCs -> gradient mapped to the free-form numbers), copies each case's results to S3 as it
  ends and terminates itself. `status.sh` / `fetch.sh` read and fetch; results are committed by hand from `fetch.sh`'s
  folder (`carhybrid/results/<batch>/`).
- GitHub-hosted runners (free, up to ~20 at once): every build, every score, the line-search candidates, and
  `fine_adjoint` only as a plumbing test on a coarse bg_div (16 GB).

## Interfaces (fixed so the workflow and the loop agree)

Request case (carhybrid/request.json):
`{"name", "type": "fine_adjoint", "layout", "u" (json list or "zero"), "k" (frame), "bg_div" (12.889 = production for
the whole car), "grid" ("5,2,3,4" default | "24,9,12,9" = the seed's own grid), "np" (MPI ranks, default all cores)}`;
on GitHub a `fine_adjoint` case is a plumbing test (coarse bg_div), on AWS (`aws/launch.sh <batch>`) the real thing

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
1. `grad` stage: one `fine_adjoint` job on AWS (`aws/launch.sh`) at u_i (frame k_i = i mod 4) -> g_i (race time per unit u).
2. `line` stage (hosted runners): candidates u_i + a * d_i for a in {0.5, 1, 2} x the trust radius, d_i the projected
   descent direction (unit box), plus u_i itself, each built and scored on two frames. A candidate the builder or the
   referee rejects counts as failed.
3. Accept the best legal candidate if its paired race-time change is <= -3 SE and <= -0.1 ms; the trust radius grows x1.5
   when the largest step wins and halves when nothing is accepted. Stop at radius < 0.02 or 20 iterations.
The no-adjoint control arm (same budget, random unit directions in place of d_i) runs alongside for the claim.

Running it with AWS: the `grad` stage's batch (`loop_adjoint_<i>_grad`, production bg_div) is skipped by the GitHub
workflow and launched by hand: `aws/launch.sh loop_adjoint_<i>_grad`, `aws/fetch.sh <run_id>`, then
`python carhybrid/loop.py advance carhybrid/results/loop_adjoint_<i>_grad <run_id> --evo <builder checkout>`, commit
and push (the `line` stage it writes runs on hosted runners, and their collect hook advances the loop as before).

### Decisions (carhybrid/loop.py; the text above does not fix these)

1. **Initial trust radius 0.25** (P1's h: 1.5 mm at 6 mm per unit), capped at **1.0** when it grows (a 2 x radius step
   then spans the whole unit box; larger ones would only be clipped). `init --radius` changes the start.
2. **The line uses the frame 0 and 1 scores of u_i re-run in the same batch** ("base"), not the grad stage's own score
   of u_i (frame i mod 4): the pairing needs the same frames for every candidate.
3. **After a rejected line, the next grad stage still runs** (at the same u, the next frame, as the text's k_i = i mod
   4) and **the line direction is the mean of every gradient taken at that u** (frames of the same derivative, as
   grad.py averages frames); the list resets when a candidate is accepted. Each grad record keeps the cosine with the
   previous gradient at that u (the fine gradient's frame-to-frame agreement, measured for free).
4. **Projection on the box**: d_i = -g / |g| with the components that push out of an active face (u_j = +-1) dropped,
   then unit length; candidates are clipped to [-1, 1]. No feasible direction stops the arm.
5. **Race time and SE**: whole car D20 = 2 drag_mN / 1e3, L20 = 2 lift_mN / 1e3 (lift + up, no sign change:
   grad.py's docstring); the block SE in race time is dT/dD20 x 2 x drag_block_se_mN (the score has no lift block
   SE; lift's race-time weight is ~1 % of drag's per N: dT/dL20 = -2.0 ms/N against dT/dD20 = +185 ms/N on the team
   car). For K frames the frame part is std(Delta_k) / sqrt(K) (= |Delta_0 - Delta_1| / 2 for two).
6. **A missing base** (either frame) cannot be judged against: the same line stage is issued once more (batch
   `loop_<arm>_<i>_line_r1`); missing twice counts as nothing accepted.
7. **Gates on every grad stage**: no `grad.json`, a non-finite or wrong-length `dT_du_s`, or a gate known to be
   missed on that job (`fine_adjoint.json` `F1_pass` / `F2_pass` / `F3_pass` false; `steady.json` `"pass": false`;
   `dafoam_sens.json` `"gmres": {f: {"converged": false}}`) stops the arm (exit 3); a gate that is unknown (None) does
   not. `loop.py reissue <arm>` re-issues a pending stage by hand (e.g. after a runner fault).
8. **The random arm's directions**: unit normal vectors from numpy `default_rng([seed, i])` (seed 20261128 by
   default), recorded per iteration in state.json; same accept, radius and stop rules, no grad stage.
9. **Both arms in one state file** (`carhybrid/loop/state.json`, keyed by arm); `advance` finds the arm by the batch
   name of the results folder. Each advance sets `launch` to that arm's next batch, so two arms running at once rely
   on the collect hook dispatching the batch `advance` wrote (as `.github/workflows/carhybrid.yml` does), not on
   `launch` alone.
