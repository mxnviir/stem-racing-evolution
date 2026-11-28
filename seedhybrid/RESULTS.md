# Seeded hybrid: results

Written by `analyze.py` after every stage. Marks: [PROTOCOL.md](PROTOCOL.md). **Inside-loop numbers are the steady model on each lineage's morphed base mesh: the loop's yardstick, not a claim** (claims: E3-E6, after the confirm stage).

State: stage `gen`, generation 1, setup run 46113354235.

## E0. Setup

| lineage | arm | status | base drag (steady, mN) | steady_pass | DAFoam at seed |
|---|---|---|---|---|---|
| team | guided (base: team) | active | 185.11 | True | converged |
| d7 | guided (base: d7) | active | 182.83 | False | converged |
| window | guided (base: window) | active | 185.33 | True | converged |
| team_noadj | no adjoint (base: team) | active | 185.11 | True | - |
| d7_noadj | no adjoint (base: d7) | active | 182.83 | False | - |
| window_noadj | no adjoint (base: window) | skipped_by_rule | 185.33 | True | - |

## Generations (inside-loop yardstick)

| lineage | gen | amp mm | parent re-run mN | children: change vs re-run (mN) | accepted |
|---|---|---|---|---|---|

## E1. Direction credit (guided lineages; paired within each generation)

| direction | generations it won | mean change per mm RMS (95 % interval over generations) | children |
|---|---|---|---|
| daf | 0 | - | 0 |
| comb | 0 | - | 0 |
| es | 0 | - | 0 |
| rnd | 0 | - | 0 |

## E2. Ablation pairs (inside-loop trajectory; the verdict needs the confirm stage)

- **team** (guided) vs **team_noadj** (no adjoint), cumulative accepted change after each generation: 
- **d7** (guided) vs **d7_noadj** (no adjoint), cumulative accepted change after each generation: 

## E7. Cost so far (CFD jobs)

| lineage | children solved | DAFoam gradients |
|---|---|---|
| team | 0 | 1 |
| d7 | 0 | 1 |
| window | 0 | 1 |
| team_noadj | 0 | 0 |
| d7_noadj | 0 | 0 |
| window_noadj | 0 | 0 |

Setup: one base state (fresh mesh, ~4.5 h) per seed, shared by an ablation arm and its guided twin.

## Claims (E2-E6)

Not yet: the confirm stage runs when every lineage has stopped.
