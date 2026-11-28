# Seeded hybrid: results

Written by `analyze.py` after every stage. Marks: [PROTOCOL.md](PROTOCOL.md). **Inside-loop numbers are the steady model on each lineage's morphed base mesh: the loop's yardstick, not a claim** (claims: E3-E6, after the confirm stage).

State: stage `gen`, generation 2, setup run 46113354235.

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
| team | 1 | 1 | +185.11 | dafh -2.03, rnd1001 +1.45, rnd1002 +1.65, rnd1000 +1.81 | dafh |
| d7 | 1 | 1 | +182.83 | dafh -1.28, rnd1001 -1.03, daf -0.45, rnd1002 +0.77, rnd1000 +1.47 | dafh |
| window | 1 | 1 | +185.36 | dafh -1.69, daf -0.45, rnd1000 +1.15, rnd1002 +1.20, rnd1001 +1.49 | dafh |
| team_noadj | 1 | 1 | +185.11 | rnd1004 +1.16, rnd1001 +1.45, rnd1002 +1.65, rnd1000 +1.81, rnd1003 +10.73 | - |
| d7_noadj | 1 | 1 | +182.83 | rnd1001 -1.03, rnd1002 +0.77, rnd1000 +1.47, rnd1004 +1.99, rnd1003 +7.21 | rnd1001 |

## E1. Direction credit (guided lineages; paired within each generation)

| direction | generations it won | mean change per mm RMS (95 % interval over generations) | children |
|---|---|---|---|
| daf | 3 | -6.56 (-12.75 to -4.31) mN/mm | 5 |
| comb | 0 | - | 0 |
| es | 0 | - | 0 |
| rnd | 0 | +2.28 (+0.79 to +3.38) mN/mm | 9 |

- daf minus random, mean change per mm RMS: -8.84 mN/mm (95 % interval -16.13 to -5.10, 3 generations)

- **E1a**: `daf` below the median random child in 2/2 generations (100 %; mark >= 70 %); interval of the difference below 0: True. Verdict: **pass**.
- **E1b**: DAFoam's predicted sign matches the measured sign on 9/13 clear children (69 %; mark >= 80 %); Pearson r per generation: +0.97, +0.41, +0.74.

**Accepted steps by direction type** (inside-loop change, mN):

- team: daf -2.03; total -2.03
- d7: daf -1.28; total -1.28
- window: daf -1.69; total -1.69
- d7_noadj: rnd -1.03; total -1.03

## E2. Ablation pairs (inside-loop trajectory; the verdict needs the confirm stage)

- **team** (guided) vs **team_noadj** (no adjoint), cumulative accepted change after each generation: g1: -2.03 vs +0.00
- **d7** (guided) vs **d7_noadj** (no adjoint), cumulative accepted change after each generation: g1: -1.28 vs -1.03

## E7. Cost so far (CFD jobs)

| lineage | children solved | DAFoam gradients |
|---|---|---|
| team | 6 | 1 |
| d7 | 6 | 1 |
| window | 6 | 1 |
| team_noadj | 6 | 0 |
| d7_noadj | 6 | 0 |
| window_noadj | 0 | 0 |

Setup: one base state (fresh mesh, ~4.5 h) per seed, shared by an ablation arm and its guided twin.

## Claims (E2-E6)

Not yet: the confirm stage runs when every lineage has stopped.
