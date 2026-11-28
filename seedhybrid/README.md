# Seeded hybrid lineages

**New here? Read [GUIDE.md](GUIDE.md) first**; the experiments and their pass marks are in
[PROTOCOL.md](PROTOCOL.md): what the workflow does, why, a diagram, a worked example, and how to
operate it. This page is the formal record (settings, pass marks, change log).

The hybrid (adjoint-led) loop of `mxnviir/stem-racing-adjoint` (`research/sidepod_opt/hybrid/`), run from several
legal seeds at once instead of only the team pod's CAD surface. Why: the drag audit (`audit/RESULTS_b1.md`) found the
hybrid's gain real on fresh meshes (-10.4 mN vs its start), three times seed evolution's, while the team pod's own
legal field rebuild already sits 7.2 mN below the CAD surface the hybrid started from. Seeds let the hybrid start from
better and different shapes.

## What runs where

- **Workflow** `.github/workflows/seedhybrid.yml`, one stage per run, chained by the collect job (`config.json`
  `"auto": true`; create `seedhybrid/STOP` to halt after the current stage).
- **CFD code**: `mxnviir/stem-racing-adjoint` at `9ab9092`, unchanged except `seedhybrid/patch_adj.py` (opt-in
  `STEP_REGION=1`: moves confined to the legal sidepod design region, `gate.region_mask`; the original loop's box
  reached x = 200 mm, 20 mm past the region's rear end).
- **Driver** `seedhybrid/driver.py`: the original loop's generation rule (`hybrid/loop.py`), per lineage.
- **State** `seedhybrid/state.json`; per lineage `seedhybrid/lineages/<name>/` (gradients, measured moves, combined
  directions, every child's result and force history).

## Stages

1. **setup** (per seed, in parallel)
   - *geom*: rebuild the seed legally (seed evolution's builder: rules, 3-axis 6.35 mm ball mill, 3 mm walls, one
     piece, cargo, clearance), remesh its surface to ~1 mm triangles, symmetry face snapped to y = 0, region taper.
   - *base*: fresh wall-resolved mesh (~2.1 M cells) of that surface, production scheme 1000 iterations, then the
     steady model (first-order upwind U) 1200 iterations: the state every child is morphed from and warm-started from.
     `steady_pass` (block span <= 0.10 mN over the last 1000 iterations) is recorded.
   - *dafoam*: coarse (105 k) steady state of the seed and its DAFoam discrete-adjoint drag gradient.
2. **generation k** (every active lineage, in parallel): children = parent + ONE move each, chained explicit morph of
   the lineage's base mesh, warm start, 1000 steady iterations, drag = mean of the last 404:
   - `zero` (parent re-run),
   - `daf` along the newest DAFoam gradient (gen 1 also at half amplitude),
   - `comb` / `es`: steepest descent in the span of generation k-1's measured moves (all / random only),
   - 3 new random smooth fields (5 if a seed has no DAFoam gradient),
   - plus a fresh DAFoam gradient wherever the parent moved.
   Accept the best child if it beats the re-run by >= 3 SE and >= 0.15 mN. No accept -> amplitude halves.
   Amplitude starts at 1.0 mm; a lineage stops at 6 generations, below 0.125 mm, or when every child folds the morph.

Seeds (`config.json`): `team` = the team pod's legal field rebuild (174.75 mN on the audit protocol), `d7` (171.44),
`window` (the windowed seed, 173.50); all from `stem-racing-adjoint` `research/sidepod_opt/seedevo/campaign_poc`.

## Marks (fixed 2026-11-28, before any result of this campaign)

- **S1 (setup)**: a lineage runs only if its base state exists; `steady_pass` and whether DAFoam converged at the seed
  are recorded per seed (the adjoint's convergence on shapes unlike the team pod is itself a result: 1 of 6 before).
- **Change log (before any generation result)**: 2026-11-28, owner: ESI removed (`use_esi: false`). The original
  loop used an ESI continuous-adjoint map as a second guide; it failed its own validation (wrong signs at 2 of 3
  wall-resolved bumps, 7/11 signs vs DAFoam's 10/11 in hybrid round 1) and computing it after each parent re-run
  added 1-2 h to every generation's slowest job. Guides are now DAFoam only, plus random and recombined moves.
- **Change log**: 2026-11-28 06:40 UTC, before any generation result: ablation arms added (PROTOCOL.md E2):
  `team_noadj` (always), `d7_noadj` (if DAFoam converged at d7), `window_noadj` (if not, and it converged at window).
  Every arm now gets the same number of children per generation (6 in generation 1, 7 after), random moves filling up.
- **G (each generation)**: the acceptance rule above. Steady-model numbers on the morphed mesh are the loop's
  internal yardstick only, never a claim.
- **F (the claim, after the lineages stop; confirm stage)**: each lineage's final design is legalised through the same
  builder (the moves are not milling-checked while the loop runs) and scored on the audit protocol (fresh meshes,
  production scheme, frames 0-3) against its own legal seed.
  - F1: a lineage's gain over its seed counts if the paired mean is <= -1.4 mN and >= 3 SE.
  - F2: the campaign's best legal design beats d7 (171.44 mN, the best legal design so far) by >= 3 SE (unpaired,
    4 frames each).
  - Reported: how much of each lineage's steady-model gain survives legalisation and fresh meshes; adjoint-led vs
    random accepted steps per lineage.

## Built since: confirm and robust stages (2026-11-28)

- **confirm** runs automatically when every lineage has stopped (PROTOCOL.md "How the stages map to the
  experiments"); `seedhybrid/legalise.py` puts a morphed design back through the legal builder; results land in
  `seedhybrid/confirm/` and `state.json` `confirm`.
- **robust** follows: the best legal final and its seed on the finer mesh and with the transition model.

## Not built yet (next)

- **re-meshing** ("phase B"): a morph reaches ~2 mm in total; when every child folds, the lineage stops. Re-meshing the
  current parent and continuing would let lineages travel further.
