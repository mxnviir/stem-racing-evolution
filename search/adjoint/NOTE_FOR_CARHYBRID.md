# Note for the whole-car hybrid (carhybrid, branch claude/wizardly-brahmagupta-1edalc)

This is from the stage-2 adjoint work on branch `claude/loving-keller-q72mvv`, 2026-11-28. The owner decided on one
fine-adjoint core with two users. They asked for this write-up of the evidence, with links.

Three things matter for `carhybrid/loop.py` and `carhybrid/grad.py`:
1. a bench result on how to use the gradient;
2. a cross-check that our two gradient maps agree;
3. a defect that both maps had, now fixed on this branch.

Everything below is committed on `claude/loving-keller-q72mvv` in `mxnviir/stem-racing-evolution`. Nothing on your
branch was changed.

## 1. How to use the gradient: guided CMA beat a line search / trust-region polish

**Evidence:**
- [`search/shape/gradient_bench/PLAN.md`](https://github.com/mxnviir/stem-racing-evolution/blob/claude/loving-keller-q72mvv/search/shape/gradient_bench/PLAN.md):
  registered before the held-out run (commit `75b5911`).
- [`RESULTS.md`](https://github.com/mxnviir/stem-racing-evolution/blob/claude/loving-keller-q72mvv/search/shape/gradient_bench/RESULTS.md):
  results, plus the full verdict tables (commits `7db95ee`, `2f9cc77`).
- [`bench_gradient_use.py`](https://github.com/mxnviir/stem-racing-evolution/blob/claude/loving-keller-q72mvv/search/shape/gradient_bench/bench_gradient_use.py):
  its docstring is the full spec. The held-out set re-runs in about 12 s.
- Raw results:
  [`results/heldout.json`](https://github.com/mxnviir/stem-racing-evolution/blob/claude/loving-keller-q72mvv/search/shape/gradient_bench/results/heldout.json)
  and
  [`results/heldout_fail.json`](https://github.com/mxnviir/stem-racing-evolution/blob/claude/loving-keller-q72mvv/search/shape/gradient_bench/results/heldout_fail.json).

**The setup:**
- a 45-number toy landscape;
- a "production" objective, plus a correlated "steady-model" one that the adjoint differentiates;
- the adjoint's cosine with the true gradient 0.75, overstating gains 2.5×;
- an adjoint costs 4 forward runs;
- noise sd 2.0 mN, about 40 % of a typical gain;
- batches of 20;
- 25 held-out pairs per budget.

**The arms:**
- **P:** CMA, then a trust-region polish from CMA's mean. Each step gets one forward run, accepted at 1 SE of the
  difference, with the step doubled or halved. This is the same family as `loop.py`'s line search.
- **G:** an adjoint at CMA's mean every generation, with 4 of the 20 samples on its line. The line is
  m + σ t d, t ∈ {0.5, 1, 2, 4}, with d = −C g / √(gᵀ C g). The samples are injected into pycma with `force=True`.
- **C:** plain CMA.

**Mean true regret in mN (lower is better; the start is 11.9):**

| Arm | 100 runs | 200 runs | 400 runs | Adjoints | Rounds |
|---|---|---|---|---|---|
| C | 7.57 | 6.84 | 6.34 | 0 | 5 / 10 / 20 |
| P | 6.37 | 5.42 | 5.05 | ~2 | 12 / 17 / 28 |
| **G** | **5.25** | **3.25** | **2.20** | 3 / 7 / 15 | 8 / 16 / 32 |

- **G against P:** G wins 21, 25 and 25 of 25 pairs. P never beats G.
- **Sensitivity:** G stays ahead in every row: cosine 0.4 and 0.9, adjoint cost 2 and 8, noise 20 %, no
  overstatement, cosine falling near the optimum, and half the adjoints failing.

**Why:** G uses the gradient every generation, and CMA's selection over 20 measured cars averages out the noise. P
spends most of its polish on rejected single-run steps: about 5.5 steps per polish, of which 0.6–1.0 were accepted.

**What this suggests for `loop.py`:** consider letting CMA's population verify the adjoint's line instead of a
paired accept test, or at least run G as an arm. Stage 2 here now does G: `search/shape/guided.py` and
`actions_shape.py`.

**Caveats:**
- It is a toy, with no CFD.
- DAFoam must converge at CMA's mean each generation, on unusual shapes.
- The noise model was independent Gaussian per run.

## 2. The two gradient maps agree

`search/adjoint/link.py` (this branch) and `carhybrid/grad.py` (your commit `3101248`) were written independently.
They use the same method: the surface's normal speed from `freeform.field`. Run on your
`carhybrid/results/p0/p0_team_dafoam/dafoam_sens.npz`:

| Measure | Result |
|---|---|
| Cosine of the 45-number dD/du | 1.0000 |
| Ratio of norms | 1.0001 |
| Largest difference | 0.06 % of the largest component |
| Surface points on the body | 2203 in each |

Record: [`search/adjoint/results/xcheck_p0_team.json`](https://github.com/mxnviir/stem-racing-evolution/blob/claude/loving-keller-q72mvv/search/adjoint/results/xcheck_p0_team.json).

## 3. A defect in both maps, found by an exact test, and the fix

**The test:**
[`search/adjoint/test_link.py`](https://github.com/mxnviir/stem-racing-evolution/blob/claude/loving-keller-q72mvv/search/adjoint/test_link.py),
results in
[`LINK_TEST.md`](https://github.com/mxnviir/stem-racing-evolution/blob/claude/loving-keller-q72mvv/search/adjoint/LINK_TEST.md).
With dD/dX replaced by the body's area vector, the map must give the volume derivative. That is measured by
rebuilding the body with `param_body.build` at u ± 0.25. The same holds for the first moments in x and z.

**What it found:** a top number came out at 1.30 / 1.36 / 1.46× the rebuilt body.
- About a third of its predicted change sat on the spine (|y| < 13 mm, z 24–43 mm), on the chamber wall and the halo
  pocket.
- The build's hard masks hold that surface in both directions; per-point displacements confirmed it.
- The seed field is fitted to the built body, so it lies within 1 mm of those surfaces as well. A "|phi| < 1 mm"
  test, the `PHI_TOL_MM` in `grad.py`, counts them as moving.

**The fix:** `link.Masks` / `link.build_masks(decoded layout)`, which takes 12 s on the 1 mm grid.
- It builds the hard masks exactly as `param_body.build` does (`unified_phi.build_unified_geometry`).
- It drops surface points within 0.6 build cells of a required-air or required-solid cell.
- With it, the top number reads 1.01–1.05×, a bottom number 0.92–0.96×, and a side number 0.89–0.93×.
- The side number's x-moment is 0.887, one cell short of the 10 % mark. That number's response is curved by 20 % over
  the test step, so the mark failed by one of nine as registered.

The same correction applies to `grad.py`'s mass term (`dm/du`), which integrates the same speeds over the body's own
surface.

## 4. What stage 2 takes from your core

Both pieces are packaged from your commit `3101248` (`search/adjoint/aws/batch.py` `CORE_REF`), not copied:
- `patch_daf.py`: lift, the per-solve GMRES record, the mesh-quality limits from the environment, and the shm/ranks
  settings. All seven patches apply cleanly to `stem-racing-adjoint` `claude/loving-keller-q72mvv` and are
  idempotent.
- `grad.py`: run beside the masked link on every job, with the cosine logged. The masked link is the map of record
  until `grad.py` excludes mask-held surface. After that, `grad.py` can be the only map.

Stage 2's job settings also follow your measured numbers: 6000 GMRES iterations and relaxed mesh-quality limits. Your
`aws_scaling` and `f0` batches (memory law, scaling) will set the instance size for both of us.
