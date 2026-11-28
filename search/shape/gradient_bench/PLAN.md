# Gradient-use bench: how stage 2 should use a fine-mesh adjoint (registered plan)

Owner 2026-11-28: "Run your tests and think about it and whether which one will be better for our use case."
Written and committed **before any held-out run**. Bench only, no CFD. Code and full spec:
`bench_gradient_use.py` (docstring). Amendments go below, dated, before the results they touch.

## The question
Stage 2 (`search/shape/simscale_shape.py`) is CMA-ES (population 20) over the 45 coarse numbers that move the body's
top / bottom / side control grids (`car/body/freeform.from_unit`, COARSE = (5, 2, 3, 4)). One evaluation = one CFD run
(2.1M cells), mesh-to-mesh sd ~2.1 mN against search gains of a few mN. A DAFoam discrete adjoint on 2.1M cells
(AWS, on a steady continuation, so not exactly the production model) can give a gradient. At **equal cost**, which use
is better?
- **P (trust-region polish):** CMA as now for most of the budget, then gradient steps from CMA's result, each step
  checked by a forward run (accept only if the measured drag falls by more than the noise bar, else halve the step;
  a new adjoint at each accepted design).
- **G (gradient-guided CMA):** an adjoint at CMA's mean every generation; 4 of the 20 samples are put on the
  adjoint's line (a fan of step sizes), all 20 forward-run, CMA's own selection judges them.
- **C (plain CMA)**, the control. Reported, not marked: **GP** (G, then P's polish), **G2** (G with an adjoint every
  other generation), **Glag** (G with a one-generation-old adjoint taken in the same round as the batch).

## The bench (frozen; details in the docstring)
- **Landscape (45-dim, x in [-2, 2]^45, x = 0 the stage-1 car, units mN):** true f = rotated, ill-conditioned quadratic
  (the grids' DCT modes mixed across grids; stiff = smooth modes, condition ~2000) + 8 smooth ripple terms (10 % of D,
  mildly nonconvex) with known optimum f* = 0 inside the box. Scale: an isotropic sigma0 = 0.5 move costs 12 mN on
  average (SPREAD); f(0) - f* = D = 12 mN (less where x* is clipped to +-1.8). Optimum in modes ~ lam^-0.25 (the seed
  is closer to optimal in the stiff, smooth directions).
- **Steady model f_s:** the Hessian rotated (QR(I + 0.10 G/sqrt n)) and re-weighted (lognormal 0.15), its own optimum
  0.02 D away (in truth), same ripple.
- **Adjoint (primary, "exact"):** cosine with the TRUE gradient exactly COS = 0.75 at every call; the error part is
  half systematic (f_s's error direction at x) and half random per call (variance). Magnitude: the linear gain it
  predicts along its own direction is KAPPA = 2.5 x the true linear gain.
- **Noise:** forward sd = NOISE x 5 mN, primary NOISE = 0.40 (2.0 mN; C's realised gain at 200 runs is ~4-5 mN, so
  noise / gain ~ 40-50 %, as on the car).
- **Cost:** forward 1 unit, adjoint A = 4 units (primary). **Rounds:** 20 forward runs in parallel = 1 round; an
  adjoint = 1 round (forward runs that do not depend on it share it).
- **Budgets** b = 100, 200, 400 units, adjoints included. Every arm keeps >= 10 % of b as the hand-over reserve.

## Arms (frozen; set before any arm was run, not tuned)
All: numpy CMA-ES (standard (mu/mu_w, 20), Hansen defaults), x0 = 0, sigma0 = 0.5, box by clip + penalty. The
gradient direction is CMA-preconditioned: d = -C g, unit length in CMA's metric.
- **C:** generations of 20.
- **P:** CMA generations while units + 20 <= 0.65 b; polish from CMA's mean (2 runs + adjoint, 1 round); step
  x0 + Delta sigma d, Delta0 = 2; one forward run per step (1 round); accept if y1 < mean(x0) - 1.0 sd sqrt(1 + 1/k0);
  on accept: one re-run of x1 + a new adjoint there (1 round), Delta x2 if measured / predicted gain >= 0.75; on reject
  Delta / 2; after 4 rejections in a row the polish stops and the left units go back to CMA generations.
- **G:** each generation: adjoint at the mean (1 round, A units), then 16 CMA samples + 4 injected points
  m + sigma t d, t = 0.5, 1, 2, 4 (1 round); injected steps clipped to Mahalanobis length sqrt(n) + 2n/(n+2).
- **GP:** G while units + 20 + A <= 0.65 b, then P's polish (and CMA generations if it stops early).
- **G2:** guided generation (with adjoint) and plain generation alternate. **Glag:** described above.
- **Hand-over (all arms):** the reserve re-runs 5 candidates (the arm's current point: CMA's mean, or P/GP's polished
  design; plus the 4 lowest-mean other evaluated designs), round-robin; the lowest mean is handed over.
- **Score:** TRUE regret f(handed) - f* (lower = better).

## Held-out set and marks
- **Held-out:** instances **41-45** (never run before this commit), seeds **1-5**: 25 pairs per budget, paired by
  (instance, seed); the random streams are keyed by (instance, seed), so P's CMA phase is exactly C's first
  generations. Development: instance 40, seeds 1-3 only (below; not evidence).
- **"X beats Y at b":** X's regret lower in **>= 18 of the 25 pairs** (one-sided sign test p = 0.022).
- **Decision rule (primary setting only):**
  1. **G** if G beats P at >= 2 of the 3 budgets, P beats G at none, and G beats C at >= 2 budgets.
  2. **P** if P beats G at >= 2 budgets, G beats P at none, P beats C at >= 1 budget and C beats P at none.
  3. **Budget-dependent** if G beats P at some budgets and P beats G at others: stated per budget; stage 2 follows the
     budget nearest its planned size.
  4. Otherwise G and P are not separated: **P (tie rule**: it reuses the existing CMA unchanged and spends fewer AWS
     adjoints), unless C beats P at some budget, then **C** (no adjoint).
- **Reported, not marked:** mean (sd) regret per arm and budget; paired wins, mean difference +- SE, Wilcoxon p for
  P-C, G-C, G-P, GP-C/G/P, G2-G/C, Glag-G/C; AWS adjoints and rounds per arm; P's steps / accepts; the measured value
  of the handed design against its truth (winner's curse). **Sensitivity rows** (same 25 pairs, every arm; C re-run
  only where the noise changes): cosine 0.4 and 0.9; A = 2 and 8; noise 0.20; kappa 1.0 (no overstatement);
  `cosmodel` (adjoint = f_s's gradient + random error to 0.75 where the model is better than that, so its cosine
  falls near the optimum where the model error dominates).
- Command: `python bench_gradient_use.py run results/heldout.json --instances 41 42 43 44 45 --seeds 1 2 3 4 5
  --settings primary cos0.4 cos0.9 A2 A8 noise0.2 kappa1 cosmodel`, then `verdict results/heldout.json`.

## Development (instance 40, seeds 1-3; `results/dev_p40.json`; not evidence)
What was set on it, and only from **C runs and the model cosine** (no arm comparison): the scale (SPREAD = 12,
D = 12, x* ~ lam^-0.25 clipped to +-1.8 so the optimum is in the box and C's gain at 200 is ~4-5 mN; a first try with
x* iid gave CMA 12 mN in 3 generations and then a stall, unlike the car); the steady model's error (with a 0.25
rotation and the model optimum 0.08 D away its cosine fell to ~0.1-0.5 after a few generations, far from the 0.82 on
record near the team pod; the primary adjoint was then defined with an exact cosine, and the falling-cosine version
kept as the `cosmodel` row). Then one run of all arms, primary (mean regret at 100 / 200 / 400; f(0) = 11.7):

| C | P | G | GP | Glag | G2 |
|---|---|---|---|---|---|
| 7.65 / 7.69 / 6.51 | 6.59 / 4.44 / 4.87 | 4.99 / 3.18 / 1.70 | 4.13 / 3.36 / 1.93 | 7.15 / 4.18 / 2.54 | 6.63 / 4.89 / 3.35 |

Nothing in the arms was changed after it (G2 was added as a reported arm after it, for its AWS cost).

## Caveats known before the run
- A toy, not aerodynamics: the real landscape's conditioning, nonconvexity and model error are unknown; the bench
  states the believed structure, and the sensitivity rows show how far the answer moves.
- The primary adjoint has the same cosine everywhere; on the car the cosine may fall near the optimum (`cosmodel`)
  and DAFoam may fail to converge on some shapes (1 of 6 converged far from the team pod): not modelled (a failed
  adjoint costs its units and gives nothing; at a failure rate q, G's adjoint cost per useful gradient is A/(1-q)).
- The noise is independent Gaussian; real mesh-to-mesh noise may be correlated between nearby shapes.

## Amendment (2026-11-28, after the held-out primary and registered rows, before any run of the row it adds)
**Post hoc, reported, not marked; the registered decision is unchanged.** On the record DAFoam converged on 1 of 6
shapes far from the team pod, and the bench gave every adjoint a gradient. Added row **`fail0.5`** (primary otherwise):
each adjoint fails with probability 0.5 (a separate random stream), its A units and its round are spent and no
gradient comes back; G / G2 / Glag then run a plain CMA generation, P stops its polish (start point kept as a
candidate) and its units go back to CMA. Same 25 pairs; written to `results/heldout_fail.json`. The code change was
checked to leave every registered result bit-identical.
