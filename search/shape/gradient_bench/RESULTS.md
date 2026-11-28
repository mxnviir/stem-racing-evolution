# Gradient-use bench: results

Marks: `PLAN.md` (committed 75b5911 before any held-out run; one dated post-hoc amendment adding the `fail0.5` row,
committed 7db95ee before that row was run). Held-out: instances 41-45 x seeds 1-5 = 25 pairs per budget
(`results/heldout.json`, `results/heldout_fail.json`; the whole set runs in ~12 s on 4 cores). Score = TRUE regret
(mN above the known optimum) of the design each arm hands over, mean (sd) over 25; lower = better. f(0) = 11.9 mN on
average (the stage-1 car). Budgets count adjoints at A units each.

## Decision as registered: **G (gradient-guided CMA)**
G beats P at all 3 budgets (21, 25, 25 of 25), P beats G at none, and G beats C at all 3 (25, 25, 25 of 25).

### Primary (cosine 0.75, adjoint overstates the gain 2.5x, A = 4, noise sd 2.0 mN = 40 % of a typical gain)
| arm | 100 | 200 | 400 | AWS adjoints (100/200/400) | rounds (100/200/400) |
|---|---|---|---|---|---|
| C plain CMA | 7.57 (1.25) | 6.84 (1.47) | 6.34 (1.08) | 0 / 0 / 0 | 5 / 10 / 20 |
| **P** polish at the end | 6.37 (1.43) | 5.42 (1.16) | 5.05 (1.38) | 1.6 / 1.8 / 2.0 | 11.9 / 17.0 / 27.5 |
| **G** guided CMA | **5.25** (0.78) | **3.25** (0.44) | **2.20** (0.41) | 3 / 7 / 15 | 8 / 16 / 32 |
| GP (G, then P) | 4.92 (0.98) | 3.47 (0.58) | 2.63 (0.41) | 3.7 / 6.6 / 11.6 | 13.7 / 21.0 / 34.7 |
| G2 (adjoint every 2nd gen.) | 6.54 (0.94) | 4.87 (0.71) | 3.73 (0.74) | 2 / 4 / 8 | 7 / 14 / 27 |
| Glag (1-gen.-old adjoint) | 6.79 (1.43) | 4.57 (0.68) | 3.38 (0.71) | 2 / 6 / 14 | 5 / 9 / 18 |

Realised true gain from the stage-1 car (f(0) - handed): C 4.3 / 5.1 / 5.6 mN, P 5.5 / 6.5 / 6.8, G 6.6 / 8.6 / 9.7.
(Noise sd 2.0 against C's gains of 4.3-5.6 = 36-46 %, as on the car.)

Paired (X wins of 25; mean X - Y +- SE; Wilcoxon p). Marked pairs in bold:
| pair | 100 | 200 | 400 |
|---|---|---|---|
| **G vs P** | **21** (-1.12 +- 0.26, 0.001) | **25** (-2.17 +- 0.21, 1e-5) | **25** (-2.85 +- 0.28, 1e-5) |
| **G vs C** | **25** (-2.32 +- 0.25) | **25** (-3.59 +- 0.31) | **25** (-4.14 +- 0.20) |
| **P vs C** | **21** (-1.20 +- 0.26) | **22** (-1.43 +- 0.31) | **20** (-1.29 +- 0.29) |
| GP vs G | 14 (-0.33 +- 0.19) | 7 (+0.22 +- 0.12) | 3 (+0.43 +- 0.10) |
| G2 vs P | 13 (+0.17) | 16 (-0.55) | 20 (-1.33) |
| Glag vs P | 10 (+0.42) | 19 (-0.85) | 23 (-1.68) |

### Sensitivity (reported; mean regret at 100 / 200 / 400 and G-vs-P wins of 25)
| row | C | P | G | G beats P | P beats C | rule applied to the row |
|---|---|---|---|---|---|---|
| primary | 7.57 / 6.84 / 6.34 | 6.37 / 5.42 / 5.05 | 5.25 / 3.25 / 2.20 | 21 / 25 / 25 | 21 / 22 / 20 | G |
| cosine 0.4 | (same) | 7.57 / 6.46 / 6.10 | 7.21 / 4.42 / 2.92 | 17 / 25 / 25 | 12 / 15 / 14 | G |
| cosine 0.9 | (same) | 5.69 / 5.05 / 4.81 | 4.67 / 3.26 / 2.16 | 19 / 24 / 25 | 23 / 23 / 20 | G |
| A = 2 | (same) | 6.30 / 5.45 / 5.05 | 4.58 / 3.09 / 2.16 | 22 / 25 / 25 | 22 / 22 / 20 | G |
| A = 8 | (same) | 6.27 / 5.39 / 5.07 | 5.29 / 3.50 / 2.37 | 21 / 24 / 25 | 21 / 22 / 19 | G |
| noise 20 % | 6.66 / 5.98 / 5.25 | 4.85 / 4.77 / 3.95 | 4.82 / 3.02 / 1.82 | 14 / 24 / 25 | 25 / 20 / 24 | G |
| kappa 1 (no overstatement) | (same) | 6.23 / 5.42 / 5.21 | 5.09 / 3.28 / 2.34 | 20 / 25 / 25 | 21 / 22 / 18 | G |
| cosmodel (cosine falls near optimum) | (same) | 6.31 / 5.44 / 5.30 | 5.04 / 3.59 / 2.50 | 16 / 25 / 25 | 19 / 19 / 20 | G |
| fail0.5 (post hoc: half the adjoints fail) | (same) | 7.09 / 6.37 / 5.68 | 6.24 / 4.44 / 3.64 | 17 / 22 / 24 | 16 / 15 / 17 | G |

G is the row's answer in every row; P never beats G anywhere (most P wins over G in any cell: 11 of 25, noise 20 % at
100). The only cells where G is not shown better than P are at 100 units (cosine 0.4, noise 20 %, cosmodel,
fail0.5), where G has had only 3 generations.

## Reading (not marked)
- **Why G wins.** In 45 dimensions a random CMA sample has cosine ~0.15 with the downhill direction; a 0.75-cosine
  adjoint is ~5x better per unit step. G spends that information EVERY generation and lets CMA's selection over 20
  measured designs judge it, which averages out the 2.0 mN noise: four injected points on the line (0.5, 1, 2, 4 sigma)
  cost 4 of the 20 runs and pull the mean along the gradient when they measure well. P spends it once, through
  single-run accept tests against a 2.4 mN bar (1 sd of a difference): over 25 runs P made 5.4-5.8 steps and accepted
  0.6-1.0 per polish, so most of P's budget is rejected steps. P still beats C (+1.2-1.4 mN, 20-22/25) but gets less
  than half of G's improvement (G: +2.3 / +3.6 / +4.1 mN over C).
- **The 2.5x overstatement does not matter** for either (kappa 1 row: same answer), because both take step sizes
  from measurements (G's fan; P's trust radius), never from the adjoint's predicted gain. It only stops P's radius
  from growing.
- **Polish after G (GP) adds nothing**: GP loses to G at 200 and 400 (7 and 3 of 25 wins); the units the polish spends
  are worth more as guided generations.
- **Cost.** At A = 8 (adjoints twice as dear) G still wins 21 / 24 / 25 against P. G's adjoint share of the budget is
  12-15 % at A = 4 (3 x 4 / 100, 7 x 4 / 200, 15 x 4 / 400) and 24 % at A = 8.
- **Fewer AWS adjoints.** G2 (an adjoint every other generation: 2 / 4 / 8) beats C everywhere (19, 23, 25 / 25) and
  P at 400 (20 / 25), ties P at 100-200; it loses to G at every budget (1, 1, 0 / 25). Halving the adjoints costs
  about 1.3-1.6 mN of the result.
- **Wall clock (rounds).** G takes 2 rounds per generation (adjoint, then batch). Per round it is still ahead of C:
  G after 8 rounds (100 units, 5.25) beats C after 20 rounds (400 units, 6.34). Glag (adjoint in the same round as
  the batch, one generation old) halves G's rounds and is the best per round: 9 rounds -> 4.57, 18 rounds -> 3.38
  (G: 8 -> 5.25, 16 -> 3.25, 32 -> 2.20); per unit it is worse than G (0-3 / 25). P is the slowest: its polish is
  sequential (one forward run per round): 17 rounds at 200 units against C's 10. If an adjoint takes 2 rounds of
  wall clock instead of 1, G's rounds become 11 / 23 / 47 (still ahead of C at matched rounds: 23 rounds at 3.25
  against C's 20 at 6.34).
- **Winner's curse.** The handed design's measured mean reads 0-1.1 mN better than its truth in every arm (C 0.0 /
  0.5 / 0.7, G 0.8 / 1.0 / 1.0 at 100 / 200 / 400), despite the reserve's re-runs; G's handed design is as good in
  truth as its CMA mean (3.25 vs 3.17 at 200).
- **Adjoint failures (post hoc row).** With half the adjoints failing (their units and rounds spent), G still beats P
  at 200 and 400 (22, 24 / 25) and C everywhere (19, 22, 25 / 25), while P's edge over C mostly disappears
  (16, 15, 17 / 25: no failed adjoint, no polish).

## Recommendation for stage 2
**Use G: an adjoint at CMA's mean every generation, with 4 of the 20 samples placed on the adjoint's line.**
- It wins at every budget, in every sensitivity row, by 1-3 mN of true drag over P (a few mN is the size of what a
  search finds), and it is robust to the two known weaknesses of the fine-mesh adjoint: the 2.5x overstatement is
  irrelevant (step sizes come from the measured fan), and a poor cosine (0.4) or a cosine that falls near the
  optimum (`cosmodel`) still leaves G ahead from 200 runs on. A wrong gradient costs G at most 4 of 20 runs per
  generation, because CMA's selection discards line points that measure badly.
- Implementation in `simscale_shape.py` (not done here): before each `ask`, take the adjoint at `es.mean`, map it to
  the 45 numbers (chain rule through `from_unit`, the amp_mm scaling and the rule masks), and replace 4 of the 20
  asked vectors with `mean + sigma * t * d`, d = -C g in CMA's metric, t = 0.5, 1, 2, 4. pycma has
  `es.inject(...)`, which does the Mahalanobis clipping the bench used. Everything else (evaluation, `tell`, the
  hand-over) is unchanged.
- **If AWS spend or DAFoam reliability binds**, G2 (an adjoint every other generation: 4 adjoints per 200 units
  instead of 7) still beats plain CMA by ~2 mN at 200 and P at 400. **If wall clock binds**, Glag (the adjoint of
  the current mean runs in the same round as the batch and steers the next one) gives about G's result in half
  G's rounds. P (the cheapest in adjoints: ~2 per run) is the fallback only if the adjoint cannot run each
  generation at all.
- Polishing at the end, after G, is not worth its budget (GP).

## Caveats
1. A toy landscape, not aerodynamics: ill-conditioning, ripple, the seed's position and the model error are
   assumptions, calibrated on development instance 40 from plain-CMA runs and the model's cosine only. The answer
   held under every row tried, but the real 45-number landscape is unmeasured.
2. The bench's adjoint has a known cosine; on the car the fine-mesh adjoint's cosine is measured once (0.82, coarse
   DAFoam near the team pod) and DAFoam converged on 1 of 6 shapes far from the team pod. G needs a gradient at
   CMA's mean every generation, i.e. on shapes unlike the team pod. **Before committing to G, check that DAFoam
   converges on 3-4 stage-2 means**; at a 50 % failure rate G still wins on the bench, but if it rarely converges
   G degenerates to plain CMA with A units wasted per generation, while P's cost stays at ~2 adjoints.
3. Noise is independent Gaussian (sd 2.0 or 1.0); correlated mesh noise between nearby shapes is not modelled. It
   would make P's single-run accept tests and G's line points both less reliable.
4. The gradient needs the chain rule through `from_unit` and the rule masks; where a mask clips a control point the
   gradient there is zero or misleading (not modelled; the bench box is the only constraint).
5. An adjoint is counted as 1 round and A units; its real wall clock (steady continuation + DAFoam on AWS) should be
   checked; the rounds above scale as stated.
6. P was implemented as specified (one forward run per step, bar 1 sd of a difference). A P that tests a fan of step
   sizes in one parallel round would be closer to G's line points; G already is that, every generation.

---
## Full tables (`python bench_gradient_use.py verdict results/heldout.json results/heldout_fail.json`)

### setting primary: {'cos': 0.75, 'A': 4, 'noise': 0.4, 'kappa': 2.5}
| arm | 100 | 200 | 400 | adjoints | rounds | mean cos |
|---|---|---|---|---|---|---|
| C | 7.57 (1.25) | 6.84 (1.47) | 6.34 (1.08) | 0.0/0.0/0.0 | 5.0/10.0/20.0 | -/-/- |
| P | 6.37 (1.43) | 5.42 (1.16) | 5.05 (1.38) | 1.6/1.8/2.0 | 11.9/17.0/27.5 | 0.75/0.75/0.75 |
| G | 5.25 (0.78) | 3.25 (0.44) | 2.20 (0.41) | 3.0/7.0/15.0 | 8.0/16.0/32.0 | 0.75/0.75/0.75 |
| GP | 4.92 (0.98) | 3.47 (0.58) | 2.63 (0.41) | 3.7/6.6/11.6 | 13.7/21.0/34.7 | 0.75/0.75/0.75 |
| Glag | 6.79 (1.43) | 4.57 (0.68) | 3.38 (0.71) | 2.0/6.0/14.0 | 5.0/9.0/18.0 | 0.75/0.75/0.75 |
| G2 | 6.54 (0.94) | 4.87 (0.71) | 3.73 (0.74) | 2.0/4.0/8.0 | 7.0/14.0/27.0 | 0.75/0.75/0.75 |

| pair (X beats Y: X wins of n) | 100 | 200 | 400 |
|---|---|---|---|
| P vs C | 21/25 (-1.20 +- 0.26, p 0.0004) | 22/25 (-1.43 +- 0.31, p 0.0002) | 20/25 (-1.29 +- 0.29, p 0.0008) |
| G vs C | 25/25 (-2.32 +- 0.25, p 1e-05) | 25/25 (-3.59 +- 0.31, p 1e-05) | 25/25 (-4.14 +- 0.20, p 1e-05) |
| G vs P | 21/25 (-1.12 +- 0.26, p 0.001) | 25/25 (-2.17 +- 0.21, p 1e-05) | 25/25 (-2.85 +- 0.28, p 1e-05) |
| GP vs C | 25/25 (-2.66 +- 0.29, p 1e-05) | 25/25 (-3.37 +- 0.32, p 1e-05) | 25/25 (-3.71 +- 0.21, p 1e-05) |
| GP vs G | 14/25 (-0.33 +- 0.19, p 0.1) | 7/25 (+0.22 +- 0.12, p 0.1) | 3/25 (+0.43 +- 0.10, p 0.0005) |
| GP vs P | 19/25 (-1.45 +- 0.33, p 0.0005) | 25/25 (-1.95 +- 0.24, p 1e-05) | 25/25 (-2.42 +- 0.30, p 1e-05) |
| Glag vs G | 3/25 (+1.54 +- 0.24, p 3e-05) | 0/25 (+1.32 +- 0.13, p 1e-05) | 1/25 (+1.18 +- 0.14, p 2e-05) |
| Glag vs C | 20/25 (-0.78 +- 0.31, p 0.01) | 22/25 (-2.28 +- 0.31, p 3e-05) | 25/25 (-2.97 +- 0.25, p 1e-05) |
| G2 vs G | 1/25 (+1.29 +- 0.19, p 2e-05) | 1/25 (+1.62 +- 0.14, p 1e-05) | 0/25 (+1.53 +- 0.14, p 1e-05) |
| G2 vs C | 19/25 (-1.03 +- 0.25, p 0.001) | 23/25 (-1.97 +- 0.27, p 2e-05) | 25/25 (-2.62 +- 0.26, p 1e-05) |

**Registered decision (primary): G [G>P at [100, 200, 400], P>G at [], G>C at [100, 200, 400], P>C at [100, 200, 400], C>P at []]**

### setting cos0.4: {'cos': 0.4, 'A': 4, 'noise': 0.4, 'kappa': 2.5}
| arm | 100 | 200 | 400 | adjoints | rounds | mean cos |
|---|---|---|---|---|---|---|
| C | 7.57 (1.25) | 6.84 (1.47) | 6.34 (1.08) | 0.0/0.0/0.0 | 5.0/10.0/20.0 | -/-/- |
| P | 7.57 (1.08) | 6.46 (1.02) | 6.10 (1.36) | 1.4/1.8/1.6 | 11.4/17.1/27.1 | 0.40/0.40/0.40 |
| G | 7.21 (1.14) | 4.42 (0.58) | 2.92 (0.56) | 3.0/7.0/15.0 | 8.0/16.0/32.0 | 0.40/0.40/0.40 |
| GP | 7.18 (1.23) | 5.02 (0.75) | 3.46 (0.49) | 3.6/6.7/11.7 | 13.6/21.2/34.9 | 0.40/0.40/0.40 |
| Glag | 7.90 (1.55) | 5.74 (0.96) | 3.92 (0.73) | 2.0/6.0/14.0 | 5.0/9.0/18.0 | 0.40/0.40/0.40 |
| G2 | 7.32 (1.55) | 5.48 (0.73) | 4.43 (0.81) | 2.0/4.0/8.0 | 7.0/14.0/27.0 | 0.40/0.40/0.40 |

| pair (X beats Y: X wins of n) | 100 | 200 | 400 |
|---|---|---|---|
| P vs C | 12/25 (-0.01 +- 0.21, p 0.8) | 15/25 (-0.38 +- 0.24, p 0.2) | 14/25 (-0.24 +- 0.31, p 0.4) |
| G vs C | 13/25 (-0.36 +- 0.28, p 0.4) | 25/25 (-2.42 +- 0.30, p 1e-05) | 25/25 (-3.42 +- 0.21, p 1e-05) |
| G vs P | 17/25 (-0.35 +- 0.21, p 0.07) | 25/25 (-2.04 +- 0.17, p 1e-05) | 25/25 (-3.18 +- 0.30, p 1e-05) |
| GP vs C | 14/25 (-0.40 +- 0.29, p 0.3) | 23/25 (-1.83 +- 0.29, p 5e-05) | 25/25 (-2.88 +- 0.21, p 1e-05) |
| GP vs G | 9/25 (-0.04 +- 0.26, p 0.6) | 3/25 (+0.60 +- 0.15, p 0.0008) | 4/25 (+0.54 +- 0.11, p 0.0003) |
| GP vs P | 15/25 (-0.39 +- 0.26, p 0.2) | 24/25 (-1.44 +- 0.21, p 2e-05) | 25/25 (-2.64 +- 0.29, p 1e-05) |
| Glag vs G | 8/25 (+0.69 +- 0.29, p 0.03) | 2/25 (+1.31 +- 0.18, p 3e-05) | 3/25 (+1.00 +- 0.19, p 0.0001) |
| Glag vs C | 10/25 (+0.33 +- 0.34, p 0.2) | 19/25 (-1.11 +- 0.33, p 0.003) | 25/25 (-2.42 +- 0.28, p 1e-05) |
| G2 vs G | 9/25 (+0.11 +- 0.25, p 0.5) | 1/25 (+1.06 +- 0.16, p 3e-05) | 0/25 (+1.51 +- 0.19, p 1e-05) |
| G2 vs C | 11/25 (-0.25 +- 0.36, p 0.5) | 22/25 (-1.37 +- 0.28, p 0.0001) | 22/25 (-1.91 +- 0.28, p 3e-05) |

### setting cos0.9: {'cos': 0.9, 'A': 4, 'noise': 0.4, 'kappa': 2.5}
| arm | 100 | 200 | 400 | adjoints | rounds | mean cos |
|---|---|---|---|---|---|---|
| C | 7.57 (1.25) | 6.84 (1.47) | 6.34 (1.08) | 0.0/0.0/0.0 | 5.0/10.0/20.0 | -/-/- |
| P | 5.69 (1.61) | 5.05 (1.08) | 4.81 (1.38) | 1.8/1.9/2.0 | 12.0/17.1/27.4 | 0.90/0.90/0.90 |
| G | 4.67 (0.74) | 3.26 (0.61) | 2.16 (0.39) | 3.0/7.0/15.0 | 8.0/16.0/32.0 | 0.90/0.90/0.90 |
| GP | 4.46 (0.86) | 3.26 (0.44) | 2.55 (0.50) | 3.6/6.6/11.6 | 13.6/20.9/34.8 | 0.90/0.90/0.90 |
| Glag | 6.73 (1.14) | 4.49 (0.89) | 3.17 (0.68) | 2.0/6.0/14.0 | 5.0/9.0/18.0 | 0.90/0.90/0.90 |
| G2 | 5.57 (0.94) | 4.35 (0.79) | 3.59 (0.75) | 2.0/4.0/8.0 | 7.0/14.0/27.0 | 0.90/0.90/0.90 |

| pair (X beats Y: X wins of n) | 100 | 200 | 400 |
|---|---|---|---|
| P vs C | 23/25 (-1.88 +- 0.28, p 4e-05) | 23/25 (-1.79 +- 0.33, p 6e-05) | 20/25 (-1.53 +- 0.30, p 0.0003) |
| G vs C | 25/25 (-2.91 +- 0.21, p 1e-05) | 25/25 (-3.59 +- 0.32, p 1e-05) | 25/25 (-4.18 +- 0.22, p 1e-05) |
| G vs P | 19/25 (-1.02 +- 0.30, p 0.005) | 24/25 (-1.79 +- 0.22, p 2e-05) | 25/25 (-2.64 +- 0.28, p 1e-05) |
| GP vs C | 25/25 (-3.11 +- 0.25, p 1e-05) | 25/25 (-3.59 +- 0.31, p 1e-05) | 25/25 (-3.80 +- 0.22, p 1e-05) |
| GP vs G | 13/25 (-0.21 +- 0.16, p 0.4) | 7/25 (+0.00 +- 0.10, p 0.4) | 5/25 (+0.38 +- 0.12, p 0.004) |
| GP vs P | 19/25 (-1.23 +- 0.33, p 0.002) | 25/25 (-1.79 +- 0.21, p 1e-05) | 25/25 (-2.26 +- 0.29, p 1e-05) |
| Glag vs G | 0/25 (+2.06 +- 0.20, p 1e-05) | 1/25 (+1.23 +- 0.14, p 1e-05) | 1/25 (+1.00 +- 0.12, p 3e-05) |
| Glag vs C | 21/25 (-0.85 +- 0.26, p 0.001) | 25/25 (-2.35 +- 0.30, p 1e-05) | 25/25 (-3.18 +- 0.24, p 1e-05) |
| G2 vs G | 3/25 (+0.90 +- 0.15, p 6e-05) | 2/25 (+1.09 +- 0.16, p 3e-05) | 0/25 (+1.43 +- 0.15, p 1e-05) |
| G2 vs C | 25/25 (-2.00 +- 0.22, p 1e-05) | 24/25 (-2.50 +- 0.33, p 2e-05) | 25/25 (-2.75 +- 0.26, p 1e-05) |

### setting A2: {'cos': 0.75, 'A': 2, 'noise': 0.4, 'kappa': 2.5}
| arm | 100 | 200 | 400 | adjoints | rounds | mean cos |
|---|---|---|---|---|---|---|
| C | 7.57 (1.25) | 6.84 (1.47) | 6.34 (1.08) | 0.0/0.0/0.0 | 5.0/10.0/20.0 | -/-/- |
| P | 6.30 (1.39) | 5.45 (1.24) | 5.05 (1.38) | 1.6/1.8/2.0 | 12.0/17.1/27.7 | 0.75/0.75/0.75 |
| G | 4.58 (0.68) | 3.09 (0.56) | 2.16 (0.45) | 4.0/8.0/16.0 | 9.0/18.0/35.0 | 0.75/0.75/0.75 |
| GP | 4.94 (0.99) | 3.46 (0.55) | 2.37 (0.48) | 3.7/6.6/12.8 | 14.1/21.6/37.2 | 0.75/0.75/0.75 |
| Glag | 6.44 (1.05) | 4.39 (0.55) | 3.27 (0.56) | 3.0/7.0/15.0 | 5.0/10.0/19.0 | 0.75/0.75/0.75 |
| G2 | 6.51 (0.95) | 4.84 (0.73) | 3.48 (0.74) | 2.0/4.0/9.0 | 7.0/14.0/29.0 | 0.75/0.75/0.75 |

| pair (X beats Y: X wins of n) | 100 | 200 | 400 |
|---|---|---|---|
| P vs C | 22/25 (-1.27 +- 0.25, p 0.0001) | 22/25 (-1.39 +- 0.30, p 0.0002) | 20/25 (-1.29 +- 0.30, p 0.0007) |
| G vs C | 25/25 (-3.00 +- 0.29, p 1e-05) | 25/25 (-3.75 +- 0.27, p 1e-05) | 25/25 (-4.18 +- 0.23, p 1e-05) |
| G vs P | 22/25 (-1.73 +- 0.28, p 6e-05) | 25/25 (-2.36 +- 0.21, p 1e-05) | 25/25 (-2.89 +- 0.30, p 1e-05) |
| GP vs C | 25/25 (-2.63 +- 0.29, p 1e-05) | 25/25 (-3.39 +- 0.29, p 1e-05) | 25/25 (-3.98 +- 0.22, p 1e-05) |
| GP vs G | 9/25 (+0.37 +- 0.22, p 0.1) | 7/25 (+0.36 +- 0.12, p 0.005) | 6/25 (+0.21 +- 0.08, p 0.02) |
| GP vs P | 18/25 (-1.36 +- 0.32, p 0.0007) | 24/25 (-2.00 +- 0.27, p 2e-05) | 25/25 (-2.69 +- 0.29, p 1e-05) |
| Glag vs G | 1/25 (+1.87 +- 0.22, p 1e-05) | 0/25 (+1.30 +- 0.11, p 1e-05) | 2/25 (+1.11 +- 0.14, p 3e-05) |
| Glag vs C | 19/25 (-1.13 +- 0.26, p 0.0005) | 22/25 (-2.45 +- 0.31, p 3e-05) | 25/25 (-3.07 +- 0.24, p 1e-05) |
| G2 vs G | 1/25 (+1.93 +- 0.21, p 1e-05) | 1/25 (+1.75 +- 0.16, p 2e-05) | 1/25 (+1.32 +- 0.17, p 2e-05) |
| G2 vs C | 19/25 (-1.07 +- 0.25, p 0.001) | 23/25 (-2.00 +- 0.29, p 2e-05) | 25/25 (-2.86 +- 0.25, p 1e-05) |

### setting A8: {'cos': 0.75, 'A': 8, 'noise': 0.4, 'kappa': 2.5}
| arm | 100 | 200 | 400 | adjoints | rounds | mean cos |
|---|---|---|---|---|---|---|
| C | 7.57 (1.25) | 6.84 (1.47) | 6.34 (1.08) | 0.0/0.0/0.0 | 5.0/10.0/20.0 | -/-/- |
| P | 6.27 (1.33) | 5.39 (1.17) | 5.07 (1.41) | 1.5/1.8/2.0 | 11.0/16.5/27.1 | 0.75/0.75/0.75 |
| G | 5.29 (0.90) | 3.50 (0.38) | 2.37 (0.34) | 3.0/6.0/12.0 | 7.0/14.0/28.0 | 0.75/0.75/0.75 |
| GP | 4.99 (1.00) | 3.69 (0.74) | 2.73 (0.49) | 3.7/5.8/10.6 | 12.3/19.1/31.8 | 0.75/0.75/0.75 |
| Glag | 7.20 (1.40) | 4.94 (0.90) | 3.68 (0.61) | 2.0/5.0/11.0 | 5.0/8.0/16.0 | 0.75/0.75/0.75 |
| G2 | 6.46 (0.98) | 4.77 (0.64) | 3.89 (0.74) | 2.0/4.0/7.0 | 7.0/13.0/25.0 | 0.75/0.75/0.75 |

| pair (X beats Y: X wins of n) | 100 | 200 | 400 |
|---|---|---|---|
| P vs C | 21/25 (-1.30 +- 0.26, p 0.0002) | 22/25 (-1.45 +- 0.31, p 0.0002) | 19/25 (-1.27 +- 0.30, p 0.001) |
| G vs C | 25/25 (-2.28 +- 0.24, p 1e-05) | 25/25 (-3.35 +- 0.31, p 1e-05) | 25/25 (-3.97 +- 0.21, p 1e-05) |
| G vs P | 21/25 (-0.98 +- 0.26, p 0.001) | 24/25 (-1.90 +- 0.24, p 2e-05) | 25/25 (-2.70 +- 0.28, p 1e-05) |
| GP vs C | 25/25 (-2.59 +- 0.29, p 1e-05) | 25/25 (-3.16 +- 0.26, p 1e-05) | 25/25 (-3.62 +- 0.21, p 1e-05) |
| GP vs G | 15/25 (-0.30 +- 0.19, p 0.1) | 9/25 (+0.19 +- 0.13, p 0.2) | 3/25 (+0.35 +- 0.08, p 0.0002) |
| GP vs P | 18/25 (-1.29 +- 0.32, p 0.002) | 24/25 (-1.71 +- 0.26, p 5e-05) | 24/25 (-2.35 +- 0.30, p 2e-05) |
| Glag vs G | 0/25 (+1.91 +- 0.21, p 1e-05) | 0/25 (+1.44 +- 0.17, p 1e-05) | 0/25 (+1.31 +- 0.16, p 1e-05) |
| Glag vs C | 15/25 (-0.38 +- 0.24, p 0.1) | 23/25 (-1.91 +- 0.30, p 3e-05) | 25/25 (-2.67 +- 0.26, p 1e-05) |
| G2 vs G | 2/25 (+1.17 +- 0.19, p 4e-05) | 0/25 (+1.28 +- 0.13, p 1e-05) | 1/25 (+1.52 +- 0.15, p 1e-05) |
| G2 vs C | 20/25 (-1.12 +- 0.24, p 0.0005) | 23/25 (-2.07 +- 0.29, p 2e-05) | 25/25 (-2.46 +- 0.25, p 1e-05) |

### setting noise0.2: {'cos': 0.75, 'A': 4, 'noise': 0.2, 'kappa': 2.5}
| arm | 100 | 200 | 400 | adjoints | rounds | mean cos |
|---|---|---|---|---|---|---|
| C | 6.66 (0.92) | 5.98 (0.95) | 5.25 (0.95) | 0.0/0.0/0.0 | 5.0/10.0/20.0 | -/-/- |
| P | 4.85 (0.97) | 4.77 (1.02) | 3.95 (0.85) | 2.6/2.2/2.3 | 13.6/18.2/28.2 | 0.75/0.75/0.75 |
| G | 4.82 (0.62) | 3.02 (0.43) | 1.82 (0.34) | 3.0/7.0/15.0 | 8.0/16.0/32.0 | 0.75/0.75/0.75 |
| GP | 4.28 (0.75) | 3.21 (0.43) | 2.20 (0.48) | 4.3/6.8/11.6 | 15.0/21.5/34.9 | 0.75/0.75/0.75 |
| Glag | 6.46 (0.94) | 4.11 (0.52) | 2.65 (0.35) | 2.0/6.0/14.0 | 5.0/9.0/18.0 | 0.75/0.75/0.75 |
| G2 | 5.93 (0.84) | 4.33 (0.60) | 3.16 (0.60) | 2.0/4.0/8.0 | 7.0/14.0/27.0 | 0.75/0.75/0.75 |

| pair (X beats Y: X wins of n) | 100 | 200 | 400 |
|---|---|---|---|
| P vs C | 25/25 (-1.81 +- 0.26, p 1e-05) | 20/25 (-1.22 +- 0.22, p 0.0001) | 24/25 (-1.30 +- 0.19, p 2e-05) |
| G vs C | 25/25 (-1.84 +- 0.21, p 1e-05) | 25/25 (-2.97 +- 0.22, p 1e-05) | 25/25 (-3.43 +- 0.23, p 1e-05) |
| G vs P | 14/25 (-0.03 +- 0.25, p 1) | 24/25 (-1.75 +- 0.21, p 1e-05) | 25/25 (-2.13 +- 0.21, p 1e-05) |
| GP vs C | 25/25 (-2.38 +- 0.25, p 1e-05) | 25/25 (-2.78 +- 0.21, p 1e-05) | 25/25 (-3.05 +- 0.24, p 1e-05) |
| GP vs G | 19/25 (-0.54 +- 0.17, p 0.009) | 7/25 (+0.19 +- 0.07, p 0.02) | 2/25 (+0.38 +- 0.07, p 8e-05) |
| GP vs P | 20/25 (-0.57 +- 0.20, p 0.006) | 24/25 (-1.56 +- 0.20, p 2e-05) | 24/25 (-1.75 +- 0.20, p 2e-05) |
| Glag vs G | 2/25 (+1.64 +- 0.22, p 2e-05) | 1/25 (+1.09 +- 0.12, p 1e-05) | 1/25 (+0.83 +- 0.07, p 1e-05) |
| Glag vs C | 14/25 (-0.20 +- 0.19, p 0.3) | 24/25 (-1.87 +- 0.20, p 2e-05) | 25/25 (-2.60 +- 0.22, p 1e-05) |
| G2 vs G | 1/25 (+1.11 +- 0.12, p 1e-05) | 0/25 (+1.31 +- 0.13, p 1e-05) | 1/25 (+1.34 +- 0.15, p 1e-05) |
| G2 vs C | 18/25 (-0.73 +- 0.26, p 0.01) | 23/25 (-1.66 +- 0.22, p 2e-05) | 25/25 (-2.09 +- 0.20, p 1e-05) |

### setting kappa1: {'cos': 0.75, 'A': 4, 'noise': 0.4, 'kappa': 1.0}
| arm | 100 | 200 | 400 | adjoints | rounds | mean cos |
|---|---|---|---|---|---|---|
| C | 7.57 (1.25) | 6.84 (1.47) | 6.34 (1.08) | 0.0/0.0/0.0 | 5.0/10.0/20.0 | -/-/- |
| P | 6.23 (1.54) | 5.42 (1.12) | 5.21 (1.50) | 1.7/1.8/2.0 | 12.0/17.0/27.5 | 0.75/0.75/0.75 |
| G | 5.09 (0.66) | 3.28 (0.43) | 2.34 (0.44) | 3.0/7.0/15.0 | 8.0/16.0/32.0 | 0.75/0.75/0.75 |
| GP | 4.93 (0.93) | 3.55 (0.58) | 2.65 (0.30) | 3.6/6.7/11.6 | 13.6/21.1/34.8 | 0.75/0.75/0.75 |
| Glag | 7.00 (1.30) | 4.67 (0.66) | 3.13 (0.58) | 2.0/6.0/14.0 | 5.0/9.0/18.0 | 0.75/0.75/0.75 |
| G2 | 6.21 (1.04) | 4.82 (0.84) | 3.74 (0.66) | 2.0/4.0/8.0 | 7.0/14.0/27.0 | 0.75/0.75/0.75 |

| pair (X beats Y: X wins of n) | 100 | 200 | 400 |
|---|---|---|---|
| P vs C | 21/25 (-1.34 +- 0.30, p 0.0003) | 22/25 (-1.42 +- 0.30, p 0.0002) | 18/25 (-1.13 +- 0.35, p 0.006) |
| G vs C | 25/25 (-2.48 +- 0.23, p 1e-05) | 25/25 (-3.56 +- 0.30, p 1e-05) | 25/25 (-4.00 +- 0.20, p 1e-05) |
| G vs P | 20/25 (-1.15 +- 0.30, p 0.004) | 25/25 (-2.14 +- 0.20, p 1e-05) | 25/25 (-2.87 +- 0.30, p 1e-05) |
| GP vs C | 25/25 (-2.64 +- 0.30, p 1e-05) | 25/25 (-3.30 +- 0.31, p 1e-05) | 25/25 (-3.70 +- 0.19, p 1e-05) |
| GP vs G | 13/25 (-0.16 +- 0.22, p 0.7) | 9/25 (+0.26 +- 0.14, p 0.08) | 5/25 (+0.30 +- 0.10, p 0.006) |
| GP vs P | 19/25 (-1.30 +- 0.34, p 0.001) | 24/25 (-1.88 +- 0.24, p 1e-05) | 25/25 (-2.57 +- 0.30, p 1e-05) |
| Glag vs G | 0/25 (+1.91 +- 0.24, p 1e-05) | 0/25 (+1.39 +- 0.14, p 1e-05) | 2/25 (+0.79 +- 0.12, p 3e-05) |
| Glag vs C | 19/25 (-0.57 +- 0.23, p 0.03) | 23/25 (-2.17 +- 0.29, p 3e-05) | 25/25 (-3.21 +- 0.23, p 1e-05) |
| G2 vs G | 1/25 (+1.13 +- 0.18, p 6e-05) | 1/25 (+1.54 +- 0.19, p 1e-05) | 1/25 (+1.40 +- 0.15, p 2e-05) |
| G2 vs C | 21/25 (-1.36 +- 0.32, p 0.0006) | 22/25 (-2.02 +- 0.36, p 5e-05) | 25/25 (-2.61 +- 0.25, p 1e-05) |

### setting cosmodel: {'cos': 0.75, 'A': 4, 'noise': 0.4, 'kappa': 2.5, 'mode': 'model'}
| arm | 100 | 200 | 400 | adjoints | rounds | mean cos |
|---|---|---|---|---|---|---|
| C | 7.57 (1.25) | 6.84 (1.47) | 6.34 (1.08) | 0.0/0.0/0.0 | 5.0/10.0/20.0 | -/-/- |
| P | 6.31 (1.57) | 5.44 (1.22) | 5.30 (1.48) | 1.6/1.9/2.0 | 11.6/17.1/27.7 | 0.74/0.71/0.69 |
| G | 5.04 (0.71) | 3.59 (0.54) | 2.50 (0.50) | 3.0/7.0/15.0 | 8.0/16.0/32.0 | 0.73/0.67/0.62 |
| GP | 4.93 (0.88) | 3.70 (0.52) | 2.92 (0.43) | 3.8/6.6/11.6 | 14.0/21.1/34.6 | 0.72/0.68/0.64 |
| Glag | 6.89 (1.34) | 4.69 (0.88) | 3.36 (0.68) | 2.0/6.0/14.0 | 5.0/9.0/18.0 | 0.75/0.73/0.70 |
| G2 | 6.26 (1.06) | 4.86 (0.87) | 3.86 (0.92) | 2.0/4.0/8.0 | 7.0/14.0/27.0 | 0.74/0.73/0.70 |

| pair (X beats Y: X wins of n) | 100 | 200 | 400 |
|---|---|---|---|
| P vs C | 19/25 (-1.26 +- 0.30, p 0.001) | 19/25 (-1.41 +- 0.34, p 0.0006) | 20/25 (-1.04 +- 0.32, p 0.005) |
| G vs C | 25/25 (-2.53 +- 0.28, p 1e-05) | 25/25 (-3.25 +- 0.30, p 1e-05) | 25/25 (-3.85 +- 0.22, p 1e-05) |
| G vs P | 16/25 (-1.27 +- 0.34, p 0.003) | 25/25 (-1.85 +- 0.23, p 1e-05) | 25/25 (-2.80 +- 0.29, p 1e-05) |
| GP vs C | 25/25 (-2.65 +- 0.29, p 1e-05) | 25/25 (-3.15 +- 0.28, p 1e-05) | 25/25 (-3.42 +- 0.22, p 1e-05) |
| GP vs G | 16/25 (-0.11 +- 0.15, p 0.5) | 8/25 (+0.10 +- 0.08, p 0.2) | 5/25 (+0.42 +- 0.10, p 0.0006) |
| GP vs P | 21/25 (-1.38 +- 0.30, p 0.0001) | 24/25 (-1.74 +- 0.23, p 1e-05) | 25/25 (-2.38 +- 0.30, p 1e-05) |
| Glag vs G | 1/25 (+1.85 +- 0.25, p 2e-05) | 1/25 (+1.10 +- 0.18, p 1e-05) | 3/25 (+0.87 +- 0.14, p 0.0001) |
| Glag vs C | 19/25 (-0.69 +- 0.28, p 0.01) | 25/25 (-2.15 +- 0.28, p 1e-05) | 25/25 (-2.98 +- 0.25, p 1e-05) |
| G2 vs G | 1/25 (+1.22 +- 0.21, p 2e-05) | 2/25 (+1.27 +- 0.16, p 2e-05) | 0/25 (+1.36 +- 0.17, p 1e-05) |
| G2 vs C | 21/25 (-1.32 +- 0.25, p 9e-05) | 23/25 (-1.98 +- 0.24, p 2e-05) | 24/25 (-2.48 +- 0.30, p 2e-05) |

### setting fail0.5: {'cos': 0.75, 'A': 4, 'noise': 0.4, 'kappa': 2.5, 'fail': 0.5}
| arm | 100 | 200 | 400 | adjoints | rounds | mean cos |
|---|---|---|---|---|---|---|
| C | 7.57 (1.25) | 6.84 (1.47) | 6.34 (1.08) | 0.0/0.0/0.0 | 5.0/10.0/20.0 | -/-/- |
| P | 7.09 (1.64) | 6.37 (1.70) | 5.68 (1.30) | 1.4/1.4/1.4 | 8.5/13.4/23.5 | 0.75/0.75/0.75 |
| G | 6.24 (1.35) | 4.44 (0.84) | 3.64 (0.78) | 3.0/7.0/15.0 | 8.0/16.0/32.0 | 0.75/0.75/0.75 |
| GP | 6.46 (1.52) | 4.83 (1.05) | 3.75 (0.78) | 3.4/6.6/11.4 | 10.7/18.0/32.2 | 0.75/0.75/0.75 |
| Glag | 7.30 (1.30) | 5.51 (0.80) | 4.45 (0.79) | 2.0/6.0/14.0 | 5.0/9.0/18.0 | 0.75/0.75/0.75 |
| G2 | 6.92 (0.93) | 5.55 (0.93) | 4.56 (0.69) | 2.0/4.0/8.0 | 7.0/14.0/27.0 | 0.75/0.75/0.75 |

| pair (X beats Y: X wins of n) | 100 | 200 | 400 |
|---|---|---|---|
| P vs C | 16/25 (-0.48 +- 0.28, p 0.09) | 15/25 (-0.47 +- 0.42, p 0.2) | 17/25 (-0.67 +- 0.29, p 0.03) |
| G vs C | 19/25 (-1.34 +- 0.31, p 0.0004) | 22/25 (-2.41 +- 0.39, p 4e-05) | 25/25 (-2.70 +- 0.25, p 1e-05) |
| G vs P | 17/25 (-0.85 +- 0.24, p 0.003) | 22/25 (-1.94 +- 0.39, p 9e-05) | 24/25 (-2.04 +- 0.26, p 2e-05) |
| GP vs C | 19/25 (-1.11 +- 0.31, p 0.002) | 21/25 (-2.01 +- 0.38, p 0.0001) | 25/25 (-2.59 +- 0.27, p 1e-05) |
| GP vs G | 10/25 (+0.22 +- 0.17, p 0.4) | 8/25 (+0.40 +- 0.17, p 0.04) | 9/25 (+0.11 +- 0.11, p 0.3) |
| GP vs P | 16/25 (-0.63 +- 0.28, p 0.04) | 20/25 (-1.54 +- 0.40, p 0.001) | 24/25 (-1.93 +- 0.27, p 2e-05) |
| Glag vs G | 6/25 (+1.06 +- 0.27, p 0.002) | 3/25 (+1.07 +- 0.22, p 0.0002) | 4/25 (+0.81 +- 0.19, p 0.0005) |
| Glag vs C | 16/25 (-0.27 +- 0.20, p 0.04) | 21/25 (-1.33 +- 0.34, p 0.0005) | 24/25 (-1.89 +- 0.23, p 1e-05) |
| G2 vs G | 7/25 (+0.68 +- 0.27, p 0.02) | 5/25 (+1.12 +- 0.26, p 0.0009) | 4/25 (+0.92 +- 0.19, p 0.0002) |
| G2 vs C | 16/25 (-0.66 +- 0.24, p 0.008) | 22/25 (-1.29 +- 0.29, p 0.0002) | 25/25 (-1.78 +- 0.22, p 1e-05) |
