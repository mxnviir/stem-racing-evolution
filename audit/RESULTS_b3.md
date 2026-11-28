# Drag audit b3: where the original hybrid's gain went (2026-11-28)

Readings R1-R3 were fixed in `RESULTS_b2.md` before these results. Audit protocol (fresh wall-resolved mesh, frames 0-3,
linearUpwind, 1000 iterations, score = mean of iterations 200-1000), half-car mN, paired by frame, SE = sd / 2.
Raw: `audit/results/b3/`.

| frame | hybrid_region | legal:hybrid_region | legal:team_rm | CAD team pod (b1) | seed evolution's legal team pod (b2) |
|---|---|---|---|---|---|
| k0 | 172.64 | 175.69 | 176.58 | 180.97 | 174.18 |
| k1 | 172.82 | 173.72 | 178.27 | 182.86 | 175.26 |
| k2 | 173.14 | 175.35 | 176.42 | 181.03 | 174.63 |
| k3 | 172.70 | 173.83 | 178.21 | 182.92 | 174.94 |
| mean | 172.83 | 174.65 | 177.37 | 181.95 | 174.75 |

## The three readings

| Reading | Comparison | Mean | SE | Verdict |
|---|---|---|---|---|
| R1 | hybrid_region - CAD team pod | **-9.12** | 0.59 | 88 % of the full hybrid's -10.39 (mark: >= 70 %). **The out-of-region moves were not the main source of the gain.** |
| R2 | legal:hybrid_region - legal:team_rm | **-2.72** | 1.01 | 2.7 SE; mark (<= -1.4 mN and >= 3 SE) **not met**, narrowly |
| R3 | legal:team_rm - seed evolution's legal team pod | **+2.62** | 0.33 | 8 SE: **the two legal builds of the same team pod are not the same car** |

## R3 corrects b2

b2 judged the legalised hybrid against seed evolution's own rebuild of the team pod (174.75) and found -0.43 +- 0.53 mN,
"no measurable gain". R3 shows that reference is 2.6 mN away from the team pod legalised by the path the hybrid went
through (`seedhybrid/legalise.py`, 177.37). Against the same-path reference:

| Comparison (post hoc: b2's design against b3's reference) | Mean | SE |
|---|---|---|
| legal:hybrid_final - legal:team_rm | **-3.05** | 0.80 (3.8 SE) |
|   body | -6.74 | 0.33 |
|   wheels | +3.70 | 0.51 |

So the legalised original hybrid **keeps about 30 % of its gain** (-3.05 of -10.39 mN), significant at 3.8 SE, and b2's
"does not survive legalisation" is withdrawn. This comparison was not pre-registered as a mark (R3 was there to check
the reference, and it found the reference wrong); it is reported as post hoc. Which reference is fair: the same-path one,
because both designs then differ only in shape; the cross-path difference is a property of the two builders.

## Two caveats the numbers carry

- **Frame pattern.** Every comparison against the team pod splits by frame: about -1.7 mN on frames 0 and 2, -4.4 mN on 1
  and 3. The team pod (CAD and legalised alike) reads ~1.9 mN higher on frames 1 and 3, all of it in wheel drag
  (151.3 vs 149.4-150.0 mN); the hybrid designs do not. The mean gain therefore leans on a frame-dependent wheel-wake state
  of the team pod. Whether avoiding that state is a real gain (robustness) or mesh luck needs more frames (k4-k7).
- **Body down, wheels up.** In-region shape lowers body drag 6.7 mN and raises wheel drag 3.7-4.0 mN once legal
  (unlegalised, the wheels went down 6.9 mN, b1). The legal shape moves flow onto the wheels; the whole-car pipeline's
  objective already counts both.

## Legalisation itself

- The legaliser moves the team pod by -4.57 +- 0.07 mN (CAD 181.95 -> 177.37), consistently on every frame.
- It costs the region-only design +1.82 +- 0.50 mN (172.83 -> 174.65).

## Consequences

- METHODS_TRIED and RESULTS_b2 are corrected (same-path reference; ~30 % of the gain kept).
- Every future comparison of legalised designs uses a reference built by the **same** path (the whole-car pipeline does
  this by construction: every car, reference included, goes through the same builder).
- If the claim matters for the paper: frames k4-k7 for legal:hybrid_final and legal:team_rm (8 runs) settle the frame
  pattern.
