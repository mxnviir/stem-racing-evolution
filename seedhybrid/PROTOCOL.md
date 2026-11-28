# Seeded hybrid: experimental protocol

Written 2026-11-28 06:30 UTC, while the setup stage (run 46113354235) was still solving and **before any generation
result existed**. Every mark below is fixed here; any later change goes into the change log at the end with its date
and reason, and the original mark is still reported.

Companion pages: [GUIDE.md](GUIDE.md) (how the workflow works), [README.md](README.md) (settings record),
`RESULTS.md` (written by `analyze.py` after each stage).

## 0. What "the best workflow, provable" means here

A workflow is a set of choices: where to start, which directions to try, how to judge a step, when to stop, how to
turn the result into a legal car. "Best" is only defensible as: **every choice was either compared against its
alternative in an experiment with a pass mark fixed beforehand, or is listed below as an untested assumption with the
experiment that would test it.** Section 1 lists every choice. Section 2 is the experiments. Section 3 is what is
still assumed, after all of them, and why.

The objective throughout is **drag** (half car, k-omega SST, the audit's protocol for claims). Race time is a later
stage by the owner's decision.

## 1. Every choice in the workflow, and how it is decided

| # | Choice | Current setting | Alternative(s) | Decided by |
|---|---|---|---|---|
| C1 | Use an adjoint (DAFoam) to propose directions | yes | no adjoint (random + recombination only) | **E2** ablation; **E1** per-generation credit |
| C2 | Which adjoint | DAFoam, coarse 105 k mesh, steady model | ESI continuous; DAFoam on a finer mesh | ESI: dropped on evidence (failed validation, change log in README). Coarse vs finer DAFoam: **E1** fidelity per lineage; finer DAFoam = gap G4 |
| C3 | Recombine measured directions (`comb`, `es`) | yes | no | **E1** per-generation credit |
| C4 | Random probes per generation | 3 (fill to 6 children in gen 1, 7 after) | more / fewer | **E1** (how often a random probe is the winner); sensitivity = gap G5 |
| C5 | Inside-loop yardstick | steady model on the lineage's morphed base mesh | production scheme on fresh meshes (the claim's yardstick) | **E3** yardstick fidelity |
| C6 | Acceptance rule | best child beats the parent's re-run by >= 3 SE and >= 0.15 mN | other thresholds | **E3** (false accepts are visible as steps that do not hold on fresh meshes) |
| C7 | Step size schedule | 1.0 mm, halve when nothing is accepted, stop below 0.125 mm | fixed size; line search | **E1** (gen 1 measures the DAFoam direction at 1.0 and 0.5 mm: linearity); otherwise gap G5 |
| C8 | Starting shape | three legal seeds (team rebuild, d7, window) | the team CAD surface (original hybrid); one seed | **E6** seeds compared; original hybrid = reference (audit T1) |
| C9 | Moves confined to the legal region, legalised at the end | yes (region taper), legalisation at confirm | legalise every accepted step | **E4** legalisation loss |
| C10 | Stop when the morph range is spent | yes (phase B re-meshing not built) | re-mesh and continue | **E8** (after phase B is built) |
| C11 | Final claim protocol | fresh meshes, production scheme, frames 0-3, against the lineage's own legal seed; finer mesh and transition model for the winner | - | **E5** (the audit's own marks, applied to the winner) |

## 2. Experiments

All drag numbers are half-car mN. "Paired" = same meshing frame (claims) or same parent and mesh (inside the loop).

### E0. Setup checks (per seed; recorded, not a gate except where stated)
- Base state exists: **gate** (a lineage without one does not run).
- `steady_pass` (the steady model settles: drag wanders <= 0.10 mN over the last 1000 iterations): recorded. If it
  fails, that lineage's SEs are reported alongside every decision, since a wandering base inflates them.
- DAFoam converged at the seed: recorded per seed (before this campaign: 1 of 6 on shapes unlike the team pod).

### E1. Direction credit (no extra runs: the loop's own children)
Every generation, all children start from the same parent on the same mesh at the same step size, so their measured
drag changes are paired comparisons of direction types.
- Per direction type (`daf`, `comb`, `es`, `rnd`), pooled over every generation of every guided lineage:
  - win rate (how often it is the generation's best child),
  - mean drag change per mm RMS of movement, with a bootstrap 95 % interval over generations,
  - share of accepted steps it produced, and of the total accepted drag change.
- DAFoam fidelity per gradient: its first-order prediction of every child's drag change (rebuilt afterwards from
  the stored gradient and move; deterministic, so not open to after-the-fact choice) against the measured change:
  sign agreement and Pearson r, per generation and pooled.
- **Mark E1a (the adjoint is a better direction than a random one):** across all guided generations with a `daf`
  child, the `daf` child's drag change is below the median random child's in >= 70 % of generations, and the pooled
  mean `daf` change per mm is below the pooled mean random change per mm with a bootstrap interval that excludes 0.
- **Mark E1b (the gradient is trustworthy on this lineage):** sign agreement >= 80 % on children whose |change| > 3 SE.
  A lineage that fails E1b is reported as "adjoint not trusted here".

### E2. Adjoint ablation (the causal test of C1)
- Arm A (guided): lineage `team` as configured.
- Arm B (no adjoint): lineage `team_noadj`, the **same seed, base mesh and base state**, the same number of children
  every generation (gen 1: 6; after: 7), the same acceptance rule, step schedule, stopping rule and generation limit.
  Its children: the parent's re-run, `comb` and `es` (from its own measured moves), and random moves filling the rest.
  Random fields use the same seeds in both arms where both draw them (common random numbers).
- Second pair, chosen by a rule fixed now, **before** any DAFoam result is known: `d7_noadj` against `d7` if DAFoam
  converged at the d7 seed; otherwise `window_noadj` against `window` if DAFoam converged at the window seed; otherwise
  no second pair (a no-adjoint copy of a lineage that has no gradient anyway would test nothing).
- Outcome per arm: the final design, legalised, scored on frames 0-3 (E5 protocol), paired by frame.
- **Mark E2 (the adjoint causes the gain):** for each pair, guided final minus no-adjoint final <= 0 with |mean| >= 3 SE.
  Both pairs pass -> "the adjoint causes the gain" is supported; one passes -> "supported on that seed"; none -> not
  supported, and the workflow's adjoint stage must be justified otherwise or dropped.
- Also reported: drag after each generation in both arms (inside-loop yardstick), and CFD runs used by each arm.

### E3. Yardstick fidelity (does the inside-loop measurement predict the claim's measurement?)
- Every accepted parent of every lineage (and each seed) is scored on fresh meshes with the production scheme,
  frames 0 and 1, **without** legalisation (it measures the yardstick, not legality).
- **Mark E3a:** Spearman rank correlation >= 0.8, pooled over all scored parents, between the cumulative inside-loop
  change (vs the seed) and the fresh-mesh change (vs the seed's own fresh score).
- **Mark E3b:** every accepted step whose inside-loop change is <= -1.0 mN has the same sign on fresh meshes in
  >= 80 % of cases.
- Fail of E3a or E3b -> the acceptance rule (C6) needs a fresh-mesh check on every accepted step; the campaign's
  inside-loop numbers may then not be cited even as a trend.

### E4. Legalisation loss
- Each lineage's final design is scored on frames 0-3 both **as the loop left it** and **after legalisation**
  (rebuilt through the legal builder: 3-axis ball-end milling, 3 mm walls, one piece, cargo, clearance).
- Reported: the drag change legalisation causes, and the volume it changes.
- **Mark E4:** legalisation keeps >= 70 % of the lineage's fresh-mesh gain. Fail -> legalising inside the loop
  (at every re-mesh, phase B) becomes part of the workflow, and E8 must include it.

### E5. Final claim (the audit's protocol, applied to each lineage's legal final)
- **F1:** a lineage's legal final beats its own legal seed: paired mean over frames 0-3 <= -1.4 mN and >= 3 SE.
- **F2:** the campaign's best legal design beats d7 (171.44 mN, the best legal design before this campaign) by >= 3 SE.
- **Robustness of the winner** (the best legal final), each against its seed on the same settings:
  - finer mesh (background divisor 15.12, ~3.6 M cells), frames 0 and 1: same sign and within 25 % of the scoring
    mesh's difference, or within 2.65 mN if the difference is under 10 mN (the audit's T2 rule);
  - transition model (kOmegaSSTLM), frame 0: same sign, within a factor 2 (the audit's T4 rule);
  - flow state: body and wheel drag split per frame; a frame whose split sits apart from the design's other frames is
    flagged (the audit's T3 rule).

### E6. Does the starting seed matter?
- The three guided lineages' legal finals, frames 0-3 each: reported with SEs; pairwise differences >= 3 SE called.
- Reported: gain from seed vs absolute final drag (a worse seed may gain more and still end worse).

### E7. Cost accounting
- CFD runs and runner-hours per lineage and per arm, split into setup, children, gradients and confirmation.
- Gain per 100 runner-hours for each arm, against seed evolution's POC (56 runs for -3.3 mN, audit T1) and the
  original hybrid (audit T1: -10.4 mN).

### E8. Continuing past the morph limit (after phase B is built)
- When a lineage's morph range is spent, its parent is legalised, re-meshed (a fresh base), given a fresh DAFoam
  gradient, and continued as a new segment.
- **Mark E8:** the continued lineage's legal final beats the legal design at the hand-over by >= 3 SE (frames 0-3).
  Fail -> stopping at the morph limit is the workflow; pass -> re-meshing is.

## 3. What remains assumed after all experiments (residual gaps)

These cannot be closed inside this campaign. Each says what would close it.

| # | Gap | Why it stays | What closes it |
|---|---|---|---|
| G1 | RANS (k-omega SST) model error | Every number is a RANS number; the Ahmed check showed +3 to +13 % drag error at a different Reynolds number | Wind-tunnel or track measurement of two designs; LES of the winner |
| G2 | One car, one region | All lineages share the trial car and the sidepod region | Repeat on the whole-car optimiser's region set |
| G3 | Race time | The owner's staging: drag first | The calibrated race model, once measured |
| G4 | DAFoam only on the coarse mesh | Memory on 16 GB runners | DAFoam on the 2.1 M mesh on a rented high-memory machine; compare E1 fidelity |
| G5 | Settings not tuned (random probes per generation, step schedule, acceptance thresholds, generations) | Taken from the original loop and fixed beforehand; tuning on results would bias them | A sensitivity campaign on a second car, reported as such |
| G6 | Few lineages | Two ablation pairs at most: the lineage-level comparison has n <= 2 | More seeds per arm; E1's within-generation pairs carry most of the statistical weight meanwhile |
| G7 | Frame-to-frame flow-state flips | Fresh-mesh scores switch between states ~1.9 mN apart (audit T3) | Label states per run (E5); URANS for the winner |

## How the stages map to the experiments

| Stage | Runs | Feeds |
|---|---|---|
| setup | per seed: geometry, base state, DAFoam at the seed | E0 |
| generations 1..6 | per lineage: 6-7 children (+ a DAFoam gradient where the parent moved) | E1, E2 trajectory, E7 |
| confirm (when every lineage has stopped) | every accepted parent unlegalised, frames 0-1; the unlegalised seed frames 0-3; each final unlegalised and legalised frames 0-3; each legal seed frames 0-3; identical designs scored once | E2, E3, E4, E5 (F1, F2), E6 |
| robust | the best legal final and its legal seed: finer mesh frames 0-1, transition model frame 0 | E5 robustness |

`analyze.py` applies every mark mechanically after each stage and writes `RESULTS.md`.

## Change log
- 2026-11-28 06:30 UTC: protocol written (before any generation result).
- 2026-11-28 07:40 UTC (before any generation result): the unlegalised seed is scored on frames 0-3 (not only 0-1), so
  E4's "gain kept" uses the same four frames on both sides. Flow-state flag made concrete: a frame whose body drag is
  more than 1.0 mN from the median body drag of that design's frames.
