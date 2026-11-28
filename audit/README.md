# Drag audit (batch b1)

Independent checks of the drag model behind the sidepod searches in `mxnviir/stem-racing-adjoint`
(pinned at `9ab9092`). Nothing in that repository is changed. The CFD code runs unchanged except for two
environment hooks added on the runner (`.github/workflows/audit_drag.yml`); with their defaults the code is
byte-for-byte the same behaviour:

- `P2_BG_DIV`: background cell = car length / `P2_BG_DIV` (default 12). Every snappy level, wake box and layer stack
  is relative to it, so changing it refines the whole mesh systematically (first-layer height fixed, as it should be).
- `PRODCHECK_TURB`: turbulence model (default `kOmegaSST`).

Every case is scored exactly as seed evolution scores (`seedevo/campaign_case.py`): fresh wall-resolved mesh,
meshing frame k (domain shifted k x 0.37 mm), production scheme, cold start, 1000 iterations, drag = mean of
iterations 200-1000 (half car, mN), SE from four 200-iteration blocks.

Designs:
- `team_rm`: the team pod surface the hybrid loop started from (`gate/original_rm`).
- `hybrid_final`: `hybrid/final_design/car_half.stl`, the loop's generation-6 design (steady model -9.94 mN on the
  morphed mesh; production scheme -6.32 mN for generation 3 on the morphed mesh; never scored on a fresh mesh).
  Not legalised (milling, mass rules off), so this measures drag only.
- `seedevo:7` / `seedevo:10`: d7 and d10 of `seedevo/campaign_poc` (frames 0-1 already run; this adds frames 2-3).

## Pass marks (written before any result of this batch exists)

Paired difference per frame: delta_k = drag(design, k) - drag(reference, k); mean over frames, SE = sd/sqrt(n).

**T1. Does the hybrid's gain hold on fresh meshes?** (`t1_*`, frames 0-3, scoring mesh)
- Reference `team_rm`. **Holds** if mean delta <= -4.0 mN and |mean| >= 3 SE.
  **Partly** if -4.0 < mean delta and |mean| >= 3 SE. **Gone** if |mean| < 3 SE.
- Seed evolution's best on the same yardstick: d7 and d10 against the seedevo team pod (`original`, frames 0-3
  already in `campaign_poc/state.json`), frames 0-3. Reported side by side; the approach with the larger confirmed
  gain against its own team-pod reference "leads". Different team-pod surfaces (CAD remesh vs rebuilt field) are why
  each design is compared with its own reference.

**T2. Is a few-mN difference resolved by the scoring mesh?** (`t2L0_*`, `t2L2_*`; L1 = T1's frames 0-1)
- Levels by background divisor: L0 9.524, L1 12, L2 15.12 (refinement ratio 1.26 per step, about 2x cells).
- **Resolved at L1** if delta(L2) has the same sign as delta(L1) and |delta(L2) - delta(L1)| <= max(1.0 mN, 0.25 |delta(L1)|).
- Reported: absolute drag per level per design, the observed order of convergence and a GCI on the difference where
  the three levels are monotone, and cells / wall time per level (whether L2 fits a runner at all).

**T3. Flow state** (no extra runs): body vs wheels drag split, lift and the block means of every run; a run whose
split or lift sits apart from its design's other frames is flagged as a possible second flow state.

**T4. Does the transition model change the answer?** (`t4_*`, frame 0, kOmegaSSTLM)
- **Same answer** if delta(SSTLM) has the same sign as T1's delta and lies within a factor 2 of it. One frame only:
  indicative, not a verdict.

Results land in `audit/results/b1/<case>/result.json` (plus force histories), committed by the workflow's collect job.

## Log of deviations (added as they happen)

- 2026-11-27 20:40 UTC: L0 frame-0 cases died in `decomposePar` with SIGFPE on the runner (the same case decomposes cleanly
  locally; L0 team_rm frame 1 ran: 183.06 mN, 1.20 M cells). Relaunched as batch b1r2 with frames 0-1 retried plus frame 2
  for both designs as a backup. If a frame-0 case fails again, L0 uses frames 1-2 (frames are interchangeable repeat
  meshes; the mark itself is unchanged).

## Batch b2: does the original hybrid's gain survive legalisation? (added 2026-11-28 07:10 UTC, before its result)

`legal:hybrid_final`: `hybrid/final_design/car_half.stl` put through the legal builder (`seedhybrid/legalise.py`:
3-axis 6.35 mm ball-end mill, 3 mm walls, one piece, cargo, clearance, outside-region reset to the team CAD). Local
build: accepted; the rules changed 1382 mm3 inside the region and reset 456 mm3 outside it; pod 67.8 cm3.
Scored on frames 0-3 as every audit case. Reference: the team pod's legal rebuild (seed evolution's `original`,
174.75 mN mean over frames 0-3, same protocol), paired by frame.
- **Mark T5:** the legal hybrid keeps >= 70 % of the hybrid's fresh-mesh gain against the team pod surface it started
  from (b1 T1: -10.39 mN), i.e. legal hybrid minus the CAD team pod (b1, 181.95) <= -7.27 mN; and separately reported
  against the legal team pod (174.75): a legal design better than the legal team pod by >= 3 SE.
