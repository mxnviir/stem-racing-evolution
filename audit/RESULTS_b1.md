# Drag audit b1: results (2026-11-27)

Marks: `audit/README.md` (fixed before any result). Raw: `audit/results/{b1,b1r,b1r2}/<case>/result.json` and force
histories. All numbers are half-car drag in mN, production scheme, fresh wall-resolved mesh, cold start, iterations
200-1000. Paired difference = design minus its reference on the same meshing frame.

## T1. Hybrid final on fresh meshes: HOLDS

| frame | team pod (CAD surface, `team_rm`) | hybrid final | difference |
|---|---|---|---|
| k0 | 180.97 | 171.42 | -9.56 |
| k1 | 182.86 | 171.23 | -11.64 |
| k2 | 181.03 | 171.82 | -9.21 |
| k3 | 182.92 | 171.76 | -11.16 |
| mean | 181.95 | 171.56 | **-10.39 +- 0.59 (SE), 17.6 SE** |

Mark: <= -4.0 mN and >= 3 SE -> **holds**. On the morphed mesh the loop had reported -6.3 (gen 3) to -8.2 (gen 5)
under the production scheme; fresh meshes give more, not less.

Seed evolution on the same yardstick (against its own team pod, the field rebuild, 174.75 mean over k0-3):

| design | k0 | k1 | k2 | k3 | mean difference |
|---|---|---|---|---|---|
| d7 | 170.11 | 171.44 | 172.22 | 172.00 | -3.31 +- 0.39 (-3.06 without its selection run k0) |
| d10 | 171.17 | 171.65 | 172.60 | 173.57 | -2.51 +- 0.50 |

**Hybrid leads**: about 3x the confirmed gain against its own reference at similar compute. In absolute drag the two
end level (hybrid 171.56, d7 171.44), because the two team pods differ by 7.2 mN (CAD surface 181.95 vs field rebuild
174.75 on the same protocol). The hybrid design is not legalised (milling, walls, mass rules off).

Where the hybrid's gain comes from (scoring mesh, means of 4 frames): **wheels -6.85 +- 0.06**, body -3.54 +- 0.59.
Two thirds of it is drag taken off the wheels by the pod's shape, not the pod's own drag.

## T2. Is a ~10 mN difference resolved by the scoring mesh? PASS (narrowly; incomplete)

| level | bg divisor | cells | team pod | hybrid final | difference (paired) |
|---|---|---|---|---|---|
| L0 | 9.524 | 1.20 M | 182.75 (k1, k2) | 171.61 (k1, k2) | -11.14 (-11.20, -11.08) |
| L1 | 12 | 2.13 M | 181.95 (k0-3) | 171.56 (k0-3) | -10.39 (k0-1: -10.60) |
| L2 | 15.12 | 3.62 M | 182.39 (k1 only) | 173.22 (k0, k1) | -8.38 (k1 pair); -9.16 unpaired |

Mark: same sign, and |delta(L2) - delta(L1)| <= max(1.0, 0.25 |delta(L1)|) = 2.65 on frames 0-1:
|-8.38 - (-10.60)| = 2.22 -> **pass**. Narrow, and on one L2 pair: `t2L2_team_rm_k0` timed out at the 320-minute
step limit (L2 runs take ~4.9 h), and L0 frame 0 dies in `decomposePar` with SIGFPE on runners for both designs
(frames 1-2 used instead, logged in README before the results).

Reading: the sign and rough size are robust across a 3x range of cell count, but the magnitude drifts by about
1-1.5 mN per refinement step (-11.1, -10.4, -8.4/-9.2) and is not converged; no meaningful observed order or GCI
can be computed (the differences between levels do not shrink). Treat the hybrid's gain as **-8 to -11 mN**, and the
scoring mesh's resolution for design differences as **about +-2-3 mN**. Gains of ~3 mN (d7, d10) are not certified by
this ladder.

Absolute drag: the team pod is flat across levels (182.75 / 181.95 / 182.39); the hybrid rises at L2 (171.6 -> 173.2).

## T3. Flow states (no extra runs)

- Team pod, scoring mesh: two groups. k0/k2 read 181.0 (body 33.1 / 32.4, lift 82.8 / 82.9); k1/k3 read 182.9
  (body 34.3 / 34.6, lift 85.9 / 84.2). Wheels are the same in both (148.3-148.6). **Flagged: a second body-wake state,
  +1.9 mN.** The paired differences inherit it (-9.6/-9.2 vs -11.6/-11.2).
- Hybrid final, scoring mesh: one group (171.2-171.8, body 29.6-30.4). Its two L2 frames are 1.6 mN apart.
- The block SEs within a run (0.03-0.6 mN) are far below these frame-to-frame jumps: one run's own SE says nothing
  about which state it sits in.

## T4. Transition model (kOmegaSSTLM, frame 0): SAME ANSWER (indicative)

team pod 182.51, hybrid 172.38, difference **-10.14** (SST on the same frame: -9.56). Same sign, within a factor 2.
Split: wheels -8.07, body -2.07.

## What this changes

1. The hybrid (adjoint-led) loop produced the larger, mesh-robust drag gain; it is the approach to build on.
2. Start it from the field-rebuilt team pod (or legalise its result through the same build): the 7 mN gap between the two
   team-pod surfaces is as large as most gains, and the hybrid design still has to pass milling and wall rules.
3. Design differences below ~3 mN need more than the 2.1 M mesh and four frames; anything near that size should be
   checked on L2 or with more frames, with the flow state of each run labelled.
4. L2 costs ~4.9 h per run on a 4-core runner: one per job, no slack; L0 frame 0 needs a fix for the decomposePar SIGFPE.
