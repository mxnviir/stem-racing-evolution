# Methods tried, and why they were rejected or demoted

A record of every optimisation, modelling and measurement approach tried in this project before the current
pipeline, with the test that judged each one, the numbers, and the reason it was dropped. Written to serve as the
"approaches considered" section and appendix of a paper. Compiled 2026-11-28 from both repositories' documents, result
files and commit history; every number carries its source.

## How to read this

**Sources.** `adj:` = `mxnviir/stem-racing-adjoint` at `9ab9092` (`sp/` = `research/sidepod_opt/`).
`evo:` = `mxnviir/stem-racing-evolution` (`main`; `lk` = branch `claude/loving-keller-q72mvv` at `e9c896b`; the commit
`457821b` holds, squashed, everything the evolution repo inherited before 2026-11-22). `audit:` / `seedhybrid:` = this
branch (`claude/wizardly-brahmagupta-1edalc`). Line numbers are file lines; `ROUTES.md` entries are dated lab-log
entries.

**Verdict labels.**
- **E**: abandoned on evidence: a test with a stated pass mark failed, or a measurement showed it wrong.
- **O**: changed by the owner's decision (a design rule, not a test result).
- **S**: superseded by a later design, without a head-to-head test.
- **N**: never finished; no verdict.

**Units.** Drag in mN: half car unless "whole car" is stated. Race time is the race model's, not a track time: no car
has been built and timed (see section 7). "ms per mN" sensitivities are the race model's.

## Summary

| # | Approach | Area | Verdict | One-line reason |
|---|---|---|---|---|
| 1.1 | Lofted parametric body + pattern search | geometry | O, S | Owner's free-form rule; tunnel was a token feature |
| 1.2 | Carving the legal block along a level-set "compass" | geometry | S (argued) | Tunes the 13-22 % body share of drag; only ~1 mm local moves |
| 1.3 | Four separate carve code paths | geometry | E | A timed body was not always the judged body |
| 1.4 | Lattice / CPPN representations | geometry | O, E | Not free-form; CPPN mutations 27/40 no-ops, rest up to 48 % change |
| 1.5 | Parts-kit recipes, feature vocabulary | geometry | O | "Parameters can't describe organic shapes" |
| 1.6 | Binary-voxel shape loop v1 | geometry | E | 1.5 mm step flipped 18 voxels of a 67 cm3 pod |
| 1.7 | Fixed 8-station loft, box sidepods, cone nose | geometry | O | Replaced by free-form surfaces and an organic nose |
| 2.1 | Low-Reynolds lattice-Boltzmann tunnels | CFD | E | Re ~2x10^3 vs ~3x10^5: said no pod beats the team pod; RANS says the opposite |
| 2.2 | Coarse mesh (105 k) for decisions | CFD | E | Adjoint design's -10.4 mN became -0.8 +- 0.6 on medium |
| 2.3 | Medium mesh (520 k) as the scoring mesh | CFD | E | Mark M1 failed: groove -3.6 on medium is +2.2 on fine |
| 2.4 | Single fine-mesh runs (3.36 M) | CFD | E | Mark M2 failed: pooled sigma 2.28 mN; two flow states |
| 2.5 | Wall-function meshes without layers | CFD | demoted | y+ up to ~180 |
| 2.6 | Steady first-order model as the judge | CFD | demoted | +26 % drag; late polish gains transfer 0.49 |
| 2.7 | Ahmed-body validation | CFD | demoted | +9.5-13.3 % drag; run pages lost |
| 2.8 | SimScale meshing/import variants | CFD | E | Several settings failed to import or mesh |
| 2.9 | OpenFOAM on Actions (evolution repo) | CFD | S | Moved to SimScale without a head-to-head |
| 3.1 | Fresh-mesh / morph noise ladder | measurement | E | 8.65 ms fresh -> 0.091 ms "method A" (precision only) |
| 3.2 | Single frozen reference state | measurement | E | Two states rank cars differently (Spearman 0.77) |
| 3.3 | Repeat meshes by nudging the design | measurement | E | Nudges moved surfaces 0.2-1.5 mm: not repeats |
| 3.4 | Plain repeats; SIMPLE early stop; 4-core runs | measurement | E | Measure nothing / hid noise / not reproducible |
| 4.1 | ESI continuous adjoint (all tests) | gradient | E | Wing r 0.42-0.70; car wrong signs 2/3; no better than random |
| 4.2 | Measured map ("compass") | gradient | demoted | Big steps raised drag; car test never judged |
| 4.3 | DAFoam beyond coarse meshes | gradient | E (resources) | ~150 kB/cell: medium ~128 GB, fine ~500 GB |
| 4.4 | Coarse adjoint shape loop | gradient | E | -33 mN steady (whole car); no established gain on fine |
| 4.5 | Porous topology optimisation | gradient | E | 18 % better in its model, 11 % worse in real CFD |
| 4.6 | In-loop ("memetic") adjoint | gradient | E | Bench T3b/T3c failed: adjoint best at the end |
| 4.7 | DAFoam on shapes unlike the team pod | gradient | E | Converged on 1 of 6 |
| 5.1 | Optimiser benchmark | search | N | Results never recorded |
| 5.2 | Seven-method LBM competition | search | S | Run on the discredited tunnel; no method beat a hand shape |
| 5.3 | 2-D organic-feasibility study | search | N | No written verdict |
| 5.4 | Sidepod CMA-ES on SimScale | search | E | 105 cars; never beat the team car (+2.6 ms) |
| 5.5 | CMA-ES as the gradient-free polish | search | E | TC2: 7 of 15 against a mark of 11 |
| 5.6 | Layout CMA-ES (49 genes) | search | S | Replaced by a 10-number Bayesian search the same day |
| 5.7 | Seed evolution, drag-only scoring | search | E (flagged) | d7 is ~10 ms slower in race time (heavier) |
| 6.x | Invented thrust curve; fixed 48.2 g sizing; "COM doesn't matter" | race model | E, O | See section 6 |
| 7.x | Track calibration | race model | not done | Every time is a model number |
| 8.x | Wing and support designs | parts | O, E | See section 8 |
| 9.x | Rule-check gaps | legality | E | Search cars unchecked on 15 rules before 2026-11-23 |

---

## 1. Geometry representation

**1.1 Lofted parametric body searched by pattern search** (O, later S).
- *What:* an 8-station superellipse loft with sidepods, a token underfloor tunnel and 56 "sculpt" modes, tuned by a
  pattern search on direct CFD (team decision 2026-09-28, `adj:part5-simulation/pattern.py:3-8`).
- *Why a pattern search:* the first adjoint check had the right sign on only 3 of 8 settings (r = 0.45), while direct
  CFD of a 2 mm sculpt mode moved drag by up to 2.1 % against 0.08 % noise (`pattern.py:4-6`).
- *Why dropped:* replaced by carving from the block on 2026-11-08/09 (`adj:part1-simulation/README.md:25-28`). The tunnel
  was a token one (`ft_h` <= 10 mm, off by default, `adj:research/APPROACH.md:36-38`); the 56 modes were frozen as "too
  many for ~110 runs" (`evo:main docs/PARAMETERS.md` L79-135). No pattern-search results survive in either repository.

**1.2 Carving the full legal block along a level-set sensitivity "compass"** (S, on argument).
- *What:* the body cut from the legal block in small legal steps along a race-time sensitivity map (Parts 1, 3, 5).
- *Judged by:* the 2026-11-19 review (`adj:research/APPROACH.md:9-38`), not an experiment:
  - it tunes the smallest share of drag: wheels are 64-75 % of the car's drag, the body 13-22 %
    (`adj:part4-simulation/assembly.py:24`);
  - each cut stays within ~1.5 mode heights (~1 mm); 25-40 rounds x 3 chains reach only the block plus ~25-40 mm of
    local moves. A tunnel or real sidepods are 10-20 mm away everywhere at once.
- Part 3's car-level tests were still "running, none judged" when it was dropped (`adj:part3-simulation/PROOF.md`).

**1.3 Four inconsistent carve code paths** (E).
One of four paths added Gaussians straight onto the field, skipping the rule zones and the milling: "a body that was
timed was not always the body that was judged" (`adj:part1-simulation/PROOF.md:60-73`). The first 360 cars of the
2,000-body gate were discarded (`PROOF.md:23`); 31 Part-3 mode runs were cancelled (`adj:part3-simulation/PROOF.md:125-129`).
Replaced by one `carve.py`. A related bug: until 2026-11-17 `build()` silently dropped the underfloor tunnel, so "Part
3's five floor variables did nothing" (`adj:part1-simulation/tests/test_min_wall.py:356-358`).

**1.4 Lattice and CPPN representations** (O + E).
- *What:* search framework `sp/methods.py` (CMA-ES, Bayesian optimisation, SAIL, NEAT, SPSA) over an 8x3x4 or 72-number
  control lattice, or a CPPN.
- *Owner:* lattice/CPPN parameters are not free-form (`sp/DECISIONS.md` D3, lines 51-57).
- *Evidence:* the lattice "bulges a smooth blob", its points 10-18 mm apart cannot make ribs or windows
  (`sp/seedevo/SEED_EVOLUTION_PLAN.md:24`). CPPN-NEAT as a genome: in a 3-D heritability probe 27 of 40 mutations changed
  nothing and the rest changed the body by up to 48 % (minimum IoU 0.517)
  (`adj:research/organic_feasibility/results3d/herit.json`). Demoted to an operator only (`SEED_EVOLUTION_PLAN.md:23`).

**1.5 Parts-kit recipes with quality-diversity search; a "feature vocabulary"** (O). Proposed 2026-11-21
(`sp/ROUTES.md:609-652`); rejected: "parameters can't describe organic shapes" (`sp/HANDOVER.md:158`). Not built.

**1.6 Shape loop v1: binary voxel field, max-normalised step, drag-only objective** (E). Dry run 2026-11-20 18:45 with
the real gradient: a 1.5 mm step flipped 18 voxels of a 67,000 mm3 pod; 0.75 mm flipped none (`sp/ROUTES.md:230-239`).
Replaced by v2 (sub-voxel signed distance, 95th-percentile step, race-time objective), in which 1.5 mm moves 44 % of
the surface band (`sp/ROUTES.md:250-268`).

**1.7 Fixed loft with box sidepods; cone nose with three optimiser numbers** (O).
- Loft: 8 stations, one squareness 2.6, constant-width box sidepods (28 mm wide, 27 mm flat top) (`evo:main
  docs/PARAMETERS.md` L33-63). Replaced by three free-form surfaces (`evo:car/body/freeform.py`, commit 20bd6c4,
  2026-11-23; owner: "no constant-width sidepods, fixed stations, squareness or cone nose"). `fit()` reproduces a car at
  IoU 0.85 on a 24x9x12 grid; the floor grid was raised to 4 points across because with 2 "the floor could only tilt
  into a V" (f98bba5).
- Cone nose: `nose_len`, `nose_k`, `nose_tip_z` in the stage-1 model; after 39 timed cars their length scales were
  3.79 / 2.96 / 1.66 (`evo:main search/layout/results/state.json`, `ask.length_scales`): the model saw an effect of the tip
  height. Replaced on `lk` by a nose grown from the body's own Ref A section (eed0ca6, 2026-11-28; owner). **No CFD
  comparison of loft vs free-form or cone vs organic nose was made.**
- Manufacturing limit found along the way (E): closed lengthwise tunnels are unbuildable under 3-axis
  top/bottom/left/right milling; a 5 mm bore was refilled, a 6 mm groove refilled, an 8 mm outboard groove buildable
  (`sp/ROUTES.md:203-208`).

## 2. CFD fidelity and meshes

**2.1 Low-Reynolds lattice-Boltzmann tunnels** (E).
- *What:* 3-D D3Q19 LBM with LES, 3 mm cells, Re ~ 2x10^3, 157x36x32 cells, with a differentiable JAX copy for
  topology optimisation (`adj:research/sidepod_region/README.md:60-65`); a 2-D version at Re ~ 290
  (`adj:research/organic_feasibility/lbm2d.py:3-4`).
- *Why rejected:* "ranked sidepods wrongly" (`sidepod_region/README.md:1-4`). The documents do not quantify it; the raw
  files do: the LBM gives the team pod 390.0 mN and no pod (cargo only) 376.6 mN, i.e. no pod is better
  (`sidepod_region/results/old_fixed.json`), while real-wall RANS gives the team pod 0.362 N and no pod 0.411 N (whole
  car, `sp/ROUTES.md:76-78`). Reynolds number ~100x too low; replaced by full-Reynolds RANS only.

**2.2 The coarse mesh (~105 k cells) for any design decision** (E; decision D6, `sp/DECISIONS.md:76-81`).
The coarse adjoint design's -10.4 mN became -1.5 mN on medium in one draw (`sp/ROUTES.md:690-741`) and -0.8 +- 0.6 mN on
clean repeats (`sp/ROUTES.md:985-991`). A pre-registered control (mark: <= -7 mN from a cold start means resolution is
the cause) gave -10.3 cold vs -10.4 warm (`sp/ROUTES.md:742-766`). Team-pod body pressure drag: 19.8 (coarse), 8.4
(medium), 2.6 mN (fine) (`sp/ROUTES.md:1021-1035`). Three claims were withdrawn ("design 8 is 4.5-5.7 % better", "early
steps transfer 81-89 %", "coarse scatter 0.32 mN resolves 1.5 mN": precision, not accuracy; `sp/ROUTES.md:755-812`).

**2.3 The medium mesh (~520 k) as the scoring mesh** (E).
Mark M1 (fixed 2026-11-21 13:00, `sp/fine_campaign/README.md:61-76`): team body drag within 3 mN of fine and every design
difference within max(2 mN, 3 SE). Result 2026-11-23: **fail**. Team body drag medium 23.95 vs fine 18.77 (-5.2 mN);
groove -3.6 +- 0.5 on medium became +2.2 +- 1.9 on fine; the no-pod gap shifted +13.9 mN against a limit of 8.3
(`sp/ROUTES.md:1240-1262, 1399-1417`; D8).

**2.4 Single fine-mesh runs (3.36 M) as a score** (E). Mark M2 (pooled sigma <= 2 mN): fail at 2.28 mN. The team pod on
three fine meshes read 176.04 / 179.86 / 173.76 mN (sd 3.08): two converged flow states (`sp/ROUTES.md:1070-1097, 1218-1254`).

**2.5 Wall-function meshes with no boundary-layer cells** (demoted). y+ up to ~180 (`sp/ROUTES.md:898-911`); on coarse,
layers formed on 0 of 5,666 wall faces (`sp/ROUTES.md:936-943`). The "resolved" mesh came out at 2,132,871 cells rather
than the planned 4-5 M, with layers on 62 % of wall faces (`sp/ROUTES.md:1099-1112`). A first wall-resolved DAFoam
boundary check "passed" on a mesh with no layers and was declared invalid (`sp/ROUTES.md:898-907`).

**2.6 The steady first-order model as the judge of designs** (demoted). Needed because the production scheme never
settles (force range 2.80 mN vs 0.02 mN for upwind), but it reads +26 % drag (228.72 vs 181.04 mN,
`sp/ROUTES.md:168-182`; D10). Late polish did not transfer: the ratio of production to steady gain fell from 0.81 to
0.49 over designs 1-8 (`sp/ROUTES.md:398-420`). Kept as "a good compass for big moves, poor judge of fine polish".

**2.7 Ahmed-body validation** (demoted). Drag above Lienhart's measurement: 25 deg fine +9.5 %, wall-resolved +13.3 %;
35 deg +5.4 %; the cube case unfinished (`adj:part2-simulation/PROOF.md:213-222`). Replaced by a fixed correction
D/1.09 +- 0.04. The run pages are on suspended accounts and cannot be opened (`PROOF.md:222`). It runs at about 15 times
the car's Reynolds number.

**2.8 SimScale meshing and import variants** (E), all `evo:`:
- Sidepod route: a seam on the region boundary with 0.05 mm overlap: "about 2 in 3 imports refused"; adopted a seam
  0.3 mm off the region: 3 of 3 imported first time (`457821b:research/sidepod_opt/EVOLUTION_SIMSCALE.md` §5).
- Automatic sizing gave 2.90 M cells on the team car but 6.12-8.20 M on evolved pods and "failed often"; manual 30 /
  0.8 mm: 3.48-4.04 M, all nine meshes first time (§3.3).
- Whole car: small-feature suppression 0.05 mm failed on 5 of 9 cars; 0.8 mm meshed first time at 5.0 M cells
  (`evo:docs/RUNBOOK.md`, "Meshing"; 321c1c5). The setting shifted the team layout's result from 0.4483 N / 0.76519 s to
  0.4009 N / 0.75642 s (`search/layout/results/state.json` refs). 0.05 mm results are not used.
- A finer 1.0 mm refinement box was rejected on cost: team car 0.3823 N (3.51 M cells, 25.8 core-hours) vs 0.3634 N
  (6.79 M, 43.7) and 0.3585 N (6.28 M, 41.7); car c00004 0.4098 vs 0.4061 N; ranking kept on these **2 cars only**
  (`search/sidepod/results/validate/results.jsonl`).
- Direct decimation replaced by a 0.3 mm shrink-wrap rebuild (needle triangles and 0.03 mm ledges rejected by the
  importer; 7531977); no failure counts recorded.

**2.9 OpenFOAM on GitHub Actions in the evolution repo** (S). Removed when the search moved to SimScale (9efa291,
a4b8f00, 2026-11-22) without a head-to-head comparison. A planned "recheck SimScale winners in OpenFOAM" stage was never
dispatched (`457821b:research/sidepod_opt/EVOLUTION_SIMSCALE.md` §7, §10).

## 3. Noise and measurement

**3.1 The stopwatch noise ladder** (Part 2, medium mesh; mark sigma <= 0.3 ms; `adj:part2-simulation/PROOF.md:77-171`).

| Set-up | sigma (95 % interval) | Verdict |
|---|---|---|
| fresh mesh, 1000 iterations | 8.65 ms (6.27-13.93) | rejected |
| fresh mesh, 3000 iterations | 5.02 ms | rejected |
| fresh mesh + warm start | 6.88 ms | rejected |
| reshaped mesh, plain projection (folded a face on every cut) | 4.09 ms | rejected |
| reshaped mesh + warm start | 0.78 ms | rejected |
| reshaped mesh only | 0.57 ms | rejected |
| reshaped from reference 3.978 | 5.52 ms (drag jumping between 0.413 / 0.428 / 0.442 / 0.454 N) | rejected |
| "method A": both frozen states, mean | **0.091 ms** (0.067-0.144), confirmed 0.087 | kept, then demoted |

Diagnosis: the wake behind the body and rear wheel has more than one steady solution (`PROOF.md:102-105`). Demoted on
2026-11-19 (`adj:research/APPROACH.md:40-53`): the low noise holds only for shapes morphed from one frozen reference; a
car with a tunnel falls back to 5-9 ms; averaging two states "hides a real physical effect". The audit later measured
precision vs accuracy directly (section 10).

**3.2 Timing from a single frozen state** (E). T(B) - T(A) over 10 cars: rms 0.50 ms, max 1.11 ms, Spearman 0.77
between states (`PROOF.md:169`). Replaced by method A.

**3.3 Repeat meshes made by nudging the design +-0.05 mm** (E). The nudges moved surfaces by 0.2-1.5 mm (no-pod +0.05:
1.33 mm at p99, 1.49 mm max); one "repeat" was a different design by -1.0 cm3 (`sp/ROUTES.md:944-962`). Medium sigma
1.52 mN, groove -3.9 and design 8 -3.5 were withdrawn (`sp/ROUTES.md:998-999`). Replaced by meshing-frame shifts of
k x 0.37 mm (D7): clean medium sigma 0.74-0.75 mN.

**3.4 Other measurement shortcuts** (E).
- Plain repeats of an identical STL return the identical value (180.85 mN to the digit): "measure nothing"
  (`sp/ROUTES.md:656-657`).
- SIMPLE early stop: the earlier 0.28 ms repeatability came from solves ending near iteration 170; without it the
  six-machine spread on fresh meshes was 0.94 ms (`PROOF.md:144-149`).
- 4-core runs (mark spread <= 0.1 ms): 0.29-0.30 ms; 0.17 with a fixed MPI order; 2 cores 0.01 ms. One core was not
  machine-independent either (1.3 ms spread, traced to the two flow states) (`PROOF.md:158-172`).
- Early wall-layer settings: sigma 1.3 / 2.2 / 1.0 ms with layers / "better" layers / none; the search mesh went
  layer-free (`adj:part5-simulation/solver.py:52-60`).
- Seed evolution's production score including the cold-start transient: block SE 5.2 / 5.6 mN vs 0.10 / 0.19 mN on
  iterations 200-1000 (`sp/seedevo/SEED_EVOLUTION_PLAN.md:431-438`).
- Evolution repo's noise reference "the same car at two face counts": both landed on 40,000 faces and the same mesh
  (identical 0.4482511 N); the replacement diverged; the pair used gives 0.54 ms, against the stage-1 model's fitted
  noise 0.905 ms (`evo:search/layout/results/state.json`; b803c8c, 098e726).

## 4. Gradient methods (adjoints)

**4.1 ESI (OpenFOAM's) continuous adjoint** (E, in every test).

| Test | Mark | Result | Source |
|---|---|---|---|
| Part 3 wing, known answer | sign >= 9/10, r >= 0.9 | official bumps r 0.53; 20 single bumps r 0.70, signs 75 %; SA setup r 0.42-0.55; 4x finer mesh r ~ -0.5; simple surface formula ~66x too large | `adj:part3-simulation/PROOF.md:5-37` |
| Car side study | physical plausibility | says growing the whole car evenly lowers drag: "physically wrong" | `PROOF.md:66-73` |
| First sidepod gate (2026-11-20) | Part 3's marks | r -0.54, signs 0/8 (matched r -0.47, 4/8); confounded: adjoint on 375,068 cells, bumps on 519,620 | `sp/gate_results/RESULT.md:85-118` |
| Coarse steady flow re-test | prediction within 0.5-2 | signs 3/3, r 0.98, but predicted/measured 0.09 / 0.02 / 0.24 | `sp/ROUTES.md:861-935` |
| Wall-resolved mesh (amended heights, 2026-11-24) | r >= 0.9, all signs, ratio 0.5-2 | measured -1.177 / -0.637 / -0.255; predicted (frozen) -0.078 / +0.031 / +0.096: 2 of 3 signs wrong, ratio <= 0.43 | `sp/ROUTES.md:1468-1497`; `sp/esi/RESOLVED_TEST.md` |
| As a step proposer (0.125 mm) | beats both random smooth steps | ESI -0.644 mN; random seed 1 -1.043, seed 2 +0.408: fail | `sp/esi/STEP_TEST.md:89-117` |

Diagnosis: ESI's normal sensitivity correlates 0.18 with DAFoam's and is 1-4 % of its size
(`sp/ROUTES.md:1499-1528`). It was later used in the hybrid loop as a CFD-checked proposal and gave one accepted step
(-1.28 mN, generation 2); removed from the seeded campaign on 2026-11-28 by the owner (`seedhybrid/README.md` change log).

**4.2 The measured map ("compass")** (demoted). 80 Gaussian modes fitted to measured bumps: r 0.977, signs 95 %; outside
its family r 0.952 (`adj:part3-simulation/README.md:50-56`). Wing step rehearsal: a 3-height step raised drag (+0.049
measured vs -0.244 predicted); steps <= 1.5 heights helped, best -0.67 % (`PROOF.md:53-64`). Car test "running, none
judged" (`PROOF.md:18`). Demoted to final polish only (`APPROACH.md:131-133`).

**4.3 DAFoam (discrete adjoint) on the production flow and on finer meshes** (E on resources). On medium the preconditioner needed
> 15 GB and the colouring > 40 GB (`sp/ROUTES.md:142-144`); on the coarse production state GMRES stalled at 1.89e-2;
the fixed-point adjoint supports only Spalart-Allmaras (`sp/ROUTES.md:142-167`). About 150 kB per cell: medium ~128 GB,
fine ~500 GB (`sp/fine_campaign/README.md:23-26`). Where it did run (coarse steady flow) it passed: predicted/measured
1.09 / 0.66 / 0.97, r 0.98 (`sp/ROUTES.md:285-296`). It is the adjoint the current hybrid uses.

**4.4 The coarse adjoint shape loop** (E, demoted). `cloud_original`: -33.2 mN whole car (-7.4 %) on the steady model
(`sp/ROUTES.md:359-386`), but -0.8 +- 0.6 on medium and -1.2 mN (z = -0.45) on fine (`sp/ROUTES.md:1284-1301`).
`cloud_lean` (from no pod): -26.0 mN steady, -2.2 mN on the production scheme, recovering 9 % of the pod's 25 mN: it
"improves a pod far better than it invents one" (`sp/ROUTES.md:665-682`). The selected route was superseded on
2026-11-24 (`sp/ROUTES.md:9`).

**4.5 Porous (Brinkman) topology optimisation, ESI** (E). From an empty region it stalled at the seed skin; its gradient
check was invalid (units bug fixed 2026-11-20; the topology adjoint itself stayed untested, `sp/ROUTES.md:73, 128-139`).
CMA-ES scored by the porous model found a design 18 % better than empty in that model (J 2.181 -> 1.782; team 1.624) but
**11 % worse in real-wall CFD**: 0.456 N vs no pod 0.411 N and team 0.362 N (`sp/topo/mill_final.json`,
`mill_team_porous.json`; `sp/ROUTES.md:74`). Replaced by real-wall CFD for every score (D4).

**4.6 The adjoint inside the evolution loop ("memetic")** (E, on the seed-evolution test bench).
T3b (adjoint inside the loop beats evolve-then-polish) failed on 2 of 3 held-out instances; T3c (a useless or
misleading adjoint costs no more than its share) failed; the adjoint was moved to the end
(`sp/seedevo/SEED_EVOLUTION_RESULTS.md:35-60`). The trust check T3e failed (trust 32 % with a good adjoint; refused
morphs counted as distrust), and T3e' after the fix by one pair (9/15 vs 10) (`RESULTS.md:66-119`). A stand-in problem,
not aerodynamics.

**4.7 DAFoam on shapes unlike the team pod** (E, for lack of a gradient). Converged on 1 of 6 attempts; the channel seed
stalled at 8.45e-4; four runners were lost to memory; a 2-rank, restart-100 retry stalled at 1.024e-2. Where it
converged it overstated the gain ~2.5x (-7.50 predicted, -2.64 measured) (`RESULTS.md:232-255, 285-288`;
`sp/seedevo/GO_NO_GO.md:20-25`). Those shapes were polished without a gradient. (The seeded hybrid campaign retests this
on three seeds: `seedhybrid/PROTOCOL.md` E0.)

**4.8 Hybrid loop details that were demoted** (E). H2 (combined directions beat evolution-only at 0.125 mm): fail,
evolution-only -1.84 vs guided -1.40 (`sp/hybrid/R2_VERDICT.txt`); at 0.5 mm, post hoc, guided -5.24, DAFoam -4.84,
evolution-only -2.62 (`sp/hybrid/HYBRID_TEST.md:192-198`). DAFoam step sizes ~2x too large (slope 0.48, `R1_VERDICT.txt`).
The loop stopped when the chained morph ran out at ~2 mm (`HYBRID_TEST.md:284-288`); its final design was never
legalised (`sp/hybrid/final_design/README.md:10`).

## 5. Optimisers and search methods

**5.1 The optimiser benchmark** (N). A stand-in problem tuned to the repo's noise numbers, comparing CMA-ES, GA,
Bayesian optimisation, TuRBO, multi-fidelity co-kriging, the measured map and an adjoint (`adj:research/optimiser_bench/bench.py`).
`sweep_out.txt` stops at "100 / 198"; the results section of `APPROACH.md` (§5) was never filled in. Its verdicts in
`APPROACH.md:80-88` were written without benchmark numbers. The audit's partial re-run (4 seeds) found the problem's
"true best" mis-computed (a negative gap observed) and that an adjoint with r = 0.55 beat r = 0.97 in one scenario
(`audit` first report).

**5.2 Seven-method competition on the LBM bench** (S). 30 tunnel runs each; best race time (old thrust sheet):
random 1.19228 / 1.19132 s, Bayesian 1.18912 / 1.18848, CMA-ES 1.18502 / 1.18497, SAIL 1.18503, SPSA 1.18339,
CPPN-NEAT 1.18176, adjoint topology 1.19732; the original car 1.18830, a hand "cargo only" shape 1.18176
(`adj:research/sidepod_region/results/old_cargo_low/runs/*.jsonl`, `old_fixed.json`). No method beat the hand shape
(NEAT's best is that seed network). The run set was discarded (cargo placed too low) and then the tunnel itself.

**5.3 2-D organic-feasibility study** (N). 200 LBM runs per method, drag (lower is better): designer-style profile
parameters 0.6126 / 0.6160, lattice CMA 0.6304 / 0.6243, CPPN-NEAT 0.6764 / 0.9585 / 0.6843, random 0.7146 / 0.7358;
references: minimal bridge 0.669, teardrop 0.779, wheels only 0.336 (`adj:research/organic_feasibility/results/`). No
written verdict; Re ~ 290.

**5.4 Sidepod CMA-ES on SimScale** (E). Lattice 6x3x4, popsize 8, sigma 0.5, amplitude 12 mm
(`evo:search/sidepod/config.simscale.json`). From the team pod: 105 cars, 14 illegal, best 0.75580 s against the team
car's 0.75316 s (+2.6 ms), plateaued by about the fifth generation, "frozen at 105 runs after plateauing"
(`evo:search/sidepod/results/state.json`). From scratch (amplitude 12 mm): best 0.76585 s, "could only grow a skin within
~12 mm of the spine"; amplitude 24 mm: 32 cars, 18 illegal, best 0.76217 s, stopped to free the keys. About 1,718
core-hours estimated over the 96 timed cars of the main run (sum of `cost_estimate`). It ran without the T4.2 virtual
cargo. Stated probable causes (untested): ~6 mm random steps, a fresh mesh per car, masks applied after sampling
(`sp/hybrid/EVOLUTION_ADJOINT_PLAN.md:56-60`).

**5.5 CMA-ES as the gradient-free polish** (E). Mark TC2 >= 11 of 15 pairs; result 7 of 15
(`sp/seedevo/SEED_EVOLUTION_RESULTS.md:159-160`). Kept in the code, off.

**5.6 Layout search by CMA-ES over 49 genes** (S). Replaced on 2026-11-22 by a Gaussian-process search over 10 numbers
(7531977): at 20-25 core-hours a run, "a few hundred runs is the whole budget" (`evo:docs/SEARCH.md`). Never run.

**5.7 Seed evolution scored on drag only** (E, flagged). Its best design d7 is 3.3 mN better than its seed in drag
(audit T1: -3.31 +- 0.39 mN over 4 frames) but weighs 49.52 g, about 10 ms slower in race time
(`sp/seedevo/SEED_EVOLUTION_RESULTS.md:295-301`).

**5.8 Smaller rejected choices in seed evolution** (E). Timid mutation operators (mark: >= 70 % in-band children;
combine 38 %, bump 60 %, cut 69 %, fill 33 %) made bolder (`RESULTS.md:14-29, 121-131`); scoring small children by
morphing showed "no consistent effect" and was replaced by a fresh mesh each (`RESULTS.md:144-147`); an exact-shape
rebuild check rejected almost every child because the build is not bit-identical (commit f0f7833); float16 field storage
broke a rebuild (34 mm3 thin foam; commit 7724f7d).

## 6. Race model and mass

- **Invented thrust curve** (E): sensitivities 0.48 ms/mN and 16.6 ms/g, replaced by the measured thrust sheet: 0.199
  ms/mN and 7.51 ms/g (`adj:part2-simulation/PROOF.md:196`). The LBM bench's 0.45 ms/mN, 16.4 ms/g and -1.5 ms per mm of
  COM height were "wrong" (`sp/README.md:50-56`). (The audit found PROOF's 0.199 predates the 1.09 drag correction: the
  current code gives 0.180-0.181 ms/mN.)
- **Sizing every car to 48.2 g with no ballast** (O): replaced on 2026-11-22 by a 48.0 g legal minimum with the race model
  choosing the lead (D2); on four re-scored cars it chose the minimum lead (`sp/ROUTES.md:1054-1069`).
- **"The centre of mass doesn't matter"** (E, corrected): true only inside a flat zone; outside it 4-22 ms cliffs from
  wheel lift and guide friction (`sp/ROUTES.md:525-549`).
- **Lift left out of the adjoint**: -0.002 vs 0.188 s/N (`sp/ROUTES.md:489-490`).
- **Mass bookkeeping** (open): non-body parts 22.52 g in one pipeline vs 19.23 g itemised in the other
  (`sp/ROUTES.md:1135-1137`).

## 7. Calibration against reality (not done)

- `TIME_COEFFICIENT` = 0.802 was set by the team lead on 2026-11-09 as "a calibration factor ... not physics"
  (`physics_contract.py:78-81`). The raw model gives 0.946 s, calibrated 0.759 s; real top cars run 1.0-1.1 s
  (`PROOF.md:196`): the factor moves the model away from real times.
- The evolution repo's "calibration against timed cars" (2026-11-27) recomputed the race time of all 39 CFD-timed cars from
  their stored inputs and found agreement to 6e-11 s (`evo:lk docs/checks/STATUS_2026-11-27.md`). That checks the code
  reproduces itself, not the physics.
- Placeholders still in the model: rolling friction 0.010, guide friction 0.25, tether drag absent, wheel inertia from
  CAD (`adj:part2-simulation/PROOF.md:224-236`). **No car has been built or track-timed.** The owner is measuring these.

## 8. Parts

- **Team v2 rear wheel-support CAD** (E): bottoms out at 1.40 mm, 0.10 mm under T3.7, and its 12.7 mm disc blocks the T7.13
  claw; replaced by a parametric beam support (`adj:part4-simulation/README.md:66-70`).
- **CFD car without wheel supports** (E): supports add 10-30 mN on coarse (180.7 without, 191.4-210.5 with; D16).
- **Beam supports in the evolution repo** (S): a fixed 5.9 g with a T7.12 waiver; replaced by per-car `support_engine.py`
  (5679660). An arm rising 6-10 mm "breaking the underside line" was replaced by a plate flush with the floor (owner,
  7f454d7; flush 0.004-0.058 mm on six layouts).
- **Rear wing on a fixed centre pylon** (O): all 39 timed layout cars had it; replaced by the rear-wing designer
  (8a53799), then restricted by owner rules (no centre-line strut; straight wings: 1131658, 7716850).
- **Front wing** (O, then E): a sequence on 2026-11-23 (58f5cd4 -> 1131658 -> 7716850 -> 62a3bb9 -> f98bba5 -> 645dd3c ->
  e4eb3f2 -> e24fe11): pylons, presets removed, lofted supports, supports moved to the wing's sides (T1.6, appendix iv),
  organic struts, free mid point, hull-union solid (E: a stitched skin folded through itself), then the "strut must not
  run back into its own span" check (E: a split wing lost 1.5 mm of T8.6.1's segment). That check made 7 of the 14
  rebuilt-illegal timed layouts fail (front wing unsupported), diagnosed on L0044 (5-8.5 mm3 in the span box), and was
  removed for the front wing (eed0ca6). Front-wing chord and strut shape then fixed (owner, b867574); underside pylons
  added as the fallback (e9c896b).
- **Rear tip-only struts** (E): two test layouts with half spans 37.5 and 38.6 mm were unsupported because "every tip
  strut left T3.4's 42.5 mm"; a pylon at 80 % of the span passes all gates (e9c896b).
- **Part 4 remesh** (E): turned 0.88 % thin triangles into 12.4 % (`adj:part4-simulation/PROOF.md:100`).

## 9. Rules and legality

- **Manufactured-file checks never ran on search cars** (E): until 2026-11-23 05:31, T4.2, T4.4.2-4, T4.6, T4.7, T5.5,
  T6.1, T7.9.1-4, T7.10, T7.11, T8.7, T8.8, T9.6, T9.7 and the block check ran only on manufacturing files, which the
  searches never wrote: "none of them was measured on a search car" (9b496b1).
- **T9.6 joint check** (E): a pylon's top buried in the wing was missed; the team car read 4.51 mm and failed; 5.51 mm
  with the fix (c4db404, `car/build/tests/test_wing_joint.py`).
- **Rebuild of the 36 timed layouts on today's builder** (E): 22 legal, 14 not (front-wing support 7, block 5, T7.9.4 2,
  T7.3 1, halo 1). The fastest that still builds legal is L0039, 0.75721 s; the previous best L0044 (0.75503 s) is
  illegal (`evo:lk docs/checks/rebuild_2026-11-27.jsonl`).
- **Thresholds chosen, not derived**: T4.1 specks < 20 mm3; block specks < 1 mm3 waived where the free-form floor grazes
  the block (919bce0, e9c896b).

## 10. Findings of the independent drag audit (2026-11-27/28, this branch)

These are not abandoned methods but tests of the surviving ones; they decided what the current pipeline keeps.
All on the audit protocol: fresh wall-resolved mesh (~2.1 M cells), production scheme, 1000 iterations, iterations
200-1000, frames shifted 0.37 mm (`audit/RESULTS_b1.md`).

| Test | Result | Consequence |
|---|---|---|
| T1: the original hybrid's final design on fresh meshes | -10.39 +- 0.59 mN vs the CAD team pod (17.6 SE); seed evolution's d7 -3.31 +- 0.39 vs its own reference | The hybrid's gain is real and ~3x seed evolution's |
| T1: the two "team pod" surfaces | CAD surface 181.95 vs legal field rebuild 174.75 mN | Always compare with a design's own reference; start from legal seeds |
| T2: mesh ladder 1.2 / 2.1 / 3.6 M cells | difference -11.1 / -10.4 / -8.4 (one pair at 3.6 M) | ~10 mN resolved in sign and rough size; resolution about +-2-3 mN |
| T3: flow states | team pod two groups 1.9 mN apart (body drag); hybrid one group | Label states per run |
| T4: transition model (one frame) | -10.14 vs -9.56 (SST) | Same answer |
| Gain breakdown | wheels -6.85 +- 0.06, body -3.54 +- 0.59 | Two thirds of the gain is wheel drag the pod shape removes |
| b2: the hybrid's final design legalised | -7.62 +- 0.85 mN vs the CAD team pod, but **-0.43 +- 0.53 vs the team pod legalised the same way**; wheel-drag gain gone (+12.7 mN vs unlegalised) | The original hybrid's shape gain does not survive legalisation (`audit/RESULTS_b2.md`); b3 separates out-of-region moves from milling loss |

The audit also found, in code: a velocity-extension bug in the level-set update (normals built from magnitudes,
`adj:part1-simulation/phi_updater.py:339-369`), seed evolution's driver dropping children and forcing `converged = True`
(`adj:sp/seedevo/campaign.py:187`, `campaign_case.py:61`), and a halved linearity metric (`adj:part3-simulation/checks/known_answer.py:431`).

## 11. Where the evidence is weakest

- No written results: the optimiser benchmark, the LBM method comparison, the 2-D feasibility study.
- Lost raw data: the Ahmed runs (suspended accounts); the SimScale CMA-ES branch's raw data are in the evolution repo,
  not the adjoint repo.
- Few runs: one fine-mesh run for several designs; three for the team pod; the 1.0 mm SimScale box decision on 2 cars;
  the supports correction on 3 pairs; the evolution repo's noise reference on one pair.
- Owner decisions with no CFD comparison: loft -> free-form, cone -> organic nose, straight wings, fixed cargo and nose,
  fixed front-wing chord and struts.
- Bench verdicts (seed evolution T3, TF, TC) come from a stand-in problem, not aerodynamics.
- Every race time is a model number; the model has never met a timed car.
