# evo_patches: mass bookkeeping fixes for the pinned evolution builder

Pinned builder: `mxnviir/stem-racing-evolution` @ `e9c896b` (branch `claude/loving-keller-q72mvv`, read only).

    python evo_patches/apply.py <evolution checkout>            # idempotent; --check only reports
    EVO_CHECKOUT=<patched checkout> [CAR_DIR=<built car>] [CAR_DIR_SHRUNK=<same layout, body shrunk>] \
        python -m pytest evo_patches/tests -q

Each patch asserts its exact anchor text (or its already-patched text) and replaces it once, with a comment
`evo_patches Pn` in the patched file. The race model (`car/race`, `search/layout/layout_car.mass`,
`simscale_layout.race_time`) is not touched, and no load, rule or margin is changed.

## What was wrong

The team car (team layout, `layout_genes.x_default()`, stage-2 numbers zero) had two masses:

| roll-up | where | what it sums | team car |
|---|---|---|---|
| referee T3.6, "as manufactured" | `car/build/legality.py:1290-1296` reads `mass_state["manufactured_mass_g"]` = `run_car.manufactured_mass_kg` (`car/build/run_car.py:107-115`) of the `--no-joints` branch (`run_car.py:83-91`) | body cells (`unified_phi.compute_mass_com`) + `parts_mass_g` (nose, front wing, rear wing, guides) + 2 x bare wheel + halo | **45.29 g** (fails by -2.71) |
| race model | `search/layout/layout_car.py:160-190` | body mesh x 0.163 + `rear_wing_mass_kg` (= the same printed parts, `car/parts/assembly.py:140-145`) + `wheels_front/rear_mass_kg` (= 2 x (wheel + support net), `assembly.py:136-138`, support net from `support_engine.build`, `assembly.py:252`) + halo; lead up to 48.2 g; cartridge 21.15 g | **53.79 g** (lead 0) |

`run_car.mass_state` (the "model's roll-up" in report.md) agrees with the race model (53.79 g).

Itemised from the team car's own files (both halves; ABS 1.04 for supports/hubcaps/wheels as the code sets
them; printed nose/wings/guides at the code's PLA 1.24 before P3):

| item | from | g (before) | g (P3: ABS 1.04) |
|---|---|---|---|
| foam | body_half.stl 160.4 cm3 x 0.163 | 26.14 | 26.14 |
| nose + front-wing pylons (solid) | nose.stl 5.397 cm3 | 6.69 | 5.61 |
| front wing | fwing.stl 1.437 cm3 | 1.78 | 1.50 |
| rear wing + struts | rwing.stl 3.557 cm3 | 4.42 | 3.70 |
| tether guides | tethers.stl 0.144 cm3 | 0.18 | 0.15 |
| halo + helmet | 6.705 cm3 LW-PLA 0.8 x 20 % | 1.07 | 1.07 |
| wheels (team_stl, 4) | 2 x 1.178 + 2 x 1.322 | 5.00 | 5.00 |
| **wheel supports + hubcaps (4)** | part + hubcap x 1.04 - slot foam x 0.163 | **8.50** | 8.50 |
| **car (T3.6 mass)** | | **53.79** | **51.68** |
| cartridge shell (separate) | | 21.15 | 21.15 |

Supports per side: front 1948 mm3 part (saddle 936, arm 410, disc + stub 602) + hubcap 500 mm3 = 2.546 g, less
936 mm3 of foam (0.153 g) = 2.393 g; rear 1389 (646 / 100 / 642) + 500 = 1.964 g - 0.105 = 1.859 g.

**Bug 1 (P1)**: the `--no-joints` manufactured mass (every search build) leaves out the supports and hubcaps,
8.50 g. The 53.79 g race-model number is the right one; T3.6's 45.29 g was wrong. With joints (finalists) the
supports are in `joints.make_all`'s manufactured files, so only the no-joints branch was affected. It also mattered
for any *sized* `--no-joints` build: the body would have been sized 8.5 g too heavy.

**Bug 2 (P3)**: `assembly.PART_MATERIAL = "PLA"` (1.24 g/cm3) weighs the search car's nose (+ pylons), wings,
struts and guides, while the team prints ABS (team 2026-11-14: `joints.DENSITY_G_CM3`, `nose.NoseCone.material`,
`beam_support.ABS_G_CM3`, and carhybrid/grad.py's nose density) and the joints path weighs the same parts at 1.04.
Both roll-ups were 2.11 g heavy (19 % on every printed gram), and a nose change weighed 19 % more than the
gradient predicts.

**Not a bug, but misleading (P2)**: "MASS NOT CONVERGED: car 45.289 g, -2.911 g" in run_car.log. Search builds
pass `--skin-mm 0.0` (layout_params `skin_mm` fixed 0.0); `param_body.build` then returns `make(skin)` at once
(`car/body/param_body.py:496-497`) and `run_car.build_car` makes one pass with no target (`run_car.py:170-171`).
No body is resized to a mass target, so the free-form body's mass changes do reach the race model. The line was
build_car comparing a fixed-skin car with the 48.2 g target anyway.

## The fixes

| patch | file | change |
|---|---|---|
| P1 | car/build/run_car.py (`parts_and_joints`, `--no-joints`) | `manufactured_g` += 2 x `info["support_mass_g_each"]` (part + hubcap - slot foam, per side: the numbers already in the race model's wheel masses); `joints["supports_g"]` records it |
| P2 | car/build/run_car.py (`build_car`) | fixed skin: `info["sized"] = False`, `mass_converged = None`, logs "fixed skin, not sized"; sized builds unchanged |
| P3 | car/parts/assembly.py | `PART_MATERIAL = "ABS"` (1.04), as `joints.DENSITY_G_CM3` |

## Before / after (real builds, 1 mm grid, carhybrid/make_car.py)

| car | T3.6 (referee) | layout_car car_g | race mass | lead | referee verdict |
|---|---|---|---|---|---|
| team, unpatched (`car_team`) | 45.29 g, FAIL -2.71 | 53.79 g | 74.94 g | 0 | 4 failed (T3.6 + 3 grid-bore checks) |
| team, patched | 51.68 g, pass +3.68 | 51.68 g (diff 0.005 g: cells vs mesh) | 72.83 g | 0 | 3 failed, all the 1 mm grid's bore (waived by layout_car); **legal** |
| team, patched, side numbers u[30:45] = -0.3 (body 1.8 mm narrower) | 48.38 g, pass +0.38 | 48.38 g | **69.53 g** | 0 | same 3 waived; **legal** |

Shrinking the body by stage-2 numbers: foam -1.98 g, nose -1.01 g, rear-wing struts -0.13 g, supports -0.18 g
(front arm longer 10.1 -> 11.9 mm so sized thicker 3.03 -> 3.36 mm, but its saddle 140 mm3 smaller; rear -0.15 g),
race mass -3.30 g (about -25 ms at 7.5 ms/g). The body was built once at its skin (no resizing); nothing undid it.
That car is 0.18 g above the 48.2 g floor: the next gram shrunk is lead.

u[2] = +0.3 (unpatched, `car_u2`): +0.384 g = foam +0.113, nose +0.274 (PLA), rear wing +0.003, supports -0.006.
The free-form field puts +0.35 g of foam there but the front-axle rule masks and milling keep only 0.11 g; the
nose at 1.24 g/cm3 was 19 % over a 1.04 prediction (P3 removes that). Note carhybrid/grad.py's mass term
samples only `body_half.stl`, which ends at Ref A, so it does not see the nose's mass change at all.

## The supports: heavy by design, not by load

Load case (beam_support.py:78-84): 3.0 N radial at the bearing, 1.2 N axial at the tyre's radius, ABS 40 MPa / SF
2.0 = 20 MPa allowed, E 2000 MPa, sag <= 0.15 mm at the bearing, 1.05 mm print minimum. On the team car the front
arm (10.1 mm from the body's side to the disc) is sized 2.00 -> 3.03 mm by the sag (0.378 -> 0.136 mm; stresses arm
4.8, saddle 7.4, disc 10.6, stub 7.8 MPa); the rear is not sized (sag 0.075 mm, stresses <= 12.5 MPa). The sizing
costs 0.24 g of the 8.50 g. The rest is set by fixed design choices that no gene reaches: hubcaps 2.08 g (12 mm
radius x 1.05 mm, `SupportDesign.hubcap`, not a gene), discs + stubs 2.59 g (`sup_disc_r` fixed 12 mm, disc at the
print minimum), saddles 3.29 g of ABS for 0.52 g of foam (`sup_saddle_h` fixed 3 mm, and the whole column printed
where less than 3 mm of foam would be left over the slot), arms 0.82 g at the design thickness. The team's v2 CAD
supports are 1.343 / 1.485 cm3 per side = 5.88 g, with no hubcaps; ours without hubcaps are 6.42 g net (6.94 g of
plastic). Not changed here (no load or design changed); making `hubcap`, `sup_saddle_h` and `sup_disc_r` stage-2
genes would let the search trade them (their load gates are already measured on every car).

Responsiveness on the team body (support_engine.build, every gate still passing): the stage-2 gene `sup_chord`
12 -> 8 mm: 8.50 -> 7.33 g; `sup_fillet` 2.5 -> 1: 8.39 g; `sup_t` 2.0 -> 1.2: 8.45 g (rear only: the front is
held at its sized 3.03 mm, so below that the gene is flat); no hubcaps: 6.46 g; `sup_saddle_h` 3 -> 1.5 (a
fixed parameter): 7.40 g. Wings respond to their numbers (rear wing 3.70 g; `rw_t` 2.3 -> 4.0: 4.82 g;
`rwsup_thick` 3.6 -> 2.0: 2.92 g; front wing 1.50 g, `fw_chord` 15.3 -> 20: 2.50 g; the team car's front wing has
no separate pylons, so `fwsup_*` do not change its mass there).

## Not fixed

- `mass_report.json` of a `--no-joints` build sums only manufactured files (there are none), wheels and halo:
  6.07 g "mass_g_without_cartridge". Nothing reads its masses (race_objective reads only the cartridge and tether
  positions); it needs the parts' COMs to fix properly.
- carhybrid/grad.py (this repo) misses the nose's mass in dm/du (see above).

## Tests (`evo_patches/tests/test_mass.py`)

Cheap (CI): patches applied; no-joints manufactured mass = race-model roll-up incl. supports; a fixed skin builds
once with no target; `param_body.build` returns at the fixed skin before any target; printed density = joints
density; layout_car.mass falls gram for gram with foam above the floor; stage-2 side numbers shrink the drawn body
and nose. With CAR_DIR: (a) T3.6 = car_g within 0.05 g and passes; supports in both roll-ups and equal to their
files; printed parts = file volume x 1.04. With CAR_DIR + CAR_DIR_SHRUNK: (b) less foam, lower race mass, no lead,
one pass at the fixed skin.
