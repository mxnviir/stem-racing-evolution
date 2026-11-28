# Status of the work on branch claude/wizardly-brahmagupta-1edalc (2026-11-28)

**The method in one page: [METHOD.md](METHOD.md).** One page per piece is linked; this page says what exists, what it showed, and what is still open.
Repositories:
- `mxnviir/stem-racing-evolution`, this branch: everything below.
- `mxnviir/stem-racing-adjoint`: the CFD and DAFoam code, pinned at `9ab9092`, patched at run time; nothing is
  committed there.

The builder is the evolution repo's branch `claude/loving-keller-q72mvv` at `e9c896b`. It is read only and patched at
run time by `evo_patches/`.

## 1. Audits and methods (paper material)

- `docs/METHODS_TRIED.md`: every method tried across both repositories, with evidence, verdicts and links.
- `audit/`: drag audits b1–b3. The original sidepod hybrid keeps −3.05 ± 0.80 mN half car against a team pod
  legalised the same way (`RESULTS_b3.md`; b2's first reading was corrected).
- **CFD spent on the sidepod before this work:**
  - SimScale CMA-ES: 133 runs (about 2,560 core-hours), no gain (best car 2.6 ms slower than the team car).
  - The adjoint hybrid: about 100 solves for the −3 mN above.

## 2. The builder: legal cars by construction (`evo_patches/`)

| Patch | What it does | Note |
|---|---|---|
| P1–P3 | mass bookkeeping | `evo_patches/README.md` |
| soft_cap v2.1 | 27.175 mm cap without flat decks; fades before the halo's rear face | `evo_patches/soft_cap.py` |
| sidepod_rules S1 | the team's 27.175 mm sidepod cap **dropped** (owner, 2026-11-28); only the halo's T4.4.2–3 limit stays; sidepod width = the block's 32.5 mm | `evo_patches/sidepod_rules.py` |
| nose_material | foam or printed nose decided per car; block front fitted to the referee | `evo_patches/NOSE.md` |
| legality_fixes L1 | rear wing and supports end before the body's back (T5.6) | `evo_patches/legality_fixes.py` docstring |
| L2 | foam body clamped to the official block as the referee places it (slot slivers, block front, bore drawn small) | same |
| L3 (+v2) | wheel supports re-tried from design variants when their own gates fail | same |
| L4 | front-wing tip moved out when a strut's flare leaves a segment under 25 mm (T8.6.1) | same |
| L5 | rear wing built further forward when no support can land ahead of the body's back | same |
| L6 | nose-joint check fails a sliver face instead of crashing the build | same |
| rear_wing_shape W1-W6, L7 | rear wing as the regulations allow with zero AoA (bend, sweep, taper, thickness, camber); T9.5.1 strict in engine and referee (own supports break the span) | `evo_patches/rear_wing_shape.py`, `docs/rules/README.md` |

Legality of the 32 test layouts:

| Sweep | Legal of 32 | Note |
|---|---|---|
| sweep2 | 17 | |
| sweep3 | 11 | regression, traced (`results/sweep3/RESULTS.md`) |
| sweep4 | 16 | |
| sweep5 | 28 | |
| sweep5 + sweep5r | **30** | |

The two left:
- **L21:** foam over the rear support's slot, 0.24 mm short.
- **S02:** halo face joined, a failure already present before any patch.

Galleries: `results/sweep3/GALLERY.png`, `results/sweep5/GALLERY.png`.

## 3. Whole-car hybrid (`carhybrid/`; design `DESIGN.md`)

| Piece | State | Doc |
|---|---|---|
| Car build + OpenFOAM adapter (`make_car.py`, `adapt.py`), shape audit (`shape_audit.py`) | working | `README.md`, `DESIGN.md` |
| P0: team car, coarse DAFoam gradient | done | `P0.md` |
| P1: does the gradient predict on the production mesh? | **both marks met**: race-time step −11.3 ms; drag step beats its opposite by 29.6 ± 2.8 mN whole car and both random steps; non-linear (slope 1.4) | `P1.md` |
| Fine-mesh adjoint (DAFoam on the 2.1 M scoring mesh; lift; gates F0–F3) | built; plumbing-tested on GitHub at 129 k cells; **not yet run at 2.1 M** | `FINE_ADJOINT.md` |
| AWS deployment (AWS CLI, spot, resume, results via S3) | built and tested locally (dry run, fake aws); **not yet launched**: needs your AWS setup | `aws/README.md` (instance choice and prices, section 5) |
| Gradient map (`grad.py`: drag, lift, mass → race time; 45 or 720 numbers; frame averaging) | working; agrees with the other agent's map (cosine 1.0000). **Known defect, not yet fixed:** counts mask-held surface as moving (their exact test: top numbers 1.30–1.46×) | `grad.py` docstring |
| Trust-region loop (`loop.py`, P2F) | built and tested (19 unit tests); not run | `FINE_ADJOINT.md` (loop, decisions) |
| P3: Bayesian layout search on OpenFOAM (`bo.py`) | **running**: batch bo_1 (16 Sobol cars) on GitHub; the collect step proposes and dispatches each next batch automatically | `BO.md` |
| Rules as constraints on the numbers (`constraints.py`) | **in progress** (see section 5) | `CONSTRAINTS.md` |

**Decision (owner, 2026-11-28 evening): no CMA in this repository.** Guided CMA is the other repository's method; ours is the hybrid trust model (adjoint steps judged by measured CFD, evolution where the adjoint cannot reach, the Bayesian layout search). **The goal is the best car** (race time on the audit protocol), not beating CMA. The earlier guided-CMA recommendation is withdrawn. Plans for wings, supports, flatness, the chamber wall and hubcaps: [PLAN_WINGS_FLATNESS.md](PLAN_WINGS_FLATNESS.md).

It waits on the constraints (section 5), the `grad.py` fix, and DAFoam converging at non-team shapes on AWS.

## 4. Seeded sidepod campaign (`seedhybrid/`)

- Protocol, marks and change log: `seedhybrid/PROTOCOL.md`.
- The setup's DAFoam gradients failed at the 2,000-iteration GMRES limit. With 6,000 iterations all three converge:
  window 1,290, d7 1,615, team 1,963.
- **Generation 1** is dispatched (run 46156806568) and queued behind the layout search for GitHub's 20 slots:
  - all three seeds adjoint-guided;
  - no-adjoint arms team and d7 (window's skipped by the pre-registered rule).

## 5. In progress

- **Rules as constraints** (`CONSTRAINTS.md`):
  - The seed body itself breaks the rules: the masks fix 263 cells of clipping air filled and 5,919 cells of
    required solid left empty. That is the origin of the flats.
  - No set of numbers meets every rule: the smallest worst violation is 12.8 mm on 45 numbers and 11.3 mm on 720.
    Some rule features are sharper than the height-map surfaces can draw.
  - The split into reachable rules (hard constraints, the search never enters them) and sharp features (the mask
    stays, the gradient holds the surface): 24,658 rows reachable together, 1,823 sharp (7 %, mostly the wheel
    zones and the halo and cargo edges).
  - Open: the projection onto the reachable set (ADMM does not converge; next, a proper QP solver), then the builder
    check of projected numbers.
- **Wing supports as a CFD-designed field** (`evo_patches/wing_field/WING_FIELD.md`):
  - Done: the support field, evolution's growth move, ball-millable pockets, 3.175 mm rolling-ball fillets, exact
    foam/plastic complements, the 3 mm chamber keep, the strength limit.
  - Prototype on L01's rear: legs grown onto the foam above the chamber wall and onto the body's side, 0.6–4.2 MPa,
    no plastic in the chamber keep.
  - Open: 2–12 mm³ of concave corner still under the ball at the 0.5 mm grid (to trace); wiring into the builder,
    `grad.py` and the loop; the front zone.

## 5b. Rear wing (built 2026-11-28 evening)

- Rules: `docs/rules/` (the owner's 2026-2027 regulations); T9.5.1 read strictly.
- Built: `evo_patches/rear_wing_shape.py`; `make_car.py --wing` takes the shape and strut numbers.
- Tested locally on L01's body: four shapes plus pylon cases. All gates met; the pylon that would break the span is
  refused; the strict referee agrees.
- On GitHub, batch `wing1` (run 46174764457; 12 random shapes plus 4 pylon cases): **10 of 16 legal**.
  - 2 broke T9.4.3/T3.5 (65 mm): a large rise (11–13 mm) on a wing set at the lowest height T9.6 allowed. This is the
    wing model inventing a shape against a rule. Fixed by **W7**: the rise is cut until the top fits, and the height is
    found again. Under local re-test.
  - 5 found no support (`rw_supported`): random strut numbers with no landing; includes the pylon at 0.7 of the half
    span (refused for the strict span as designed, and then no fallback landed). Not a rule break (the car is not
    built), but each costs a Bayesian-search slot. Open: trace why L5's move forward did not land them.
- Not yet:
  - re-judging sweep5's 32 cars under the strict span;
  - a second vane;
  - the wing and strut numbers in the Bayesian search and in our hybrid (gene vector, gradient through
    DAFoam's surface sensitivity);
  - the front wing (owner: unchanged for now).

## 6. Waiting on you

- **AWS:** run `carhybrid/aws/aws_setup.sh`, check the spot quota (128+ vCPU for m7a.32xlarge in us-east-2), then
  `launch.sh aws_smoke`. Recommended instance: m7a.32xlarge spot in Ohio; about $85–180 for the whole programme on
  spot (estimate).
- **Go-ahead** to: fix `grad.py` (held surface); build guided CMA; wire the support field into the builder.
