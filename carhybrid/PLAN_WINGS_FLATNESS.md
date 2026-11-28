# Decisions and plans of 2026-11-28 (evening): our method, wings and supports, flatness, the chamber wall, hubcaps

## 0. Decision: no CMA in this repository

The owner, 2026-11-28: CMA (plain or gradient-guided) is the other repository's method. This branch is our own
hybrid, the trust model: adjoint-driven steps judged by measured CFD, evolution where the adjoint cannot reach, and
the Bayesian layout search. The guided-CMA recommendation (STATUS.md, section 3) is withdrawn.

**Goal restated (owner, later the same evening):** the goal is the best car, not beating CMA. A win over CMA is not
evidence of a good car if CMA is a weak method. So we work with our own method and judge it by the car: race time on
the audit protocol, against the team car and the best car found so far. A comparison with the other repository is
at most a side note.

## 1. Rear-wing supports: what "searched" means, and taking theirs

The builder already has organic strut generators, so this needs no new geometry code:
- **The rear strut:** `wing_strut.py`, 13 numbers per strut (aim, S-bend mid point, swell, section, flare, where on
  the span it starts).
- **The front:** the same, through `front_wing_engine.py`.
- **The other agent's stage 2** searches the rear ones because they are marked `shape` in
  `search/layout/layout_params.json`.
- **Ours never did.** `carhybrid` only moves the body numbers, so all 32 sweep cars carry the default struts. The
  front ones are `fixed` everywhere.

Taking them: our pipeline reads those numbers through `layout_genes.decode`. Making them searchable needs:
- a gene vector for them in `make_car.py` / `bo.py` (no builder change; they are decoded already);
- the front ones moved from `fixed` to searchable in our copy of the parameters, through an evo_patches entry;
- a gradient for them: DAFoam already gives dD/dX on the strut surfaces. dD/dθ is the sum of dD/dX · dX/dθ, with
  dX/dθ by finite differences of the strut generator, which is cheap (geometry only, no CFD).

So the struts become part of our hybrid (adjoint plus measured steps), not random-only as in theirs.

## 2. Wings: a model that is legal by the rule's own definition, with nothing hard-coded

**Problem (owner).** A strut that touches or passes through a wing can split it into elements. For the rear wing,
the result may then not count as one legal rear wing.

**The rule text, now in the repository** (`docs/rules/`, the owner's PDF; readings in `docs/rules/README.md`).
T9.5.1 requires a single, unbroken span of at least 50 mm, measured along the shortest of the leading edge, trailing
edge, top and bottom surface. Every figure shows the wing held at its ENDS. **Strict reading, adopted:** any part
joined to or passing through the rear wing breaks its span there.

**What the code does today:**
- The referee (`legality.py` `_wing_split`, `_wing_measures`) never counts the rear wing's own supports as breaking
  the span. Only the body and other parts do.
- That reading ("its own pylon under it does not break it") is **not supported by T9.5.1's text or figures**. Cars
  whose rear struts fell back to under-wing pylons (the engine's root fallback 0.8 / 0.6 / 0.4 of the half span:
  sw_L08's render shows struts joining the wing well inboard) may be illegal under the strict reading. The referee
  must be fixed first: own supports count as span breaks for the rear wing, and the measure is the shortest of
  LE, TE, top and bottom. Then sweep5 is re-judged.

**Plan:**
1. **Rule definitions first.** Take T8.5–T8.8 and T9.4–T9.7 verbatim from the rules PDF into `docs/RULES_TEXT.md`.
   For each, record the measured quantity: what is an element, what breaks a span, whether a pylon under the wing
   or a strut at its tip breaks it, what the vanes and overlap rules are. Write referee tests just inside and just
   outside each (the builder already has this pattern in `test_rules_edges.py`), including "strut through the
   wing", "pylon under the wing" and "strut at the tip".
2. **One wing model per rule set, legal by construction:**
   - each wing is 1–2 elements (vanes, T9.5.2) with chord, thickness, span, height and position;
   - its supports are joints of declared types: under-pylon, tip, side, into the nose;
   - each joint type is allowed only where the verified rule says it does not break the span, and the generator
     never lets a support pass through an element unless that type is legal;
   - legal-by-construction means each measured quantity is computed by the same function the referee uses, inside
     the generator, before the part is accepted.
3. **Nothing hard-coded.**
   - Element count, joint types and the support topology become discrete choices for the Bayesian layout search.
   - Chord, thickness, span, height, position and the strut numbers become continuous numbers for our hybrid
     (adjoint where DAFoam's surface sensitivity reaches them, measured steps elsewhere).
   - **Owner, 2026-11-28 evening:** the rear wing may take any shape the regulations allow ("have it bend"), but
     with **no angle of attack**. This replaces the earlier "straight, level, no camber, no AoA" choice for the
     rear wing. What that opens, each as searched numbers with nothing hard-coded:
     - **spanwise bend:** the wing's height may vary along the span (gull, arch, droop, dihedral), within T9.5.4's
       15.0 mm top-surface deviation over the minimum span;
     - **planform:** sweep and a chord that varies along the span, each x-z section's chord 15.0–25.0 mm (T9.5.2);
     - **thickness** varying along the span, 2.0–6.0 mm (T9.5.3);
     - **a second vane** (T9.5.2: up to two elements, overlapping).

     Every section's chord line stays level (AoA 0 everywhere, so no twist either). Open questions for the owner:
     (i) camber: is a cambered section with a level chord line allowed (zero AoA but not symmetric)?
     (ii) does the same freedom apply to the front wing (T8.6)? Until answered: rear symmetric sections; front
     unchanged.
4. **Tests before it is relied on:**
   - the referee's edge tests;
   - a sweep of random wing and support choices through builder and referee (target: every built car legal);
   - one CFD sanity check that the drag response to wing size has the expected sign.

## 3. Flatness: why, and the plan

Measured on the team layout (and on sw_L03, identical):
- the seed's sidepod top is drawn at 23.4–30.2 mm, median 26.9;
- 38 % of it is above 27.175 mm;
- the soft cap (`evo_patches/soft_cap.py`) squashes every height above 21.2 mm into the 6 mm below 27.175;
- so the sidepods come out at 23.3–26.6 mm: a 7 mm range pressed into 3.3 mm.

That is the flat look. Causes:
1. **The 27.175 mm sidepod cap is a TEAM rule, not a competition rule** (`docs/PARAMETERS.md`: "team: nothing above
   27.175 in the sidepods"). The competition rules only keep the halo visible: T4.4.2 front and side views, T4.4.3
   top. That limits the body round the halo, not the sidepods ahead of and behind it.
2. **Rules after the numbers.** The seed is fitted to a car whose pods rise above the cap, and the numbers are 0
   (no optimisation has shaped any car yet). So the squash band sets the shape.
3. sweep5's extra ~200 mm² of flat per car is underneath or internal: the official block's slot walls and the bore
   end, from the block clamp L2. It is mandated geometry (the official block's slot is flat by specification), not
   the sidepods.

**Plan:**
- **(a)** ~~Owner decision: keep or drop the team sidepod cap.~~ **Dropped (owner, 2026-11-28):**
  `evo_patches/sidepod_rules.py` S1. Only T4.4.2–3 limit the body near the halo; the sidepods' width limit is the
  official block's (32.5 mm each side, T1.14/T3.1.2), which the envelope already used.
- **(b)** Either way, no squashing: the limits become constraints on the numbers (`constraints.py`).
  - For the cap and halo-view limits this is easy: they bound the top surface only, a simple bound the projection
    handles without the full set.
  - The seed is projected once onto them (re-seeded), so the body starts below the limits instead of being squashed.
  - The soft cap is then switched off for projected numbers.
- **(c)** What is left flat is then either a rule (halo rim flush, pocket floor, block slot) or the CFD's choice.

## 4. The 3 mm chamber wall (T5.5)

- **Measured:** `T5.5_wall` passes on all 32 cars of sweep4 and of sweep5 (margin 0.0).
- **In the builder:** the safety zone is required solid that overrides every air mask (`unified_phi`: "the 3 mm wall
  round the chamber outranks every design-choice void").
- **In the support field** (`evo_patches/wing_field`): the bore + 3 mm is a keep-out, and no plastic or pocket
  enters it (0 mm³ in every prototype run).
- **To add:** the same keep-out for every plastic joint the wing model in section 2 makes. The `T5.5_wall` check
  stays a hard gate in every batch.

## 5. Hubcaps

**Fact:** the CFD wheel is a solid flat-faced cylinder (L01: volume 8,300 mm³ against 8,303 for a solid cylinder;
the outer face is fully closed). The builder's hubcap is a flat 12 mm disc pressed on the stub, flush with the wheel's
outer face (`beam_support.py`). It is not a CFD surface. So every drag number so far assumes a flat closed wheel face.

**Two defensible positions for the paper:**
- **A. Fixed flat, justified by a controlled test, which has not been run.** A pre-registered comparison on the team
  car, audit protocol, frames 0–3:
  - flat flush (today);
  - flat recessed 1 mm;
  - domed (2–3 mm);
  - (if the rules allow) an open rim.

  About 16 runs, free on GitHub. If flat is best or within noise, the paper can state: "the hubcap was fixed as a
  flat flush cover: lowest mass (the printer's minimum wall), and no measurable drag penalty against the
  alternatives tested; it was excluded from the search to keep the design space within the CFD budget." Without the
  test, this statement is an assumption, not a result.
- **B. Optimise it with the car.** It adds 2–3 numbers (dome height, recess, cover radius). It needs the CFD wheel
  built with the hubcap's shape (adapt.py), and a rule check (wheel rules T7.x). The wheels are 65 % of the car's
  drag (P0: 151 of 233 mN, half car), so the wheel faces are not a negligible area.

**Recommendation: run A's test first.** If any variant beats flat by ≥ 3 SE, take B. Otherwise fix it flat with the
test as the paper's justification.
