# The nose: foam where it can be, printed where it must be, flush with the body

`evo_patches/nose_material.py` patches the pinned evolution builder (mxnviir/stem-racing-evolution @ e9c896b) at run
time, like `soft_cap.py`:

    python evo_patches/nose_material.py <evolution checkout>        # idempotent; asserts every anchor
    EVO_NOSE_SPLIT=0 ...                                              # turns the split off at run time (builder as before)

It writes one new module, `car/body/nose_split.py` (the decision and the geometry), and patches 9 files at 28
anchors (listed at the end). Applying it before or after `evo_patches/apply.py` gives byte-identical files.
Tests: `evo_patches/tests/test_nose_material.py` (`EVO_CHECKOUT=<checkout>`, optional `CAR_DIR=<out>/evo`).

## What was wrong

1. **The whole nose was printed.** `freeform.build_nose` made one printed part from the free-form field ahead of
   Ref A; `assembly.build` weighed it solid at `PART_MATERIAL` once the front wing's pylons were united with it
   (the pinned builder: PLA 1.24 g/cm3; `apply.py` P3 makes it ABS 1.04). Team car: 6.69 g, half of its 13.07 g of
   print. Foam is 0.163 g/cm3: the same volume in foam weighs 6.4-7.6 times less.
2. **The lip at Ref A** (team car, centre plane: print 24.3 mm, body 22.2 mm just behind). Not the skin (0 on every
   search car) and not the mesh smoothing: `machining.make_machinable` relieves "the zones' own sharp corners, for
   every body" (`corner = Hm & ~reachable_air(H, r, any_side=True)`, dilated by the 3.175 mm ball). T8.5.1's
   ceilings over the nose (the envelope's 25 mm nose box, `unified_phi.build_rule_masks`' 24.5 mm cap running
   1 mm + a cell aft of Ref A) end in a square corner that no ball fits, so a ball-sized dip was cut into the top
   of the body at the corner, whatever the body: 21.7 mm flat on the centre plane round Ref A. The printed nose,
   drawn from the same field without that dip, stood 2.0 mm proud (measured: top step 2.04 mm at Ref A + 1, side
   0.25 mm at the widest height). The cap's own end behind Ref A then made a second, 0.8 mm step 2-3 mm aft.

## The rules relied on

Regulation text is not in the checkout (`docs/RULES.md` cites `regulations_text_2026-27.txt`, absent); the rule
map in `docs/RULES.md` and the referee's code were used.

| rule | what it says (docs/RULES.md) | how the split keeps it |
|---|---|---|
| T1.3 / T4.1 | body = one piece of block material behind Ref A | behind Ref A nothing changes; the foam nose is the same field and mesh (`T4.1_single_body`: 1 body) |
| T1.7 | nose cone only forward of Ref A | no material rule in the map. The team: "Nose cone: may mix foam and ABS. Ahead of Ref A plastic is allowed, not required" (stem-racing-adjoint `research/sidepod_opt/PARTS_PLAN.md`, answers of 2026-11-19). `unified_phi.py`'s older note (nose printed, owner 2026-07-20) predates it |
| T1.14, App. ii | made from the official block (223 x 65 x 50, 6 x 6 slot) | foam only from 1.0 mm behind the block's front (rear face - 223); the slot and underside are the field's hard air already; the referee's `block_body_inside_official_block` measures the whole milled body, foam nose included |
| T8.2 | nose overhang <= 40 mm ahead of Ref A | unchanged shape; referee: front of every non-wheel part incl. the body |
| T8.5.1 | nose and wing support ahead of Ref A <= 25 high, <= 15 each side | the field's cap (24.5 / 14.5) on the foam, the same cap on the print; referee: body vertices ahead of Ref A - 0.5 and the printed nose (`check()`, already there) |
| T8.7 / T8.8 | 5 mm clear air from the front wing; nothing ahead of it | wing designed against the whole drawn nose as before; the referee's "within 5 mm (sideways) of the nose" now takes the nose's widest point over print AND foam (`_foam_nose_half_width`) |
| T3.7 | 1.5 mm clearance | print clipped at 1.7 mm (as before); foam is above the block's underside (2.0 mm) |
| machining | 3-axis, 6.35 mm ball, 4 setups, 3.0 mm foam wall (`machining.TOOL_R_MM`, `rules_constants.MACHINED_MIN_WALL_MM`) | the foam nose is milled with the body (same masks, `make_machinable`, skin) |
| `nose_has_tenon` | the nose is held by a glued tenon, not only a surface bond | the print carries `joints.nose_part`'s tenon (8 mm aft, 8 mm wide, underside to the section's centre) into a pocket milled from the underside; as-drawn referee check `nose_has_tenon` (>= 4 mm in), and the finalists' joints path reports it as before |

No rule forbids foam ahead of Ref A: the referee already measured body material there (T8.2, T8.5.1).

## The assignment (car/body/nose_split.py `plan`)

The nose stays ONE drawn shape: the free-form field ahead of Ref A. The front wing is designed against it exactly as
before (`front_wing_engine.build` on `freeform.build_nose`: the plan calls the same function first, 7-11 s, to know
where the wing holds on). Only the material inside is split, by one plane `x_cut`, printed ahead, foam behind:

    x_cut = max( block front + 1.0 mm                                   (b) outside the block
                 last x where the wing touches the nose or its struts /
                   underside pylons end, + 1.2 mm ABS collar             (a) the wing's attachment
                 last x of a section with > 1 mm2 under the 3 mm wall    (c) thin (opening by a 1.5 mm ball)
                 Ref A - 40 + 0.2 )                                     T8.2
    then the first section from there whose glued face carries the load case (d), never behind Ref A.

(d) is the strut generator's case (`wing_strut.LOAD_N`: 1 N drag, 0.5 N side, 2 N down per half wing, both halves,
side loads the same way) at the wing's centre, on the glued face (the block-limited section at `x_cut`): axial +
both bendings + shear + torsion against 0.5 MPa on the foam (`support_engine.FOAM_BEARING_MPA`, an assumption) and
20 MPa in the ABS (`wing_strut.ALLOW_MPA`). It never bound on the cars built (0.12-0.19 MPa); it moves `x_cut` aft
only on a small section far ahead of the wing. `x_cut` at Ref A means "all printed" (long layouts, a wing set close
to Ref A on a long body).

Why one plane: the printed piece must hold the wing and be glued to foam that is one piece with the body, and the
foam must be millable. Everything ahead of the wing's attachment is then printed, so the search's own genes decide
how much foam a car gets: a shorter body (the block reaches further forward) and a wing set further forward
(`fw_te_gap`) both move `x_cut` forward. The layout search trades it through race time with no new gene.

Geometry, so the junction is flush:
- the body's field runs `LAP_MM` = 5 mm ahead of `x_cut` (`param_body.build(nose_split=)`), labelled main body
  (foam density, not `DENSITY_NOSE_KGM3` 1000); `run_car.export_body` cuts the mesh exactly at `x_cut` and takes the
  tenon's pocket out (`nose_split.cut_body`), so the glued face is sharp;
- the T8.5.1 cap on the field stops at Ref A (where the rule stops), and the milling does not relieve the cap's
  corner (`machining.make_machinable`, split cars only): the dip is gone; a body reaching a ceiling would be cut flat;
- the print is the same field ahead of `x_cut` + 0.3 mm (into the foam: one closed solid in the CFD union), on the
  0.25 mm grid, with the tenon.

Mass and COM: foam nose = body (`unified_phi.compute_mass_com` weights the field's cells ahead of `x_cut` and in the
pocket 0; `layout_car.mass` weighs the cut mesh); print = `parts_mass_g["nose"]` (`nose_split.printed_mass_com`: the
piece ahead of `x_cut` and the tenon, not the 0.3 mm overlap; solid with pylons, else `nose_mass_com`'s shell rule, as
before), at the printed density, in `fixed_hardware_kwargs` and so in the race model's "printed". The referee's mass
(T3.6, from the field) and the race model's (from the mesh) agree within 0.03 g on both cars.
`layout_car.mass()["nose_split"]` and `assembly.json` `info.nose_split` carry the split (foam_g, printed_g, COMs).

## Before / after (team layout and a short layout; soft_cap applied in all; printed density PLA 1.24 as pinned)

| | team before | team after | short before | short after |
|---|---|---|---|---|
| genes | x_default | same | `[0,0.18,0.5,0.1304,0.5,1,0.6562,0.5172,0,0,0.175,0]` | same |
| layout | W 138, rear 63.8 | | W 120, rear 45.2, wing 20 mm ahead of Ref A | |
| block ahead of Ref A | 5.2 mm | | 41.8 mm | |
| x_cut (decided by) | Ref A | 37.00 (attachment: struts to 35.68) | Ref A | 23.75 (attachment: struts to 22.54) |
| foam nose | 0 | 2.4 mm, 0.107 g | 0 | 15.65 mm, 0.629 g |
| printed nose (+ pylons, tenon) | 6.69 g | 5.75 g | 9.98 g | 4.58 g |
| car without cartridge (race model) | 53.79 g* | 52.85 g | 50.65 g | 45.91 g + 2.29 g lead to 48.2 |
| race mass | 74.94 g* | 74.00 g | 71.80 g | 69.35 g (-2.45 g, ~-18 ms) |
| junction step, top / side / bottom | 2.04 / 0.25 / 0.03 mm | 0.03 / 0.12 / 0.02 mm | 0.52 / 0.27 / 0.10 mm | 0.003 / 0.18 / 0.50 mm |
| glued face, stress / 0.5 MPa | - | 259 mm2, 0.17 MPa | - | 191 mm2, 0.29 MPa |
| legal (waived: T3.6, T5.1, grid bore) | yes | yes | yes | yes |
| CFD one closed solid | yes | yes | yes | yes |

\* team "before": `$S/car_team` (pinned builder without soft_cap; its soft_cap rebuild is in TEAM_BEFORE below).
In ABS (apply.py P3) every printed figure is x 0.839: team print 5.61 -> 4.82 g, short 8.37 -> 3.84 g.
Steps: the largest jump of the union's outer surface across the seam, sampled every 0.05 mm (centre plane top and
bottom, side at the nose's widest height). The short car's 0.5 mm bottom step is the block's 6 x 6 mm slot: its nose
dips 0.5 mm into the slot's band, the foam cannot (the body behind Ref A has the same groove).
Shape audit (`carhybrid/shape_audit.py`): rectangularity unchanged (short 0.756 / 0.873 both); new flat patches are
the glued face at `x_cut` (96-130 mm2 a half, covered by the print) and the tenon pocket's walls (internal).
Renders: `evo_patches/renders/{team,short}_{before,after}.png` (purple: print; grey: foam).

## Limits

- The split is one plane: ahead of the wing's attachment everything is printed even inside the block (a foam tip
  would be a second foam piece). A printed sleeve or saddle round the attachment with foam core through it would save
  more on cars whose wing sits far back; not done.
- Only fixed-skin builds are split (`--skin-mm`, every search build); a mass-sized build keeps the all-printed nose
  (its skin moves the foam; the print follows the plan's skin).
- 0.5 MPa on the foam side of a glued face is the repo's assumption (`FOAM_BEARING_MPA`), not a test.
- The finalists' joints path (`joints.make_all`): the split nose is kept as printed with its tenon (no second pocket),
  the body milled from `x_cut`, `nose_has_tenon` reported; see the test result in the report (JOINTS below).
- The thin test is an opening by a 1.5 mm ball with a 1 mm2 per-section tolerance (the blended edges read 0.1-0.8 mm2).
- T8.7 in the referee now exempts wing points within 5 mm sideways of the foam nose too (the whole nose), as it did
  for the printed nose.
- The tenon's pocket is a box in the as-drawn model (as `joints.nose_part`); the mill rounds its plan corners.

## Patched anchors

`car/body/nose_split.py` (new); `param_body.py` (build(nose_split=), field front, labels, cap switch); `unified_phi.py`
(T851_CAP_AFT_MM, compute_mass_com weights); `exact_zones.py` (t851_cap); `machining.py` (no cap-corner relief);
`run_car.py` (_nose_split plan, build_body, export_body cut, part4_kwargs); `assembly.py` (printed piece, wing on the
drawn nose, masses, info); `legality.py` (_foam_nose_half_width for T8.7/T8.8, nose_joint_checks); `joints.py`
(make_all); `layout_car.py` (mass breakdown).
