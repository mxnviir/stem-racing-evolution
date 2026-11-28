# Wing supports as a designed field, and exact plastic-foam joints

Owner, 2026-11-28:
- the front and rear wing supports design themselves from the actual CFD (evolution and adjoint), not from a proxy;
- the rear wing may attach anywhere on the foam body where legal: rear pods, the area above the cartridge chamber;
- 3 mm round the chamber is always foam, and only beyond that may there be air or plastic;
- pods may be milled foam merging into the body, and a rear-wing support may land on a pod;
- wherever plastic meets foam, every fillet is at least 3.175 mm in radius, and the two parts are exact
  positives and negatives of each other.

Code: `evo_patches/wing_field/appendage.py`. The builder patch (`evo_patches/wing_field.py`) writes it into the
pinned builder as `car/parts/appendage.py`.

## 1. What decides the shape: CFD only

A wing's supports are the solid where a field is negative:

    P(p) = seed(p) - AMP_MM * sum_j u_j * exp(-|p - c_j|^2 / (2 sigma^2))

- `seed` is the signed distance of a starting support: today's loft or strut engine output, or nothing.
- `c_j` are control points on a 4 mm lattice over the wing's legal support zone.
- `u_j` in [-1, 1] are the numbers the search moves.

A positive number grows material round its point; a negative one removes it. So any topology the numbers reach is
possible, none of them drawn by hand: one strut or two, a single fairing, a bridge, two legs that merge, a support
with a hole, a pod-like body under the wing.

The numbers are moved by the same machinery as the body:
- **Evolution:** random and recombined moves, scored by CFD. These are needed for topology. Material that appears
  where there is no surface yet has no adjoint sensitivity to grow from.
- **The adjoint:** DAFoam gives dD/dX and dL/dX on every surface of the car, wings and supports included. The
  support's surface moves along its normal at V_j = AMP_MM * g_j / |grad P| per unit u_j. So dD/du_j is the sum
  of (dD/dX . n) V_j over the support's surface points, the same map `carhybrid/grad.py` makes for the body's 720
  numbers.

There is no drag proxy, no frontal-area term and no topology optimisation in the objective. Two things are hard
limits that make a design legal or not; neither is ever traded against drag:
- **Strength:** a linear-elastic voxel finite-element check of the printed solid (wing + supports, ABS, E 2 GPa,
  limit 20 MPa as `wing_strut` / `support_engine`), glued where it meets foam or the nose, under `wing_strut`'s
  load case at the wing.
- **The print wall:** material thinner than 1.2 mm is opened away.

## 2. Where a support may be

Rear (T9.4: the wing and its support between Ref B and Ref B + 40 mm):
- x from Ref B to Ref B + 40 mm, but ahead of the body's back face (behind it, T5.6 wants the protruding cartridge
  visible all round);
- 2.5 mm or more off the centre line (team rule: nothing on a rear wing's centre line);
- inside T3.4's width;
- under T9.4.3's 65 mm;
- clear of the rear wheel (1 mm) and of T7.9's zones;
- never in the bore or the 3 mm wall round it (below).

T9.6 (5 mm of clear air round the wing except at a support's joint) and T9.7 are measured by the referee.

So a support can land anywhere on the foam in that band:
- on top of the body above the chamber's wall (the "above canister" area);
- on its sides;
- on a rear pod.

**Pods are foam.** They are the body's own free-form field (the stage-2 numbers) growing out beside or above the
chamber, milled with the body. A support landing on a pod is an ordinary plastic-foam joint (section 3).

Front (T8.5.1: a support ahead of Ref A, no higher than 25 mm, within 15 mm each side), between the wing and the
nose. Ahead of the foam nose's cut it lands on the printed nose (plastic to plastic: one print). Where it reaches
the foam behind the cut, it is a plastic-foam joint like any other.

## 3. The cartridge chamber: 3 mm of foam, always

`chamber_keep()` marks the bore (measured on the body: rim, depth, diameter, axis) and every point within 3.0 mm of
it, over the bore's depth. Nothing is plastic there and no pocket is milled there; the joint's ball may not touch
it. Outside that wall (above it, beside it) a cell may be foam, air or plastic.

So "above the canister" means above axis + bore radius + 3 mm: about z > 43 mm on a car whose chamber is at
z = 31 mm.

## 4. Plastic-foam joints: 3.175 mm radii, exact positive and negative (`joint()`, `mesh_joint()`)

R = 3.175 mm (the 6.35 mm ball end mill), built 0.15 mm rounder so the measured radius is never under 3.175.

1. **Pocket (the foam the plastic takes).**
   - The pocket is only what the ball can cut, 3-axis from the side the support comes from.
   - The ball's centres are the cells where a whole ball fits inside (the wanted pocket + the open air) without
     touching kept foam (the erosion by R).
   - A centre counts only if its whole vertical line out to the air is also a centre: the tool's shank, no
     undercuts.
   - The pocket is the union of those balls inside the body: the morphological opening.
   - Every concave corner the foam keeps is therefore at least R round. Wanted material the ball cannot reach stays
     foam, and the support's foot is that much shallower.
2. **The outside corner.** Where the plastic's outside meets the foam's outside, the corner is filled with plastic
   as a rolling-ball fillet of radius R. This is the closing of their union by the ball, kept only within R of both
   (the joint line).
3. **Exact complements.**
   - foam = body - pocket
   - plastic = (support + fillet - body) + pocket
   - Both come from Boolean operations on the same two meshes (the body as built, and the pocket's mesh). Every face
     the parts share is literally the same face, so the printed part is the exact positive of the milled cavity, to
     floating point. Nothing is modelled twice.
4. **Measured, not assumed.**
   - `concave_radius_check` gives the mm³ of air a ball of 3.175 mm cannot reach in a concave corner, round the
     joint, on the foam alone and on foam + plastic. It must be 0 to within the grid.
   - The chamber keep is checked for any plastic (0 mm³).
   - The FE check gives the stress.

The same joint is to replace today's joints:
- the rear-wing loft feet, butt-glued with no fillet: through `joint()`;
- the nose tenon, a square-edged 8 x 8 mm box the ball cannot cut the corners of: its pocket becomes the ball's
  opening and the tenon its exact positive;
- the wheel-support slots, which already have a ball-opened ceiling: to be measured with the same checks.

## 5. Making new topology: evolution's growth move

Scattered positive numbers almost never join the wing to the body. Two runs of 40 random growth points each gave no
support. So evolution gets a connected move, `grow_path`. It raises the numbers along a random smooth path from a
point on the wing to a random legal landing point on the foam (above the chamber wall, the body's side, a pod). A
new leg appears, and the CFD then decides whether it stays. Further moves and the adjoint reshape it: thicken, thin,
sweep, merge two legs, fair it.

## 6. Prototype results

Car v3_L01's rear (rear-wing zone, 0.5 mm grid, 1,388 numbers), `prototype_rear.png`:

| support | support mm3 | pocket mm3 | fillet mm3 | plastic in the chamber keep | stress MPa (limit 20) |
|---|---|---|---|---|---|
| today's loft strut as the seed, u = 0 | 2,174 | 0 | 77 | 0 | 4.2 |
| seed + random u | 2,887 | 1.6 | 72 | 0 | 4.2 |
| no seed: one leg grown onto the foam above the chamber wall | 1,031 | 6.6 | 154 | 0 | 1.5 |
| no seed: one leg grown onto the body's side | 1,726 | 94 | 199 | 0 | 1.4 |
| no seed: two legs | 2,680 | 51 | 361 | 0 | 0.6 |

Volumes are for the half car.

**Open: the 3.175 mm radius is not yet met everywhere.** 0.25-3.25 mm3 of the foam's concave corners, and 4-12 mm3
round the whole joint, still do not admit the ball on the 0.5 mm grid. Each such corner must be traced to a grid
artifact or a real gap (where the fillet ends at the zone or keep boundary, or a pocket against the kept chamber
wall), fixed, and the check re-run on a 0.25 mm grid. Until then the requirement is not claimed as met.

## 7. Still to do (in order)

1. The open radius item above.
2. Builder integration:
   - `rear_wing_engine` / `front_wing_engine` gain a support mode "field": today's lofts are the seed, so
     u = 0 rebuilds today's car plus the fillet and pocket;
   - the support numbers join the stage-2 numbers;
   - the mass counts the pocket foam out;
   - the manufactured foam file = body - pocket;
   - the referee gets the field gates: connected, stress, keep, radius.
3. `carhybrid/grad.py`: dD/du_j for the support numbers (the support's surface points; section 1).
4. Evolution's support moves (grow_path, thicken / thin, merge) in the loop and the layout search.
5. The nose tenon and the loft feet through `joint()` (section 4).
6. The front zone (T8.5.1) run as the rear has been.
