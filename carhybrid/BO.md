# P3: the stage-1 layout search on OpenFOAM (carhybrid/bo.py)

What it is: the evolution builder's own Bayesian layout search (`search/layout/simscale_layout.Campaign.ask`), run
on this pipeline's OpenFOAM scores instead of SimScale's:
- Sobol until 16 cars are timed;
- then a Gaussian process (Matern 2.5, one length scale per number, a noise term) on the race times;
- every legal layout on the 1 mm grid x 3 draws of the other numbers (plus 4,000 round the best car) is scored by
  expected improvement, and the batch of 8 is picked one by one with the kriging believer;
- candidates within 0.12 (unit cube) of a car that could not be built are dropped.
`bo.py` only supplies the data (`OFCampaign.data`), so the search is the same code the SimScale campaign ran.

Each car: `make_car.py` (the patched builder: evo_patches incl. legality_fixes, the stage-1 car's organic free-form
body, u = 0) -> the audit protocol on ONE frame (k = car number mod 4: the GP's noise term absorbs the frame scatter,
1.75 mN half car) -> race time = the evolution race model at the car's own mass and layout with that frame's drag and
lift (`loop.frame_score`, the same race time the loop uses). An illegal car is kept as "could not be built".

Budget (DESIGN.md P3): 48 cars = the 16-car Sobol start + 4 model batches of 8. Then `bo_confirm`: the three fastest
layouts and two references (L0039, P2's start; the team car) on frames 0-3. **P3's mark:** the best layout beats
L0039 by >= 3 SE in race time, paired by frame; then P2 (the loop) is re-run from it.

Running: `python carhybrid/bo.py init` writes batch `bo_1` into request.json and sets `launch`; a push of that
launches it. The workflow's collect hook runs `bo.py advance` on each finished `bo_*` batch and dispatches the next,
until `bo_confirm` is ingested (exit 3). State: `carhybrid/bo/state.json`, one row per car in `designs.jsonl`.

Tested offline (simulated scores, a known optimum): Sobol start, model batches that close in on the optimum
(distance of gene 0 from it 0.29 -> 0.01 after one model batch), budget stop, confirm with the paired mark.

## Paused (2026-11-28)

The owner: the Bayesian search will run on the newer car later. `state.json` carries `"paused"`, so when bo_2
finishes, `advance` ingests it and stops (exit 3) instead of dispatching bo_3. bo_1 and bo_2 were built before
`rear_wing_shape` and the sidepod cap drop. Their race times describe the older car and should not be mixed with the
newer car's in one GP. Restart: `bo.py init` on a fresh state, after the fixes before compute tests.
