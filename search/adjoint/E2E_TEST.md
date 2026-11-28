# The adjoint chain, end to end on a coarse mesh (written 2026-11-28, before any run)

Owner 2026-11-28: link the adjoint to the free-form body; the real gradients come from DAFoam on the 2.1 M-cell
wall-resolved mesh on AWS Batch; before AWS, "prove the whole chain end to end on a coarse mesh on GitHub Actions"
("Yes, run it"). The coarse gradient is only the general direction; this test is about the MECHANISM: that every
piece runs and joins up, inside the very container AWS will run.

## What runs (`.github/workflows/adjoint_e2e.yml`)
1. **image**: `search/adjoint/aws/Dockerfile` built on a runner (DAFoam v5.1.1 + OpenFOAM v2412 + the pipeline's
   Python) and pushed to GHCR, tagged with the commit.
2. **adjoint** (one runner, 4 cores, 16 GB + 16 GB swap): `batch.py package` makes the job archive exactly as for
   AWS; the container runs `entrypoint.sh` on it (IN_TGZ / OUT_DIR instead of S3). The car: the team layout
   (`layout_genes.x_default()`), its own free-form body as the seed, CMA vector v = 0. The mesh: the pipeline's
   coarse mesh without wall layers (the setting on which DAFoam's coarse gradient was validated, r 0.98 on bumps;
   a wall-resolved mesh's adjoint does not fit a runner). Production 1000 iterations from cold; steady continuation
   2000 iterations (linearUpwind; upwind if not steady; coarse steadiness mark 0.02 mN); DAFoam Krylov, PC fill 0,
   level 1, restart 200 (the coarse settings that converged on runners); the link -> 45 coarse gradients.
3. **fd** (six runners, same container): the three body numbers with the largest |dD/dv| (fixed rule), each at
   v = +delta e_k and -delta e_k, delta_k = clip(3 mN / |dD/dv_k (half car)|, 0.1, 0.6): each car built, meshed fresh,
   production + the same steady continuation (`stop_after: steady`). Measured slope s_k = (D+ - D-) / (2 delta_k) on
   the steady model's half-car drag (mean of its last 1000 iterations).
4. **verdict** (`e2e_verdict.py`): the marks below; results committed to `docs/checks/adjoint_e2e/<run id>/`.

## Marks (fixed now)
- **E1 (the chain):** the image builds; the adjoint job ends with DAFoam `ok` (no adjoint failure), a steady base
  (either rung), DAFoam's primal within 2 % of OpenFOAM's steady drag, and a gradient for all 45 body numbers.
- **E2 (the gradient means something):** for the three tested numbers, the predicted slope p_k and the measured
  s_k have the same sign for all three, and p_k / s_k lies in [0.5, 4] for at least two of three (DAFoam's coarse
  steady gradient was seen to overstate realised gains about 2.5x; a wrong sign convention, frame, unit or chain
  would miss by far more, or flip signs). A measured slope smaller than 2 SE of itself (fresh-mesh noise, SE from
  the coarse mesh-to-mesh sd 0.32 mN: SE(s_k) = 0.32 sqrt(2) / (2 delta_k)) counts as no evidence either way and is
  reported; if two or more are like that, E2 is INCONCLUSIVE (re-run with larger delta).
- **PASS = E1 and E2.**
- Reported, not marked: each step's time and memory, the steady rung used, the active points and their share of
  the sensitivity, DAFoam's residual history, the sum of dD/dX (translation check), the link's exact geometric
  test (`test_link.py`, run separately on a built car: volume and first moments).

## What it does NOT show
Anything about the 2.1 M-cell wall-resolved adjoint: its memory, its convergence (it stalled on some designs
before), the wall-resolved boundary conditions in DAFoam. The first AWS job is that test; its result is read the
same way (E1, and E2 with fresh-mesh checks on the fine mesh if wanted).

## Amendment 1 (2026-11-28 18:10, after e2e-4; before any probe)
**e2e-4 (run 46163110187): E1 FAIL, the adjoint did not converge.**
- **What worked in the container:** the car built; the production run was coarse (93,932 cells, 229.20 mN); the
  linearUpwind continuation was steady (span 0.0002 mN, so the upwind rung was not used); DAFoam's primal reproduced
  OpenFOAM's drag to 0.0001 %.
- **What failed:** both Krylov adjoints (D and L) stalled. D went 3.15e-2 to 2.12e-2 by iteration 200, then stayed
  flat to 2000 (PetscConvergedReason -3). L stalled at 5.2e-2. The settings were the registered coarse ones (PC fill 0,
  level 1, restart 200).
- **Why that is not new:** the same settings converged on the coarse mesh's first-order upwind steady flow before
  (adjoint repo ROUTES.md 2026-11-20) and stalled on linearUpwind states.

**Probes (one run, no forward checks).** The same car and chain, three variants:
- `pc`: linearUpwind base, stronger preconditioner (fill 1, level 2, restart 400, 3000 iterations);
- `upwind`: an upwind-only steady base, the registered preconditioner;
- `upwind_pc`: both changes.

**Read-out:** whether DAFoam's drag adjoint converges (`converged_D`), its iterations and time. The full test (E1, E2
with the six forward checks) is then re-run with the first converging variant in this order: `pc` (it keeps the
owner's linearUpwind-first choice), then `upwind`, then `upwind_pc`. Taking the adjoint on an upwind base for the AWS
jobs changes the owner's choice, so it goes to the owner first. Marks unchanged.
