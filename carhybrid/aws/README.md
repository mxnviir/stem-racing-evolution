# Fine-mesh DAFoam on AWS, driven by the AWS CLI

carhybrid/FINE_ADJOINT.md runs DAFoam on the production mesh (2.1 M cells; ~150 kB per cell for the Krylov adjoint,
so ~315 GB plus the GMRES basis). No GitHub-hosted runner has that memory, so those cases run on an EC2 instance that
you start from your computer with the AWS CLI. The instance runs the same job script as the GitHub workflow
(`carhybrid/fine_adjoint.sh`), on the same pinned code, and terminates itself when the batch is done.

```
your computer                         S3 bucket                      EC2 instance (fresh, one per batch)
launch.sh <batch>  -- code.tgz, -->   runs/<run_id>/code.tgz   -->   user data: aws-cli, fetch code, start
                      manifest.json   runs/<run_id>/manifest.json    instance_run.sh:
                                                                       setup (as the workflow: Python 3.11 + pinned
                                                                       deps, OpenFOAM v2412, DAFoam image, patches)
                                                                       for each case: fine_adjoint.sh
status.sh          <-- every minute   runs/<run_id>/status.json  <--   status, log, running case's log
                                      runs/<run_id>/instance.log
fetch.sh           <-- per case       runs/<run_id>/out/<case>/  <--   results the moment a case ends
                                      runs/<run_id>/DONE         <--   end; then shutdown = terminate
```

## 1. Once: AWS setup

Requirements on your computer: AWS CLI v2 configured with an account that can create S3 buckets, IAM roles and EC2
instances; git; python3; a clone of this repository on branch `claude/wizardly-brahmagupta-1edalc` and a clone of
`mxnviir/stem-racing-adjoint` next to it (`../stem-racing-adjoint`, or set `ADJ_REPO`).

```bash
REGION=us-east-1 BUCKET=<a globally unique bucket name> bash carhybrid/aws/aws_setup.sh
```

It creates (idempotent):
- the S3 bucket (public access blocked, encrypted);
- an IAM role and instance profile `carhybrid-dafoam`. Its only data permission is read and write on that one
  bucket. With the default `SSM=1` it also gets Session Manager, so you can open a shell on a running instance with
  `aws ssm start-session --target <instance id>` and no SSH;
- a security group `carhybrid-dafoam` with **no inbound rules**. The instance only makes outbound connections: apt,
  PyPI, Docker Hub, SourceForge for OpenFOAM, S3.

It writes the names to `~/.config/carhybrid/aws.env`. No keys or secrets are stored anywhere: the instance gets its
permissions from the role.

**Quota:** a new account's limit for on-demand R instances is often below 64 vCPU. Check Service Quotas -> EC2 ->
"Running On-Demand Standard (A, C, D, H, I, M, R, T, Z) instances". The r7i.16xlarge needs 64 vCPU, and the smoke
test's r7i.4xlarge needs 16. Request an increase before the first launch if needed.

## 2. The tests, in order

Each step is one command. Wait for its result, and launch the next only if its gates pass. The batches are defined in
`carhybrid/request.json`.

| # | command | instance | what it proves | about |
|---|---|---|---|---|
| 1 | `bash carhybrid/aws/launch.sh aws_smoke` | r7i.4xlarge (16 vCPU, 128 GiB) | the whole chain on AWS on the 129 k-cell mesh of the hosted plumbing test (team car, bg_div 4.5): setup, score with kept state, DAFoam D and L, gates F1-F3, gradient map, S3 hand-back | 2-4 h |
| 2 | `bash carhybrid/aws/launch.sh f0` | r7i.16xlarge (64 vCPU, 512 GiB) | **F0** memory law: team car at bg_div 8.12 / 10.23 / 12.889 (~0.74 / 1.24 / 2.12 M cells), peak memory of each; **F1-F3** at each size; the 12.889 case is the first production-mesh gradient (frame 0) | 15-30 h |
| 3 | `bash carhybrid/aws/launch.sh p1f` | r7i.16xlarge | **P1F** gradients: frames 1-3 of the team car on the production mesh (frame 0 is in f0) | 15-25 h |

Prices change, and the times above are estimates from the coarse runs, which F0 replaces with measurements. For
on-demand prices in your region, see the EC2 pricing page. In us-east-1 at the time of writing, an r7i.16xlarge was
roughly 4 USD per hour and an r7i.4xlarge roughly 1 USD per hour. `SPOT=1` is usually much cheaper, but AWS can
reclaim the machine mid-case. Finished cases are already in S3, and the interrupted one must be re-run.

After step 2, check F0: the three peaks should fit a line in cells, and the line must give less than 450 GB at
2.1 M cells (`fetch.sh` prints cells and peak GB per case). If the 12.889 case already finished, it measured the
production peak directly.

## 3. Watching, fetching, stopping

```bash
bash carhybrid/aws/status.sh                 # the latest run: instance state, status.json, the running case's log tail
bash carhybrid/aws/status.sh <run_id>
TAIL=200 bash carhybrid/aws/status.sh <run_id>
```

```bash
bash carhybrid/aws/fetch.sh <run_id>         # results into carhybrid/results/<batch>/ (+ _aws/<run_id>/ logs)
git add carhybrid/results/<batch> && git commit -m "carhybrid <batch>: results (AWS run <run_id>)" && git push
```

```bash
aws ec2 terminate-instances --instance-ids <id>    # stop at once (launch.sh prints the id; status.sh shows it)
```

You can run `fetch.sh` while a batch is still running: it fetches the cases that have finished.

Safety stops, so the instance cannot run on unnoticed:
- the instance **terminates itself** when the batch ends (`--instance-initiated-shutdown-behavior terminate`;
  shutdown at the end of `instance_run.sh`);
- a **hard stop** after `MAX_HOURS` (scheduled with `shutdown` in the user data at boot; default per batch:
  aws_smoke 6 h, f0 36 h, p1f 36 h; override with `MAX_HOURS=`);
- a case that fails, or misses gate F1 or F3, **stops the batch** (nothing downstream runs on a failed gate;
  `STOP_ON_FAIL=0` runs every case anyway);
- a setup step that fails (Python, OpenFOAM, the DAFoam image) ends the run with its reason in `status.json` and
  terminates the instance.

Termination deletes the instance's disk. Everything worth keeping is in S3 already.

## 4. What is run, exactly

- **Code:** this repository at `HEAD`. `launch.sh` warns that uncommitted changes are not shipped, so commit first.
  The builder (`EVO_SHA`) and the adjoint code (`ADJ_SHA`) are the commits pinned in
  `.github/workflows/carhybrid.yml`, so the instance runs what the workflow would. `carhybrid/results` and
  `audit/results` are left out of the archive.
- **Setup:** the same as the workflow's case job: Python 3.11 (uv's standalone build, standing in for
  actions/setup-python) with the workflow's pinned `PY_DEPS`; OpenFOAM v2412 from the adjoint repo's checksummed
  `install_openfoam.sh`; `dafoam/opt-packages:v5.1.1`; `seedhybrid/patch_adj.py`, `evo_patches/apply.py`, and
  `carhybrid/patch_daf.py` (applied by `fine_adjoint.sh`).
- **Each case:** `carhybrid/fine_adjoint.sh` with the same environment the workflow's fine_adjoint step sets, built
  from the case's fields:

  | field | default |
  |---|---|
  | k | 0 |
  | bg_div | request.json bg_div |
  | grid | |
  | np | all cores |
  | score_np | 4, the audit protocol |
  | funcs | D,L |
  | gmres_iters | 6000 |
  | gmres_restart | 200 |
  | gmres_tol | 1e-6 |
  | pc_fill | 0 |
  | pc_level | 1 |
  | primal_only | |
  | daf_primal_iters | 300 |

- **Results:** per case, the files the workflow's collect step keeps: json, logs, npz, md, png. `dafoam_sens.npz` is
  surface sensitivities only. The meshes and fields stay on the instance and are deleted with it.

## 5. Budget for the 2.1 M-cell adjoint (from the coarse 105 k-cell runs; F0 replaces these with measurements)

- adjoint states: ~10 per cell (coarse: 1.04 M states for 105 k cells), so ~21 M;
- Jacobian colouring + preconditioner: ~150 kB per cell, so ~315 GB;
- GMRES Krylov basis: (restart + 1) x states x 8 B = ~34 GB at restart 200 (the default) and ~51 GB at 300. It is
  part of the measured peak, and `resources.json` prints its estimate (`estimate_GB_krylov_basis`);
- time: coarse GMRES ran 0.77 s per iteration on 4 ranks. At 20 x the states on 64 ranks, expect ~1 s per iteration.
  Coarse solves needed 1,300 to more than 2,000 iterations; the cap is 6,000. Two functions then take up to ~3.5 h,
  plus the colouring and the ~2 h score.
