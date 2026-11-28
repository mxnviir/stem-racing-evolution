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

**Quota:** spot and on-demand have separate vCPU limits per region, and a new account's are often small. Check
Service Quotas -> EC2 -> "All Standard (A, C, D, H, I, M, R, T, Z) Spot Instance Requests" (spot) and "Running On-Demand
Standard (A, C, D, H, I, M, R, T, Z) instances" (on-demand). Each must cover the vCPUs of the largest type you will
launch: 128 for m7a.32xlarge, 192 for m7a.48xlarge. hpc7a has its own quota ("Running On-Demand HPC instances").
Request any increase a day ahead.

## 2. The tests, in order

Each step is one command. Wait for its result, and launch the next only if its gates pass. The batches are defined in
`carhybrid/request.json`. The default is spot, trying the region's list of types (section 5) until one starts.

| # | command | what it proves | time (model) | spot cost |
|---|---|---|---|---|
| 1 | `bash carhybrid/aws/launch.sh aws_smoke` | the whole chain on AWS on the 129 k-cell mesh of the hosted plumbing test (team car, bg_div 4.5): setup, score with kept state, DAFoam D and L, gates F1-F3, gradient map, S3 hand-back. Runs on m7a.4xlarge or similar. | 1-2 h | < 1 USD |
| 2 | `TYPE=m7a.32xlarge STOP_ON_FAIL=0 bash carhybrid/aws/launch.sh aws_scaling` | **how DAFoam scales with cores**: the same 1.24 M-cell state solved on 128, 64 and 32 ranks (wall time per phase, s per GMRES iteration, iteration count). This decides the instance size with data. | 3-5 h | 5-10 USD |
| 3 | `bash carhybrid/aws/launch.sh f0` | **F0** memory law: team car at bg_div 8.12 / 10.23 / 12.889 (~0.74 / 1.24 / 2.12 M cells), peak memory of each; **F1-F3** at each size. The 12.889 case is the first production-mesh gradient (frame 0). | 3-6 h | 5-15 USD |
| 4 | `bash carhybrid/aws/launch.sh p1f` | **P1F** gradients: frames 1-3 of the team car on the production mesh (frame 0 is in f0) | 4-8 h | 6-20 USD |

After step 3, check F0: the three peaks should fit a line in cells, and the line must give less than 450 GB at
2.1 M cells (`fetch.sh` prints cells and peak GB per case). If the 12.889 case already finished, it measured the
production peak directly.

**Spot interruptions.** AWS can reclaim a spot instance with two minutes' notice. Results already in S3 stay. The case
that was running is lost and must be re-run with `RESUME=<run_id> bash carhybrid/aws/launch.sh <batch>`, which
starts a new instance on the same run and skips the finished cases. Each case is short (1-3 h for a 2.1 M-cell
gradient on the recommended types), and those types are in AWS's lowest interruption band (< 5 % per month), so an
interruption should be rare and cost at most one case.

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
  | np | one per physical core |
  | score_np | one per physical core (GitHub keeps the audit protocol's 4: the AWS score only provides the converged state DAFoam starts from) |
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

## 5. Which instance (prices checked 2026-11-28)

### What the work needs

- **Memory is fixed by the mesh, not by the cores.** Estimated from the coarse runs (F0 measures it):
  - ~150 kB per cell for DAFoam's Jacobian colouring and preconditioner;
  - plus the GMRES basis, (restart + 1) x states x 8 B, ~34 GB at 2.1 M cells;
  - plus ~0.5 GB per MPI rank.
  - So 2.1 M cells needs ~370-420 GB. **512 GiB is the safe size.** 384 GiB is possible only if F0 measures under
    ~350 GB. 768 GiB holds ~4 M cells, and 1 TiB ~5.5 M.
- **Speed scales with physical cores and memory bandwidth, not clock.** The score (OpenFOAM) and DAFoam's GMRES are
  sparse linear algebra, so they are limited by memory bandwidth:
  - one MPI rank per physical core (`instance_run.sh` now does this; the score uses them all too);
  - a hyperthread adds ~10-15 %;
  - a higher clock helps only the serial parts.
  - On AMD 7th/8th gen (m7a, r7a, c7a, m8a, r8a, hpc7a) one vCPU is one physical core. On Intel and on AMD 6th gen
    (m6a, r6a) two vCPUs share a core, so an r7i.16xlarge's 64 vCPU are 32 cores.
- **A serial floor of ~0.8 h per gradient** that no core count shortens: the build, snappyHexMesh (it meshes
  serially before decomposing), and DAFoam's renumberMesh / ASCII conversion of the 2.1 M-cell case.
- **Diminishing returns past ~128 cores at 2.1 M cells.** That is ~16 k cells per core. DAFoam's block preconditioner
  also needs more GMRES iterations as the subdomains get smaller. Step 2 (`aws_scaling`) measures where this bites.

So the cost of one gradient is about (price per hour) / (physical cores x speed per core x efficiency) + the serial
floor. Within one family the price per core is constant, so a bigger size is faster at nearly the same cost, until
efficiency drops.

### Time and cost per 2.1 M-cell gradient

Model, from the hosted runs:
- the 2.12 M-cell score took 1.9 h on 4 vCPU (EPYC Milan);
- coarse DAFoam: colouring and preconditioner 30 min, GMRES 0.5-0.8 s per iteration on 4 ranks at 105 k cells;
- scaled by cells (x20) and by physical cores, with 1,500-6,000 GMRES iterations per function, two functions (D, L).

Treat the times as +-2x until `aws_scaling` and F0 measure them.

Prices: on-demand from the AWS Price List API (published 2026-11-24). Spot from AWS's public spot feed where it lists
the type ("feed"), otherwise on-demand x (1 - the Spot Instance Advisor's savings) ("est."). Interruption bands are
from the Spot Instance Advisor (frequency over the last month). Spot prices move by zone and hour:
`carhybrid/aws/spot_prices.sh` gives the live ones.

**us-east-2 (Ohio), at least 512 GiB:**

| type | cores | RAM | CPU | h / gradient | on-demand $/h | per gradient | spot $/h | per gradient | interruption |
|---|---|---|---|---|---|---|---|---|---|
| **m7a.32xlarge** | 128 | 512 | EPYC Genoa | **1.4** (1.2-2.0) | 7.42 | 10.6 | **2.04** (feed) | **2.9** | **< 5 %** |
| r7a.16xlarge | 64 | 512 | EPYC Genoa | 1.8 (1.4-2.9) | 4.87 | 9.0 | 1.52 (feed) | 2.8 | 15-20 % |
| m7a.48xlarge | 192 | 768 | EPYC Genoa | 1.3 (1.1-1.8) | 11.13 | 14.3 | 2.36 (feed) | 3.0 | 5-10 % |
| m8a.48xlarge | 192 | 768 | EPYC Turin | 1.2 (1.0-1.6) | 11.69 | 14.2 | 2.80 (est.) | 3.4 | < 5 % |
| r8a.16xlarge | 64 | 512 | EPYC Turin | 1.7 (1.3-2.6) | 5.11 | 8.7 | 2.35 (est.) | 4.0 | < 5 % |
| m6a.32xlarge | 64 | 512 | EPYC Milan | 2.0 (1.5-3.2) | 5.53 | 11.2 | 1.88 (feed) | 3.8 | 5-10 % |
| r7i.16xlarge | 32 | 512 | Sapphire Rapids | 2.7 (1.9-4.5) | 4.23 | 11.5 | 1.53 (feed) | 4.2 | > 20 % |
| hpc7a.96xlarge | 192 | 768 | EPYC Genoa | 1.3 (1.1-1.8) | **7.20** | **9.2** | no spot | - | - |

**ap-south-1 (Mumbai), at least 512 GiB.** Mumbai has no AMD Genoa (7a) types in the price list.

| type | cores | RAM | CPU | h / gradient | on-demand $/h | per gradient | spot $/h | per gradient | interruption |
|---|---|---|---|---|---|---|---|---|---|
| **m8i.32xlarge** | 64 | 512 | Granite Rapids | **1.6** (1.3-2.4) | 7.13 | 11.5 | **1.21** (est.) | **2.0** | **< 5 %** |
| m6a.32xlarge | 64 | 512 | EPYC Milan | 2.0 (1.5-3.2) | 3.56 | 7.2 | 1.13 (feed) | 2.3 | < 5 % |
| r8i.16xlarge | 32 | 512 | Granite Rapids | 2.4 (1.8-4.0) | 4.59 | 11.1 | 1.05 (est.) | 2.6 | < 5 % |
| m8i.48xlarge | 96 | 768 | Granite Rapids | 1.4 (1.2-2.0) | 10.69 | 15.1 | 1.50 (est.) | 2.1 | 5-10 % |
| r6a.16xlarge | 32 | 512 | EPYC Milan | 3.2 (2.2-5.7) | 2.29 | 7.4 | 0.89 (feed) | 2.9 | < 5 % |
| r7i.16xlarge | 32 | 512 | Sapphire Rapids | 2.7 (1.9-4.5) | 4.37 | 11.7 | 1.42 (feed) | 3.8 | < 5 % |

### Recommendation

- **Spot works for this.** A gradient is a 1-3 h case, results go to S3 per case, and `RESUME` re-runs only the
  interrupted one. On the types above, spot is 3-6x cheaper than on-demand.
- **Ohio: m7a.32xlarge on spot.** It has 128 real cores, the fastest per dollar, and the < 5 % band. Fallbacks, in
  `launch.sh`'s default order: r7a.16xlarge, m7a.48xlarge (also the size for 3-4 M cells), m8a.48xlarge.
  **Mumbai: m8i.32xlarge on spot**, then m6a.32xlarge. The Mumbai spot prices marked "est." come from the advisor's
  savings, so check them with `spot_prices.sh` first. Spot prices are often lower in Mumbai, but the cores are fewer
  or slower.
- **If spot capacity is not there: hpc7a.96xlarge in Ohio on-demand**, about 9 USD per gradient. It is the cheapest
  on-demand per core by far, needs the HPC quota, and has no spot.
- **More cores than ~128 do not pay at 2.1 M cells.** The serial ~0.8 h floor and the falling efficiency cap the
  gain: 192 cores is ~10 % faster than 128 at ~40 % more per hour. They do pay at 3-5 M cells, where the 768 GiB
  sizes are needed anyway.
- **Bucket region = instance region.** Run `aws_setup.sh` with the region you launch in, so code and results
  never cross regions.

### Whole programme, on the recommended spot type

These are estimates, gradients of 2.1 M cells or their equivalent.

| part | gradients (2.1 M-cell equivalents) | machine time | spot cost |
|---|---|---|---|
| aws_smoke + aws_scaling | ~0.1 + ~2 | ~5 h | ~10 USD |
| f0 | ~2 | ~4 h | ~6 USD |
| p1f | 3 | ~4-6 h | ~9 USD |
| loop, 20 iterations (adjoint arm) | 20 | ~30 h of AWS. The line search runs on GitHub's free runners, ~2-3 h per iteration, so ~3-4 days elapsed. | ~60 USD |
| **total** | ~27 | ~45 h | **~85 USD** (about 300 USD all on-demand; about 250 USD on hpc7a on-demand) |

Not included: S3 storage (a few GB, cents) and data transfer to your computer (`fetch.sh` downloads tens to hundreds
of MB per batch: cents).
