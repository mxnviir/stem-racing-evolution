# Fine adjoint on AWS Batch

Owner, 2026-11-28: the free-form body's drag gradient comes from DAFoam on the 2.1 M-cell wall-resolved mesh, on AWS.
Coarse meshes give only the general direction.

## What one job does

`search/adjoint/fine_adjoint.py`, inside this folder's container. The same chain runs on a GitHub runner in the coarse
end-to-end test (`../E2E_TEST.md`).

1. **The car.** Built from its genes and its free-form body, then judged by the referee.
2. **Production run.** The route's own CFD: a fresh snappy mesh, k-ω SST, linearUpwind, 1000 iterations from cold.
   Its run folder is kept, so the adjoint uses exactly the mesh the score was computed on.
3. **Steady.** The production state is continued until the drag stops wandering (`../steady.py`):
   - first with linearUpwind;
   - if that doesn't settle, with first-order upwind;
   - the mark: 100-iteration block means span ≤ 0.10 mN over the last 1000 iterations.
4. **DAFoam.** The primal from the steady flow, then the discrete adjoint of drag and lift. The shared core's patch
   adds lift, a per-solve GMRES record, mesh-quality limits from the environment, and up to 6000 GMRES iterations.
5. **The link** (`../link.py`). DAFoam's dD/dX per surface point becomes dD/d(number):
   - 45 numbers for CMA's coarse grid;
   - 720 control points for the fine grid.

   Surface that the build's rule masks hold doesn't move with the numbers, so it is left out.

   The shared core's `grad.py` runs beside it: drag, lift, mass and race time from the same file, with the cosine
   between the two maps logged.

Out: `result.json` (each step's numbers and health, `usable`) and `gradient.json`.

### The shared core

The owner decided on one core for both users (2026-11-28). The whole-car hybrid
(`carhybrid/`, branch `claude/wizardly-brahmagupta-1edalc`) and this stage 2 use one gradient engine:
- `carhybrid/patch_daf.py` and `carhybrid/grad.py` are packaged from a pinned commit (`batch.py` `CORE_REF`), not
  copied;
- `link.py` reproduces `grad.py`'s drag map (cosine 1.0000 on the p0 team-car adjoint);
- `link.py` adds the rule-mask exclusion that `grad.py` lacks (see `../NOTE_FOR_CARHYBRID.md`).

## 1. One-time setup

You need the AWS CLI v2 and an account that can create the stack's resources. CloudFormation creates:
- an S3 bucket;
- an ECR repository;
- a Batch compute environment of 512–768 GB instances (on-demand);
- a job queue and a job definition;
- optionally, a role that GitHub Actions can assume.

```bash
aws cloudformation deploy --stack-name stem-adjoint --template-file search/adjoint/aws/cloudformation.yaml \
  --capabilities CAPABILITY_IAM --region us-east-2 \
  --parameter-overrides VpcSubnetIds=<subnet-a>,<subnet-b> SecurityGroupIds=<sg-id> \
    GitHubRepo=mxnviir/stem-racing-evolution
```

```bash
aws cloudformation describe-stacks --stack-name stem-adjoint --region us-east-2 --query "Stacks[0].Outputs"
```

**Quota.** Check Service Quotas → EC2 → "Running On-Demand Standard (A, C, D, H, I, M, R, T, Z) instances". It must
cover the largest instance type you use: 64 vCPUs for `r7i.16xlarge`, 96 for `r7i.24xlarge`. `MaxVcpus` defaults to
192. Request any increase a day ahead.

## 2. The container, into ECR

Build it from the evolution repository's root:

```bash
docker build -t stem-adjoint:v1 -f search/adjoint/aws/Dockerfile .
```

```bash
aws ecr get-login-password --region us-east-2 | docker login --username AWS --password-stdin <account>.dkr.ecr.us-east-2.amazonaws.com
```

```bash
docker tag stem-adjoint:v1 <RepositoryUri>:v1 && docker push <RepositoryUri>:v1
```

The end-to-end workflow also pushes every build to `ghcr.io/mxnviir/stem-adjoint:<commit>`. To use one of those,
pull it and tag it as above.

## 3. Jobs from your computer

Copy `aws_config.example.json` to `aws_config.json` and fill in the stack's outputs. Then:

```bash
python search/adjoint/e2e_job.py base job_team.json --res resolved --wall 1
```

```bash
python search/adjoint/aws/batch.py submit job_team.json --config search/adjoint/aws/aws_config.json
```

```bash
python search/adjoint/aws/batch.py wait <job id> --config search/adjoint/aws/aws_config.json
```

```bash
python search/adjoint/aws/batch.py fetch <job id> results/team --config search/adjoint/aws/aws_config.json
```

`submit` packages both repositories at their committed HEAD plus the pinned core. It refuses an uncommitted tree.
While a job runs, `progress/out.tgz` in its S3 folder holds its partial results every 10 minutes. The logs go to
CloudWatch (`/aws/batch/fine-adjoint`).

**First jobs, in order.** Launch each only after the one before has passed:
1. The team car on a mid-size mesh, as the smoke test on AWS. Use `--res medium` and set `memory_mib` and `vcpus`
   lower in a copy of the config.
2. The team car on the production mesh, as above. This is the first real gradient: check `result.json`'s `steady`,
   `dafoam.converged_D`, `dafoam.primal_vs_openfoam_pct` (within 2 %) and `gradient.core.cos_dD_link_vs_core`.
3. Only then, a guided campaign (section 4).

The whole-car hybrid's batches `aws_scaling` and `f0` measure DAFoam's memory per cell and its scaling with cores on
this same mesh. Use their numbers for the instance size when they land.

## 4. A guided stage-2 campaign (GitHub Actions + AWS)

**Repository settings:**
- secret `AWS_ROLE_ARN`: the stack's `GitHubRoleArn`;
- variables `AWS_REGION`, `AWS_BUCKET`, `AWS_JOB_QUEUE`, `AWS_JOB_DEFINITION`.

Then:

```bash
python search/shape/actions_shape.py init runs/shape1 --layout-run <stage-1 run> --budget 200 --popsize 20 --adjoint guided
```

```bash
python search/shape/actions_shape.py next runs/shape1
```

Commit `runs/shape1` and `.github/shape_openfoam_request.json` (`{"dir": "runs/shape1", "auto": true}`). If `next`
wrote `adjoint_request.json`, start `adjoint_aws.yml` for the directory instead.

Each generation:
1. **`adjoint_aws.yml`:**
   - submits the job at CMA's mean and waits;
   - a job longer than one Actions run is picked up by a fresh run;
   - fetches the result and stores the gradient;
   - asks CMA for the guided generation.
2. **`shape_openfoam.yml`:** the 20 cars on OpenFOAM (4 on the adjoint's line), then ingest, tell CMA and ask the next
   step.

A generation whose adjoint failed, or whose primal wasn't steady, runs as plain CMA.

## Sizing and cost (estimates, not measured)

| | |
|---|---|
| **Memory** | DAFoam's Krylov adjoint was measured at more than 150 kB per cell on the coarse mesh. For 2.1 M cells that is about 315 GB, plus the GMRES basis. Hence 512 GB instances, with `memory_mib` 480000. |
| **Time per job** | Car about 10 min; production run with 64 ranks about 15–30 min; steady continuation, 2000–4000 iterations, about 30–90 min; DAFoam primal plus two adjoints about 1–4 h. Total about 2–6 h. |
| **Cost per gradient** | `r7i.16xlarge` on demand, about 4.2 USD/h (us-east-2, 2026), times about 2–6 h: about 10–25 USD. A guided campaign of 10 generations takes about 10 gradients. |

The whole-car hybrid's `pricing_2026-11-28/` folder has measured spot and on-demand prices for its instance types.
