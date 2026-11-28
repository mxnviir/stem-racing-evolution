# The big-memory AWS machine: a self-hosted runner for the fine adjoint

carhybrid/FINE_ADJOINT.md runs DAFoam on the production mesh (2.1 M cells, ~150 kB per cell for the Krylov adjoint,
so ~315 GB). No GitHub-hosted runner has that memory, so one AWS machine is registered as a **self-hosted runner** of
this repository with the labels `self-hosted` and `bigmem`. The workflow (`.github/workflows/carhybrid.yml`) sends
every `fine_adjoint` case there (unless the case says `"runner": "hosted"`); everything else stays on hosted runners.

> **Security: use a self-hosted runner only while this repository is private.** On a public repository anyone can
> open a pull request whose workflow runs arbitrary code on the machine (and the runner user has sudo and Docker).
> If the repository is ever made public, remove the runner first (Settings -> Actions -> Runners -> Remove).

## 1. The instance

| | choice | why |
|---|---|---|
| type | **r7i.16xlarge** (64 vCPU, 512 GiB) or the same class: r6i.16xlarge, r7a.16xlarge (AMD), r6a.16xlarge | memory: 2.1 M cells x 150 kB = 315 GB plus OS and the primal; 64 ranks for DAFoam |
| price | **check the current price** in the EC2 console for your region (on-demand and spot); at the time of writing an r7i.16xlarge was a few USD per hour on demand. Spot is cheaper but can be reclaimed mid-job (the job then fails and must be re-run) | |
| OS | **Ubuntu Server 24.04 LTS** (x86_64), the Canonical AMI | the workflow and the OpenFOAM packages are for noble |
| disk | **500 GB gp3** root volume (3000 IOPS / 125 MB/s are enough) | Docker image (~10 GB), OpenFOAM, python deps, the 2.1 M-cell state and DAFoam case (ASCII fields, processor dirs: tens of GB), headroom |
| network | a security group with **no inbound rules** except SSH from your own IP (the runner only makes outbound HTTPS to GitHub) | |
| shutdown behaviour | **Stop** (the default for EBS-backed instances) | the idle timer powers the machine off; Stop keeps the disk |

Memory headroom: gate F0 measures the real law first (team car, three mesh sizes, batch `f0`); the line it fits
must extrapolate to < 450 GB at 2.1 M cells before anything relies on the production-mesh adjoint.

## 2. Setup (once)

```bash
# on your computer
scp -i key.pem carhybrid/aws/setup_runner.sh ubuntu@<instance>:~
scp -i key.pem <adjoint repo>/.github/scripts/install_openfoam.sh ubuntu@<instance>:~    # optional, see below
ssh -i key.pem ubuntu@<instance>

# on the instance: paste the registration token when asked (repo Settings -> Actions -> Runners ->
# New self-hosted runner -> Linux; the token in the "./config.sh --token ..." line, valid for one hour)
sudo bash setup_runner.sh ~/install_openfoam.sh
```

`setup_runner.sh` (idempotent; re-run it after a failure):
1. installs git, Python 3, the system libraries the car builder needs, Docker (the runner user `ubuntu` joins the
   `docker` group and gets passwordless sudo: the workflow's setup step may call `sudo apt-get`);
2. pulls the DAFoam image `dafoam/opt-packages:v5.1.1`;
3. installs OpenFOAM v2412 with the adjoint repo's `.github/scripts/install_openfoam.sh` (checksummed packages);
   without that file the first job installs it (the workflow runs the same script whenever OpenFOAM is missing);
4. downloads the GitHub Actions runner, registers it with `--labels bigmem` (GitHub adds `self-hosted`, `Linux`,
   `X64`) and installs it as a systemd service (`~ubuntu/actions-runner/svc.sh status`);
5. installs a timer that **powers the machine off after 60 idle minutes** (no `Runner.Worker` process;
   `IDLE_MIN=0` disables it, `IDLE_MIN=30` shortens it).

The token is never written to the repository or to disk by the script; it is used once by `config.sh`. Python
itself comes from `actions/setup-python` in each job (it caches 3.11 in the runner's tool folder).

## 3. Running and stopping

- **Start** the instance (EC2 console, or `aws ec2 start-instances --instance-ids <id>`) before dispatching a batch
  with `fine_adjoint` cases; the runner service starts with it and picks up queued jobs (a queued job waits up to
  24 h for a runner).
- Jobs run **one at a time** on this machine (one runner). A batch with several `fine_adjoint` cases runs them in
  sequence; the hosted cases of the same batch run in parallel as usual.
- **Stop** it when idle: the timer does it after `IDLE_MIN`; or by hand (`aws ec2 stop-instances --instance-ids <id>`).
  A stopped instance costs only its disk (500 GB gp3: check the per-GB-month price).
- Logs: the job log on GitHub (everything the job prints); on the machine `journalctl -u actions.runner.*`.
- Long jobs: the workflow allows `fine_adjoint` 1400 minutes on self-hosted (the job token lasts 24 h).

## 4. Verify

1. Settings -> Actions -> Runners shows `bigmem-1` **Idle** with labels `self-hosted, Linux, X64, bigmem`.
2. A tiny test job. Add this as `.github/workflows/bigmem_check.yml` on the branch (or paste it into any branch),
   run it from the Actions tab, then delete it:

   ```yaml
   name: bigmem check
   on: workflow_dispatch
   jobs:
     check:
       runs-on: [self-hosted, bigmem]
       timeout-minutes: 15
       steps:
         - run: |
             nproc; free -g; df -h /
             test -f /usr/lib/openfoam/openfoam2412/etc/bashrc && echo "OpenFOAM ok" || echo "OpenFOAM missing (first job installs it)"
             docker run --rm dafoam/opt-packages:v5.1.1 bash -c '. /home/dafoamuser/dafoam/loadDAFoam.sh >/dev/null && python -c "import dafoam; print(\"DAFoam ok\")"'
   ```
3. Then the real thing, cheapest first: batch `f0` (the F0 memory ladder: team car, frame 0, bg_div 8.12 / 10.23 /
   12.889) with `gh workflow run carhybrid.yml --ref claude/wizardly-brahmagupta-1edalc -f batch=f0`. Each job prints
   cells, ranks, peak memory and wall times, and writes them to `resources.json`.

## 5. What the fine_adjoint job does on this machine

`carhybrid/fine_adjoint.sh`: build the car -> the audit-protocol score of frame k (`audit/run_case.py`, 4 ranks as
every hosted score; `"score_np"` in the case overrides) with the converged state kept -> steadiness (gate F2) ->
`dafoam/make_case.py` on that state -> DAFoam D and L with all cores (`"np"` overrides), memory sampled every 10 s
-> `carhybrid/grad.py`. Outputs: `carhybrid/results/<batch>/<case>/` (result, steady, dafoam_sens, resources, grad,
fine_adjoint.json with the gates F1-F3).
