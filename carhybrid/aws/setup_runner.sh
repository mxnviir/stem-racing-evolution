#!/bin/bash
# One-time setup of the big-memory AWS machine as this repository's self-hosted GitHub Actions runner (labels
# self-hosted, bigmem) for the fine_adjoint jobs (carhybrid/FINE_ADJOINT.md). Read carhybrid/aws/README.md first.
#
#   sudo RUNNER_TOKEN=<registration token> bash setup_runner.sh [path/to/install_openfoam.sh]
#
# RUNNER_TOKEN: from the repo's Settings -> Actions -> Runners -> New self-hosted runner (valid one hour). Never commit
#   it; if it is not in the environment the script asks for it (input hidden).
# Optional env: RUNNER_USER (default ubuntu: the account the runner service runs as), RUNNER_NAME (default bigmem-1),
#   REPO_URL (default https://github.com/mxnviir/stem-racing-evolution), RUNNER_VERSION (default: the latest release),
#   IDLE_MIN (default 60: the machine powers off after that many minutes without a job; 0 = never),
#   SWAP_GB (default 0: no swap file).
# install_openfoam.sh: the adjoint repo's .github/scripts/install_openfoam.sh (copy it to the machine). Without it,
#   the first job installs OpenFOAM itself (the workflow runs that script whenever OpenFOAM is missing; the runner
#   user gets passwordless sudo for it below).
# Idempotent: re-running it skips what is already done.
set -euo pipefail
[ "$(id -u)" = 0 ] || { echo "run as root (sudo)"; exit 1; }
RUNNER_USER=${RUNNER_USER:-ubuntu}
RUNNER_NAME=${RUNNER_NAME:-bigmem-1}
REPO_URL=${REPO_URL:-https://github.com/mxnviir/stem-racing-evolution}
IDLE_MIN=${IDLE_MIN:-60}
OF_SH=${1:-}
IMAGE=dafoam/opt-packages:v5.1.1
RDIR=/home/$RUNNER_USER/actions-runner

echo "== 1. system packages (Ubuntu 24.04)"
. /etc/os-release; [ "$VERSION_ID" = "24.04" ] || echo "WARNING: tested for Ubuntu 24.04, this is $PRETTY_NAME"
export DEBIAN_FRONTEND=noninteractive
apt-get update -q
apt-get install -y -q ca-certificates curl git jq unzip rsync build-essential python3 python3-pip python3-venv \
  libopengl0 libglu1-mesa libegl1 libgl1 libgomp1 libxkbcommon0 libicu74 docker.io
id "$RUNNER_USER" >/dev/null 2>&1 || useradd -m -s /bin/bash "$RUNNER_USER"
usermod -aG docker "$RUNNER_USER"
# the workflow's setup step may call sudo apt-get / install_openfoam.sh; a dedicated machine, so passwordless sudo
echo "$RUNNER_USER ALL=(ALL) NOPASSWD:ALL" > /etc/sudoers.d/90-actions-runner && chmod 440 /etc/sudoers.d/90-actions-runner
systemctl enable --now docker

if [ "${SWAP_GB:-0}" -gt 0 ] && ! swapon --show | grep -q /swapfile_bigmem; then
  echo "== swap ${SWAP_GB} GB"
  fallocate -l "${SWAP_GB}G" /swapfile_bigmem && chmod 600 /swapfile_bigmem && mkswap /swapfile_bigmem && swapon /swapfile_bigmem
  grep -q swapfile_bigmem /etc/fstab || echo "/swapfile_bigmem none swap sw 0 0" >> /etc/fstab
fi

echo "== 2. the DAFoam image"
for i in 1 2 3; do docker pull -q "$IMAGE" && break; sleep 30; done
docker image inspect "$IMAGE" >/dev/null

echo "== 3. OpenFOAM v2412"
if [ -f /usr/lib/openfoam/openfoam2412/etc/bashrc ]; then
  echo "already installed"
elif [ -n "$OF_SH" ] && [ -f "$OF_SH" ]; then
  sudo -u "$RUNNER_USER" OF_DEBS="/home/$RUNNER_USER/ofdebs" bash "$OF_SH"
else
  echo "no install_openfoam.sh given: the first job installs OpenFOAM (adj/.github/scripts/install_openfoam.sh)"
fi

echo "== 4. the runner (labels self-hosted, bigmem) as a service"
if [ -f "$RDIR/.runner" ]; then
  echo "already registered: $(jq -r .agentName "$RDIR/.runner" 2>/dev/null)"
else
  if [ -z "${RUNNER_VERSION:-}" ]; then
    RUNNER_VERSION=$(curl -sf https://api.github.com/repos/actions/runner/releases/latest | jq -r .tag_name | sed 's/^v//')
  fi
  echo "runner version $RUNNER_VERSION"
  mkdir -p "$RDIR"
  curl -sfL -o /tmp/actions-runner.tgz \
    "https://github.com/actions/runner/releases/download/v${RUNNER_VERSION}/actions-runner-linux-x64-${RUNNER_VERSION}.tar.gz"
  tar -xzf /tmp/actions-runner.tgz -C "$RDIR" && rm -f /tmp/actions-runner.tgz
  chown -R "$RUNNER_USER:$RUNNER_USER" "$RDIR"
  "$RDIR/bin/installdependencies.sh" || true
  if [ -z "${RUNNER_TOKEN:-}" ]; then read -r -s -p "registration token (Settings -> Actions -> Runners): " RUNNER_TOKEN; echo; fi
  sudo -u "$RUNNER_USER" bash -c "cd '$RDIR' && ./config.sh --unattended --replace --url '$REPO_URL' --token '$RUNNER_TOKEN' \
      --name '$RUNNER_NAME' --labels bigmem --work _work"
  unset RUNNER_TOKEN
  (cd "$RDIR" && ./svc.sh install "$RUNNER_USER" && ./svc.sh start)
fi

echo "== 5. power off when idle (IDLE_MIN=$IDLE_MIN)"
if [ "$IDLE_MIN" -gt 0 ]; then
  cat > /usr/local/bin/bigmem-idle-stop.sh <<EOF
#!/bin/bash
# powers the machine off after $IDLE_MIN minutes without a running job (Runner.Worker is the process of a job).
# On EBS-backed EC2, "shutdown" from inside = instance STOPPED (disk kept, no compute billed) if the instance's
# shutdown behaviour is "stop" (the default).
S=/var/lib/bigmem-idle-since
if pgrep -f Runner.Worker >/dev/null; then rm -f \$S; exit 0; fi
[ -f \$S ] || date +%s > \$S
if [ \$(( \$(date +%s) - \$(cat \$S) )) -ge $((IDLE_MIN * 60)) ]; then rm -f \$S; logger "bigmem idle ${IDLE_MIN} min: power off"; shutdown -h now; fi
EOF
  chmod 755 /usr/local/bin/bigmem-idle-stop.sh
  cat > /etc/systemd/system/bigmem-idle-stop.service <<'EOF'
[Unit]
Description=Power off the bigmem runner when idle
[Service]
Type=oneshot
ExecStart=/usr/local/bin/bigmem-idle-stop.sh
EOF
  cat > /etc/systemd/system/bigmem-idle-stop.timer <<'EOF'
[Unit]
Description=Check every 5 min whether the bigmem runner is idle
[Timer]
OnBootSec=15min
OnUnitActiveSec=5min
[Install]
WantedBy=timers.target
EOF
  rm -f /var/lib/bigmem-idle-since
  systemctl daemon-reload && systemctl enable --now bigmem-idle-stop.timer
fi

echo "== check"
echo "cores $(nproc), RAM $(awk '/MemTotal/{printf "%.0f GB", $2/1048576}' /proc/meminfo), disk free $(df -h / | awk 'NR==2{print $4}')"
sudo -u "$RUNNER_USER" docker run --rm "$IMAGE" bash -c 'echo docker ok as runner user'
(cd "$RDIR" && ./svc.sh status | head -5) || true
echo "done: the runner should show as Idle under Settings -> Actions -> Runners (labels self-hosted, Linux, X64, bigmem)"
