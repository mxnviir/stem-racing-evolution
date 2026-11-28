#!/bin/bash
# Install OpenFOAM v2412 (openfoam.com) on an Ubuntu 24.04 runner without depending on one download mirror.
# dl.openfoam.com redirects to SourceForge mirrors, and one of them has served truncated files ("File has unexpected
# size ... Mirror sync in progress?"). So: take the packages from $OF_DEBS (an actions/cache directory) if they are
# there, otherwise fetch each one from a list of named SourceForge mirrors, keep it only if its SHA-256 matches the
# repo's index (dists/noble/main/binary-*/Packages, release 2412.260127-1), and install the set with apt, which takes
# the Ubuntu dependencies from Ubuntu's own mirrors. Stops with an error if OpenFOAM is still missing.
set -u
DEBS=${OF_DEBS:-$HOME/ofdebs}
POOL=project/openfoam/repos/deb/dists/noble/main/pool/2412_260127
MIRRORS="netix phoenixnap cfhcable altushost-swe sinalbr deac-riga freefr master"
PKGS="
binary-all/openfoam2412-common_2412.260127-1_all.deb 8d734715802a207de0875ecfc2a2b44399c28cbd4a3f59471a682c56ddaf0310
binary-all/openfoam2412-source_2412.260127-1_all.deb 48fe314aa58984506e082c76ade986ac4b5e95f7472c177903cda72be25a41e0
binary-all/openfoam2412-tutorials_2412.260127-1_all.deb a11835f716b90d0aeda93d78f2d4dfe3fd7c7b74ff59b4892ca2eaa6d8dbca80
binary-amd64/openfoam2412-default_2412.260127-1_amd64.deb 80ddf15bac5fcd80f1d504c94ac32abcc3cf8358cd7c6357845f438be386123f
binary-amd64/openfoam2412-dev_2412.260127-1_amd64.deb 976aa19da55836c92657f787c42a2af27b39579c6944d1fc8cf2d964178fec83
binary-amd64/openfoam2412-tools_2412.260127-1_amd64.deb f19b8fc2ec2a0b25f8489bf9f3f3a7da8fb3451c6b76fab90cc6550608f52ed3
binary-amd64/openfoam2412_2412.260127-1_amd64.deb 7022b3d4d20fbda9198e7e85510785ce094d91424713233f7571e97cb0d1c80c
"
mkdir -p "$DEBS"
ok=1
while read -r p sum; do
  [ -z "$p" ] && continue
  f=$DEBS/$(basename "$p")
  if [ -f "$f" ] && echo "$sum  $f" | sha256sum -c --quiet - 2>/dev/null; then echo "cached $(basename "$p")"; continue; fi
  got=0
  for round in 1 2 3; do
    for m in $MIRRORS; do
      curl -sfL -m 600 -o "$f.part" "https://$m.dl.sourceforge.net/$POOL/$p" \
        && echo "$sum  $f.part" | sha256sum -c --quiet - 2>/dev/null && { mv "$f.part" "$f"; echo "ok $m $(basename "$p")"; got=1; break 2; }
      echo "bad $m $(basename "$p")"
    done
    sleep 20
  done
  rm -f "$f.part"
  [ $got = 1 ] || ok=0
done <<< "$PKGS"
if [ $ok = 1 ]; then
  for i in 1 2 3 4; do sudo apt-get install -y -q "$DEBS"/*.deb && break; sleep 30; sudo apt-get update -q || true; done
fi
if [ ! -f /usr/lib/openfoam/openfoam2412/etc/bashrc ]; then
  echo "pinned packages did not install; trying the openfoam.com apt repo"
  for i in 1 2 3; do (curl -sf https://dl.openfoam.com/add-debian-repo.sh | sudo bash) && sudo apt-get install -y -q openfoam2412-default && break; sudo apt-get clean; sleep 60; done
fi
test -f /usr/lib/openfoam/openfoam2412/etc/bashrc || { echo "OpenFOAM install failed"; exit 1; }
echo "OpenFOAM v2412 installed"
