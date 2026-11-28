# The fine adjoint's container (owner 2026-11-28: DAFoam adjoint on the 2.1 M-cell mesh, on AWS Batch).
#
# One image runs the whole chain of search/adjoint/fine_adjoint.py: the car build, the production run (OpenFOAM
# v2412, the pipeline's own solver), the steady continuation (v2412), DAFoam's primal + discrete adjoint (DAFoam's
# own OpenFOAM v2506 AD build, loaded only for that step) and the link. The CODE is not baked in: every job brings
# a snapshot of both repositories at a git commit (batch.py package), so the image changes only when the software
# does.
#
#   docker build -t stem-adjoint:v1 -f search/adjoint/aws/Dockerfile .          (from the evolution repo root)
#
# Base: DAFoam v5.1.1 (Ubuntu 24.04, the same OS as the GitHub runners, so the checksum-pinned OpenFOAM v2412
# packages of .github/scripts/install_openfoam.sh install unchanged). The base image is an import with no ENV, so
# PATH etc. are set here.
FROM dafoam/opt-packages:v5.1.1

USER root
ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin \
    HOME=/root \
    DEBIAN_FRONTEND=noninteractive \
    LANG=C.UTF-8

RUN for i in 1 2 3 4 5; do apt-get update -q && break; sleep 20; done \
 && apt-get install -y -q --no-install-recommends sudo curl ca-certificates python3 python3-venv python3-pip \
      libopengl0 libglu1-mesa libegl1 libgl1 libgomp1 libxkbcommon0 \
 && rm -rf /var/lib/apt/lists/*

# OpenFOAM v2412 at /usr/lib/openfoam/openfoam2412 (checksum-pinned packages, several mirrors)
COPY .github/scripts/install_openfoam.sh /tmp/install_openfoam.sh
RUN for i in 1 2 3; do apt-get update -q && break; sleep 20; done \
 && OF_DEBS=/tmp/ofdebs bash /tmp/install_openfoam.sh \
 && test -f /usr/lib/openfoam/openfoam2412/etc/bashrc \
 && rm -rf /tmp/ofdebs /var/lib/apt/lists/*

# The pipeline's Python (car builder, meshing scripts, link) in its own venv; DAFoam keeps its conda Python
COPY requirements.txt /tmp/requirements.txt
RUN python3 -m venv /opt/venv \
 && /opt/venv/bin/pip install -q --no-cache-dir -r /tmp/requirements.txt boto3

COPY search/adjoint/aws/entrypoint.sh /usr/local/bin/entrypoint.sh
COPY search/adjoint/aws/s3io.py /usr/local/bin/s3io.py
RUN chmod +x /usr/local/bin/entrypoint.sh

ENV PATH=/opt/venv/bin:$PATH \
    DAFOAM_RUNNER=native \
    DAFOAM_LOAD=/home/dafoamuser/dafoam/loadDAFoam.sh \
    OF_BASHRC=/usr/lib/openfoam/openfoam2412/etc/bashrc \
    QT_QPA_PLATFORM=offscreen \
    OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
    OMPI_MCA_rmaps_base_oversubscribe=1 OMPI_MCA_hwloc_base_use_hwthreads_as_cpus=1 \
    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

WORKDIR /scratch
ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
