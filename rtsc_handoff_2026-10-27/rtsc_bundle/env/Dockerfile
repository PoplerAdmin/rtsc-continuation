# syntax=docker/dockerfile:1
# rtsc image: Ubuntu 24.04 + conda-forge Quantum ESPRESSO 7.x (pw.x, ph.x, q2r.x, matdyn.x,
# epw.x; OpenMPI) + the rtsc Python environment (env/environment.yml).
#
# Build from the BUNDLE ROOT (the directory that contains rtsc/, env/, verify_all.sh):
#   docker build -f env/Dockerfile -t rtsc:0.1.0 .
# Run:
#   docker run --rm -it -v "$PWD/work:/work" -v /data/SSSP:/pseudo:ro rtsc:0.1.0 rtsc demo
#   docker run --rm -it rtsc:0.1.0 /opt/rtsc_bundle/verify_all.sh --skip-lean --workdir /work/verify
#
# Why conda-forge QE and not the Ubuntu package: Ubuntu 24.04's `quantum-espresso` (6.7) is built
# with _FORTIFY_SOURCE, and its md5 helper (snprintf with maxlen = sizeof(char*) = 8 on a 33-byte
# buffer) makes pw.x abort at start-up with "buffer overflow detected". The upstream line is
# unchanged up to QE 7.6, so a from-source build with Ubuntu's default gcc flags can hit it too.
# The conda-forge qe 7.6 binaries do not import __snprintf_chk (checked with nm -D), so the abort
# cannot occur there. The Ubuntu-package alternative with the LD_PRELOAD shim is sketched at the
# end of this file.
#
# NOT TESTED: this image has not been built (no Docker build was run while preparing the bundle).
# The environment file was checked with a micromamba dry-run solve (qe 7.6 mpi_openmpi, python
# 3.12). Running QE inside containers on a cluster usually needs the host MPI (Apptainer with a
# matching OpenMPI, or the site's QE build); for multi-node runs prefer the conda environment or
# site modules directly.

FROM ubuntu:24.04

ARG MAMBA_VERSION=2.9.0
ENV DEBIAN_FRONTEND=noninteractive \
    LANG=C.UTF-8 LC_ALL=C.UTF-8 \
    MAMBA_ROOT_PREFIX=/opt/conda \
    PATH=/opt/conda/envs/rtsc/bin:/opt/conda/bin:/usr/local/bin:/usr/bin:/bin \
    OMP_NUM_THREADS=1 \
    PYTHONDONTWRITEBYTECODE=1

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

# hadolint ignore=DL3008
RUN apt-get update \
 && apt-get install -y --no-install-recommends ca-certificates curl bzip2 tini \
 && rm -rf /var/lib/apt/lists/*

# micromamba: a single static binary, from the conda-forge channel itself
RUN curl -fsSL "https://conda.anaconda.org/conda-forge/linux-64/micromamba-${MAMBA_VERSION}-0.tar.bz2" \
      | tar -xj -C /usr/local bin/micromamba \
 && micromamba --version

COPY env/environment.yml /tmp/environment.yml
RUN micromamba create -y -n rtsc -f /tmp/environment.yml \
 && micromamba clean -a -y \
 && for x in pw.x ph.x q2r.x matdyn.x epw.x mpirun; do command -v "$x" || exit 1; done

# the bundle (Lean sources included; the Lean toolchain itself is not installed in this image)
COPY . /opt/rtsc_bundle
RUN python -m pip install --no-deps --no-build-isolation --no-cache-dir /opt/rtsc_bundle/rtsc \
 && rm -rf /opt/rtsc_bundle/rtsc/build /opt/rtsc_bundle/rtsc/*.egg-info \
 && rtsc --help > /dev/null \
 && python -c "import flint, numpy, scipy, pymatgen, ase, spglib, mp_api, rtsc.certify"

# unprivileged user (uid 1000 is taken by "ubuntu" in ubuntu:24.04); OpenMPI refuses to run
# as root without extra flags. The bundle is made writable for verify_all.sh/pytest.
RUN useradd --create-home --uid 10001 --user-group rtsc \
 && mkdir -p /work /pseudo \
 && chown -R 10001:10001 /work /opt/rtsc_bundle
USER 10001:10001
WORKDIR /work

ENTRYPOINT ["/usr/bin/tini", "--"]
CMD ["rtsc", "--help"]

# --- Alternative (not recommended): Ubuntu's QE 6.7 package + FORTIFY shim --------------------
# Replace the micromamba/QE part with:
#   RUN apt-get update && apt-get install -y --no-install-recommends quantum-espresso \
#         openmpi-bin gcc libc6-dev python3-venv && rm -rf /var/lib/apt/lists/*
#   RUN gcc -O2 -shared -fPIC -o /usr/local/lib/libsnprintf_chk_shim.so \
#         /opt/rtsc_bundle/rtsc/rtsc/dfpt/contrib/snprintf_chk_shim.c
# and in config.toml (never set LD_PRELOAD globally: it weakens FORTIFY for every program):
#   [dfpt.launcher]  env = { LD_PRELOAD = "/usr/local/lib/libsnprintf_chk_shim.so" }
#                    mpi_args = ["-x", "LD_PRELOAD"]
# This route (QE 6.7 + shim, mpirun -np 2) is the one the bundled real-QE test data were made
# with (rtsc/tests/data_dfpt/*/PROVENANCE.txt). EPW 5.3 in that package worked the same way.
