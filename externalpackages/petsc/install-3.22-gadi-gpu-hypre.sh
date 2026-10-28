#!/bin/bash
# PETSc 3.22 with CUDA + Hypre BoomerAMG for NCI Gadi GPU nodes.
#
# Identical to install-3.22-gadi-gpu.sh but adds --download-hypre=1.
# Installs to install-gpu-hypre (separate from the base GPU install)
# so the working install is not disrupted.
#
# Runtime: point LD_LIBRARY_PATH at install-gpu-hypre/lib before
# install-gpu/lib to use this build without rebuilding ISSM.
#
# Before running:
#   module load cuda/12.3.2 openmpi/4.1.3
#   export ISSM_DIR=/g/data/au88/jh7060/ISSM
#   bash install-3.22-gadi-gpu-hypre.sh
#
# CUDA_ARCH (default 80 = A100) sets the GPU architecture for both PETSc and hypre:
# --with-cuda-gencode does not reach hypre's own build, which otherwise defaults
# to sm_70 when configured without a GPU present. Use 90 for H100 / Grace Hopper.
set -eu

## Constants
VER="3.22.3"
CUDA_ARCH="${CUDA_ARCH:-80}"

PETSC_DIR="${ISSM_DIR}/externalpackages/petsc/src-gpu-hypre"
PREFIX="${ISSM_DIR}/externalpackages/petsc/install-gpu-hypre"

# Require CUDA_HOME / CUDA_ROOT to be set by the module system
CUDA_DIR="${CUDA_HOME:-${CUDA_ROOT:-/usr/local/cuda}}"
if [ ! -d "${CUDA_DIR}" ]; then
    echo "ERROR: CUDA directory not found at ${CUDA_DIR}."
    echo "       Load the CUDA module first: module load cuda/12.3.2"
    exit 1
fi

# Build PETSc against the SAME MPI that ISSM uses (system OpenMPI), rather than
# letting PETSc download its own MPICH.  A mismatched MPI causes issm.exe to be
# linked against two MPI runtimes at once, which breaks any np>1 (multi-rank)
# run: each rank initialises as its own standalone MPI_COMM_WORLD rank-0.
# The OpenMPI compiler wrappers must be on PATH (module load openmpi/4.1.3).
for w in mpicc mpicxx mpif90; do
    if ! command -v ${w} &>/dev/null; then
        echo "ERROR: ${w} not found. Load OpenMPI first: module load openmpi/4.1.3"
        exit 1
    fi
done
if ! mpicc --showme:version 2>/dev/null | grep -qi "open mpi"; then
    echo "WARNING: mpicc does not appear to be OpenMPI; check your modules."
fi

# Environment
if [ -z ${LDFLAGS+x} ]; then
    LDFLAGS=""
fi

# Download source (skip if tarball already present)
if [ ! -f "petsc-${VER}.tar.gz" ]; then
    ${ISSM_DIR}/scripts/DownloadExternalPackage.sh \
        "https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-${VER}.tar.gz" \
        "petsc-${VER}.tar.gz"
fi

# Unpack source
tar -zxvf petsc-${VER}.tar.gz

# Cleanup
rm -rf ${PREFIX} ${PETSC_DIR}
mkdir -p ${PETSC_DIR}

# Move source to $PETSC_DIR
mv petsc-${VER}/* ${PETSC_DIR}
rm -rf petsc-${VER}

# Configure
cd ${PETSC_DIR}
./configure \
    --prefix="${PREFIX}" \
    --PETSC_DIR="${PETSC_DIR}" \
    --LDFLAGS="${LDFLAGS}" \
    --CFLAGS="-g -O2" --CXXFLAGS="-g -O2" --FFLAGS="-g -O2" \
    --with-debugging=0 \
    --with-valgrind=0 \
    --with-x=0 \
    --with-ssl=0 \
    --with-pic=1 \
    --with-cc=mpicc \
    --with-cxx=mpicxx \
    --with-fc=mpif90 \
    --download-fblaslapack=1 \
    --download-metis=1 \
    --download-mumps=1 \
    --download-parmetis=1 \
    --download-scalapack=1 \
    --download-zlib=1 \
    --with-cuda=1 \
    --with-cuda-dir="${CUDA_DIR}" \
    --with-cudac="${CUDA_DIR}/bin/nvcc" \
    --with-cuda-gencode="arch=compute_${CUDA_ARCH},code=sm_${CUDA_ARCH}" \
    --with-hypre-gpu-arch=${CUDA_ARCH} \
    --download-hypre=1

# Compile and install (parallel on PBS node; fall back to 4 locally)
J="${PBS_NCPUS:-4}"
make -j "${J}"
make install
