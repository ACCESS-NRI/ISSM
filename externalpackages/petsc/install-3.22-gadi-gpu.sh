#!/bin/bash
# PETSc 3.22 with CUDA support for NCI Gadi GPU nodes (V100 / A100).
#
# Before running, load the CUDA module and set ISSM_DIR, e.g.:
#   module load cuda/12.3.2
#   export ISSM_DIR=/g/data1b/au88/jh7060/ISSM
#   bash install-3.22-gadi-gpu.sh
#
# The GPU-enabled PETSc is installed to a separate prefix so it can
# coexist with the CPU-only install.
set -eu

## Constants
VER="3.22.3"

PETSC_DIR="${ISSM_DIR}/externalpackages/petsc/src-gpu"  # DO NOT CHANGE THIS
PREFIX="${ISSM_DIR}/externalpackages/petsc/install-gpu"  # separate from CPU install

# Require CUDA_HOME / CUDA_ROOT to be set by the module system
CUDA_DIR="${CUDA_HOME:-${CUDA_ROOT:-/usr/local/cuda}}"
if [ ! -d "${CUDA_DIR}" ]; then
    echo "ERROR: CUDA directory not found at ${CUDA_DIR}."
    echo "       Load the CUDA module first: module load cuda/12.3.2"
    exit 1
fi

# Environment
if [ -z ${LDFLAGS+x} ]; then
    LDFLAGS=""
fi

# Download source
${ISSM_DIR}/scripts/DownloadExternalPackage.sh \
    "https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-${VER}.tar.gz" \
    "petsc-${VER}.tar.gz"

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
    --with-debugging=0 \
    --with-valgrind=0 \
    --with-x=0 \
    --with-ssl=0 \
    --with-pic=1 \
    --download-fblaslapack=1 \
    --download-metis=1 \
    --download-mpich=1 \
    --download-mumps=1 \
    --download-parmetis=1 \
    --download-scalapack=1 \
    --download-zlib=1 \
    --with-cuda=1 \
    --with-cuda-dir="${CUDA_DIR}" \
    --with-cudac="${CUDA_DIR}/bin/nvcc"

# Compile and install
make
make install
