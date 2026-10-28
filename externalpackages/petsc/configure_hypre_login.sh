#!/bin/bash
# Configure PETSc 3.22 + Hypre on a login node (where internet is available).
# Run this BEFORE submitting run_gpu_hypre_a100.pbs.
#
# IMPORTANT: This script purges all modules and reloads ONLY what the
# compute node uses (openmpi/4.1.3 + cuda/12.3.2), so configure picks
# up gfortran via mpif90 and links against openmpi 4.1.3 — matching the
# PBS environment exactly.  Do NOT run with pre-loaded Intel modules.
#
# Usage:
#   export ISSM_DIR=/g/data/au88/jh7060/ISSM
#   bash $ISSM_DIR/externalpackages/petsc/configure_hypre_login.sh
#
# CUDA_ARCH (default 80 = A100) sets the GPU architecture for PETSc and for hypre's
# own build (--with-hypre-gpu-arch; --with-cuda-gencode does not reach it).
# Use 90 for H100 / Grace Hopper.
#
# On success: "Configure stage complete".  Then:
#   qsub $ISSM_DIR/test/run_gpu_hypre_a100.pbs

set -eu

# Save ISSM_DIR before module purge can clear env vars
_ISSM_DIR="${ISSM_DIR:-/g/data/au88/jh7060/ISSM}"
CUDA_ARCH="${CUDA_ARCH:-80}"

# ── enforce a clean module environment matching the PBS job ────────────────
# module purge may unset variables (including PETSC_DIR, ISSM_DIR) set via
# modules — so we run it first, then re-derive all our variables.
# shellcheck disable=SC1091
if command -v module &>/dev/null; then
    module purge 2>/dev/null || true
    module load openmpi/4.1.3 2>/dev/null || { echo "ERROR: module load openmpi/4.1.3 failed"; exit 1; }
    module load cuda/12.3.2   2>/dev/null || { echo "ERROR: module load cuda/12.3.2 failed";   exit 1; }
else
    echo "WARNING: module command not found — make sure openmpi/4.1.3 and cuda/12.3.2 are in PATH"
fi

# Re-derive all paths AFTER module purge (which may have cleared them)
export ISSM_DIR="${_ISSM_DIR}"
PETSC_DIR="${ISSM_DIR}/externalpackages/petsc/src-gpu-hypre"
PREFIX="${ISSM_DIR}/externalpackages/petsc/install-gpu-hypre"

# Paranoia: strip any Intel compiler overrides that might shadow mpif90 → gfortran
unset OMPI_FC OMPI_F77 OMPI_F90
unset I_MPI_FC I_MPI_F90 I_MPI_F77
unset F77 F90 FC VT_FC VT_F77 VT_F90 VT_F
unset MPIFC MPIF77 MPIF90

# Verify the compiler chain we expect
CUDA_DIR="${CUDA_HOME:-${CUDA_ROOT:-/usr/local/cuda}}"
if [ ! -d "${CUDA_DIR}" ]; then
    echo "ERROR: CUDA dir not found at ${CUDA_DIR} — check: module load cuda/12.3.2"
    exit 1
fi
for w in mpicc mpicxx mpif90; do
    if ! command -v ${w} &>/dev/null; then
        echo "ERROR: ${w} not found — check: module load openmpi/4.1.3"
        exit 1
    fi
done

# Confirm gfortran is driving mpif90 (not ifort)
FC_IMPL=$(mpif90 --showme:version 2>&1 || mpif90 -show 2>/dev/null | head -1)
echo "mpif90 backend: ${FC_IMPL}"
if echo "${FC_IMPL}" | grep -qi "ifort"; then
    echo "ERROR: mpif90 is wrapping ifort — cannot proceed. Unset OMPI_FC and Intel modules."
    exit 1
fi
MPI_VER=$(mpicc --showme:version 2>&1 || true)
echo "MPI version:    ${MPI_VER}"

if [ ! -f "${PETSC_DIR}/configure" ]; then
    echo "ERROR: ${PETSC_DIR}/configure not found."
    echo "  Extract the source first:"
    echo "    cd ${ISSM_DIR}/externalpackages/petsc"
    echo "    tar -zxf petsc-3.22.3.tar.gz"
    echo "    mkdir -p src-gpu-hypre && mv petsc-3.22.3/* src-gpu-hypre/"
    exit 1
fi

# Wipe any previous tainted build artifacts (keep downloaded package tarballs)
# Also wipe git.hypre build dir so --with-hypre-gpu-arch change takes effect
echo "Cleaning previous build artifacts..."
rm -rf "${PETSC_DIR}/arch-linux-c-opt/obj" \
       "${PETSC_DIR}/arch-linux-c-opt/include" \
       "${PETSC_DIR}/arch-linux-c-opt/lib" \
       "${PETSC_DIR}/arch-linux-c-opt/externalpackages/git.hypre" \
       "${PREFIX}"

echo "========================================"
echo " PETSc 3.22 + Hypre: configure"
echo " PETSC_DIR: ${PETSC_DIR}"
echo " PREFIX:    ${PREFIX}"
echo " CUDA_DIR:  ${CUDA_DIR}"
echo " date:      $(date)"
echo "========================================"

if [ -z ${LDFLAGS+x} ]; then LDFLAGS=""; fi

cd "${PETSC_DIR}"
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

echo "========================================"
echo " Configure stage complete: $(date)"
echo " mpif90 backend confirmed: ${FC_IMPL}"
echo " Now submit:  qsub ${ISSM_DIR}/test/run_gpu_hypre_a100.pbs"
echo "========================================"
