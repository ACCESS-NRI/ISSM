# ISSM GPU Solver — Change Notes

Work to make ISSM run its linear solves on NVIDIA V100 GPUs (NCI Gadi
`gpuvolta` queue) via a CUDA-enabled PETSc. Chronological summary of what was
changed and why.

## Files changed / added

| File | Change |
|------|--------|
| `externalpackages/petsc/install-3.22-gadi-gpu.sh` | Build PETSc against **system OpenMPI** instead of downloading MPICH |
| `configure.sh` (gpu step) | Added `--with-python=...` so Python wrappers build |
| `src/m/classes/gpuoptions.py` | New solver toolkit: CUDA vec/mat + GMRES/GAMG, with block-size and host-PtAP tuning |
| `src/c/Makefile.am` | Added missing `AndersonAccelerator.cpp` to the build |
| `test/gpu_kspsolve_test.py` | Resolution-sweep CPU-vs-GPU benchmark |
| `test/run_gpu_kspsolve_test.pbs` | PBS job script for the gpuvolta queue |
| `test/gpu_solver_report.tex` | LaTeX write-up of the work |

## Problems hit and how they were fixed (in order)

1. **`ModuleNotFoundError: IssmConfig_python`** — GPU `configure` was missing
   `--with-python`, so the compiled Python wrappers (`*.so`) were never built.
   Added `--with-python=/apps/intel-python3/.../python3`. Wrappers install to
   `$ISSM_DIR/lib`, which must be on `PYTHONPATH`.

2. **`libk5crypto.so.3: undefined symbol EVP_KDF_ctrl`** — Intel Python ships an
   old `libcrypto` lacking `EVP_KDF_ctrl`. Fixed at *runtime* in the PBS script
   by `LD_PRELOAD`-ing system OpenSSL 1.1.1k (`/lib64/libcrypto.so.1.1`,
   `libssl.so.1.1`). Not a build change.

3. **Tiny problem, GPU slower than CPU** — expected. At 50 km (680 DOFs) MUMPS
   is sub-second; GPU overhead (CUDA init, CUSPARSE assembly, H2D transfer)
   dominates. Built the resolution-sweep benchmark to find the crossover.

4. **GMRES + bjacobi/ILU stalls at 5 km+** — ILU preconditioning quality decays
   with mesh refinement. Switched default `pc_type` from `bjacobi` to **`gamg`**
   (algebraic multigrid).

5. **`np=4` multi-rank → 4 independent rank-0 jobs** — root cause: `issm.exe`
   was linked against **two MPI runtimes at once** (`libmpi.so.12` MPICH from
   GPU PETSc + `libmpi.so.40` OpenMPI from `mpicxx`). Harmless at `np=1`, breaks
   any `np>1`. **Fixed by rebuilding GPU PETSc against system OpenMPI**: removed
   `--download-mpich=1`, added `--with-cc=mpicc --with-cxx=mpicxx
   --with-fc=mpif90` + an MPI sanity check in the install script. Requires
   `rm -rf install-gpu src-gpu` before re-running `source configure.sh gpu`.

6. **Link error `undefined reference to AndersonAccelerator::...`** — after the
   clean rebuild, `AndersonAccelerator.cpp` was missing from `src/c/Makefile.am`
   (the file existed and was referenced by `solutionsequence_nonlinear.cpp` but
   never compiled). Added it to the source list.

7. **GAMG won't converge at 2 km / 1 km** — GAMG treated the SSA system as a
   scalar matrix. SSA has **2 coupled DOFs/node (vx, vy)**. Added
   `mat_block_size = 2` so GAMG aggregates by node. Fixed convergence for all
   resolutions that fit in memory.

8. **`np=4` on a single GPU is ~4× slower** — 4 ranks contend for one V100 and
   add MPI overhead inside GAMG with no parallelism gain. Reverted benchmark to
   `NP = 1`. (Multi-rank only helps with one GPU per rank; the MPI fix is still
   needed for that future case.)

9. **2 km / 1 km → `cudaErrorMemoryAllocation` in `MatProductSymbolic`** —
   GAMG's coarse-grid Galerkin product (`MatPtAP`) uses CUSPARSE SpGEMM, whose
   scratch buffers overflow 32 GB at 397k / 1.59M DOFs. Fix is to run the matrix
   products on the host (96 GB) while the Krylov iteration stays on the GPU —
   but **the option name matters**: GAMG calls `MatPtAP()` directly
   (`api_user=true`), so PETSc reads `-matptap_backend_cpu`, *not* the generic
   `-mat_product_algorithm_backend_cpu` (which it silently ignored —
   `Option left: ... value: true`). Set `matptap_backend_cpu = true` plus
   `matmatmult_backend_cpu = true` (for the inner A·B products of the decomposed
   PtAP). This cleared the OOM: **2 km then passed at 2.2×**.

10. **1 km → true residual stuck at 1.51e-6 (> 1e-6 gate)** — *not* OOM and *not*
    iteration-starved: the residual was **bit-identical** across `ksp_max_it`
    2000→5000 and `ksp_gmres_restart` 30→100, proving GMRES was converging to its
    relative tolerance, not exhausting iterations. `ksp_rtol = 1e-10` (relative
    reduction of the *preconditioned* residual) left the *true* residual
    `||KU-F||/||F||` at 1.5e-6. Tightened `ksp_rtol = 1e-12` → true residual
    ~1.5e-8. **1 km then passed at 3.2×.** (Longer `ksp_gmres_restart = 100` kept
    so the extra iterations to reach the tighter tol don't stagnate.)

## Multi-GPU (one rank per GPU) — the parallel saga

With single-GPU complete, the next goal was multi-GPU scaling (1 MPI rank per
V100 on a 4-GPU gpuvolta node). PETSc maps rank `r` -> GPU `r % ndev`
automatically (`PETSC_DECIDE`, verified in `cupmdevice.cxx:259`), so no binding
wrapper is needed. Getting it to actually converge took a four-step hunt:

11. **CPU MUMPS scaled, every GPU solve failed.** First confirmation the OpenMPI
    rebuild works at np>1: 4-rank MUMPS sped up (1km: 865→395 s). GPU solves
    crashed in `PCGAMGOptimizeProlongator_AGG`: *"Computed maximum singular value
    as zero."* (smoothed-aggregation eigenvalue estimate).

12. **`pc_gamg_agg_nsmooths` silently ignored.** Setting it to integer `0` hit the
    `toolkits.py` marshaller bug (`if not optionvalue:` writes a *valueless flag*
    for falsy values), so PETSc kept its default. **Pass option values as
    STRINGS** (`'0'`) to survive marshalling. (Latent foot-gun for any option set
    to integer 0.)

13. **Unsmoothed aggregation then *diverged*** (residual 3–7, worse than the RHS).
    A `-use_gpu_aware_mpi 0` guess changed nothing (bit-identical residuals), so
    that was a red herring — ruling it out mattered.

14. **Isolation diagnostic pinned the real cause** (`gpu_mgpu_diagnose.py`: one
    small mesh, np=4, several PCs back-to-back in one job). Results:
    `bjacobi+ilu` CONVERGED (so the **parallel GPU matvec is correct**); GAMG with
    **host-PtAP FAILED** but the *identical* GAMG with **GPU-PtAP CONVERGED
    (1e-11)**. So **`matptap_backend_cpu=true` builds a corrupt coarse operator in
    parallel** — the very flag that fixed the single-GPU OOM. The rule:
    **host-PtAP for single-GPU (memory); GPU-PtAP for multi-GPU (correctness, and
    it fits since each rank holds ~1/N of the matrix).** Encoded as an NP-based
    switch in `gpu_kspsolve_test.py`.

### Multi-GPU benchmark (4× V100, unsmoothed aggregation + GPU-PtAP) — all PASS

| Resolution | DOFs | 4-GPU (s) | 1-GPU (s) | 1→4 GPU scaling |
|-----------:|-----:|----------:|----------:|:---:|
|  5 km |    64,362 |   9.58 | 10.42 | 1.09× |
|  2 km |   396,752 |  40.61 | 63.30 | 1.56× |
|  1 km | 1,585,078 | 149.85 | 289.49 | **1.93×** |

Scales sublinearly (~48% efficiency at 1.6M DOFs, growing with size — small
problems are communication-bound). NB the 1→4 comparison mixes a preconditioner
change (multi-GPU used unsmoothed, single-GPU smoothed); a clean strong-scaling
study needs the same config at np=1,2,4.

## A100 benchmark and extended investigation

### 15. Three-column benchmark on A100 (July 2026)

Submitted the same three-solver sweep to the `dgxa100` queue (single A100 80GB,
`ISSM_GPU_PTAP=1`). With 80GB HBM2e, the CUSPARSE SpGEMM scratch for
`MatPtAP` fits without OOM (unlike 32GB V100), so GPU-PtAP can run on the device.

Result: GPU-PtAP converges on all resolutions up to 1km — **but only with
`pc_gamg_agg_nsmooths=0`** (unsmoothed aggregation). Smoothed aggregation
(`nsmooths=1`) diverges: `solver residue too high: 5.38 > 1e-6`. Root cause
identified later (see #20).

| Resolution | DOFs | MUMPS (s) | CPU-GMRES (s) | GPU-GMRES (s) | vs MUMPS | vs CPU-GMRES |
|-----------:|-----:|----------:|--------------:|--------------:|---------:|-------------:|
|  50 km |       680 |   2.638 |   2.171 |   2.697 | 0.98× | 0.80× |
|  20 km |     4,172 |   2.525 |   2.296 |   3.350 | 0.75× | 0.69× |
|  10 km |    16,304 |   4.741 |   4.494 |   4.883 | 0.97× | 0.92× |
|   5 km |    64,362 |  16.055 |  12.815 |  12.397 | 1.30× | 1.03× |
|   2 km |   396,752 | 135.062 |  81.540 |  67.680 | 2.00× | 1.20× |
|   1 km | 1,585,078 | 821.372 | 374.803 | 267.222 | **3.07×** | **1.40×** |

Hardware speedup improved to **1.40×** (vs 1.35× on V100) because GPU-PtAP
eliminates the H2D transfer tax that dominated on V100 (confirmed by nsys: only
~2% of V100 walltime was on GPU with host-PtAP). PBS reported 16% GPU
utilisation with GPU-PtAP vs ~2% with host-PtAP.

### 16. Solver variants sweep (Chebyshev, BiCGStab, pipelining)

Tested five options combinations on A100 NP=1 at 2km and 1km
(`gpu_speedup_variants_test.py`):

| Variant | 1km (s) | vs MUMPS | vs CPU-GMRES |
|---------|---------|---------|-------------|
| baseline (GMRES, SOR) | 300.4 | 2.88× | 1.37× |
| Chebyshev+Jacobi smoother | 303.6 | 2.84× | 1.35× |
| BiCGStab | 296.5 | 2.91× | 1.38× |
| BiCGStab + Chebyshev | **294.2** | **2.94×** | **1.40×** |
| nsmooths=1 + CG eigenestimator | FAIL | — | — |

All within ~6% of each other. No variant breaks through the ~1.40× hardware
ceiling. The CG-eigenestimator fix did not rescue nsmooths=1 (root cause is
different — see #20). **cheb+bcgs is marginally best at 1km** but the gain is
too small to justify changing the default.

### 17. Multi-GPU A100 (NP=4)

Ran the benchmark with NP=4 (one rank per A100) on the dgxa100 node.
`NP > 1` triggers GPU-PtAP + nsmooths=0 automatically in `gpu_kspsolve_test.py`.

| Config | 1km GPU-GMRES | 1km CPU-GMRES | vs CPU-GMRES | vs MUMPS |
|--------|:---:|:---:|:---:|:---:|
| A100 NP=1 | 267s | 375s | **1.40×** | 3.07× |
| A100 NP=4 | 182s | 177s | **0.97×** | 2.14× |

Going from NP=1 to NP=4: CPU-GMRES scales 2.1× (375→177s) but GPU-GMRES only
scales 1.5× (267→182s). At NP=4 the GPU barely matches CPU. Reason: with 4
ranks each holding ~400k DOFs, the GMRES global AllReduce (one per iteration)
and halo exchanges dominate. GPU is ~64% utilised but spending most of its time
waiting on MPI rather than computing. PBS: GPU utilisation 64%, GPU memory only
6.35 GB (1/4 of the single-GPU problem).

**Single-GPU is the sweet spot** for this problem size. Multi-GPU is
communication-limited, not compute-limited.

### 18. GPU-aware MPI (NP=4)

OpenMPI 4.1.3 on Gadi is compiled with `--with-cuda` (`smcuda` BTL,
`MCA coll: cuda`), so GPU-aware MPI can be enabled by removing
`-use_gpu_aware_mpi 0`. Tried this on NP=4 to eliminate the
GPU→CPU→MPI→CPU→GPU round-trip per iteration.

Result: **GPU GMRES fails at 2km and 1km** (residual 42026 at 2km, 1.43 at 1km).
Small resolutions pass (GPU-PtAP not triggered). The GPU-aware MPI transfers
corrupt the parallel coarse operator assembly during GPU-PtAP. UCX's CUDA IPC
path cannot safely share device pointers across CUDA contexts on different GPUs.

**GPU-aware MPI is incompatible with GPU-PtAP at NP>1.** Keep
`-use_gpu_aware_mpi 0` for all multi-GPU runs.

### 19. Pipelined Krylov (NP=4, GPU-aware MPI off)

Pipelined variants (`pipefgmres`, `pipebcgs`) restructure the Krylov recurrence
to overlap AllReduce(k) with SpMV(k+1), hiding the MPI latency.

| Variant | 2km (s) | 1km (s) | vs CPU-GMRES |
|---------|--------|--------|-------------|
| gmres (baseline) | 52.2 | 183.8 | 2.29× |
| bcgs | 52.9 | 184.0 | 2.29× |
| pipefgmres | 56.6 | 193.0 | 2.18× |
| pipebcgs | 54.6 | 191.6 | 2.20× |

Pipelining is **slightly slower** than standard GMRES at NP=4. On a single DGX
node the 4 A100s are connected by NVLink (~600 GB/s), so the AllReduce latency
is already negligible. The pipeline overhead (extra synchronisation checkpoints
in the recurrence) outweighs any latency hiding. Pipelining would only help for
inter-node communication (high-latency Infiniband), not intra-node NVLink.

Note: the CPU-GMRES baseline in this test was inadvertently set to nsmooths=0
(same as GPU), making it appear slower (421s at 1km). The apparent 2.29×
hardware speedup is misleading — it measures GPU vs CPU at the same (weaker)
preconditioner quality, not best-vs-best. The true hardware speedup at NP=4
remains ~0.97× (GPU-nsmooths=0 vs CPU-nsmooths=1).

### 20. nsmooths=1 GPU-PtAP root cause (CUSPARSE_STATUS_INSUFFICIENT_RESOURCES)

Definitive diagnosis via `gpu_nsmooths_diag.py` — five variants at A100 NP=1,
sweeping GPU/CPU assignment for MatMatMult (prolongator smoothing: A×P₀) and
MatPtAP (coarse operator: P^T×A×P) independently.

The PETSc error log:
```
CUSPARSE_STATUS_INSUFFICIENT_RESOURCES (error 11)
MatProductSymbolic_SeqAIJCUSPARSE_SeqAIJCUSPARSE()  ← GPU SpGEMM
MatProductSymbolic_ABC_Basic()
MatPtAP()
PCGAMGCreateLevel_GAMG()
```

**Root cause:** with nsmooths=1, the smoothed prolongator P has more non-zeros
than the tentative P₀ (nsmooths=0). The CUSPARSE SpGEMM symbolic phase for
P^T×A×P then exceeds CUSPARSE's internal hardware resource limit (not device
memory — a workspace/occupancy limit in the SpGEMM kernel). This is a hard
CUSPARSE constraint that cannot be worked around via command-line options.

Isolation results:
- `n1-all-gpu` (GPU-MatMul + GPU-PtAP, n=1): FAIL (CUSPARSE crash → residual 5.38)
- `n1-mul-cpu` (CPU-MatMul + GPU-PtAP, n=1): FAIL (same CUSPARSE crash in MatPtAP)
- `n1-ptap-cpu` (GPU-MatMul + CPU-PtAP, n=1): ~926s at 2km (extremely slow but passes)
- `n1-all-cpu` (CPU-MatMul + CPU-PtAP, n=1): **49s at 2km (PASS, only 7% slower than n0-gpu)**

The fix for n1 is to put **both** MatMatMul and MatPtAP on CPU. At 2km this adds
only 7% overhead vs n0-gpu because the better preconditioner (nsmooths=1) needs
fewer GMRES iterations. This is essentially the same trade-off as the original
V100 host-PtAP approach: CPU setup overhead offset by fewer Krylov iterations.

**The GPU-PtAP improvement and the nsmooths degradation cancel each other out:**
both V100-host-PtAP+nsmooths=1 and A100-GPU-PtAP+nsmooths=0 land near 3.07–3.13×
vs MUMPS at 1km. This is a fundamental result: the CUSPARSE SpGEMM workspace
limit prevents smoothed aggregation from running on GPU with large matrices,
and the preconditioner quality loss from nsmooths=0 offsets the GPU-PtAP gain.

## Hypre BoomerAMG (A100, July 2026)

Following the CUSPARSE SpGEMM ceiling diagnosis (#20), we pursued Hypre BoomerAMG
as an alternative AMG preconditioner — one that runs its own GPU-native AMG setup
entirely on-device, bypassing PETSc's GAMG and CUSPARSE's SpGEMM.

New files:
- `externalpackages/petsc/configure_hypre_login.sh` — login-node configure script
- `externalpackages/petsc/src-gpu-hypre/` — PETSc source tree with Hypre
- `externalpackages/petsc/install-gpu-hypre/` — build prefix
- `test/gpu_hypre_test.py` — 5-solver benchmark (MUMPS / CPU-GAMG / GPU-GAMG / Hypre-CPU / Hypre-GPU)
- `test/run_gpu_hypre_a100.pbs` — PBS job (Phase 1: compile, Phase 2: benchmark)
- `test/gpu_hypre_native_test.py` — isolation test for mat_type=hypre crash

### 21. Hypre SM 80 architecture mismatch

`--with-cuda-gencode="arch=compute_80,code=sm_80"` does **not** propagate to
Hypre's sub-build. Hypre's autoconf configure runs independently inside PETSc's
external-package system and defaults to SM 70 when no GPU is present at configure
time (login node). `cuobjdump -arch libHYPRE.so` showed only `sm_70`.

**Fix:** `--with-hypre-gpu-arch=80` (PETSc's `config/BuildSystem/config/packages/hypre.py`
line 26, `formGNUConfigureArgs()`). This flag passes `--with-gpu-arch=80` directly
to Hypre's own configure, overriding the SM 70 default. Confirmed SM 80 via
`cuobjdump` after rebuild.

**Lesson:** Each external package's GPU arch must be set with its dedicated option,
not the top-level `--with-cuda-gencode`.

### 22. Build guard fires too early

PBS job guard checked `[ -d lib ]` — but `install-gpu-hypre/lib/` exists as soon
as configure downloads packages, long before `libpetsc.so` is compiled. Job always
skipped the build.

**Fix:** changed guard to `[ -f "${PETSC_HYPRE}/lib/libpetsc.so" ]`.

### 23. BoomerAMG scalar divergence (residue 5.37)

With default options, BoomerAMG treated the SSA system as a scalar matrix and
diverged completely (residue 5.37 — worse than the initial guess).

**Fix:** `mat_block_size=2`. PETSc maps this to
`HYPRE_BoomerAMGSetNumFunctions(2)`, telling BoomerAMG the system has 2 DOFs
per node (vx, vy) and to coarsen by node rather than by scalar DOF. This is the
same near-null-space argument as for GAMG (item #7).

**Note:** `pc_hypre_boomeramg_nodal_coarsen=6` (Frobenius-norm coarsening) is a
stronger version of the same correction but triggers SIGSEGV in Hypre 2.31.0 when
combined with `mat_block_size=2` and CUDA. `mat_block_size=2` alone is sufficient.

### 24. VecSetValues "locked for read-only access" (Hypre-CPU)

After GPU-GAMG runs first and initialises CUDA globally, `PCCreate_HYPRE` (even
for CPU AMG with `mat_type=aij`) triggers further CUDA init that locks existing
PETSc vector objects. ISSM then fails in `VecSetValues` when writing the RHS.

**Fix:** `vec_type=cuda` for all Hypre variants, including Hypre-CPU. ISSM's
GPU-aware assembly path is used for vectors; BoomerAMG still runs on CPU (the
matrix remains `aij` in CPU memory). The result is a hybrid: CPU AMG, GPU Krylov.

### 25. ISSM recovery mode blocks benchmark

When `solver_residue_threshold` (default 1e-6) is exceeded, ISSM enters a C++
recovery sub-loop that never raises a Python exception — runs for 5-12+ min per
failed solve, blocking the benchmark loop silently.

**Fix:** `md.settings.solver_residue_threshold = 1e6` to disable the gate
entirely. Failed solves return immediately with a wrong `max_vel` value.

### 26. mat_type=hypre crash (open bug in Hypre 2.31.0)

`mat_type=hypre` (native Hypre ParCSR format on GPU) immediately crashes with:
```
CUDA ERROR (code=1, invalid argument) at
  .../git.hypre/src/utilities/memory.c:777
  cudaMemcpy(dst, src, size, cudaMemcpyDeviceToDevice)
```
`cudaErrorInvalidValue` means src or dst is a stale/uninitialised device pointer.
Hypre does not abort — it prints this in a tight loop for the entire walltime,
generating a 253 MB output file.

An isolation test (`gpu_hypre_native_test.py`) ran `mat_type=hypre` *before* any
`aijcusparse` solve. Same crash. Contamination hypothesis disproved: the bug is
intrinsic to Hypre 2.31.0's PARCSR GPU matrix setup on SM 80.

**Status:** excluded from production. `mat_type=aijcusparse` is the working GPU
Hypre path.

### Hypre benchmark results (A100, NP=1, July 2026) — all PASS

Five-solver sweep extended to 500 m (6.3M DOFs):

| Res | DOFs | MUMPS | CPU-GAMG | GPU-GAMG | Hypre-CPU | Hypre-GPU |
|----:|-----:|------:|---------:|---------:|----------:|----------:|
| 2 km | 396,752 | 157.5 s | 81.8 s (1.93×) | 66.2 s (2.38×) | 59.3 s (2.66×) | 56.9 s (2.77×) |
| 1 km | 1,585,078 | 959.3 s | 372.2 s (2.58×) | 268.8 s (3.57×) | 240.5 s (3.99×) | 230.3 s (4.17×) |
| 500 m | 6,332,712 | 7423.7 s | 1810.3 s (4.10×) | 1228.5 s (6.04×) | 1008.4 s (7.36×) | **973.2 s (7.63×)** |

All solutions correct — `max_vel` matches MUMPS to floating-point at each resolution.

**Scaling exponents (power-law fit, 3 points):**
- MUMPS: α ≈ 1.39 (theory: O(N^{3/2}) for 2D)
- GPU-GAMG: α ≈ 1.05 (theory: O(N))
- Hypre-GPU: α ≈ 1.03 (theory: O(N))

The speedup over MUMPS grows as N^{0.36}: a 4× increase in DOFs multiplies the
advantage by ~1.7×. At 30M DOFs (Antarctic scale), the projected speedup is ~14×.

## Continental AIS3 SSA (A100, October 2026)

ACCESS-AIS3's Antarctic SSA model: 1.59M vertices, 2,948,476 unknowns, Budd friction. Baseline is
what production runs: MUMPS on 48 CPU cores (one Gadi `normal` node). Note that the SquareShelf
tables above compare against MUMPS on **one** core. GPU runs use one A100 GPU + 16 cores (1/8 of a
DGX node, `ngpus=1 ncpus=16`).

### 27. Recovery mode was sticky (`Solverx.cpp`)

When a solve misses `md.settings.solver_residue_threshold`, `Solverx` re-solves with the
`RecoveryAnalysis` toolkit (MUMPS) but never switches back: the options are only reset by the next
`SetCurrentConfiguration`, so every later nonlinear iteration of that core ran on MUMPS. One hard
solve pinned the whole stress balance to the direct solver, and the run's timings were MUMPS's.

**Fix:** after a successful recovery solve, re-apply the current analysis' options
(`AnalysisTypeEnum`); prints `Recovery Mode succeeded, back to <analysis> toolkit options`.

### 28. Recovery solve on GPU objects returned zero (residue exactly 1)

The direct solve's solution vector is created by `NewVec` → `VecSetFromOptions`, i.e. from the
*recovery* options. With no `vec_type` there it was a host vector next to the CUSPARSE matrix and
CUDA right-hand side; the residual check then read 1 and ISSM aborted.

**Fix:** give `RecoveryAnalysis` `vec_type=cuda` (MUMPS accepts `mpiaijcusparse`: PETSc matches
solver types by type-name prefix).

### 29. First Picard matrix near-singular from `md.initialization = vel_obs`

Every AMG setting (GAMG, GAMG + rigid-body modes, BoomerAMG) broke down on Picard iteration 1 at
residual ~1e-2, while all later iterations converged. Dumping the matrix (`-ksp_view_mat binary:`)
showed its diagonal spanning 16.8 orders of magnitude: the observed velocity's gap fills are
piecewise constant, so ~9.5k vertices touch elements with strain rate ~1e-25 /s and Glen viscosity
blows up. Even MUMPS only reached residue 1e-6 there (later iterations 1e-10–1e-15). With #27 and
#28, BoomerAMG fails over to MUMPS for those first iterations and returns. Cap `ksp_max_it` (300) so
the hand-over is quick.

### 30. Several ranks on one GPU need CUDA MPS

Assembly stays on CPU ranks, so one rank per GPU is slow (724 s vs 137 s end to end). Without MPS,
ranks sharing the GPU time-slice it: BoomerAMG setup + solve went from 3.2 s (1 rank) to 25 s
(16 ranks). Starting `nvidia-cuda-mps-control -d` in the job (pipe/log dirs in `$PBS_JOBFS`) keeps
it at 4.2 s with 16 ranks.

### 31. GPU architecture as a variable; `gpuhypreoptions`

`install-3.22-gadi-gpu-hypre.sh` and `configure_hypre_login.sh` take `CUDA_ARCH` (default 80) for
both `--with-cuda-gencode` and `--with-hypre-gpu-arch` (#21); use 90 for H100 / Grace Hopper.
`src/m/solvers/gpuhypreoptions.{py,m}` returns the options below, for
`md.toolkits.StressbalanceAnalysis` only (`mat_block_size=2` breaks scalar systems).

### Continental results (one SSA stress balance, 20 Picard iterations)

| | MUMPS, 48 cores | BoomerAMG, 1 A100 GPU + 16 ranks + MPS |
|---|---:|---:|
| linear solves | 125 s | 66 s |
| assembly etc. | 12 s | 43 s |
| total core | 137 s | **109 s** |
| SU | 4.21 | **2.84** |

Velocity vs MUMPS: relative difference 3.1e-7 (99th percentile), max 0.17 m/yr. Includes two
MUMPS fallbacks (#29). On the converged-viscosity matrix alone (outside ISSM): BoomerAMG on one A100
3.1 s (2.3 s setup + 0.8 s solve, 110 its) vs MUMPS on 16 cores 15.7 s.

GPU options: `vec_type cuda, mat_type aijcusparse, mat_block_size 2, ksp_type gmres,
ksp_pc_side right, ksp_rtol 1e-8, ksp_max_it 300, ksp_gmres_restart 100, pc_type hypre,
pc_hypre_type boomeramg, pc_hypre_boomeramg_strong_threshold 0.25`; recovery and other analyses
MUMPS with `vec_type cuda`.

## Final `gpuoptions` defaults

```
toolkit              = petsc
vec_type             = cuda
mat_type             = aijcusparse
ksp_type             = gmres
pc_type              = gamg
mat_block_size       = 2        # SSA: 2 DOFs/node (GAMG aggregation)
matptap_backend_cpu  = true     # GAMG PtAP on host (avoid GPU OOM); api_user name
matmatmult_backend_cpu = true   # inner A*B of decomposed PtAP, also on host
mat_product_algorithm_backend_cpu = true  # belt-and-braces (MatProduct path)
pc_gamg_coarse_eq_limit = 2000
pc_gamg_threshold    = 0.08
ksp_rtol             = 1e-12     # tight enough that TRUE residual clears 1e-6 gate
ksp_max_it           = 5000
ksp_gmres_restart    = 100
```

These defaults target **single-GPU**. For **multi-GPU** the driver overrides PtAP
to run on the GPU (host-PtAP is broken in parallel):
`matptap_backend_cpu=false`, `matmatmult_backend_cpu=false`,
`mat_product_algorithm_backend_cpu=false`. Pass all option values as **strings**
(toolkits.py marshalling drops falsy/zero values written as bare flags).

## Benchmark status (single V100, SquareShelf SSA, np=1) — COMPLETE, all PASS

Three-solver sweep isolates the hardware and algorithm contributions separately:

| Resolution | DOFs | CPU MUMPS | CPU GMRES+GAMG | GPU GMRES+GAMG | vs MUMPS | vs CPU-GMRES |
|-----------:|-----:|----------:|---------------:|---------------:|---------:|-------------:|
| 50 km |       680 |   1.44 s |   0.84 s |   2.25 s | 0.64× | 0.37× |
| 20 km |     4,172 |   1.14 s |   1.21 s |   1.85 s | 0.62× | 0.65× |
| 10 km |    16,304 |   3.24 s |   3.04 s |   3.16 s | 1.03× | 0.96× |
|  5 km |    64,362 |  14.14 s |  12.54 s |   9.99 s | 1.42× | 1.25× |
|  2 km |   396,752 | 132.36 s |  78.16 s |  59.87 s | 2.21× | 1.31× |
|  1 km | 1,585,078 | 854.94 s | 368.82 s | 272.73 s | **3.13×** | **1.35×** |

The CPU GMRES+GAMG column uses identical GAMG tuning (`mat_block_size=2`,
`ksp_rtol=1e-12`, etc.) with only `vec_type=standard` and `mat_type=aij` (no
CUDA) — so the two rightmost columns change exactly one variable: the device.

**Reading the numbers:**
- **Algorithm effect (MUMPS → GMRES+GAMG on CPU):** ~2.3× at 1.6M DOFs. The
  iterative O(N) solver beats direct O(N^{3/2}) on the CPU alone. This is the
  dominant effect and grows with N.
- **Hardware effect (CPU → GPU, same GMRES+GAMG):** ~1.35× at 1.6M DOFs. The
  true GPU contribution — and it will grow further with the next resolution level
  as the bandwidth advantage of HBM2 over DRAM compounds.
- **Combined (MUMPS → GPU-GMRES, what a user experiences):** 3.13×.

Correctness: GPU vs CPU-MUMPS `Vx` relative diff ≤ ~4e-11 throughout.

## Rebuild / run cheatsheet

```bash
# interactive GPU node
qsub -I -q gpuvolta -l ncpus=12,ngpus=1,mem=96GB,walltime=2:00:00,storage=gdata/au88+scratch/ek9
module load cuda/12.3.2
export ISSM_DIR=/g/data/au88/jh7060/ISSM   # /g/data1b no longer exists (Oct 2026)
cd $ISSM_DIR

# only needed if PETSc itself must be rebuilt (e.g. MPI/CUDA change):
rm -rf externalpackages/petsc/install-gpu externalpackages/petsc/src-gpu
source configure.sh gpu

# sanity check: must show ONLY libmpi.so.40 (OpenMPI), no libmpi.so.12
ldd $ISSM_DIR/bin/issm.exe | grep -i mpi

# run the benchmark (from a login node)
qsub $ISSM_DIR/test/run_gpu_kspsolve_test.pbs
```

Note: edits to `gpuoptions.py` or the test script are **Python-only** — no
rebuild needed, just resubmit the PBS job. Only C/PETSc/Makefile.am changes
require `source configure.sh gpu`.
