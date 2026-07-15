from collections import OrderedDict
from pairoptions import pairoptions


def gpuoptions(*args):
    """
    GPUOPTIONS - PETSc options for GPU-accelerated (CUDA) solves.

    Uses aijcusparse matrices and CUDA vectors so that KSPSolve runs on GPU.
    Requires a CUDA-enabled PETSc build (see
    externalpackages/petsc/install-3.22-gadi-gpu.sh).

    Defaults to GMRES preconditioned with smoothed-aggregation algebraic
    multigrid (GAMG), whose convergence rate is nearly independent of mesh size
    -- essential for the iterative GPU solve to scale, since the simpler
    block-Jacobi/ILU preconditioner degrades badly under mesh refinement.  The
    matrix block size and the host-side coarse-grid products below are required
    for GAMG to converge and to fit in GPU memory on large meshes.

    Usage:
        options = gpuoptions()
        options = gpuoptions('ksp_type', 'cg', 'pc_type', 'bjacobi')
    """

    # Retrieve any overrides passed by the caller
    options = pairoptions(*args)
    gpu = OrderedDict()

    gpu['toolkit']   = 'petsc'
    gpu['vec_type']  = options.getfieldvalue('vec_type',  'cuda')
    gpu['mat_type']  = options.getfieldvalue('mat_type',  'aijcusparse')
    gpu['ksp_type']  = options.getfieldvalue('ksp_type',  'gmres')
    gpu['pc_type']   = options.getfieldvalue('pc_type',   'gamg')
    # SSA stress-balance has 2 coupled DOFs per node (vx, vy).  Telling PETSc the
    # block size lets GAMG aggregate by node rather than by scalar DOF, which is
    # essential for the multigrid hierarchy to converge on large vector-elliptic
    # systems.  Read by MatSetFromOptions when the stiffness matrix is created.
    gpu['mat_block_size'] = options.getfieldvalue('mat_block_size', 2)
    # GAMG builds each coarse level with a sparse matrix-matrix product (PtAP).
    # CUSPARSE's SpGEMM needs huge scratch buffers and runs the GPU out of memory
    # on large meshes (cudaErrorMemoryAllocation in MatProductSymbolic).  Compute
    # these products on the host (96 GB) instead; the resulting coarse operators
    # are small and the Krylov iteration itself stays on the GPU.
    #
    # GAMG calls MatPtAP() directly (api_user=true), so PETSc reads the
    # api-specific option name '-matptap_backend_cpu' -- NOT the generic
    # '-mat_product_algorithm_backend_cpu', which is only read on the MatProduct
    # interface path.  PtAP decomposes into inner A*B products, so cover those too.
    gpu['matptap_backend_cpu']    = options.getfieldvalue('matptap_backend_cpu', 'true')
    gpu['matmatmult_backend_cpu'] = options.getfieldvalue('matmatmult_backend_cpu', 'true')
    gpu['mat_product_algorithm_backend_cpu'] = options.getfieldvalue('mat_product_algorithm_backend_cpu', 'true')
    # Stop coarsening when grid has <= 2000 equations (direct solve on coarse grid)
    gpu['pc_gamg_coarse_eq_limit'] = options.getfieldvalue('pc_gamg_coarse_eq_limit', 2000)
    # Aggressive dropping of weak connections improves coarse grid quality
    gpu['pc_gamg_threshold']  = options.getfieldvalue('pc_gamg_threshold',  0.08)
    # Prolongator smoothing.  Default to smoothed aggregation (nsmooths=1), which
    # is the stronger preconditioner and works on single-GPU.  In PARALLEL the
    # multi-GPU driver overrides this to '0' (unsmoothed) for robustness.
    # NB: values are STRINGS -- toolkits.py marshals with `if not optionvalue:`,
    # so an integer 0 would be written as a bare valueless flag (PETSc keeps its
    # default of 1).  A non-empty string serialises correctly.
    gpu['pc_gamg_agg_nsmooths'] = options.getfieldvalue('pc_gamg_agg_nsmooths', '1')
    # ksp_rtol is the relative reduction of the *preconditioned* residual at which
    # GMRES stops.  Empirically a 1e-10 reduction leaves the *true* residual
    # ||KU-F||/||F|| at ~1.5e-6 on the 1.6M-DOF case -- just above ISSM's 1e-6
    # acceptance gate (md.settings.solver_residue_threshold).  1e-12 drives the
    # true residual to ~1.5e-8, clearing the gate with margin.  (Confirmed not an
    # iteration/restart limit: result was bit-identical across max_it and restart
    # changes, i.e. GMRES converges to rtol well before exhausting iterations.)
    gpu['ksp_rtol']  = options.getfieldvalue('ksp_rtol',  1e-12)
    gpu['ksp_max_it'] = options.getfieldvalue('ksp_max_it', 5000)
    # Longer restart than the GMRES(30) default, so the extra iterations needed to
    # reach the tighter rtol on large systems don't stagnate.
    gpu['ksp_gmres_restart'] = options.getfieldvalue('ksp_gmres_restart', 100)

    return gpu
