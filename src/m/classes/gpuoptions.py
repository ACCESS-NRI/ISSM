from collections import OrderedDict
from pairoptions import pairoptions


def gpuoptions(*args):
    """
    GPUOPTIONS - PETSc options for GPU-accelerated (CUDA) solves.

    Uses aijcusparse matrices and CUDA vectors so that KSPSolve runs on GPU.
    Requires a CUDA-enabled PETSc build (see
    externalpackages/petsc/install-3.22-gadi-gpu.sh).

    The preconditioner defaults to block-Jacobi with ILU sub-solves, which
    PETSc can run on GPU via CUDA.  For larger problems consider switching to
    'gamg' (algebraic multigrid) which also has CUDA support.

    Usage:
        options = gpuoptions()
        options = gpuoptions('ksp_type', 'cg', 'pc_type', 'gamg')
    """

    # Retrieve any overrides passed by the caller
    options = pairoptions(*args)
    gpu = OrderedDict()

    gpu['toolkit']   = 'petsc'
    gpu['vec_type']  = options.getfieldvalue('vec_type',  'cuda')
    gpu['mat_type']  = options.getfieldvalue('mat_type',  'aijcusparse')
    gpu['ksp_type']  = options.getfieldvalue('ksp_type',  'gmres')
    gpu['pc_type']   = options.getfieldvalue('pc_type',   'bjacobi')
    gpu['sub_pc_type'] = options.getfieldvalue('sub_pc_type', 'ilu')
    gpu['ksp_rtol']  = options.getfieldvalue('ksp_rtol',  1e-10)
    gpu['ksp_max_it'] = options.getfieldvalue('ksp_max_it', 500)

    return gpu
