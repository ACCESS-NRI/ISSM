from pairoptions import pairoptions
from collections import OrderedDict


def gpuhypreoptions(*args):
    """
    GPUHYPREOPTIONS - GMRES + hypre BoomerAMG with CUDA vectors and matrices.

    For the stress balance on a CUDA + hypre PETSc build
    (externalpackages/petsc/install-3.22-gadi-gpu-hypre.sh). Set it on that
    analysis only, e.g.

        md.toolkits.StressbalanceAnalysis = gpuhypreoptions()

    mat_block_size=2 groups each node's two velocity components (SSA/HO), which
    BoomerAMG needs to converge; it breaks scalar systems (mass transport, level
    sets), so keep those on their defaults. Give RecoveryAnalysis vec_type=cuda
    too, or its direct solve on the GPU matrices returns zeros. Right
    preconditioning makes ksp_rtol bound the true residual that
    md.settings.solver_residue_threshold checks.

    Values are strings: an integer 0 would be written as a valueless flag.
    """

    #retrieve options provided in *args
    options = pairoptions(*args)
    solverOptions = OrderedDict()
    solverOptions['toolkit'] = 'petsc'
    solverOptions['vec_type'] = options.getfieldvalue('vec_type', 'cuda')
    solverOptions['mat_type'] = options.getfieldvalue('mat_type', 'aijcusparse')
    solverOptions['mat_block_size'] = options.getfieldvalue('mat_block_size', '2')
    solverOptions['ksp_type'] = options.getfieldvalue('ksp_type', 'gmres')
    solverOptions['ksp_pc_side'] = options.getfieldvalue('ksp_pc_side', 'right')
    solverOptions['ksp_rtol'] = options.getfieldvalue('ksp_rtol', '1e-8')
    solverOptions['ksp_max_it'] = options.getfieldvalue('ksp_max_it', '300')
    solverOptions['ksp_gmres_restart'] = options.getfieldvalue('ksp_gmres_restart', '100')
    solverOptions['pc_type'] = options.getfieldvalue('pc_type', 'hypre')
    solverOptions['pc_hypre_type'] = options.getfieldvalue('pc_hypre_type', 'boomeramg')
    solverOptions['pc_hypre_boomeramg_strong_threshold'] = options.getfieldvalue('pc_hypre_boomeramg_strong_threshold', '0.25')
    return solverOptions
