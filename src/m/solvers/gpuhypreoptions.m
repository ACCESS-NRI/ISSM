function solverOptions=gpuhypreoptions(varargin)
%GPUHYPREOPTIONS - GMRES + hypre BoomerAMG with CUDA vectors and matrices
%
%   For the stress balance on a CUDA + hypre PETSc build
%   (externalpackages/petsc/install-3.22-gadi-gpu-hypre.sh). Set it on that
%   analysis only:
%      md.toolkits.StressbalanceAnalysis=gpuhypreoptions();
%
%   mat_block_size=2 groups each node's two velocity components (SSA/HO), which
%   BoomerAMG needs to converge; it breaks scalar systems (mass transport, level
%   sets), so keep those on their defaults. Give RecoveryAnalysis vec_type=cuda
%   too, or its direct solve on the GPU matrices returns zeros. Right
%   preconditioning makes ksp_rtol bound the true residual that
%   md.settings.solver_residue_threshold checks.

%retrieve options provided in varargin
options=pairoptions(varargin{:});
solverOptions=struct();
solverOptions.toolkit='petsc';
solverOptions.vec_type=getfieldvalue(options,'vec_type','cuda');
solverOptions.mat_type=getfieldvalue(options,'mat_type','aijcusparse');
solverOptions.mat_block_size=getfieldvalue(options,'mat_block_size','2');
solverOptions.ksp_type=getfieldvalue(options,'ksp_type','gmres');
solverOptions.ksp_pc_side=getfieldvalue(options,'ksp_pc_side','right');
solverOptions.ksp_rtol=getfieldvalue(options,'ksp_rtol','1e-8');
solverOptions.ksp_max_it=getfieldvalue(options,'ksp_max_it','300');
solverOptions.ksp_gmres_restart=getfieldvalue(options,'ksp_gmres_restart','100');
solverOptions.pc_type=getfieldvalue(options,'pc_type','hypre');
solverOptions.pc_hypre_type=getfieldvalue(options,'pc_hypre_type','boomeramg');
solverOptions.pc_hypre_boomeramg_strong_threshold=getfieldvalue(options,'pc_hypre_boomeramg_strong_threshold','0.25');
