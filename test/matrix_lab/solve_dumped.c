static char help[] = "Solve a linear system dumped from ISSM (-ksp_view_mat/-ksp_view_rhs binary:...).\n"
                     "  -f A.bin -rhs b.bin   matrix and right-hand side (required)\n"
                     "  -ref x.bin            reference solution (e.g. MUMPS) to report the error against\n"
                     "  -coords coords.bin    interleaved node coordinates per DOF block: attach rigid-body\n"
                     "                        modes as the near-null space (used by GAMG)\n"
                     "  -bs 2                 block size (DOFs per node), default 2\n"
                     "All solver choices come from the usual PETSc options (-ksp_type, -pc_type, -mat_type, ...).\n";

#include <petscksp.h>

/* Ignore the dump's .info file: it records ISSM's block size at dump time (1 unless the dumping
   run set -mat_block_size), which would override -bs. */
static PetscErrorCode OpenNoInfo(const char *path, PetscViewer *vw)
{
  PetscFunctionBeginUser;
  PetscCall(PetscViewerCreate(PETSC_COMM_WORLD, vw));
  PetscCall(PetscViewerSetType(*vw, PETSCVIEWERBINARY));
  PetscCall(PetscViewerBinarySetSkipInfo(*vw, PETSC_TRUE));
  PetscCall(PetscViewerFileSetMode(*vw, FILE_MODE_READ));
  PetscCall(PetscViewerFileSetName(*vw, path));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode LoadVec(const char *path, Vec v)
{
  PetscViewer vw;

  PetscFunctionBeginUser;
  PetscCall(OpenNoInfo(path, &vw));
  PetscCall(VecLoad(v, vw));
  PetscCall(PetscViewerDestroy(&vw));
  PetscFunctionReturn(PETSC_SUCCESS);
}

int main(int argc, char **argv)
{
  char               fA[PETSC_MAX_PATH_LEN], fb[PETSC_MAX_PATH_LEN], fref[PETSC_MAX_PATH_LEN], fc[PETSC_MAX_PATH_LEN];
  PetscBool          flg, hasref, hascoords;
  PetscInt           bs = 2, n, its;
  Mat                A;
  Vec                b, x, r, ref = NULL, coords = NULL;
  KSP                ksp;
  KSPConvergedReason reason;
  PetscViewer        vw;
  PetscLogDouble     t0, t1, t2, t3;
  PetscReal          bnorm, rnorm, err = -1, refnorm;

  PetscCall(PetscInitialize(&argc, &argv, NULL, help));
  PetscCall(PetscOptionsGetString(NULL, NULL, "-f", fA, sizeof(fA), &flg));
  PetscCheck(flg, PETSC_COMM_WORLD, PETSC_ERR_USER, "need -f A.bin");
  PetscCall(PetscOptionsGetString(NULL, NULL, "-rhs", fb, sizeof(fb), &flg));
  PetscCheck(flg, PETSC_COMM_WORLD, PETSC_ERR_USER, "need -rhs b.bin");
  PetscCall(PetscOptionsGetString(NULL, NULL, "-ref", fref, sizeof(fref), &hasref));
  PetscCall(PetscOptionsGetString(NULL, NULL, "-coords", fc, sizeof(fc), &hascoords));
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-bs", &bs, NULL));

  /* Matrix: type from -mat_type; block size set before the load so GAMG/BoomerAMG see bs nodes */
  PetscCall(PetscTime(&t0));
  PetscCall(MatCreate(PETSC_COMM_WORLD, &A));
  PetscCall(MatSetBlockSize(A, bs));
  PetscCall(MatSetFromOptions(A));
  PetscCall(OpenNoInfo(fA, &vw));
  PetscCall(MatLoad(A, vw));
  PetscCall(PetscViewerDestroy(&vw));
  PetscCall(MatGetBlockSize(A, &bs));
  PetscCall(MatGetSize(A, &n, NULL));

  PetscCall(MatCreateVecs(A, &x, &b));
  PetscCall(VecDuplicate(b, &r));
  PetscCall(LoadVec(fb, b));
  if (hasref) {
    PetscCall(VecDuplicate(b, &ref));
    PetscCall(LoadVec(fref, ref));
  }
  if (hascoords) {
    MatNullSpace sp;
    PetscCall(VecDuplicate(b, &coords));
    PetscCall(LoadVec(fc, coords));
    PetscCall(VecSetBlockSize(coords, bs));
    PetscCall(MatNullSpaceCreateRigidBody(coords, &sp));
    PetscCall(MatSetNearNullSpace(A, sp));
    PetscCall(MatNullSpaceDestroy(&sp));
  }
  PetscCall(PetscTime(&t1));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "loaded n=%" PetscInt_FMT " bs=%" PetscInt_FMT " rigid-body modes %s, load %.2f s\n", n, bs,
                        hascoords ? "yes" : "no", t1 - t0));

  /* Solve: setup (preconditioner build, incl. host->device copy) timed separately from the iteration */
  PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));
  PetscCall(KSPSetOperators(ksp, A, A));
  PetscCall(KSPSetFromOptions(ksp));
  PetscCall(VecSet(x, 0.0));
  PetscCall(PetscTime(&t1));
  PetscCall(KSPSetUp(ksp));
  PetscCall(PetscTime(&t2));
  PetscCall(KSPSolve(ksp, b, x));
  PetscCall(PetscTime(&t3));
  PetscCall(KSPGetIterationNumber(ksp, &its));
  PetscCall(KSPGetConvergedReason(ksp, &reason));

  /* True residual ||b - A x|| / ||b|| (what ISSM's solver_residue_threshold checks) */
  PetscCall(MatMult(A, x, r));
  PetscCall(VecAYPX(r, -1.0, b));
  PetscCall(VecNorm(r, NORM_2, &rnorm));
  PetscCall(VecNorm(b, NORM_2, &bnorm));
  if (hasref) {
    PetscCall(VecNorm(ref, NORM_2, &refnorm));
    PetscCall(VecWAXPY(r, -1.0, ref, x)); /* r reused for x - ref */
    PetscCall(VecNorm(r, NORM_2, &err));
    err /= refnorm;
  }
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "RESULT setup %.3f s  solve %.3f s  total %.3f s  its %" PetscInt_FMT "  reason %s  relres %.3e  relerr_vs_ref %.3e\n",
                        t2 - t1, t3 - t2, t3 - t1, its, KSPConvergedReasons[reason], (double)(rnorm / bnorm), (double)err));

  PetscCall(KSPDestroy(&ksp));
  PetscCall(VecDestroy(&coords));
  PetscCall(VecDestroy(&ref));
  PetscCall(VecDestroy(&r));
  PetscCall(VecDestroy(&x));
  PetscCall(VecDestroy(&b));
  PetscCall(MatDestroy(&A));
  PetscCall(PetscFinalize());
  return 0;
}
