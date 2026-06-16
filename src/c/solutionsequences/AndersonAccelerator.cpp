/*!\file:  AndersonAccelerator.cpp
 * \brief: Anderson acceleration for Picard fixed-point iterations.
 *
 * See AndersonAccelerator.h for theory and usage.
 */

#include "AndersonAccelerator.h"

/* ──────────────────────────────────────────────────────────────────────────
 * Construction / destruction
 * ────────────────────────────────────────────────────────────────────────── */

AndersonAccelerator::AndersonAccelerator(int m) : depth(m), iter(0) {
    /*reserve to avoid repeated allocation for typical depths <= 10*/
    if(m > 0){
        f_hist.reserve(m+1);
        G_hist.reserve(m+1);
    }
}

AndersonAccelerator::~AndersonAccelerator(){
    for(int i=0; i<(int)f_hist.size(); i++) delete f_hist[i];
    for(int i=0; i<(int)G_hist.size(); i++) delete G_hist[i];
}

/* ──────────────────────────────────────────────────────────────────────────
 * SmallSolve
 *
 * Gaussian elimination with partial pivoting for a dense mk x mk system.
 * On entry:  A[i*mk+j], b[i].
 * On exit:   theta[i] holds the solution.
 * Returns false and leaves theta zeroed if A is numerically singular.
 * ────────────────────────────────────────────────────────────────────────── */
bool AndersonAccelerator::SmallSolve(IssmDouble* A, IssmDouble* b,
                                     IssmDouble* theta, int mk){

    const IssmDouble tol = 1.0e-14;

    /*copy b into theta; we will overwrite it as we go*/
    for(int i=0; i<mk; i++) theta[i] = b[i];

    /*forward elimination*/
    for(int col=0; col<mk; col++){

        /*partial pivot: find row with largest absolute value in this column*/
        int pivot = col;
        IssmDouble maxval = fabs(A[col*mk+col]);
        for(int row=col+1; row<mk; row++){
            IssmDouble val = fabs(A[row*mk+col]);
            if(val > maxval){ maxval = val; pivot = row; }
        }

        if(maxval < tol){
            /*singular or near-singular — fall back to Picard (theta = 0)*/
            for(int i=0; i<mk; i++) theta[i] = 0.0;
            return false;
        }

        /*swap rows col and pivot in A and theta*/
        if(pivot != col){
            for(int j=0; j<mk; j++){
                IssmDouble tmp = A[col*mk+j];
                A[col*mk+j]    = A[pivot*mk+j];
                A[pivot*mk+j]  = tmp;
            }
            IssmDouble tmp  = theta[col];
            theta[col]      = theta[pivot];
            theta[pivot]    = tmp;
        }

        /*eliminate column below the pivot*/
        IssmDouble inv_diag = 1.0 / A[col*mk+col];
        for(int row=col+1; row<mk; row++){
            IssmDouble factor = A[row*mk+col] * inv_diag;
            for(int j=col; j<mk; j++)
                A[row*mk+j] -= factor * A[col*mk+j];
            theta[row] -= factor * theta[col];
        }
    }

    /*back substitution*/
    for(int row=mk-1; row>=0; row--){
        for(int j=row+1; j<mk; j++)
            theta[row] -= A[row*mk+j] * theta[j];
        theta[row] /= A[row*mk+row];
    }

    return true;
}

/* ──────────────────────────────────────────────────────────────────────────
 * Apply
 *
 * Core of Anderson acceleration.  Called once per nonlinear iteration
 * after the linear solve has produced a new Picard update in *puf.
 *
 * Algorithm (Walker & Ni 2011, "Type 1"):
 *
 *   f_k  = G(u_k) - u_k                    Picard residual
 *   mk   = min(iter, depth)                 effective window
 *
 *   if mk == 0:  u_{k+1} = G(u_k)          plain Picard
 *   else:
 *     dF[:,j] = f_hist[j+1] - f_hist[j]    mk difference columns
 *     dG[:,j] = G_hist[j+1] - G_hist[j]
 *     theta   = argmin_t || f_k - dF t ||  (mk x mk normal equations)
 *     u_{k+1} = G(u_k) - dG theta          accelerated update
 * ────────────────────────────────────────────────────────────────────────── */
void AndersonAccelerator::Apply(Vector<IssmDouble>** puf,
                                Vector<IssmDouble>* old_uf){

    if(depth == 0) return;   /* disabled: leave *puf untouched */

    Vector<IssmDouble>* G_k = *puf;   /* Picard update G(u_k) from linear solve */

    /* ── 1. Compute Picard residual  f_k = G_k - u_k  ───────────────────── */
    Vector<IssmDouble>* f_k = G_k->Duplicate();
    G_k->Copy(f_k);
    f_k->AXPY(old_uf, -1.0);   /* f_k = G_k - old_uf */

    /* ── 2. Store G_k and f_k in history before any modification ─────────── */
    Vector<IssmDouble>* G_store = G_k->Duplicate();
    G_k->Copy(G_store);

    f_hist.push_back(f_k);
    G_hist.push_back(G_store);

    /* trim history to at most depth+1 entries (depth+1 gives depth differences)*/
    while((int)f_hist.size() > depth+1){
        delete f_hist.front();  f_hist.erase(f_hist.begin());
        delete G_hist.front();  G_hist.erase(G_hist.begin());
    }

    /* ── 3. Determine effective window size ──────────────────────────────── */
    int mk = (int)f_hist.size() - 1;   /* number of difference columns */

    if(mk == 0){
        /* first iteration or depth==1 first step: keep plain Picard update */
        iter++;
        return;
    }

    /* ── 4. Build difference vectors dF[j] and dG[j]  ───────────────────── */
    /*  dF[j] = f_hist[j+1] - f_hist[j]   (j = 0 ... mk-1)
     *  dG[j] = G_hist[j+1] - G_hist[j]                                      */
    std::vector<Vector<IssmDouble>*> dF(mk, NULL);
    std::vector<Vector<IssmDouble>*> dG(mk, NULL);

    for(int j=0; j<mk; j++){
        dF[j] = f_hist[j+1]->Duplicate();
        f_hist[j+1]->Copy(dF[j]);
        dF[j]->AXPY(f_hist[j], -1.0);      /* dF[j] = f_hist[j+1] - f_hist[j] */

        dG[j] = G_hist[j+1]->Duplicate();
        G_hist[j+1]->Copy(dG[j]);
        dG[j]->AXPY(G_hist[j], -1.0);      /* dG[j] = G_hist[j+1] - G_hist[j] */
    }

    /* ── 5. Form normal equations  A theta = b  where
     *       A[i][j] = dF[i] . dF[j]
     *       b[i]    = dF[i] . f_k                                            */
    IssmDouble* A     = new IssmDouble[mk*mk];
    IssmDouble* b_rhs = new IssmDouble[mk];
    IssmDouble* theta = new IssmDouble[mk]();

    Vector<IssmDouble>* f_k_cur = f_hist.back();   /* = f_k we just pushed */

    for(int i=0; i<mk; i++){
        b_rhs[i] = dF[i]->Dot(f_k_cur);
        for(int j=0; j<=i; j++){
            A[i*mk+j] = A[j*mk+i] = dF[i]->Dot(dF[j]);
        }
    }

    /* ── 6. Solve the mk x mk system ────────────────────────────────────── */
    bool ok = SmallSolve(A, b_rhs, theta, mk);

    /* ── 7. Anderson update:  G_k -= sum_j theta[j] * dG[j]  ─────────────
     *   If SmallSolve failed (singular) theta is zeroed so G_k is unchanged,
     *   effectively falling back to Picard for this iteration.               */
    if(ok){
        for(int j=0; j<mk; j++){
            G_k->AXPY(dG[j], -theta[j]);   /* G_k -= theta[j] * dG[j] */
        }
        if(VerboseConvergence())
            _printf0_("   Anderson(" << depth << "): mk=" << mk << " applied\n");
    }
    else{
        if(VerboseConvergence())
            _printf0_("   Anderson(" << depth << "): mk=" << mk
                      << " singular, falling back to Picard\n");
    }

    /* ── 8. Clean up ─────────────────────────────────────────────────────── */
    for(int j=0; j<mk; j++){ delete dF[j]; delete dG[j]; }
    delete[] A;
    delete[] b_rhs;
    delete[] theta;

    /* *puf still points to G_k which has been modified in-place */
    iter++;
}
