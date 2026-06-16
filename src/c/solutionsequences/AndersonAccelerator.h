/*!\file:  AndersonAccelerator.h
 * \brief: Anderson acceleration (depth-m mixing) for Picard fixed-point iterations.
 *
 * Theory: Walker & Ni (2011), SIAM J. Numer. Anal. 49(4), 1715-1735.
 *
 * Given the Picard map  G : u -> u_new = K(u)^{-1} f(u),
 * standard Picard sets  u_{k+1} = G(u_k).
 * Anderson(m) instead solves a small m x m least-squares problem each
 * iteration to find the best linear combination of the last m Picard updates,
 * replacing the plain update at negligible extra cost.
 *
 * Usage:
 *   AndersonAccelerator aa(m);          // m=0 -> pure Picard (no-op)
 *   for(;;){
 *     ... linear solve -> uf ...
 *     aa.Apply(&uf, old_uf);            // uf replaced with accelerated update
 *     ... UpdateInputFromSolution(uf) ...
 *   }
 */

#ifndef _ANDERSONACCELERATOR_H_
#define _ANDERSONACCELERATOR_H_

#include <vector>
#include "../toolkits/toolkits.h"
#include "../shared/shared.h"

class AndersonAccelerator {

  public:

    int  depth;  /* window size m; 0 = disabled (pure Picard) */
    int  iter;   /* iteration counter (reset on construction)  */

    /* history of Picard residuals  f_k = G(u_k) - u_k  */
    std::vector<Vector<IssmDouble>*> f_hist;
    /* history of raw Picard updates  G(u_k)  */
    std::vector<Vector<IssmDouble>*> G_hist;

    AndersonAccelerator(int m);
    ~AndersonAccelerator();

    /*! Replace *puf (= G(u_k) from linear solve) with the Anderson-accelerated
     *  update, given old_uf = u_k (solution vector from the previous iteration).
     *  On return *puf holds u_{k+1}. */
    void Apply(Vector<IssmDouble>** puf, Vector<IssmDouble>* old_uf);

  private:

    /*! Solve the small mk x mk system  A theta = b  in-place via
     *  Gaussian elimination with partial pivoting.  A and b are
     *  overwritten.  Returns false if the system is singular. */
    bool SmallSolve(IssmDouble* A, IssmDouble* b, IssmDouble* theta, int mk);
};

#endif  /* _ANDERSONACCELERATOR_H_ */
