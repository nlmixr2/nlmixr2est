// NPB -- Nonparametric Bayes: a truncated stick-breaking Dirichlet-process
// mixture sampled by a blocked Metropolis-within-Gibbs sampler (Tatarinova 2013).
// Shares the Psi conditional-likelihood primitive with NPAG (npag.cpp) and reuses
// the FOCEI inner machinery in inner.cpp.  Called from foceiFitCpp_ in place of
// foceiOuter when est=="npb".
#include <RcppArmadillo.h>
#include "np.h"

using namespace Rcpp;

// Module M0: scaffold only.  Later milestones implement the blocked Gibbs sampler
// (cluster assignments -> stick weights -> MH location moves), Gelman-Rubin
// convergence, and posterior / credible-interval summaries.
void npbOuter(Environment e) {
  stop("NPB estimation ('est=\"npb\"') is not yet implemented");
}
