// NPAG -- Nonparametric Adaptive Grid (Yamada 2021).  The estimation loop lives
// entirely here and reuses the FOCEI inner likelihood machinery in inner.cpp
// (via the np.h / imp.h numeric interfaces) to fill the Psi matrix of per-subject
// conditional likelihoods at each support point.  Called from foceiFitCpp_ in
// place of foceiOuter when est=="npag".
#include <RcppArmadillo.h>
#include "np.h"

using namespace Rcpp;

// Module M0: scaffold only.  Later milestones assemble the adaptive-grid cycle
// (Psi -> Burke IPM -> condense -> gamma optimization -> expand -> converge).
void npagOuter(Environment e) {
  stop("NPAG estimation ('est=\"npag\"') is not yet implemented");
}
