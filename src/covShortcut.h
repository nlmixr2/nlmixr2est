#ifndef __COVSHORTCUT_H__
#define __COVSHORTCUT_H__
#if defined(__cplusplus)

#include "armahead.h"
#include "fdHess.h"

// The verified shortcut of the covariance step's finite-difference Hessian
// (covShortcut): after the measured diagonal, a full Hessian predicted from it and
// a hint's correlations is checked along a few fixed directions instead of the
// off-diagonals being measured pair by pair.

// What covShortcutVerify() found: status 1 when the prediction explains every
// check direction, 0 when one does not (the caller measures the off-diagonals),
// -1 when a probe failed; checks has one row per direction: the measured and
// predicted second differences and the allowance between them.
struct CovShortcutResult {
  int status = -1;
  arma::mat checks;
};

CovShortcutResult covShortcutVerify(FdHessObj &obj, double *x, int n, double f0,
                                    const double *h, const arma::mat &H, double fac);

#endif
#endif
