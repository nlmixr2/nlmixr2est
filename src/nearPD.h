#ifndef __NEARPD_H__
#define __NEARPD_H__
#if defined(__cplusplus)

using namespace arma;

// The defaults are Matrix::nearPD()'s, the same as the R nmNearPD()
bool nmNearPD(mat &ret, mat x, bool keepDiag = false,
             bool do2eigen = true, bool doDykstra = true, bool only_values = false,
             double eig_tol   = 1e-6, double conv_tol  = 1e-7, double posd_tol  = 1e-8,
             int maxit    = 100, bool trace = false // set to TRUE (or 1 ..) to trace iterations
             );

// A repair that keeps the diagonal (the variances, or the per-eta curvature)
// and only moves the off-diagonals, allowing more Dykstra iterations for it
static inline bool nmNearPDKeepDiag(mat &ret, mat x) {
  return nmNearPD(ret, x, true, true, true, false, 1e-6, 1e-7, 1e-8, 1000);
}

bool chol_sym(mat &Hout, mat& Hin);

#endif
#endif
