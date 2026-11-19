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
// and only moves the off-diagonals, allowing more Dykstra iterations for it.
// It cannot fix a zero or negative diagonal, so when its result is not
// positive definite the diagonal is let go (keepDiag = false).  lotri's
// nearPD gives up when its iterations do not converge, which a zero or
// negative diagonal causes, so the last resort is the eigenvalue floor
// Matrix::nearPD() ends with.
// Returns 0 when nothing works, 1 for the kept diagonal, 2 when the diagonal
// had to move.
static inline int nmNearPDKeepDiag(mat &ret, mat x) {
  if (nmNearPD(ret, x, true, true, true, false, 1e-6, 1e-7, 1e-8, 1000) && ret.is_sympd()) {
    return 1;
  }
  if (nmNearPD(ret, x, false) && ret.is_sympd()) {
    return 2;
  }
  mat s = 0.5 * (x + x.t());
  vec d;
  mat Q;
  if (!s.is_finite() || !eig_sym(d, Q, s) || d.max() <= 0) return 0;
  d = arma::max(d, arma::vec(d.n_elem, arma::fill::value(1e-6 * d.max())));
  ret = Q * diagmat(d) * Q.t();
  ret = 0.5 * (ret + ret.t());
  return ret.is_sympd() ? 2 : 0;
}

bool chol_sym(mat &Hout, mat& Hin);

#endif
#endif
