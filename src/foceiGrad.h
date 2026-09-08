#ifndef __FOCEIGRAD_H__
#define __FOCEIGRAD_H__
#if defined(__cplusplus)
#include <functional>

struct FoceiHessianNode {
  arma::vec f, R;
  arma::mat a, aR;
  arma::cube A, AR;
};

struct FoceiHessianQuadrature {
  arma::mat points, weights;
  std::function<bool(const arma::vec&, FoceiHessianNode&)> evaluate;
};

arma::mat foceiRSubjectFR_(const arma::mat& a, const arma::cube& A, const arma::cube& Ath,
                         const arma::mat& aR, const arma::cube& AR, const arma::cube& AthR,
                         const arma::mat& dvSens, const arma::mat& dvSens2,
                         const arma::ivec& censv, const arma::vec& limv,
                         const arma::vec& fv, const arma::vec& yv, const arma::vec& Rv,
                         const arma::vec& ehat, const arma::mat& Oi,
                         const arma::cube& dOi, const arma::cube& d2Oi, const arma::mat& d2LD,
                         int neta, int ndir, int ndirP, int nom, const arma::ivec& dirP,
                         const FoceiHessianQuadrature *quadrature = nullptr);

arma::mat foceiRSubjectFoceFR_(const arma::mat& a, const arma::cube& A, const arma::cube& Ath,
                             const arma::mat& aRe, const arma::mat& aRc,
                             const arma::cube& ARe, const arma::cube& ARc,
                             const arma::mat& dvSens, const arma::mat& dvSens2,
                             const arma::ivec& censv, const arma::vec& limv,
                             const arma::vec& fv, const arma::vec& yv, const arma::vec& Rv,
                             const arma::vec& ehat, const arma::mat& Oi,
                             const arma::cube& dOi, const arma::cube& d2Oi, const arma::mat& d2LD,
                             int neta, int ndir, int ndirP, int nom, const arma::ivec& dirP);

// Per-subject FOCEI (f,R) outer-gradient kernel.  Declared here so the pooled driver in
// inner.cpp can call it directly: that driver already holds the per-subject sensitivity
// structures in C++, and routing them out to R only to have R stack them and hand them
// back was the round trip Phase 8E removes.
void foceiGradSubjectFR_(const arma::mat& a, const arma::cube& A,
                         const arma::mat& aR, const arma::cube& AR,
                         const arma::mat& Rsig, const arma::cube& RsigDir,
                         const arma::mat& dvSens,
                         const arma::ivec& censv, const arma::vec& limv, int censOpt,
                         const arma::vec& fv, const arma::vec& yv, const arma::vec& Rv,
                         const arma::vec& ehat, const arma::mat& Oi,
                         const arma::cube& dOiEst, const arma::vec& tr28,
                         int neta, int nth, int nsg, int nom,
                         const arma::ivec& dirTh, const arma::ivec& sigCol,
                         arma::vec& g_out, arma::mat& etaP_out);

// Per-subject FOCE (frozen-R0) kernel.  aRe drives the eta block (0 for nonmem, live
// E$aR for foce+), aRc the parameter columns; no AR cube.  fp = foce+ (1) / nonmem (0).
void foceiGradSubjectFoceFR_(const arma::mat& a, const arma::cube& A,
                             const arma::mat& aRe, const arma::mat& aRc,
                             const arma::mat& R0sig, const arma::mat& dvSens,
                             const arma::ivec& censv, const arma::vec& limv,
                             const arma::vec& fv, const arma::vec& yv, const arma::vec& R0v,
                             const arma::vec& ehat, const arma::mat& Oi,
                             const arma::cube& dOiEst, const arma::vec& tr28,
                             int neta, int nth, int nsg, int nom,
                             const arma::ivec& dirTh, const arma::ivec& sigCol, int fp,
                             arma::vec& g_out, arma::mat& etaP_out, int censOpt = 0);

// Per-subject AGQ kernel.  The node arrays are node-major: node k occupies rows
// k*nobs .. k*nobs+nobs-1, and y is node-invariant.  Reports failure through ok_out
// (non-throwing chol/inv) rather than throwing, unlike the two kernels above.
void foceiGradSubjectAgqFR_(const arma::mat& a, const arma::cube& A,
                            const arma::mat& aR, const arma::cube& AR,
                            const arma::mat& Rsig, const arma::cube& RsigDir,
                            const arma::vec& fv, const arma::vec& yv, const arma::vec& Rv,
                            const arma::mat& aN, const arma::mat& aRN, const arma::mat& RsigN,
                            const arma::vec& fN, const arma::vec& RN,
                            const arma::mat& qx, const arma::mat& qw,
                            const arma::vec& ehat, const arma::mat& Oi,
                            const arma::cube& dOiEst, const arma::vec& tr28,
                            int neta, int nth, int nsg, int nom,
                            const arma::ivec& dirTh, const arma::ivec& sigCol,
                            arma::vec& g_out, arma::mat& etaP_out, bool& ok_out,
                            const arma::ivec& censv = arma::ivec(),
                            const arma::vec& limv = arma::vec(), int censOpt = 0);
#endif
#endif
