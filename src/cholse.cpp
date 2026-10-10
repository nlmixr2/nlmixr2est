#define STRICT_R_HEADER
#include "armahead.h"


using namespace Rcpp;
using namespace arma;

arma::mat gershNested(arma::mat A, int j, int n) {
  arma::mat g(n, 1, fill::zeros);
  double sumToI, sumAfterI;
  for (int ii = j; ii < n; ++ii){
    if (ii == 0){
      sumToI=0.0;
    } else if (j == ii){
      sumToI=arma::sum(arma::abs(A(ii, span(ii-1, j))));
    } else {
      sumToI=arma::sum(arma::abs(A(ii, span(j, ii-1))));
    }
    if (ii == n-1){
      sumAfterI = 0;
    } else {
      sumAfterI = arma::sum(arma::abs(A(span(ii+1, n-1), ii)));
    }
    g(ii, 0) = sumToI+sumAfterI-A(ii,ii);
  }
  return g;
}
// Suggested from https://gking.harvard.edu/files/help.pdf
// Translated from
//http://www.dynare.org/dynare-matlab-m2html/matlab/chol_SE.html
// Use tau1=sqrt(eps) instead of eps^1/3; In my tests eps^1/3 produces NaNs
bool cholSE0(arma::mat &Ao, arma::mat &E, arma::mat A, double tol) {
  int n = A.n_rows;
  double tau1 = tol;//
  double tau2 = tol;//tau1;
  bool phase1 = true;
  double delta = 0;
  int j;
  arma::mat P(n,1);
  // for (j = n; j--;) p(j,0) = j+1;
  arma::mat g(n,1, fill::zeros);
  // arma::mat E(n,1, fill::zeros);
  E = mat(n, 1, fill::zeros);
  double gamma = A(n-1,n-1);
  if (gamma < 0) phase1 = false;
  for (j = 0; j < n-1; j++){
    if (A(j, j) < 0) phase1 = false;
    if (A(j, j) > gamma) gamma = A(j,j);
  }
  double taugam = tau1*gamma;
  if (!phase1) g = gershNested(A, 0, n);
  // N=1 case
  if (n == 1){
    delta = tau2*std::fabs(A(0,0)) - A(0,0);
    if (delta > 0) E(0,0) = delta;
    if (A(0,0) == 0) E(0,0) = tau2;
    A(0,0)=_safe_sqrt(A(0,0)+E(0,0));
    Ao = A;
    return true;
  }
  int jp1, ii, k;
  double tempjj, temp=1., normj, tmp;
  for (j = 0; j < n-1;  j++){
    // Pivoting not included
    if (phase1){
      jp1 = j+1;
      if (A(j,j)>0){
        arma::mat tmp = (A(span(jp1,n-1),span(jp1,n-1))).diag() - A(span(jp1, n-1),j)%A(span(jp1, n-1),j)/A(j,j);
        double mintmp = tmp[0];
        for (ii = 1; ii < (int)tmp.size(); ii++) mintmp = (mintmp < tmp[ii]) ? mintmp : tmp[ii];
        if (mintmp < taugam) phase1=false;
      } else phase1 = false;

      if (phase1){
        // Do the normal cholesky update if still in phase 1
        A(j,j) = _safe_sqrt(A(j,j));
        tempjj = A(j,j);
        for (ii = jp1; ii < n; ii++){
          A(ii,j) = A(ii,j)/tempjj;
        }
        for (ii=jp1; ii <n; ii++){
          temp=A(ii,j);
          for (k = jp1; k < ii+1; k++){
            A(ii,k) = A(ii,k)-(temp * A(k,j));
          }
        }
        if (j == n-2){
          A(n-1,n-1)=_safe_sqrt(A(n-1,n-1));
        }
      } else {
        // Calculate the negatives of the lower Gershgorin bounds
        g=gershNested(A,j,n);
      }
    }

    if (!phase1){
      if (j != n-2){
        // Calculate delta and add to the diagonal. delta=max{0,-A(j,j) + max{normj,taugam},delta_previous}
        // where normj=sum of |A(i,j)|,for i=1,n, delta_previous is the delta computed at the previous iter and taugam is tau1*gamma.
        normj=arma::sum(arma::abs(A(span(j+1, n-1),j)));
        if (delta < 0) delta = 0;
        tmp  = -A(j,j)+normj;
        if (delta < tmp) delta = tmp;
        tmp  = -A(j,j)+taugam;
        if (delta < tmp) delta = tmp;
        // get adjustment based on formula on bottom of p. 309 of Eskow/Schnabel (1991)
        E(j,0) =  delta;
        A(j,j) = A(j,j) + E(j,0);
        // Update the Gershgorin bound estimates (note: g(i) is the negative of the Gershgorin lower bound.)
        if (A(j,j) != normj){
          temp = (normj/A(j,j)) - 1;
          for (ii = j+1; ii < n; ii++){
            g(ii) = g(ii) + std::fabs(A(ii,j)) * temp;
          }
        }
        for (int ii = j+1; ii < n; ii++){
          g(ii,0) = g(ii,0) + std::fabs(A(ii,j)) * temp;
        }
        // Do the cholesky update
        A(j,j) = _safe_sqrt(A(j,j));
        tempjj = A(j,j);
        for (ii = j+1; ii < n; ii++){
          A(ii,j) = A(ii,j) / tempjj;
        }
        for (ii = j+1; ii < n; ii++){
          temp = A(ii,j);
          for (k = j+1; k < ii+1; k++){
            A(ii,k) = A(ii,k) - (temp * A(k,j));
          }
        }
      } else {
        // Find eigenvalues of final 2 by 2 submatrix
        // Find delta such that:
        // 1.  the l2 condition number of the final 2X2 submatrix + delta*I <= tau2
        // 2. delta >= previous delta,
        // 3. min(eigvals) + delta >= tau2 * gamma, where min(eigvals) is the smallest eigenvalue of the final 2X2 submatrix
        // A(n-2,n-1)=A(n-1,n-2);
        //set value above diagonal for computation of eigenvalues
        A(n-2,n-1)=A(n-1,n-2); //set value above diagonal for computation of eigenvalues
        mat Ain = A(span(n-2, n-1),span(n-2, n-1));
        Ain = 0.5*(Ain+Ain.t());
        vec eigvals  = eig_sym(Ain);
        // Formula 5.3.2 of Schnabel/Eskow (1990)
        if (delta < 0) delta = 0;
        tmp= (max(eigvals)-min(eigvals))/(1-tau1);
        if (tmp < gamma) tmp = gamma;
        tmp=tau2*tmp;
        tmp =tmp - min(eigvals);
        if (delta < tmp) delta = tmp;
        if (delta > 0){
          A(n-2, n-2) = A(n-2,n-2) + delta;
          A(n-1, n-1) = A(n-1,n-1) + delta;
          E(n-2, 0) = delta;
          E(n-1, 0) = delta;
        }
        // Final update
        A(n-2,n-2) = _safe_sqrt(A(n-2,n-2));
        A(n-1,n-2) = A(n-1,n-2)/A(n-2,n-2);
        A(n-1,n-1) = A(n-1,n-1) - A(n-1,n-2)*A(n-1,n-2);
        A(n-1,n-1) = _safe_sqrt(A(n-1,n-1));
      }
    }
  }
  Ao = (trimatl(A)).t();
  return phase1;
}

arma::mat cholSE__(arma::mat A, double tol) {
  arma::mat Ao, E;
  cholSE0(Ao, E, A, tol);
  return Ao;
}
//[[Rcpp::export]]
NumericMatrix cholSE_(NumericMatrix A, double tol){
  arma::mat Ao, E;
  cholSE0(Ao, E, as<arma::mat>(A), tol);
  return wrap(Ao);
}

// Whether M0's eigenvalues reach down to rounding level (or cannot be computed)
static bool covRankDeficient(const arma::mat &M0) {
  arma::vec ev;
  if (!arma::eig_sym(ev, arma::symmatu(M0))) return true;
  ev = arma::abs(ev);
  return ev.min() <= ev.max() * M0.n_rows * arma::datum::eps;
}

// The "|.|" rung: U = chol(sqrtm(M0 %*% M0)), and that matrix in *Mabs when given
static bool covAbsRung(const arma::mat &M0, arma::mat &U, arma::mat *Mabs) {
  arma::cx_mat H1;
  arma::mat ch;
  if (!arma::sqrtmat(H1, M0*M0) || arma::any(arma::any(arma::imag(H1), 0)) ||
      !arma::chol(ch, arma::real(H1))) return false;
  U = ch;
  if (Mabs != nullptr) *Mabs = arma::real(H1);
  return true;
}

// Whether an information matrix (R) or score cross-product (S) M0 can be used, after
// cholSE0 (pd, E, U): 1 as it is; 2 corrected, cholSE0's factor of M0 + diag(E), when
// every added diagonal is within cholAccept ("+"); 3 chol(sqrtm(M0 %*% M0)) ("|.|", U
// replaced, and the matrix itself in *Mabs when given); 0 not usable.  The "+" rung also
// needs a positive largest diagonal (cholSE0 scales E by it) and a finite E and factor.
// cholSE0 calls every 1x1 matrix positive definite, so a 1x1 M0 is judged by its value.
// A numerically rank-deficient M0 may be "+" but never "|.|" (sqrtm lifts rounding-level
// eigenvalues to about sqrt(eps)).
int covAcceptRule(const arma::mat &M0, bool pd, const arma::vec &E, arma::mat &U,
                  double cholAccept, arma::mat *Mabs) {
  // cholSE0 calls a matrix with NaN positive definite
  if (M0.n_elem == 0 || !M0.is_finite()) return 0;
  if (pd && (M0.n_elem != 1 || M0(0, 0) > 0)) return 1;
  if (M0.diag().max() > 0 && E.is_finite() && U.is_finite() && !arma::any(E > cholAccept)) {
    return 2;
  }
  if (covRankDeficient(M0)) return 0;
  return covAbsRung(M0, U, Mabs) ? 3 : 0;
}

// covAcceptRule() for R callers: type "" (as it is), "+", "|" or "failed", the factor U
// and the matrix it factors (sqrtm(A %*% A) for "|", else A).
//[[Rcpp::export]]
List covAccept_(NumericMatrix A, double cholSEtol, double cholAccept) {
  arma::mat M0 = as<arma::mat>(A), U, E, Mabs;
  if (M0.n_elem == 0 || !M0.is_finite()) {
    return List::create(_["type"] = "failed", _["U"] = U, _["M"] = M0);
  }
  bool pd = cholSE0(U, E, M0, cholSEtol);
  arma::vec Ev = arma::vectorise(E);
  int rc = covAcceptRule(M0, pd, Ev, U, cholAccept, &Mabs);
  static const char *types[4] = {"failed", "", "+", "|"};
  arma::mat M = (rc == 3) ? Mabs : M0;
  return List::create(_["type"] = types[rc], _["U"] = U, _["M"] = M);
}

// cholSE0's factor U (U'U = A + diag(E)), E and whether A was factored without
// adding anything (pd), for callers that apply foceiCovUsable()'s rules in R.
//[[Rcpp::export]]
List cholSEpd_(NumericMatrix A, double tol) {
  arma::mat Ao, E;
  bool pd = cholSE0(Ao, E, as<arma::mat>(A), tol);
  // the objects themselves: create() wraps each once its result is protected
  return List::create(_["U"] = Ao, _["E"] = arma::vec(arma::vectorise(E)),
                      _["pd"] = pd);
}
