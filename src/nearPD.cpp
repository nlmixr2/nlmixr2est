// https://github.com/cran/Matrix/blob/master/R/nearPD.R
// Contributors Martin Maechler, Jens Oehlschlägel
#define ARMA_WARN_LEVEL 1
#define STRICT_R_HEADER
#include "armahead.h"
#include "nearPD.h"

using namespace arma;
using namespace Rcpp;
#include <lotri.h>
extern "C" {
#define iniLotriPtr _nlmixr2est_iniLotriPtr
  iniLotri
}

lotriNearPDarmaSetup

bool nmNearPD(mat &ret, mat x
              , bool keepDiag// = false
              , bool do2eigen// = true  // if TRUE do a sfsmisc::posdefify() eigen step
              , bool doDykstra// = true // do use Dykstra's correction
              , bool only_values// = false // if TRUE simply return lambda[j].
              , double eig_tol//   = 1e-6 // defines relative positiveness of eigenvalues compared to largest
              , double conv_tol//  = 1e-7 // convergence tolerance for algorithm
              , double posd_tol//  = 1e-8 // tolerance for enforcing positive definiteness
              , int maxit//    = 100 // maximum number of iterations allowed
              , bool trace// = false // set to TRUE (or 1 ..) to trace iterations
              ){
  return lotriNearPDarma(ret, x, keepDiag, do2eigen, doDykstra, only_values, eig_tol, conv_tol, posd_tol, maxit, trace);
}

//[[Rcpp::export]]
RObject nmNearPD_(RObject x
                  , bool keepDiag = false
                  , bool do2eigen = true  // if TRUE do a sfsmisc::posdefify() eigen step
                  , bool doDykstra = true // do use Dykstra's correction
                  , bool only_values = false // if TRUE simply return lambda[j].
                  , double eig_tol   = 1e-6 // defines relative positiveness of eigenvalues compared to largest
                  , double conv_tol  = 1e-7 // convergence tolerance for algorithm
                  , double posd_tol  = 1e-8 // tolerance for enforcing positive definiteness
                  , int maxit    = 100 // maximum number of iterations allowed
                  , bool trace = false // set to TRUE (or 1 ..) to trace iterations
                  ){
  arma::mat ret;
  if (nmNearPD(ret, as<arma::mat>(x), keepDiag, do2eigen,
               doDykstra, only_values, eig_tol, conv_tol, posd_tol, maxit, trace)) {
    return wrap(ret);
  } else {
    stop("nearest PD calculation failed");
  }
  return R_NilValue;
}

bool chol_sym(mat &Hout, mat &Hin) {
  mat H = 0.5*(Hin+Hin.t());
  if (!H.is_symmetric()) return false;
  return chol(Hout, H);
}
