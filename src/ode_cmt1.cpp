#define STRICT_R_HEADER
#include <stan/math/prim/mat/fun/Eigen.hpp> // must come before #include <RcppEigen.h>
#include "../inst/include/nlmixr2est_types.h"
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;

#include <vector>
#include <stan/math/rev/core.hpp>
#include <stan/math.hpp>
#include "PKPDLib_WW.h"

// Change to the Rcpp interface and let Rcpp handle the BEGIN_RCPP interface.
// It may have changed in the new RcppEigen...

//[[Rcpp::export]]
SEXP lin_cmt_stan(Eigen::Map<Eigen::VectorXd> obs_time,
		  Eigen::Map<Eigen::VectorXd> dose_time,
		  Eigen::Map<Eigen::VectorXd> dose,
		  Eigen::Map<Eigen::VectorXd> Tinf,
		  Eigen::Map<Eigen::VectorXd> params,
		  SEXP oralSEXP,
		  SEXP infusionSEXP,
		  SEXP ncmtSEXP,
		  SEXP parameterizationSEXP ) {
  const int oral = as<int>(oralSEXP);
  const int infusion = as<int>(infusionSEXP);
  const int ncmt = as<int>(ncmtSEXP);
  const int parameterization = as<int>(parameterizationSEXP);
  stan::math::lin_cmt_fun f(obs_time, dose_time, dose, Tinf, ncmt, oral, infusion, parameterization);
  Eigen::VectorXd fx;
  Eigen::Matrix<double, -1, -1> J;
  // ncol = npar
  /// nrow = nobs
  stan::math::jacobian(f, params, fx, J);

  return Rcpp::List::create(Rcpp::Named("fx") = wrap(fx),
			    Rcpp::Named("J") = wrap(J));
}

extern void lin_cmt_stanC(double *obs_timeD, const int nobs, double *dose_timeD, const int ndose, double *doseD, double *TinfD,
			  double *paramsD, const int oral, const int infusion, const int ncmt, const int parameterization,
			  const int neta, double *fxD, double *dvdxD, double *fpD){
  Eigen::Map<Eigen::VectorXd> obs_time(obs_timeD, nobs);
  Eigen::Map<Eigen::VectorXd> dose_time(dose_timeD, ndose);
  Eigen::Map<Eigen::VectorXd> dose(doseD, ndose);
  Eigen::Map<Eigen::VectorXd> Tinf(TinfD, ndose);
  Eigen::Map<Eigen::VectorXd> params(paramsD, (int)(2*ncmt+2));
  stan::math::lin_cmt_fun f(obs_time, dose_time, dose, Tinf, ncmt, oral, infusion, parameterization);
  Eigen::VectorXd fx;
  Eigen::Matrix<double, -1, -1> J;
  // Jacobian ncols=npars
  // nrows=nobs
  
  stan::math::jacobian(f, params, fx, J);
  std::copy(&fx[0],&fx[0]+nobs, fpD);
  // dvdx
  // ncol = netas
  // nrow = npars
  Eigen::Map<Eigen::MatrixXd> dvdx(dvdxD, 2*ncmt+2, neta);
  Eigen::Map<Eigen::MatrixXd> fp(fpD, nobs, neta);
  fp = J.transpose() * dvdx;
}




// I'm not sure why the below is here...

#if 0

require(Rcpp)
dyn.load("nlmixr2est.so")
lin_cmt <- function(obs_time,dose_time,dose,Tinf,params,oral,infusion,ncmt,parameterization)
   .Call('lin_cmt_stan', obs_time,dose_time,dose,Tinf,params,oral,infusion,ncmt,parameterization)
a = lin_cmt(0:72*1.0, 0:1*24*1.0, rep(10.0,2), 0, c(.1, 1, .2, 0), 1, 0, 1, 1)
b = matrix(scan("1"),  nrow=73)
range(b[,1] - a$fx)
range(as.vector(b[,-1]) - a$J)


ode_sol <- function(inits,time,evid,amt,params,absolute_tolerance,relative_tolerance,nobs,wh)
   .Call('ode_sol', inits,time,evid,amt,params,absolute_tolerance,relative_tolerance,nobs,wh)

require(rxode2)
e = eventTable()
e$add.dosing(10, 2, 24)
e$add.sampling(0:72)
x = e$get.EventTable()
x$amt[is.na(x$amt)] = 0

ode_sol(c(0,0), x$time, x$evid, x$amt, c(.2, .1), 1e-8, 1e-8, 73, 2)

#endif

