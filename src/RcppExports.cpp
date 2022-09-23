// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/nlmixr2est.h"
#include "../inst/include/nlmixr2est_types.h"
#include <RcppEigen.h>
#include <Rcpp.h>
#include <string>
#include <set>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// cholSE_
NumericMatrix cholSE_(NumericMatrix A, double tol);
RcppExport SEXP _nlmixr2est_cholSE_(SEXP ASEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type A(ASEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(cholSE_(A, tol));
    return rcpp_result_gen;
END_RCPP
}
// nlmixrExpandFdParNlme_
List nlmixrExpandFdParNlme_(CharacterVector state, CharacterVector vars);
static SEXP _nlmixr2est_nlmixrExpandFdParNlme__try(SEXP stateSEXP, SEXP varsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type state(stateSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type vars(varsSEXP);
    rcpp_result_gen = Rcpp::wrap(nlmixrExpandFdParNlme_(state, vars));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _nlmixr2est_nlmixrExpandFdParNlme_(SEXP stateSEXP, SEXP varsSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_nlmixr2est_nlmixrExpandFdParNlme__try(stateSEXP, varsSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// filterNormalLikeAndDoses
List filterNormalLikeAndDoses(IntegerVector& inCmt, IntegerVector& inDistribution, IntegerVector& inDistCmt);
RcppExport SEXP _nlmixr2est_filterNormalLikeAndDoses(SEXP inCmtSEXP, SEXP inDistributionSEXP, SEXP inDistCmtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector& >::type inCmt(inCmtSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type inDistribution(inDistributionSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type inDistCmt(inDistCmtSEXP);
    rcpp_result_gen = Rcpp::wrap(filterNormalLikeAndDoses(inCmt, inDistribution, inDistCmt));
    return rcpp_result_gen;
END_RCPP
}
// freeFocei
void freeFocei();
RcppExport SEXP _nlmixr2est_freeFocei() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    freeFocei();
    return R_NilValue;
END_RCPP
}
// foceiInnerLp
NumericVector foceiInnerLp(NumericVector eta, int id);
RcppExport SEXP _nlmixr2est_foceiInnerLp(SEXP etaSEXP, SEXP idSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< int >::type id(idSEXP);
    rcpp_result_gen = Rcpp::wrap(foceiInnerLp(eta, id));
    return rcpp_result_gen;
END_RCPP
}
// likInner
double likInner(NumericVector eta, int id);
RcppExport SEXP _nlmixr2est_likInner(SEXP etaSEXP, SEXP idSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< int >::type id(idSEXP);
    rcpp_result_gen = Rcpp::wrap(likInner(eta, id));
    return rcpp_result_gen;
END_RCPP
}
// foceiLik
double foceiLik(NumericVector theta);
RcppExport SEXP _nlmixr2est_foceiLik(SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(foceiLik(theta));
    return rcpp_result_gen;
END_RCPP
}
// foceiOfv
double foceiOfv(NumericVector theta);
RcppExport SEXP _nlmixr2est_foceiOfv(SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(foceiOfv(theta));
    return rcpp_result_gen;
END_RCPP
}
// foceiNumericGrad
NumericVector foceiNumericGrad(NumericVector theta);
RcppExport SEXP _nlmixr2est_foceiNumericGrad(SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(foceiNumericGrad(theta));
    return rcpp_result_gen;
END_RCPP
}
// foceiSetup_
NumericVector foceiSetup_(const RObject& obj, const RObject& data, NumericVector theta, Nullable<LogicalVector> thetaFixed, Nullable<LogicalVector> skipCov, RObject rxInv, Nullable<NumericVector> lower, Nullable<NumericVector> upper, Nullable<NumericMatrix> etaMat, Nullable<List> control);
RcppExport SEXP _nlmixr2est_foceiSetup_(SEXP objSEXP, SEXP dataSEXP, SEXP thetaSEXP, SEXP thetaFixedSEXP, SEXP skipCovSEXP, SEXP rxInvSEXP, SEXP lowerSEXP, SEXP upperSEXP, SEXP etaMatSEXP, SEXP controlSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const RObject& >::type obj(objSEXP);
    Rcpp::traits::input_parameter< const RObject& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< Nullable<LogicalVector> >::type thetaFixed(thetaFixedSEXP);
    Rcpp::traits::input_parameter< Nullable<LogicalVector> >::type skipCov(skipCovSEXP);
    Rcpp::traits::input_parameter< RObject >::type rxInv(rxInvSEXP);
    Rcpp::traits::input_parameter< Nullable<NumericVector> >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< Nullable<NumericVector> >::type upper(upperSEXP);
    Rcpp::traits::input_parameter< Nullable<NumericMatrix> >::type etaMat(etaMatSEXP);
    Rcpp::traits::input_parameter< Nullable<List> >::type control(controlSEXP);
    rcpp_result_gen = Rcpp::wrap(foceiSetup_(obj, data, theta, thetaFixed, skipCov, rxInv, lower, upper, etaMat, control));
    return rcpp_result_gen;
END_RCPP
}
// foceiOuterF
double foceiOuterF(NumericVector& theta);
RcppExport SEXP _nlmixr2est_foceiOuterF(SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(foceiOuterF(theta));
    return rcpp_result_gen;
END_RCPP
}
// foceiOuterG
NumericVector foceiOuterG(NumericVector& theta);
RcppExport SEXP _nlmixr2est_foceiOuterG(SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(foceiOuterG(theta));
    return rcpp_result_gen;
END_RCPP
}
// foceiOuter
Environment foceiOuter(Environment e);
RcppExport SEXP _nlmixr2est_foceiOuter(SEXP eSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Environment >::type e(eSEXP);
    rcpp_result_gen = Rcpp::wrap(foceiOuter(e));
    return rcpp_result_gen;
END_RCPP
}
// nlmixr2Gill83_
List nlmixr2Gill83_(Function what, NumericVector args, Environment envir, LogicalVector which, double gillRtol, int gillK, double gillStep, double gillFtol, bool optGillF);
RcppExport SEXP _nlmixr2est_nlmixr2Gill83_(SEXP whatSEXP, SEXP argsSEXP, SEXP envirSEXP, SEXP whichSEXP, SEXP gillRtolSEXP, SEXP gillKSEXP, SEXP gillStepSEXP, SEXP gillFtolSEXP, SEXP optGillFSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Function >::type what(whatSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type args(argsSEXP);
    Rcpp::traits::input_parameter< Environment >::type envir(envirSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type which(whichSEXP);
    Rcpp::traits::input_parameter< double >::type gillRtol(gillRtolSEXP);
    Rcpp::traits::input_parameter< int >::type gillK(gillKSEXP);
    Rcpp::traits::input_parameter< double >::type gillStep(gillStepSEXP);
    Rcpp::traits::input_parameter< double >::type gillFtol(gillFtolSEXP);
    Rcpp::traits::input_parameter< bool >::type optGillF(optGillFSEXP);
    rcpp_result_gen = Rcpp::wrap(nlmixr2Gill83_(what, args, envir, which, gillRtol, gillK, gillStep, gillFtol, optGillF));
    return rcpp_result_gen;
END_RCPP
}
// nlmixr2Eval_
double nlmixr2Eval_(NumericVector theta, std::string md5);
RcppExport SEXP _nlmixr2est_nlmixr2Eval_(SEXP thetaSEXP, SEXP md5SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< std::string >::type md5(md5SEXP);
    rcpp_result_gen = Rcpp::wrap(nlmixr2Eval_(theta, md5));
    return rcpp_result_gen;
END_RCPP
}
// nlmixr2Unscaled_
RObject nlmixr2Unscaled_(NumericVector theta, std::string md5);
RcppExport SEXP _nlmixr2est_nlmixr2Unscaled_(SEXP thetaSEXP, SEXP md5SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< std::string >::type md5(md5SEXP);
    rcpp_result_gen = Rcpp::wrap(nlmixr2Unscaled_(theta, md5));
    return rcpp_result_gen;
END_RCPP
}
// nlmixr2Grad_
NumericVector nlmixr2Grad_(NumericVector theta, std::string md5);
RcppExport SEXP _nlmixr2est_nlmixr2Grad_(SEXP thetaSEXP, SEXP md5SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< std::string >::type md5(md5SEXP);
    rcpp_result_gen = Rcpp::wrap(nlmixr2Grad_(theta, md5));
    return rcpp_result_gen;
END_RCPP
}
// nlmixr2ParHist_
RObject nlmixr2ParHist_(std::string md5);
RcppExport SEXP _nlmixr2est_nlmixr2ParHist_(SEXP md5SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type md5(md5SEXP);
    rcpp_result_gen = Rcpp::wrap(nlmixr2ParHist_(md5));
    return rcpp_result_gen;
END_RCPP
}
// nlmixr2Hess_
RObject nlmixr2Hess_(RObject thetaT, RObject fT, RObject e, RObject gillInfoT);
RcppExport SEXP _nlmixr2est_nlmixr2Hess_(SEXP thetaTSEXP, SEXP fTSEXP, SEXP eSEXP, SEXP gillInfoTSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< RObject >::type thetaT(thetaTSEXP);
    Rcpp::traits::input_parameter< RObject >::type fT(fTSEXP);
    Rcpp::traits::input_parameter< RObject >::type e(eSEXP);
    Rcpp::traits::input_parameter< RObject >::type gillInfoT(gillInfoTSEXP);
    rcpp_result_gen = Rcpp::wrap(nlmixr2Hess_(thetaT, fT, e, gillInfoT));
    return rcpp_result_gen;
END_RCPP
}
// sqrtm
NumericMatrix sqrtm(NumericMatrix m);
RcppExport SEXP _nlmixr2est_sqrtm(SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(sqrtm(m));
    return rcpp_result_gen;
END_RCPP
}
// foceiCalcCov
NumericMatrix foceiCalcCov(Environment e);
RcppExport SEXP _nlmixr2est_foceiCalcCov(SEXP eSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Environment >::type e(eSEXP);
    rcpp_result_gen = Rcpp::wrap(foceiCalcCov(e));
    return rcpp_result_gen;
END_RCPP
}
// foceiFitCpp_
Environment foceiFitCpp_(Environment e);
RcppExport SEXP _nlmixr2est_foceiFitCpp_(SEXP eSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Environment >::type e(eSEXP);
    rcpp_result_gen = Rcpp::wrap(foceiFitCpp_(e));
    return rcpp_result_gen;
END_RCPP
}
// boxCox_
NumericVector boxCox_(NumericVector x, double lambda, int yj);
RcppExport SEXP _nlmixr2est_boxCox_(SEXP xSEXP, SEXP lambdaSEXP, SEXP yjSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type yj(yjSEXP);
    rcpp_result_gen = Rcpp::wrap(boxCox_(x, lambda, yj));
    return rcpp_result_gen;
END_RCPP
}
// iBoxCox_
NumericVector iBoxCox_(NumericVector x, double lambda, int yj);
RcppExport SEXP _nlmixr2est_iBoxCox_(SEXP xSEXP, SEXP lambdaSEXP, SEXP yjSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type yj(yjSEXP);
    rcpp_result_gen = Rcpp::wrap(iBoxCox_(x, lambda, yj));
    return rcpp_result_gen;
END_RCPP
}
// nmNearPD_
RObject nmNearPD_(RObject x, bool keepDiag, bool do2eigen, bool doDykstra, bool only_values, double eig_tol, double conv_tol, double posd_tol, int maxit, bool trace);
RcppExport SEXP _nlmixr2est_nmNearPD_(SEXP xSEXP, SEXP keepDiagSEXP, SEXP do2eigenSEXP, SEXP doDykstraSEXP, SEXP only_valuesSEXP, SEXP eig_tolSEXP, SEXP conv_tolSEXP, SEXP posd_tolSEXP, SEXP maxitSEXP, SEXP traceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< RObject >::type x(xSEXP);
    Rcpp::traits::input_parameter< bool >::type keepDiag(keepDiagSEXP);
    Rcpp::traits::input_parameter< bool >::type do2eigen(do2eigenSEXP);
    Rcpp::traits::input_parameter< bool >::type doDykstra(doDykstraSEXP);
    Rcpp::traits::input_parameter< bool >::type only_values(only_valuesSEXP);
    Rcpp::traits::input_parameter< double >::type eig_tol(eig_tolSEXP);
    Rcpp::traits::input_parameter< double >::type conv_tol(conv_tolSEXP);
    Rcpp::traits::input_parameter< double >::type posd_tol(posd_tolSEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< bool >::type trace(traceSEXP);
    rcpp_result_gen = Rcpp::wrap(nmNearPD_(x, keepDiag, do2eigen, doDykstra, only_values, eig_tol, conv_tol, posd_tol, maxit, trace));
    return rcpp_result_gen;
END_RCPP
}
// lin_cmt_stan
SEXP lin_cmt_stan(Eigen::Map<Eigen::VectorXd> obs_time, Eigen::Map<Eigen::VectorXd> dose_time, Eigen::Map<Eigen::VectorXd> dose, Eigen::Map<Eigen::VectorXd> Tinf, Eigen::Map<Eigen::VectorXd> params, SEXP oralSEXP, SEXP infusionSEXP, SEXP ncmtSEXP, SEXP parameterizationSEXP);
RcppExport SEXP _nlmixr2est_lin_cmt_stan(SEXP obs_timeSEXP, SEXP dose_timeSEXP, SEXP doseSEXP, SEXP TinfSEXP, SEXP paramsSEXP, SEXP oralSEXPSEXP, SEXP infusionSEXPSEXP, SEXP ncmtSEXPSEXP, SEXP parameterizationSEXPSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type obs_time(obs_timeSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type dose_time(dose_timeSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type dose(doseSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type Tinf(TinfSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type oralSEXP(oralSEXPSEXP);
    Rcpp::traits::input_parameter< SEXP >::type infusionSEXP(infusionSEXPSEXP);
    Rcpp::traits::input_parameter< SEXP >::type ncmtSEXP(ncmtSEXPSEXP);
    Rcpp::traits::input_parameter< SEXP >::type parameterizationSEXP(parameterizationSEXPSEXP);
    rcpp_result_gen = Rcpp::wrap(lin_cmt_stan(obs_time, dose_time, dose, Tinf, params, oralSEXP, infusionSEXP, ncmtSEXP, parameterizationSEXP));
    return rcpp_result_gen;
END_RCPP
}
// llik_binomial_c
SEXP llik_binomial_c(Eigen::Map<Eigen::VectorXd> y, Eigen::Map<Eigen::VectorXd> N, Eigen::Map<Eigen::VectorXd> params);
RcppExport SEXP _nlmixr2est_llik_binomial_c(SEXP ySEXP, SEXP NSEXP, SEXP paramsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type N(NSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type params(paramsSEXP);
    rcpp_result_gen = Rcpp::wrap(llik_binomial_c(y, N, params));
    return rcpp_result_gen;
END_RCPP
}
// llik_poisson
SEXP llik_poisson(Eigen::Map<Eigen::VectorXd> y, Eigen::Map<Eigen::VectorXd> params);
RcppExport SEXP _nlmixr2est_llik_poisson(SEXP ySEXP, SEXP paramsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type params(paramsSEXP);
    rcpp_result_gen = Rcpp::wrap(llik_poisson(y, params));
    return rcpp_result_gen;
END_RCPP
}
// llik_normal
SEXP llik_normal(Eigen::Map<Eigen::VectorXd> y, Eigen::Map<Eigen::VectorXd> params);
RcppExport SEXP _nlmixr2est_llik_normal(SEXP ySEXP, SEXP paramsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type params(paramsSEXP);
    rcpp_result_gen = Rcpp::wrap(llik_normal(y, params));
    return rcpp_result_gen;
END_RCPP
}
// llik_betabinomial
SEXP llik_betabinomial(Eigen::Map<Eigen::VectorXd> y, Eigen::Map<Eigen::VectorXd> N, Eigen::Map<Eigen::VectorXd> params);
RcppExport SEXP _nlmixr2est_llik_betabinomial(SEXP ySEXP, SEXP NSEXP, SEXP paramsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type N(NSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type params(paramsSEXP);
    rcpp_result_gen = Rcpp::wrap(llik_betabinomial(y, N, params));
    return rcpp_result_gen;
END_RCPP
}
// llik_student_t
SEXP llik_student_t(Eigen::Map<Eigen::VectorXd> y, Eigen::Map<Eigen::VectorXd> params);
RcppExport SEXP _nlmixr2est_llik_student_t(SEXP ySEXP, SEXP paramsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type params(paramsSEXP);
    rcpp_result_gen = Rcpp::wrap(llik_student_t(y, params));
    return rcpp_result_gen;
END_RCPP
}
// llik_beta
SEXP llik_beta(Eigen::Map<Eigen::VectorXd> y, Eigen::Map<Eigen::VectorXd> params);
RcppExport SEXP _nlmixr2est_llik_beta(SEXP ySEXP, SEXP paramsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type params(paramsSEXP);
    rcpp_result_gen = Rcpp::wrap(llik_beta(y, params));
    return rcpp_result_gen;
END_RCPP
}
// llik_neg_binomial
SEXP llik_neg_binomial(Eigen::Map<Eigen::VectorXd> y, Eigen::Map<Eigen::VectorXd> params);
RcppExport SEXP _nlmixr2est_llik_neg_binomial(SEXP ySEXP, SEXP paramsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type params(paramsSEXP);
    rcpp_result_gen = Rcpp::wrap(llik_neg_binomial(y, params));
    return rcpp_result_gen;
END_RCPP
}
// augPredTrans
RObject augPredTrans(NumericVector& pred, NumericVector& ipred, NumericVector& lambda, RObject& yjIn, NumericVector& low, NumericVector& hi);
RcppExport SEXP _nlmixr2est_augPredTrans(SEXP predSEXP, SEXP ipredSEXP, SEXP lambdaSEXP, SEXP yjInSEXP, SEXP lowSEXP, SEXP hiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type pred(predSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type ipred(ipredSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< RObject& >::type yjIn(yjInSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type low(lowSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type hi(hiSEXP);
    rcpp_result_gen = Rcpp::wrap(augPredTrans(pred, ipred, lambda, yjIn, low, hi));
    return rcpp_result_gen;
END_RCPP
}
// saem_user_opt_ll_fun
double saem_user_opt_ll_fun(NumericVector& resParsNV);
RcppExport SEXP _nlmixr2est_saem_user_opt_ll_fun(SEXP resParsNVSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type resParsNV(resParsNVSEXP);
    rcpp_result_gen = Rcpp::wrap(saem_user_opt_ll_fun(resParsNV));
    return rcpp_result_gen;
END_RCPP
}
// saem_do_pred
SEXP saem_do_pred(SEXP in_phi, SEXP in_evt, SEXP in_opt);
RcppExport SEXP _nlmixr2est_saem_do_pred(SEXP in_phiSEXP, SEXP in_evtSEXP, SEXP in_optSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type in_phi(in_phiSEXP);
    Rcpp::traits::input_parameter< SEXP >::type in_evt(in_evtSEXP);
    Rcpp::traits::input_parameter< SEXP >::type in_opt(in_optSEXP);
    rcpp_result_gen = Rcpp::wrap(saem_do_pred(in_phi, in_evt, in_opt));
    return rcpp_result_gen;
END_RCPP
}
// saem_fit
SEXP saem_fit(SEXP xSEXP);
RcppExport SEXP _nlmixr2est_saem_fit(SEXP xSEXPSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type xSEXP(xSEXPSEXP);
    rcpp_result_gen = Rcpp::wrap(saem_fit(xSEXP));
    return rcpp_result_gen;
END_RCPP
}
// nlmixr2Parameters
List nlmixr2Parameters(NumericVector theta, DataFrame eta);
RcppExport SEXP _nlmixr2est_nlmixr2Parameters(SEXP thetaSEXP, SEXP etaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type eta(etaSEXP);
    rcpp_result_gen = Rcpp::wrap(nlmixr2Parameters(theta, eta));
    return rcpp_result_gen;
END_RCPP
}

// validate (ensure exported C++ functions exist before calling them)
static int _nlmixr2est_RcppExport_validate(const char* sig) { 
    static std::set<std::string> signatures;
    if (signatures.empty()) {
        signatures.insert("List(*nlmixrExpandFdParNlme_)(CharacterVector,CharacterVector)");
    }
    return signatures.find(sig) != signatures.end();
}

// registerCCallable (register entry points for exported C++ functions)
RcppExport SEXP _nlmixr2est_RcppExport_registerCCallable() { 
    R_RegisterCCallable("nlmixr2est", "_nlmixr2est_nlmixrExpandFdParNlme_", (DL_FUNC)_nlmixr2est_nlmixrExpandFdParNlme__try);
    R_RegisterCCallable("nlmixr2est", "_nlmixr2est_RcppExport_validate", (DL_FUNC)_nlmixr2est_RcppExport_validate);
    return R_NilValue;
}
