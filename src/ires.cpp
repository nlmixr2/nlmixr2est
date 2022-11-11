#define STRICT_R_HEADER
#include "ires.h"

extern "C" SEXP _nlmixr2est_iresCalc(SEXP ipredDfLstSXP, SEXP dvIn, SEXP evidIn, SEXP censIn, SEXP limitIn,
				 SEXP relevantLHSSEXP, SEXP stateSXP, SEXP covSXP, SEXP IDlabelSEXP,
				 SEXP iresOpt) {
BEGIN_RCPP

  List ipredL = as<List>(ipredDfLstSXP);
  int ncalc = Rf_length(ipredL[0]);

  int npred = getPredIndex(ipredL);

  arma::vec ipredt(REAL(ipredL[npred]), ncalc, false, true);
  arma::vec ipred(ipredt.size());

  arma::vec dv(REAL(dvIn), ncalc, false, true);
  arma::vec dvt(ncalc);

  arma::vec riv(REAL(ipredL[npred+1]), ncalc, false, true);


  arma::Col<int> cens;
  if (Rf_isNull(censIn)) {
    cens = arma::Col<int>(ncalc, fill::zeros);
  } else {
    cens = as<arma::Col<int>>(censIn);
  }
  arma::Col<int> evid;
  if (Rf_isNull(evidIn)) {
    evid = arma::Col<int>(ncalc, fill::zeros);
  } else {
    evid = as<arma::Col<int>>(evidIn);
  }

  arma::vec limit;
  int hasLimit=0;
  getLimitFromInput(limitIn, ncalc, limit, hasLimit);

  arma::vec     hi(REAL(ipredL[ipredL.size()-1]), ncalc, false, true);
  arma::vec    low(REAL(ipredL[ipredL.size()-2]), ncalc, false, true);
  arma::vec     yj(REAL(ipredL[ipredL.size()-3]), ncalc, false, true);
  arma::vec lambda(REAL(ipredL[ipredL.size()-4]), ncalc, false, true);
  arma::vec lowerLim(ncalc);
  arma::vec upperLim(ncalc);

  bool doSim = true;
  List opt = as<List>(iresOpt);
  if (opt.containsElementNamed("doSim")) {
    RObject tmp = opt["doSim"];
    if (TYPEOF(tmp) == LGLSXP) {
      doSim = as<bool>(tmp);
    }
  }
  int censMethod = CENS_TNORM;
  if (opt.containsElementNamed("censMethod")) {
    RObject tmp = opt["censMethod"];
    if (TYPEOF(tmp) == INTSXP) {
      censMethod = as<int>(opt["censMethod"]);
    }
  }

  arma::vec pred(ipred.size());
  arma::vec predt = ipredt;
  
  arma::uvec normRelated(dv.size());
  arma::uvec normIdx;
  arma::uvec nonNormIdx;

  bool interestingLimits =
    censTruncatedMvnReturnInterestingLimits(dv, dvt, ipred, ipredt, pred, predt, cens, limit,
                                            lambda, yj, low, hi, lowerLim, upperLim,
                                            riv, doSim, censMethod,
                                            normRelated, normIdx, nonNormIdx);

  int ncalc2 = sum(normRelated); // This is the true cwres calculations

  arma::Col<int> ID(INTEGER(ipredL[0]), ncalc, false, true);

  arma::vec iwres=(dvt-ipredt);
  uvec riv0 = find(riv!=0);
  iwres.elem(riv0) /= sqrt(riv.elem(riv0));
  iwres.elem(nonNormIdx).fill(NA_REAL);

  arma::vec ires = dv - ipred;
  ires.elem(nonNormIdx).fill(NA_REAL);

  for (unsigned int j = ires.size(); j--; ) {
    if (censMethod == CENS_OMIT && cens[j] != 0) {
      dv[j]	= NA_REAL;
      ipred[j]	= NA_REAL;
      ires[j]	= NA_REAL;
      iwres[j]	= NA_REAL;
    } else if (evid[j] != 0) {
      dv[j]	= NA_REAL;
      ires[j]	= NA_REAL;
      iwres[j]	= NA_REAL;
    }
  }
  int ncol;
  if (ncalc2 == 0) {
    ncol = 1;
  } else {
    ncol = 3;
  }
  if (interestingLimits) {
    ncol += 3 + hasLimit;
  }
  List retDF(ncol);
  CharacterVector nm(ncol);
  int i=0;
  nm[i] = "IPRED"; retDF[i++] = wrap(ipred);
  if (ncalc2 != 0) {
    nm[i] = "IRES"; retDF[i++] = wrap(ires);
    nm[i] = "IWRES"; retDF[i++] = wrap(iwres);
  }
  if (interestingLimits) {
    nm[i] = "CENS"; retDF[i++] = wrap(cens);
    if (hasLimit){
      nm[i] = "LIMIT"; retDF[i++] = wrap(limit);
    }
    nm[i] = "lowerLim"; retDF[i++] = wrap(lowerLim);
    nm[i] = "upperLim"; retDF[i++] = wrap(upperLim);
  }
  retDF.names() = nm;
  retDF.attr("row.names") = IntegerVector::create(NA_INTEGER,-ncalc);
  retDF.attr("class") = "data.frame";
  List retC = List::create(retDF, R_NilValue,
			   getDfSubsetVars(ipredL, stateSXP),
			   getDfSubsetVars(ipredL, relevantLHSSEXP),
			   getDfSubsetVars(ipredL, covSXP));
  dfSetStateLhsOps(retC, opt);
  retC = dfCbindList(wrap(retC));
  List ret(3);
  ret[0] = getDfIdentifierCols(ipredL, npred, stateSXP, IDlabelSEXP);
  ret[1] = List::create(_["DV"]=wrap(dv));
  ret[2] = retC;
  return dfCbindList(wrap(ret));
END_RCPP
}
