#define STRICT_R_HEADER
#include "res.h"
#include "cwres.h"
static inline void calculateCwresDerr(arma::mat& fppm, arma::mat& fpim,
                                      arma::Col<int>& ID, arma::mat &etas,
                                      arma::vec &dErr_dEta_i, arma::vec &dErr_dEta_p,
                                      int &nid) {
  int lastId = ID[ID.size()-1], lastCol = nid-1, lastIndex=ID.size()-1;
  for (unsigned int j = fppm.n_rows; j--; ){
    if (lastId != ID[j]){
      // FIXME do it without copy?
      arma::vec tmp = fppm.rows(j+1, lastIndex) * trans(etas.row(lastCol));
      std::copy(tmp.begin(),tmp.end(),dErr_dEta_p.begin()+j+1);
      tmp = fpim.rows(j+1, lastIndex) * trans(etas.row(lastCol));
      std::copy(tmp.begin(),tmp.end(),dErr_dEta_i.begin()+j+1);
      lastId=ID[j];
      lastIndex=j;
      lastCol--;
      if (lastCol == 0){
        // Finalize dErr_dEta
        arma::vec tmp = fppm.rows(0, lastIndex) * trans(etas.row(lastCol));
        std::copy(tmp.begin(),tmp.end(),dErr_dEta_p.begin());
        tmp = fpim.rows(0, lastIndex) * trans(etas.row(lastCol));
        std::copy(tmp.begin(),tmp.end(),dErr_dEta_i.begin());
        break;
      }
    }
  }
}

extern "C" SEXP _nlmixr2est_cwresCalc(SEXP ipredPredListSEXP, SEXP omegaMatSEXP,
                                      SEXP etasDfSEXP, SEXP dvIn, SEXP evidIn, SEXP censIn, SEXP limitIn,
                                      SEXP relevantLHSSEXP, SEXP stateSXP, SEXP covSXP, SEXP IDlabelSEXP,
                                      SEXP cwresOpt) {
  BEGIN_RCPP
    List ipredPredList = as<List>(ipredPredListSEXP);
  if (ipredPredList.size() != 4) return R_NilValue; //Rcpp::stop("malformed cwres calc");
  List ipredL = ipredPredList[0];
  List predL  = ipredPredList[1];
  List etaLst = ipredPredList[2];
  List ebeL   = ipredPredList[3];

  int ncalc = Rf_length(ipredL[0]);
  List etasDf = as<List>(etasDfSEXP);
  int nid = Rf_length(etasDf[0]);
  int npred = getPredIndex(ipredL);

  arma::vec ipredt(REAL(ipredL[npred]), ncalc, false, true);
  arma::vec ipred(ipredt.size());

  arma::vec predt(REAL(predL[npred]), ncalc, false, true);
  arma::vec pred(predt.size());

  arma::vec dv(REAL(dvIn), ncalc, false, true);
  arma::vec dvt(ncalc);

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

  arma::mat omegaMat = as<arma::mat>(omegaMatSEXP);
  unsigned int neta = omegaMat.n_rows;

  arma::vec rpv(REAL(predL[npred+1+neta]), ncalc, false, true);
  arma::vec riv(REAL(ipredL[npred+1+neta]), ncalc, false, true);

  bool doSim = true;
  List opt = as<List>(cwresOpt);
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
  arma::uvec normRelated(dv.size());
  arma::uvec normIdx;
  arma::uvec nonNormIdx;
  bool interestingLimits = censTruncatedMvnReturnInterestingLimits(dv, dvt, ipred, ipredt,
                                                                   pred, predt, cens, limit,
                                                                   lambda, yj, low, hi,
                                                                   lowerLim, upperLim,
                                                                   riv, doSim, censMethod,
                                                                   normRelated, normIdx, nonNormIdx);
  int ncalc2 = sum(normRelated); // This is the true cwres calculations
  arma::Col<int> ID(INTEGER(predL[0]), ncalc, false, true);

  arma::mat fppm(ncalc,neta);
  arma::mat fpim(ncalc,neta);

  arma::mat etas(nid, neta);
  List etasDfFull(neta);
  //etasDfFull.names()=etasDf1.names();
  CharacterVector etaN1 = etasDf.names();
  CharacterVector etaN2(neta);
  for (unsigned int j = neta; j--;) {
    fppm.col(j) = arma::vec(REAL(predL[j + 1 + npred]), ncalc, false, true);
    fpim.col(j) = arma::vec(REAL(ipredL[j + 1 + npred]), ncalc, false, true);
    etas.col(j) = arma::vec(REAL(etasDf[j+1]), nid, false, true);
    etaN2[j] = etaN1[j+1];
    etasDfFull[j] = NumericVector(ncalc);
  }
  etasDfFull.names() = etaN2;
  etasDfFull.attr("row.names")=IntegerVector::create(NA_INTEGER,-ncalc);
  etasDfFull.attr("class") = "data.frame";

  arma::vec resFinal(dv.size());
  arma::vec wresFinal(dv.size());
  arma::vec iresFinal(dv.size());
  arma::vec iwresFinal(dv.size());
  arma::vec cpredFinal(dv.size());
  arma::vec cresFinal(dv.size());
  arma::vec cwresFinal(dv.size());

  if (ncalc2 > 0) {
    arma::vec hiNorm = hi.elem(normIdx);
    arma::vec lowNorm = low.elem(normIdx);
    arma::vec yjNorm = yj.elem(normIdx);
    arma::vec lambdaNorm = lambda.elem(normIdx);
    arma::vec dvNorm = dv.elem(normIdx);
    arma::vec predtNorm = predt.elem(normIdx);
    arma::Col<int> IDnorm = ID.elem(normIdx);

    arma::mat fppm2 = fppm.rows(normIdx);
    arma::mat fpim2 = fpim.rows(normIdx);

    arma::vec rpvNorm = rpv.elem(normIdx);
    arma::vec rivNorm = riv.elem(normIdx);

    arma::mat V_fo_p = (fppm2 * omegaMat * fppm2.t()); // From Mentre 2006 p. 352
    arma::mat V_fo_i = (fpim2 * omegaMat * fpim2.t()); // From Mentre 2006 p. 352
    // There seems to be a difference between how NONMEM and R/S types
    // of software calculate WRES.  Mentre 2006 states that the
    // Variance under the FO condition should only be diag(Vfo_full) + Sigma,
    // but Hooker 2007 claims there is a
    // diag(Vfo_full)+diag(dh/deta*Sigma*dh/deta).
    // h = the additional error from the predicted function.
    //
    // In the nlmixr2/FOCEi implemented here, the variance of the err
    // term is 1, or Sigma is a 1 by 1 matrix with one element (1)
    //
    // The dh/deta term would be the sd term, or sqrt(r), which means
    // sqrt(r)*sqrt(r)=|r|.  Since r is positive, this would be r.
    //
    // Also according to Hooker, WRES is calculated under the FO
    // assumption, where eta=0, eps=0 for this r term and Vfo term.
    // However, conditional weighted residuals are calculated under
    // the FOCE condition for the Vfo and the FO conditions for
    // dh/deta
    //

    arma::vec Vfop = V_fo_p.diag();
    arma::vec Vfoi = V_fo_i.diag();

    arma::vec dErr_dEta_i(ncalc2);
    arma::vec dErr_dEta_p(ncalc2);
    calculateCwresDerr(fppm2, fpim2, IDnorm, etas, dErr_dEta_i, dErr_dEta_p, nid);
    arma::vec rest = dvt.elem(normIdx) - predtNorm;
    arma::vec vsum = abs(Vfop+rpvNorm);
    arma::vec wres = rest;
    arma::uvec vsum0 = find(vsum != 0);
    wres.elem(vsum0) /= sqrt(vsum.elem(vsum0));

    arma::vec cpredt = ipredt.elem(normIdx) - dErr_dEta_i;
    arma::vec crest = dvt.elem(normIdx) - cpredt;

    arma::vec res = dvNorm - pred.elem(normIdx);
    vsum = Vfoi+rivNorm;
    vsum0 = find(vsum != 0);
    arma::vec cwres = crest;
    cwres.elem(vsum0) /= sqrt(vsum.elem(vsum0));
    arma::uvec riv0 = find(rivNorm!=0);
    arma::vec iwres=(dvt.elem(normIdx)-ipredt.elem(normIdx));
    iwres.elem(riv0)/=sqrt(rivNorm.elem(riv0));
    arma::vec ires = dv.elem(normIdx) - ipred.elem(normIdx);
    arma::vec cpred(ires.size());
    arma::vec cres(ires.size());
    arma::vec predNorm(ires.size());
    for (unsigned int i = 0; i < cres.size(); ++i) {
      cpred[i] = _powerDi(cpredt[i], lambdaNorm[i], yjNorm[i], lowNorm[i], hiNorm[i]);
      cres[i]  = dvNorm[i] - cpred[i];
      predNorm[i] = _powerDi(predtNorm[i], lambdaNorm[i], yjNorm[i], lowNorm[i], hiNorm[i]);
    }
    resFinal.elem(normIdx)   = res;
    wresFinal.elem(normIdx)  = wres;
    iresFinal.elem(normIdx)  = ires;
    iwresFinal.elem(normIdx) = iwres;
    cpredFinal.elem(normIdx) = cpred;
    cresFinal.elem(normIdx)  = cres;
    cwresFinal.elem(normIdx) = cwres;
    pred.elem(normIdx)  = predNorm;
    // fill rest with na
    resFinal.elem(nonNormIdx).fill(NA_REAL);
    wresFinal.elem(nonNormIdx).fill(NA_REAL);
    iresFinal.elem(nonNormIdx).fill(NA_REAL);
    iwresFinal.elem(nonNormIdx).fill(NA_REAL);
    cpredFinal.elem(nonNormIdx).fill(NA_REAL);
    cresFinal.elem(nonNormIdx).fill(NA_REAL);
    cwresFinal.elem(nonNormIdx).fill(NA_REAL);
  }
  calculateDfFull(ID, etas, etasDfFull, nid, neta);


  for (unsigned int j = dv.size(); j--; ) {
    if (censMethod == CENS_OMIT && cens[j] != 0) {
      dv[j]         = NA_REAL;
      pred[j]       = NA_REAL;
      resFinal[j]	= NA_REAL;
      wresFinal[j]	= NA_REAL;
      ipred[j]      = NA_REAL;
      iresFinal[j]	= NA_REAL;
      iwresFinal[j]	= NA_REAL;
      cpredFinal[j]	= NA_REAL;
      cresFinal[j]	= NA_REAL;
      cwresFinal[j]	= NA_REAL;
    } else if (evid[j] != 0) {
      dv[j]	= NA_REAL;
      resFinal[j]	= NA_REAL;
      wresFinal[j]	= NA_REAL;
      iresFinal[j]	= NA_REAL;
      iwresFinal[j]	= NA_REAL;
      cresFinal[j]	= NA_REAL;
      cwresFinal[j]	= NA_REAL;
    }
  }
  int ncol;
  if (ncalc2 != 0) {
    ncol = 9;
  } else {
    ncol = 2;
  }
  if (interestingLimits) {
    ncol += 3 + hasLimit;
  }
  List retDF(ncol);
  CharacterVector nm(ncol);
  int i=0;
  //nm[i] = "DV"; retDF[i++] = wrap(dv);
  nm[i] = "PRED"; retDF[i++] = wrap(pred);
  if (ncalc2 != 0) {
    nm[i] = "RES"; retDF[i++] = wrap(resFinal);
    nm[i] = "WRES"; retDF[i++] = wrap(wresFinal);
  }
  nm[i] = "IPRED"; retDF[i++] = wrap(ipred);
  if (ncalc2 != 0) {
    nm[i] = "IRES"; retDF[i++] = wrap(iresFinal);
    nm[i] = "IWRES"; retDF[i++] = wrap(iwresFinal);
    nm[i] = "CPRED"; retDF[i++] = wrap(cpredFinal);
    nm[i] = "CRES"; retDF[i++] = wrap(cresFinal);
    nm[i] = "CWRES"; retDF[i++] = wrap(cwresFinal);
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
  calcShrinkFinalize(omegaMat, nid, etaLst, iwresFinal, evid, etaN2, 1);
  List retC = List::create(retDF, etasDfFull, getDfSubsetVars(ipredL, stateSXP),
                           getDfSubsetVars(ebeL, relevantLHSSEXP),
                           getDfSubsetVars(ebeL, covSXP));
  dfSetStateLhsOps(retC, opt);
  retC = dfCbindList(wrap(retC));
  List ret(4);
  ret[0] = wrap(dv);
  ret[1] = getDfIdentifierCols(ebeL, npred, stateSXP, IDlabelSEXP);
  ret[2] = retC;
  ret[3] = etaLst;
  return wrap(ret);
  END_RCPP
    }
