#define STRICT_R_HEADER
#include "res.h"
#include <boost/algorithm/string.hpp>
#include <string>

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("rxode2", String)
/* replace pkg as appropriate */
#else
#define _(String) (String)
#endif


void calculateDfFull(arma::Col<int>& ID, arma::mat &etas,
                     List &etasDfFull, int &nid, unsigned int &neta) {
  int lastId = ID[ID.size()-1], lastCol = nid-1, lastIndex=ID.size()-1;
  int etaFulli = nid-1;
  double curEta=0.0;
  for (unsigned int j = ID.size(); j--; ){
    if (lastId != ID[j]){
      // Fill in full eta data frame
      for (unsigned int i = neta; i--;){
        curEta = etas(etaFulli, i);//(as<NumericVector>(etasDf1[i]))[etaFulli];
        NumericVector cur = etasDfFull[i];
        std::fill_n(cur.begin()+j+1,lastIndex-j,curEta);
      }
      etaFulli--;
      lastId=ID[j];
      lastIndex=j;
      lastCol--;
      if (lastCol == 0){
        // Finalize ETA
        for (unsigned int i = neta; i--;){
          curEta = etas(0, i);//(as<NumericVector>(etasDf1[i]))[0];
          NumericVector cur = etasDfFull[i];
          std::fill_n(cur.begin(),lastIndex+1,curEta);
        }
        break;
      }
    }
  }
}

int getPredIndex(List &ipredL) {
  CharacterVector names= ipredL.attr("names");
  for (int i = 0; i < names.size(); ++i) {
    if (names[i] == "time") return (i+1);
  }
  return -1;
}

void getLimitFromInput(SEXP limitIn, int& ncalc, arma::vec& limit, int &hasLimit) {
  hasLimit=0;
  if (Rf_isNull(limitIn)) {
    limit = arma::vec(ncalc);
    std::fill_n(limit.begin(), ncalc, R_NegInf);
  } else {
    limit = as<arma::vec>(limitIn);
    hasLimit=1;
  }
}

List getDfIdentifierCols(List &ipred, int &npred, SEXP cmtNames, SEXP idLabels) {
  SEXP cmtVar = PROTECT(getDfSubsetVars(ipred,wrap(CharacterVector::create("CMT","cmt","Cmt"))));
  int extra = 0;
  IntegerVector cmt;
  if (TYPEOF(cmtVar) == VECSXP) {
    extra  = 1;
    cmt = as<IntegerVector>(VECTOR_ELT(cmtVar, 0));
    cmt.attr("class") = "factor";
    cmt.attr("levels") = cmtNames;
  }
  List ret(npred+extra);
  CharacterVector nm(npred+extra);
  CharacterVector nmIn = ipred.names();
  for (int i = 0; i < npred; ++i) {
    nm[i] = boost::to_upper_copy<std::string>(as<std::string>(nmIn[i]));
    ret[i] = ipred[i];
  }
  if (extra) {
    nm[npred] = "CMT";
    ret[npred] = cmt;
  }
  if (TYPEOF(idLabels) == STRSXP) {
    RObject tmp = ret[0];
    if (tmp.sexp_type() == INTSXP) {
      tmp.attr("levels") = idLabels;
      tmp.attr("class") = "factor";
      ret[0] = tmp;
    }
  }
  ret.names()=nm;
  ret.attr("row.names") = ipred.attr("row.names");
  ret.attr("class") = "data.frame";
  UNPROTECT(1);
  return ret;
}

static inline SEXP dfProtectedNames(SEXP inS, std::string what) {
  if (TYPEOF(inS) != VECSXP) return R_NilValue;
  SEXP nmS = PROTECT(Rf_getAttrib(inS, R_NamesSymbol));
  if (Rf_isNull(nmS)) {
    UNPROTECT(1);
    return R_NilValue;
  }
  CharacterVector nm = as<CharacterVector>(nmS);
  const char *badNames[28] = {"IPRED", "IRES", "IWRES", "CENS", "LIMIT",
    "lowerLim", "upperLim","PRED", "RES", "CPRED",
    "CRES", "CWRES", "EPRED", "ERES", "NPDE",
    "ID", "RESETNO", "EVID", "CMT", "SS",
    "RATE", "DUR", "II", "TIME", "rxLambda",
    "rxYj", "rxLow", "rxHi"};
  for (unsigned int i = 0; i < nm.size(); ++i) {
    for (unsigned int j = 0; j < 28; ++j) {
      if (!strcmp(badNames[j], CHAR(nm[i]))) {
        std::string cur = as<std::string>(nm[i]);
        cur += "." + what;
        Rf_warning("change model defined '%s' to '%s' in table (conflicts with reserved names)", CHAR(nm[i]), cur.c_str());
        nm[i] = cur;
      }
    }
  }
  //in.names() = nm;
  Rf_setAttrib(inS, R_NamesSymbol, wrap(nmS));
  UNPROTECT(1);
  return inS;
}

void dfSetStateLhsOps(List& in, List& opt) {
  bool doState=true;
  if (opt.containsElementNamed("state")) {
    RObject tmp = opt["state"];
    if (TYPEOF(tmp) == LGLSXP) {
      doState = as<bool>(tmp);
    }
  }
  bool doLhs = true;
  if (opt.containsElementNamed("lhs")) {
    RObject tmp = opt["lhs"];
    if (TYPEOF(tmp) == LGLSXP) {
      doLhs = as<bool>(tmp);
    }
  }
  bool doEtas = true;
  if (opt.containsElementNamed("eta")) {
    RObject tmp = opt["eta"];
    if (TYPEOF(tmp) == LGLSXP) {
      doEtas = as<bool>(tmp);
    }
  }
  bool doCov = true;
  if (opt.containsElementNamed("covariates")) {
    RObject tmp = opt["covariates"];
    if (TYPEOF(tmp) == LGLSXP) {
      doCov = as<bool>(tmp);
    }
  }
  if (!doEtas) {
    in[1] = R_NilValue;
  } else {
    in[1] = dfProtectedNames(in[1], "etas");
  }
  if (!doState) {
    in[2] = R_NilValue;
  } else {
    in[2] = dfProtectedNames(in[2], "state");
  }
  if (!doLhs) {
    in[3] = R_NilValue;
  } else {
    in[3] = dfProtectedNames(in[3], "lhs");
  }
  if (!doCov) {
    in[4] = R_NilValue;
  } else {
    in[4] = dfProtectedNames(in[4], "cov");
  }
}

extern "C" SEXP _nlmixr2est_resCalc(SEXP ipredPredListSEXP, SEXP omegaMatSEXP,
                                    SEXP etasDfSEXP, SEXP dvIn, SEXP evidIn, SEXP censIn, SEXP limitIn,
                                    SEXP relevantLHSSEXP,  SEXP stateSXP, SEXP covSEXP, SEXP IDlabelSEXP,
                                    SEXP resOpt) {
  BEGIN_RCPP
    List ipredPredList = as<List>(ipredPredListSEXP);
  if (ipredPredList.size() !=3) return R_NilValue;
  List ipredL = ipredPredList[0];
  List predL = ipredPredList[1];
  List etaLst = ipredPredList[2];
  int ncalc = Rf_length(ipredL[0]);
  List etasDf = as<List>(etasDfSEXP);
  int nid = Rf_length(etasDf[0]);
  int npred = getPredIndex(ipredL);
  if (npred == -1) {
    stop(_("malformed dataframes, no time present in ipred data.frame"));
  }

  arma::vec ipredt(REAL(ipredL[npred]), ncalc, false, true);
  arma::vec ipred(ipredt.size());

  arma::vec predt(REAL(predL[npred]), ncalc, false, true);
  arma::vec pred(predt.size());

  arma::vec dv(REAL(dvIn), ncalc, false, true);
  arma::vec dvt(ncalc);

  arma::vec rpv(REAL(predL[npred+1]), ncalc, false, true);
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

  arma::mat omegaMat = as<arma::mat>(omegaMatSEXP);
  unsigned int neta = omegaMat.n_rows;

  bool doSim = true;
  List opt = as<List>(resOpt);
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

  bool interestingLimits = censTruncatedMvnReturnInterestingLimits(dv, dvt, ipred, ipredt, pred, predt, cens, limit,
                                                                   lambda, yj, low, hi, lowerLim, upperLim,
                                                                   riv, doSim, censMethod,
                                                                   normRelated, normIdx, nonNormIdx);


  arma::Col<int> ID(INTEGER(predL[0]), ncalc, false, true);

  arma::mat etas(nid, neta);
  List etasDfFull(neta);
  //etasDfFull.names()=etasDf1.names();
  CharacterVector etaN1 = etasDf.names();
  CharacterVector etaN2(neta);
  for (unsigned int j = neta; j--;) {
    etas.col(j) = arma::vec(REAL(etasDf[j+1]), nid, false, true);
    etaN2[j] = etaN1[j+1];
    etasDfFull[j] = NumericVector(ncalc);
  }
  etasDfFull.names() = etaN2;
  etasDfFull.attr("row.names")=IntegerVector::create(NA_INTEGER,-ncalc);
  etasDfFull.attr("class") = "data.frame";

  calculateDfFull(ID, etas, etasDfFull, nid, neta);

  for (unsigned int i = pred.size(); i--;){
    int cyj, dist, yj0 = (int)yj[i];
    _splitYj(&yj0, &dist,  &cyj);
    if (dist == rxDistributionNorm ||
        dist == rxDistributionT ||
        dist == rxDistributionCauchy ||
        dist == rxDistributionDnorm) {
      pred[i] = _powerDi(predt[i], lambda[i], yj0, low[i], hi[i]);
    }
  }
  arma::vec res = dv - pred;
  res.elem(nonNormIdx).fill(NA_REAL);
  arma::uvec riv0 = find(riv != 0);
  arma::vec iwres=(dvt-ipredt);
  iwres.elem(riv0) /= sqrt(riv.elem(riv0));
  iwres.elem(nonNormIdx).fill(NA_REAL);
  
  arma::vec ires = dv - ipred;
  ires.elem(nonNormIdx).fill(NA_REAL);

  for (unsigned int j = ires.size(); j--; ) {
    if (censMethod == CENS_OMIT && cens[j] != 0) {
      dv[j]	= NA_REAL;
      pred[j]	= NA_REAL;
      res[j]	= NA_REAL;
      ipred[j]	= NA_REAL;
      ires[j]	= NA_REAL;
      iwres[j]	= NA_REAL;
    } else if (evid[j] != 0) {
      dv[j]	= NA_REAL;
      res[j]	= NA_REAL;
      ires[j]	= NA_REAL;
      iwres[j]	= NA_REAL;
    }
  }
  int ncol = 5;
  if (interestingLimits) {
    ncol += 3 + hasLimit;
  }
  List retDF(ncol);
  CharacterVector nm(ncol);
  int i=0;
  nm[i] = "PRED"; retDF[i++] = wrap(pred);
  nm[i] = "RES"; retDF[i++] = wrap(res);
  nm[i] = "IPRED"; retDF[i++] = wrap(ipred);
  nm[i] = "IRES"; retDF[i++] = wrap(ires);
  nm[i] = "IWRES"; retDF[i++] = wrap(iwres);
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
  calcShrinkFinalize(omegaMat, nid, etaLst, iwres, evid, etaN2, 1);

  List retC = List::create(retDF, etasDfFull,
                           getDfSubsetVars(ipredL, stateSXP),
                           getDfSubsetVars(ipredL, relevantLHSSEXP),
                           getDfSubsetVars(ipredL, covSEXP));
  dfSetStateLhsOps(retC, opt);
  retC = dfCbindList(wrap(retC));
  List ret(4);
  ret[0] = wrap(dv);
  ret[1] = getDfIdentifierCols(ipredL, npred, stateSXP, IDlabelSEXP);
  ret[2] = retC;
  ret[3] = etaLst;
  return wrap(ret);
  END_RCPP
    }



extern "C" SEXP _nlmixr2est_popResFinal(SEXP inList) {
  BEGIN_RCPP
    List l = as<List>(inList);
  if (l.size() != 2) return R_NilValue;
  if (Rf_isNull(l[1])) {
    // Only resid in 1
    List l1 = l[0];
    if (l1.size() != 4) return R_NilValue;
    List retC = List::create(l1[1],
                             List::create(_["DV"] = l1[0]),
                             l1[2]);
    return(List::create(_["resid"]=dfCbindList(wrap(retC)),
                        _["shrink"]=l1[3]));
  }
  List l1 = l[0];
  List l2 = l[1];
  List shrinkage;
  List finalLst(3);
  // Regardless of whcih method of CENS imputation, updated DV is in the first element
  NumericVector dv = l1[0];
  List l4;
  if (l1.size() == 2 && l2.size() == 4) {
    l4 = l2;
    l2 = l1;
  } else if (l1.size() == 4 && l2.size() == 2) {
    l4 = l1; // l2 remains the same
  } else {
    return R_NilValue;
  }
  List retC = List::create(l4[1],
                           List::create(_["DV"] = dv),
                           l2[1],
                           l4[2]);
  return List::create(_["resid"]=dfCbindList(wrap(retC)),
                      _["shrink"]=l4[3]);
  END_RCPP  
    }
