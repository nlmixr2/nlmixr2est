#define STRICT_R_HEADER
#include "res.h"
#include "cwres.h"
// Eta sensitivities of the rows `idx` as a sparse neta x nrow matrix: column r
// is row r's gradient.  IOV occasion and mixture-owned etas are zero on most
// rows, and nothing here needs a dense nobs x nobs product.
static inline arma::sp_mat cwresSensT(List &L, int first, unsigned int neta,
                                      const arma::uvec &idx) {
  std::vector<arma::uword> loc;
  std::vector<double> val;
  for (arma::uword r = 0; r < idx.n_elem; ++r) {
    for (unsigned int j = 0; j < neta; ++j) {
      double v = REAL(L[first + j])[idx[r]];
      if (v != 0.0) {
        loc.push_back(j);
        loc.push_back(r);
        val.push_back(v);
      }
    }
  }
  if (val.empty()) return arma::sp_mat(neta, idx.n_elem);
  arma::umat locations(loc.data(), 2, val.size(), false, true);
  arma::vec values(val.data(), val.size(), false, true);
  return arma::sp_mat(locations, values, neta, idx.n_elem, false, false);
}

// diag(F Omega F') from the sparse transposed sensitivities, row by row
static inline arma::vec cwresVfoDiag(const arma::sp_mat &Ft, const arma::mat &omega) {
  arma::vec ret(Ft.n_cols, arma::fill::zeros);
  std::vector<arma::uword> r;
  std::vector<double> v;
  for (arma::uword c = 0; c < Ft.n_cols; ++c) {
    r.clear();
    v.clear();
    for (arma::sp_mat::const_col_iterator it = Ft.begin_col(c); it != Ft.end_col(c); ++it) {
      r.push_back(it.row());
      v.push_back(*it);
    }
    double acc = 0.0;
    for (size_t k = 0; k < r.size(); ++k) {
      double t = 0.0;
      for (size_t j = 0; j < r.size(); ++j) t += v[j] * omega(r[j], r[k]);
      acc += t * v[k];
    }
    ret[c] = acc;
  }
  return ret;
}

// dErr/dEta of each row at its subject's etas, under the population (fppt)
// and individual (fpit) models
static inline void calculateCwresDerr(const arma::sp_mat& fppt, const arma::sp_mat& fpit,
                                      arma::Col<int>& ID, arma::mat &etas,
                                      arma::vec &dErr_dEta_i, arma::vec &dErr_dEta_p,
                                      SEXP etaId) {
  arma::uvec rows = etaRowsById(ID, etaId);
  for (unsigned int c = fppt.n_cols; c--; ) {
    double p = 0.0, i = 0.0;
    for (arma::sp_mat::const_col_iterator it = fppt.begin_col(c); it != fppt.end_col(c); ++it)
      p += (*it) * etas(rows[c], it.row());
    for (arma::sp_mat::const_col_iterator it = fpit.begin_col(c); it != fpit.end_col(c); ++it)
      i += (*it) * etas(rows[c], it.row());
    dErr_dEta_p[c] = p;
    dErr_dEta_i[c] = i;
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
  int nidCol = getPredIndex(ipredL);
  // the eta sensitivities of rx_pred_ follow it, then rx_r_
  int npred = getDfColIndex(ipredL, "rx_pred_");
  int npredP = getDfColIndex(predL, "rx_pred_");

  arma::vec ipredt(REAL(ipredL[npred]), ncalc, false, true);
  arma::vec ipred(ipredt.size());

  arma::vec predt(REAL(predL[npredP]), ncalc, false, true);
  arma::vec pred(predt.size());

  // a copy: the rows below write NA and imputed values into it, and the
  // caller's DV is shared (.calcTables() hands the same predictions to NPDE)
  arma::vec dv(REAL(dvIn), ncalc);
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

  arma::vec rpv(REAL(predL[getDfColIndex(predL, "rx_r_")]), ncalc, false, true);
  arma::vec riv(REAL(ipredL[getDfColIndex(ipredL, "rx_r_")]), ncalc, false, true);

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

    arma::sp_mat fppt = cwresSensT(predL, npredP + 1, neta, normIdx);
    arma::sp_mat fpit = cwresSensT(ipredL, npred + 1, neta, normIdx);

    arma::vec rpvNorm = rpv.elem(normIdx);
    arma::vec rivNorm = riv.elem(normIdx);

    // Only the diagonal of V_fo = F Omega F' (Mentre 2006 p. 352) is used.
    // NONMEM and R/S software differ on WRES: Mentre 2006 uses diag(V_fo) + Sigma,
    // Hooker 2007 diag(V_fo) + diag(dh/deta Sigma dh/deta).  Here Sigma is 1 and
    // dh/deta is sqrt(r), giving r; WRES is under FO (eta=0) and CWRES uses the
    // FOCE V_fo with the FO dh/deta.
    arma::vec Vfop = cwresVfoDiag(fppt, omegaMat);
    arma::vec Vfoi = cwresVfoDiag(fpit, omegaMat);

    arma::vec dErr_dEta_i(ncalc2);
    arma::vec dErr_dEta_p(ncalc2);
    calculateCwresDerr(fppt, fpit, IDnorm, etas, dErr_dEta_i, dErr_dEta_p, etasDf[0]);
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
  calculateDfFull(ID, etas, etasDfFull, neta, etasDf[0]);


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
  List retC = dfTableParts(retDF, etasDfFull, ipredL, ebeL,
                           stateSXP, relevantLHSSEXP, covSXP);
  dfSetStateLhsOps(retC, opt);
  retC = dfCbindList(wrap(retC));
  List ret(4);
  ret[0] = wrap(dv);
  ret[1] = getDfIdentifierCols(ebeL, nidCol, stateSXP, IDlabelSEXP);
  ret[2] = retC;
  ret[3] = etaLst;
  return wrap(ret);
  END_RCPP
    }
