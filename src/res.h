#ifndef __RES_H__
#define __RES_H__
#if defined(__cplusplus)
#include "armahead.h"
#include "censResid.h"
#include "shrink.h"
#include "utilc.h"
void calculateDfFull(arma::Col<int>& ID, arma::mat &etas,
		     List &etasDfFull, int &nid, unsigned int &neta);

int getPredIndex(List &ipredL);

void getLimitFromInput(SEXP limitIn, int& ncalc, arma::vec& limit, int &hasLimit);

List getDfIdentifierCols(List &ipred, int &npred, SEXP cmtNames, SEXP IDlabel);

void dfSetStateLhsOps(List& in, List& opt);

extern "C" {
#endif
  SEXP _nlmixr2est_resCalc(SEXP ipredPredListSEXP, SEXP omegaMatSEXP,
		       SEXP etasDfSEXP, SEXP dvIn, SEXP evidIn, SEXP censIn, SEXP limitIn,
		       SEXP relevantLHSSEXP, SEXP stateSXP, SEXP covSXP, SEXP IDlabelSEXP,
		       SEXP resOpt);

  SEXP _nlmixr2est_popResFinal(SEXP inList);
#if defined(__cplusplus)
}
#endif
#endif
