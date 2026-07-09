// Importance-sampling EM (IMPMAP / IMP) kernel.
//
// This is the peer of saem.cpp / inner.cpp that holds the NONMEM-style Monte
// Carlo importance-sampling EM.  It reuses the FOCEI inner machinery (mode +
// Hessian + joint-density evaluation) through the thin numeric interface in
// imp.h, so all the model-evaluation heavy lifting stays in inner.cpp.
//
// Module M1: impOuter() performs a single MAP pass and collects each subject's
// mode, individual likelihood, and eta Hessian.  The proposal/sampling/weight
// and EM-update steps are added in later modules.
#include <RcppArmadillo.h>
#include "imp.h"

using namespace Rcpp;

void impOuter(Environment e) {
  // Single MAP pass (reuses the FOCEI posthoc path); leaves each subject's
  // conditional mode in the live inner state.
  impMapPass(e);

  int nsub = impNsub();
  int neta = impNeta();

  arma::mat modes(nsub, neta, arma::fill::zeros);
  arma::vec indLik(nsub, arma::fill::zeros);
  List hess(nsub);

  arma::vec mode(neta);
  arma::mat H(neta, neta, arma::fill::zeros);
  for (int id = 0; id < nsub; ++id) {
    impGetMode(id, mode);
    modes.row(id) = mode.t();
    indLik[id] = impGetIndLik(id);
    if (impGetHessian(id, H)) {
      hess[id] = wrap(H);
    } else {
      hess[id] = R_NilValue;
    }
  }

  e["impEtaMode"] = wrap(modes);
  e["impIndLik"]  = wrap(indLik);
  e["impEtaHess"] = hess;
}
