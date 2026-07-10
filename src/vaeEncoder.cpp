// vaeEncoder.cpp -- native LSTM encoder for the VAE-NLME estimation method.
// Implements plans/vae-encoder-spec.md: single-layer unidirectional LSTM
// (PyTorch i,f,g,o gate convention) + FC head -> (mu, logSigma, L), with
// reparameterization z = mu + L eps. Forward caches all gate/state activations
// so the analytic backward (LSTM BPTT + FC + reparam) is exact.
//
// The backward consumes the two upstream signals the ELBO provides:
//   gZ             = dLoss/dz                    [N, zDim]
//   gLogSigmaDirect = dLoss/d logSigma (direct)  [N, zDim]  (entropy term q(z|x))
// mu and the strictly-lower entries of L receive gradient only through the
// reparameterization (dz/dmu = I, dz/dL = eps^T), so no separate upstream is
// needed for them.
//
// Validated against the torch autograd oracle (tests/testthat/fixtures/vae/
// encoder_golden.rds) and finite differences to ~1e-6.
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

static inline arma::vec vaeSigmoid(const arma::vec& x) {
  return 1.0 / (1.0 + arma::exp(-x));
}

// [[Rcpp::export]]
List vaeEncoderFwdBwd(arma::cube dataIn,        // [N, Tmax, xDim] standardized inputs
                      IntegerVector lengths,    // [N] valid steps per subject
                      arma::mat covIn,          // [N, nCov] FC-head covariates
                      arma::mat eps,            // [N, zDim] reparam noise
                      arma::mat Wih,            // [4h, xDim]  (rows i,f,g,o)
                      arma::mat Whh,            // [4h, h]
                      arma::vec bih,            // [4h]
                      arma::vec bhh,            // [4h]
                      arma::mat fcW,            // [outDim, h + nCov]
                      arma::vec fcB,            // [outDim]
                      int zDim,
                      arma::mat gZ,             // [N, zDim] dLoss/dz
                      arma::mat gLogSigmaDirect // [N, zDim] direct dLoss/dlogSigma
                      ) {
  const int N = dataIn.n_rows;
  const int h = Whh.n_cols;
  const int xDim = dataIn.n_slices;
  const int nCov = covIn.n_cols;
  const int nOff = zDim * (zDim - 1) / 2;
  const int outDim = 2 * zDim + nOff;

  // strictly-lower (row-major, offset -1) index pairs, matching the spec
  arma::uvec offR(nOff), offC(nOff);
  {
    int idx = 0;
    for (int r = 1; r < zDim; ++r)
      for (int c = 0; c < r; ++c) { offR[idx] = r; offC[idx] = c; ++idx; }
  }

  arma::mat mu(N, zDim), logSigma(N, zDim), zOut(N, zDim);
  arma::cube Lout(zDim, zDim, N, arma::fill::zeros);

  // parameter gradient accumulators
  arma::mat gWih(4 * h, xDim, arma::fill::zeros);
  arma::mat gWhh(4 * h, h, arma::fill::zeros);
  arma::vec gbih(4 * h, arma::fill::zeros);
  arma::vec gbhh(4 * h, arma::fill::zeros);
  arma::mat gFcW(outDim, h + nCov, arma::fill::zeros);
  arma::vec gFcB(outDim, arma::fill::zeros);

  for (int i = 0; i < N; ++i) {
    const int Ti = lengths[i];
    // forward caches for this subject
    arma::mat I(h, Ti), F(h, Ti), G(h, Ti), O(h, Ti), C(h, Ti), H(h, Ti), X(xDim, Ti);
    arma::vec hPrev(h, arma::fill::zeros), cPrev(h, arma::fill::zeros);

    for (int t = 0; t < Ti; ++t) {
      arma::vec xt(xDim);
      for (int d = 0; d < xDim; ++d) xt[d] = dataIn(i, t, d);
      arma::vec pre = Wih * xt + bih + Whh * hPrev + bhh;
      arma::vec it = vaeSigmoid(pre.subvec(0, h - 1));
      arma::vec ft = vaeSigmoid(pre.subvec(h, 2 * h - 1));
      arma::vec gt = arma::tanh(pre.subvec(2 * h, 3 * h - 1));
      arma::vec ot = vaeSigmoid(pre.subvec(3 * h, 4 * h - 1));
      arma::vec ct = ft % cPrev + it % gt;
      arma::vec ht = ot % arma::tanh(ct);
      I.col(t) = it; F.col(t) = ft; G.col(t) = gt; O.col(t) = ot;
      C.col(t) = ct; H.col(t) = ht; X.col(t) = xt;
      hPrev = ht; cPrev = ct;
    }

    arma::vec hLast = H.col(Ti - 1);
    arma::vec combined(h + nCov);
    combined.subvec(0, h - 1) = hLast;
    if (nCov > 0) combined.subvec(h, h + nCov - 1) = covIn.row(i).t();
    arma::vec out = fcW * combined + fcB;

    arma::vec muI = out.subvec(0, zDim - 1);
    arma::vec lsI = out.subvec(zDim, 2 * zDim - 1);
    arma::vec lmask = (nOff > 0) ? out.subvec(2 * zDim, outDim - 1) : arma::vec();

    arma::mat Li(zDim, zDim, arma::fill::zeros);
    Li.diag() = arma::exp(lsI);
    for (int k = 0; k < nOff; ++k) Li(offR[k], offC[k]) = lmask[k];
    arma::vec epsI = eps.row(i).t();
    arma::vec zI = muI + Li * epsI;

    mu.row(i) = muI.t();
    logSigma.row(i) = lsI.t();
    Lout.slice(i) = Li;
    zOut.row(i) = zI.t();

    // ---- backward ----
    arma::vec gZi = gZ.row(i).t();
    arma::vec gMu = gZi;                          // dz/dmu = I
    arma::mat gL = gZi * epsI.t();                // dLoss/dL via z = mu + L eps
    arma::vec gLogSigma = gLogSigmaDirect.row(i).t() + gL.diag() % arma::exp(lsI);
    arma::vec gLmask(nOff);
    for (int k = 0; k < nOff; ++k) gLmask[k] = gL(offR[k], offC[k]);

    arma::vec gOut(outDim);
    gOut.subvec(0, zDim - 1) = gMu;
    gOut.subvec(zDim, 2 * zDim - 1) = gLogSigma;
    if (nOff > 0) gOut.subvec(2 * zDim, outDim - 1) = gLmask;

    gFcB += gOut;
    gFcW += gOut * combined.t();
    arma::vec gCombined = fcW.t() * gOut;
    arma::vec dh = gCombined.subvec(0, h - 1);    // grad on final hidden state
    arma::vec dc(h, arma::fill::zeros);

    for (int t = Ti - 1; t >= 0; --t) {
      arma::vec tanhC = arma::tanh(C.col(t));
      arma::vec doo = dh % tanhC;
      dc += dh % O.col(t) % (1.0 - tanhC % tanhC);
      arma::vec cPrevT = (t > 0) ? C.col(t - 1) : arma::vec(h, arma::fill::zeros);
      arma::vec df = dc % cPrevT;
      arma::vec di = dc % G.col(t);
      arma::vec dg = dc % I.col(t);
      arma::vec dcPrev = dc % F.col(t);

      arma::vec diPre = di % I.col(t) % (1.0 - I.col(t));
      arma::vec dfPre = df % F.col(t) % (1.0 - F.col(t));
      arma::vec dgPre = dg % (1.0 - G.col(t) % G.col(t));
      arma::vec doPre = doo % O.col(t) % (1.0 - O.col(t));

      arma::vec dpre(4 * h);
      dpre.subvec(0, h - 1) = diPre;
      dpre.subvec(h, 2 * h - 1) = dfPre;
      dpre.subvec(2 * h, 3 * h - 1) = dgPre;
      dpre.subvec(3 * h, 4 * h - 1) = doPre;

      arma::vec hPrevT = (t > 0) ? H.col(t - 1) : arma::vec(h, arma::fill::zeros);
      gWih += dpre * X.col(t).t();
      gWhh += dpre * hPrevT.t();
      gbih += dpre;
      gbhh += dpre;

      dh = Whh.t() * dpre;   // grad to h_{t-1}
      dc = dcPrev;
    }
  }

  return List::create(_["mu"] = mu, _["logSigma"] = logSigma, _["L"] = Lout,
                      _["z"] = zOut,
                      _["gWih"] = gWih, _["gWhh"] = gWhh,
                      _["gbih"] = gbih, _["gbhh"] = gbhh,
                      _["gFcW"] = gFcW, _["gFcB"] = gFcB);
}
