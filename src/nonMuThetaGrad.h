#ifndef __NONMUTHETAGRAD_H__
#define __NONMUTHETAGRAD_H__
// Exact-gradient refinement for thetas that are otherwise moved by a
// derivative-free search alone.
//
// saem's refinePhi0Lik (newuoa/nelderMead/optimize), npag and npb's bounded
// bobyqa on the extended-least-squares objective, and vae's regressor step all
// move their non-mu thetas by search with NO derivative information, spending a
// full-population solve on every objective evaluation -- saem's default budget
// is 25 of them per refinement.
//
// nlmixr2 already emits exact symbolic sensitivities.  One solve of the
// theta-sensitivity model (odeSlotThetaSens) yields d(f)/d(theta) for every
// theta at once, through the same pooled, threaded solve path the estimators
// already use, so a Gauss-Newton step is both far cheaper than the search and
// better directed than anything a search recovers from finite budget.
//
// The step is a HEDGE, not a replacement: it moves along the local quadratic
// model, then the existing search runs warm-started from there and corrects
// wherever that model was poor.  Linearity and non-linearity each covered by
// the method suited to it.
//
// This header is the shared arithmetic -- per-observation score and
// information, and the damped step off them.  Each estimator supplies its own
// solve loop and says which objective it is minimizing; none of them owns a
// copy of this.
//
// ---- the schedule each caller owes ----------------------------------------
//
// Because the two halves cost so differently, each estimator needs a schedule
// for BOTH, not one switch covering both:
//
//   * the gradient step is ONE pooled, threaded solve giving every theta's
//     derivative at once -- affordable every iteration/cycle;
//   * the search is many full-population objective evaluations -- saem's
//     nonMuThetaMaxEval is 25, and npag/npb's bobyqa RE-DERIVES the
//     posterior-mean etas for every candidate under muExpand=FALSE, so its
//     evaluations cost more than a solve apiece.
//
// So the cadence inverts what made sense when the search worked alone: run the
// cheap directed step often, the expensive undirected one rarely, and let the
// search start warm from the gradient step rather than from wherever the
// stochastic update left things.
//
//   saem      nonMuThetaGrad / nonMuThetaGradEvery  (gradient)
//             nonMuThetaStart / nonMuThetaEvery / nonMuThetaMaxEval  (search)
//   npag/npb  the gradient step every cycle; residOptimize
//             ("alternate"/"final"/"none") keeps governing the bobyqa search
//   vae       the same pair over its own regress-mode theta set
//
// The objective differs per caller (Gauss for saem, ELS for npag/npb) but the
// SCHEDULE SHAPE does not, which is the other half of why this lives in one
// place.
#include <RcppArmadillo.h>

// Which objective the caller is minimizing.  All three are sums over
// observations of a function of the prediction f (and, for the first two, a
// residual scale that may itself depend on f).
typedef enum {
  // 0.5*((y-f)/g)^2 + log(g), g = add + prop*|f|  -- saem's phi0NormalSSR
  nonMuObjGauss = 0,
  // (y-f)^2/g^2 + log(g^2) -- npag/npb's extended least squares.  Twice the
  // Gaussian, so it shares its derivatives up to that factor.
  nonMuObjEls,
  // f IS the per-observation log-likelihood (saem distribution==4, general
  // likelihood models): the score is just -d(f)/d(theta) and the information
  // is the BHHH outer product.
  nonMuObjLl
} nonMuObjKind;

// Accumulate ONE observation's contribution to the score and information.
//
// `dfdth` is d(f)/d(theta) for the nth thetas being refined -- the exact
// sensitivity outputs, not a finite difference.  `gsd` is the residual SD at
// this observation and `dgsdf` its derivative with respect to f (prop*sign(f)
// for a proportional or combined error model, 0 for additive), which is what
// makes the scale move with the prediction.
//
// `w` weights the contribution.  npag and npb need it: their objective is a
// sum over (subject, SUPPORT POINT) with the support and its weights held fixed
// while the residual/regressor thetas move, so the chain rule carries each
// point's posterior weight -- d(ELS)/d(theta) = sum_ij w_ij * dELS/df * df/dth.
// imp's importance weights enter the same way.  Pass 1.0 for an unweighted
// objective (saem, vae).
//
// Ignores an observation whose scale is not positive or whose sensitivities are
// not finite: one bad row must not poison the whole step.
void nonMuGradAccumObs(nonMuObjKind kind, double y, double f,
                       double gsd, double dgsdf,
                       const double *dfdth, int nth, double w,
                       arma::vec &score, arma::mat &info);

// Damped Newton step from an accumulated score and information.
//
// Gated on CONDITIONING, not on step size.  arma::solve SUCCEEDS on a
// near-singular information matrix and still returns a FINITE step -- just an
// astronomically large one -- so "solved and finite" is not a guard.  imp's
// structural-theta M-step learned this the hard way: one such step took a log
// relative variance from -2.46 to +60.3 in a single iteration.  An earlier fix
// there that gated on the step's magnitude was wrong in the other direction --
// it damped large-but-legitimate steps off perfectly good Hessians and
// perturbed every ordinary fit.
//
// So: take the exact undamped step whenever rcond says the information is
// trustworthy, and otherwise escalate Levenberg-Marquardt damping until the
// step comes inside a trust region relative to each theta's own magnitude.
// Returns false when even heavy damping cannot produce a usable step, in which
// case the caller should leave its thetas alone and let the search do the work.
bool nonMuGradStep(const arma::vec &score, const arma::mat &info,
                   const arma::vec &cur, double trust, arma::vec &step);

#endif // __NONMUTHETAGRAD_H__
