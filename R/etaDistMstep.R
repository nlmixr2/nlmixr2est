## ODE-free M-step for the parameters of a declared random-effect distribution.
##
## In the (y, eta) augmentation the complete-data likelihood factors as
##
##     log p(y | eta)  +  log p(eta | theta_dist)
##
## and theta_dist appears ONLY in the second term.  Its M-step is therefore a
## pure distribution fit to the sampled etas -- no data term, no ODE solve --
## exactly as saem's residual-error step is a fit to accumulated residuals
## rather than a re-solve.
##
## The `rxEtaDistExpand()` rewrite is what breaks that: writing
## `eta = Q(phi(z))` with `z ~ N(0,1)` moves theta_dist out of the prior and
## into the data likelihood, so it becomes a structural parameter needing an
## ODE solve per objective evaluation and lands in `refinePhi0Lik()`'s
## derivative-free search.  Measured on Bauer's gamma data that search converges
## to a STABLE wrong point (CL 2.49 against a truth of 5.03, unchanged across
## 20/60/200 iterations) and a 12-fold larger evaluation budget does not move it
## -- i.e. the formulation, not the optimizer.
##
## EM lets the augmentation be chosen freely: sample in z-space (good MCMC
## geometry, which is the whole point of Bauer's technique) but take the M-step
## in eta-space.  Both are valid EM algorithms converging to the same MLE.
##
## It is NOT circular.  If z were PRIOR draws then `Q(phi(z))` would be exactly
## `family(theta_old)` and the MLE would return theta_old.  The sampled z are
## POSTERIOR draws -- which is precisely why their spread is measured at ~0.94
## rather than 1.0 -- so the implied etas carry data information.

## ===========================================================================
## DESIGN: the general declared-distribution M-step
## ===========================================================================
##
## This is the complete specification.  It is written down before the work so
## the pieces that are NOT yet settled are visible as open questions rather
## than discovered halfway through.
##
## ---------------------------------------------------------------------------
## 1.  Why the step exists
## ---------------------------------------------------------------------------
##
## In the (y, eta) augmentation the complete-data likelihood factors as
##
##     log p(y | eta)  +  log p(eta | theta_dist)
##
## and theta_dist appears ONLY in the second term.  Its M-step is therefore a
## distribution fit to the random effects -- no data term, no ODE solve.
##
## rxEtaDistExpand() breaks that factorization: writing eta = Q(phiU(z)) with
## z ~ N(0,1) moves theta_dist out of the prior and into the data likelihood,
## where it becomes a structural parameter needing a solve per objective
## evaluation, and lands in refinePhi0Lik()'s derivative-free search.  Measured
## on Bauer's gamma data that search converges to a stable WRONG point.
##
## EM lets the augmentation be chosen freely: sample in z-space (good MCMC
## geometry, which is the point of the technique) and take the M-step in
## eta-space.  Both are valid EM algorithms for the same MLE.
##
## ---------------------------------------------------------------------------
## 2.  What is implemented today, and exactly where it stops
## ---------------------------------------------------------------------------
##
##   (a) pool the sampled latents            w_k   (copula-combined if paired)
##   (b) implied random effects              eta = Q(phiU(w); a_old)
##   (c) MLE the family to them              a_new = (shape, rate, ...)
##   (d) invert the argument expressions     a_new -> the user's thetas
##
## Step (d) is a numeric inversion (.etaDistArgsToThetas) that binds only theta
## names.  It assumes ONE population-level native parameter set.  A covariate
## breaks that assumption outright -- every subject has its own args, and with
## a time-varying covariate every observation does -- so the inversion cannot
## represent the model and returns NULL.  Steps (b) and (c) are equally wrong
## in that case: there is no single population `a` to imply etas from or to fit.
##
## Measured on a covariate model: rvCL collapsed 0.135 -> 0.020 and rvV1 to
## 0.005, CL landed at 3.36 against a truth near 5.1.
##
## STATUS.  The peer that section 4 describes was BUILT AND THEN REMOVED, and
## sections 3 and 4 are kept only as the record of why.  It scored the declared
## FAMILY against the sampled etas, which section 3a shows is the wrong
## question; and the covariate problem it existed to solve dissolves once the
## objective is the observation likelihood, because the model recomputes eta
## per record on its own.  `git log` has it if the family term is ever wanted
## as a PENALTY rather than as the objective.
##
## What replaced it is much smaller and lives in R/etaDistPeer.R:
## rxUiGet.etaDistThetaSens() builds the theta sensitivities THROUGH THE ETAS,
## emitting the same lhs names the ordinary construction does so the existing
## gradient step consumes it unchanged.  See section 3b.
##
## ---------------------------------------------------------------------------
## 3a.  THE OBJECTIVE, CORRECTED: the observation likelihood
## ---------------------------------------------------------------------------
##
## Section 3 below states the objective this file was built around -- fit the
## declared family to the sampled etas.  That is a HEURISTIC, not the M-step
## the construction implies, and section 3 is kept because the heuristic is
## what the default route still runs and measures well.  But it is not the
## right answer, and here is why.
##
## The complete data is (y, z), with z the latent standard normal:
##
##     log p(y, z | theta) = log p(y | z, theta) + log p(z)
##
## and log p(z) is theta-FREE.  eta = Q(phiU(z); args(theta)) is a
## deterministic transform of z, not observed data, so the family density
## NEVER APPEARS in the Q-function.  The declared parameters are structural
## parameters of the mean function and belong to the observation likelihood,
## exactly like any other non-mu theta:
##
##     maximize over theta:
##       sum over (subject i, observation j) of
##          log p( y_ij | f(eta_i(theta), ...), sigma )
##       with   eta_i(theta) = Q( phiU(z_i) ; args_i(theta) ),  z_i FIXED
##
## The difference from section 3 is where theta enters.  There the etas are
## frozen and the family's parameters move under them; here the LATENTS are
## frozen and the etas are recomputed at every candidate theta, which is what
## the EM actually holds fixed.
##
## Three consequences:
##
##   * A covariate on a distribution parameter needs NO machinery.  The model
##     already recomputes eta per record from the candidate thetas, so fixed
##     and time-varying covariates are handled by the ordinary solve.  This is
##     what section 4's peer was built to provide, and it obsoletes most of it.
##
##   * There is nothing to invert.  Step (d) of section 2 disappears, for the
##     same reason it does under the peer: thetas are optimized directly.
##
##   * NO SPREAD GUARD IS NEEDED.  A theta that makes the PREDICTIONS worse is
##     rejected whatever the latents look like, so an over-dispersed latent
##     cannot drive a runaway.  That failure mode belongs to fitting the
##     family to the eta sample -- measured: the peer-density route, being
##     unguarded, widened the family to cover etas spanning 1.3-9.9 against a
##     gamma of sd ~1.0 and drove the copula to 0.95 (MARE 29.0% against the
##     guarded MLE route's 18.3%).
##
## WHAT THIS ACTUALLY CHANGES -- corrected after measuring.
##
## The paragraph below said refinePhi0Lik had to be let in because a normal
## model's gate barred it.  That was wrong.  `nonMuTheta` DEFAULTS to
## "regress", so `nonMuThetaRegress` is 1 on essentially every saem fit and
## this refinement ALREADY runs, against exactly this objective.  Traced on
## Bauer's g1:
##
##   [phi0] nphi0=5 nFree=5 free={0 1 2 3 4} obsLikRoute=1 regress=1 dist4=0
##          doFreeze=0 optType=2 thetaSensActive=0
##
## nphi0 is 5 -- the four declared thetas plus rxCor -- because a residual
## parameter lives in ares/bres, not phi0.  So the free set was already the
## declared thetas, the local trust region was already on, and BOTH gate
## changes made no difference: the fits came back byte identical.
##
## So `etaDistLoglik = TRUE` is not a new estimation method.  It is a
## SCHEDULING and OWNERSHIP change to machinery that already runs: it turns the
## family M-step off and hands those columns to the regression unconditionally.
## What is genuinely restricted about that regression is
##
##   nonMuThetaStart    barred until half of (nBurn + nEm) -- the whole
##                      exploratory phase
##   nonMuThetaEvery    thins it further
##   nonMuThetaMaxEval  25 evaluations per firing
##   nonMuThetaOpt      "newuoa", DERIVATIVE-FREE; the exact gradient needs
##                      nonMuThetaOpt="n1qn1" AND nonMuThetaGrad=TRUE
##
## and any claim for this control has to be measured against moving those --
## a full optimization every X from iteration 0 -- or the schedule and the
## heuristic-removal are conflated.  The conflated numbers are not even
## uniformly good: g1 14.5 against 18.3 and g2 30.3 against 34.5, but g3 17.6
## against 5.5.
##
## NOT solve-free, either.  phi0Objective() calls user_fn -- a full population
## solve per evaluation.  `doFreeze` skips re-solving only when
## phi0AffectsOde() is false, and a declared eta drives the structural model.
## The peer of section 4 WAS solve-free; that is what this gives up.
##
## IMPLEMENTATION.  Almost none of this is new code, because saem already has
## the objective: `phi0Objective()` is the observation -log-likelihood at
## candidate phi0 values with the phi1 samples held fixed, and
## `refinePhi0Lik()` optimizes it.  The declared thetas are already IN the
## non-mu set that machinery serves -- `.impmapEstTheta()` returns them
## (struct 1, 2, 5, 6, 8 on Bauer's model; only tq/tv2 are mu-referenced), and
## the theta-sensitivity peer already differentiates the observation
## likelihood through `gammapInv`/`phiU` to give the exact gradient.
##
## (An earlier commit message here claimed .impmapEstTheta() excluded all four
## of Bauer's declaration thetas.  It does not; that was inferred from an
## empty result whose real cause was the derivative failing.)
##
## So `saemControl(etaDistLoglik = TRUE)` now means: hand the declared thetas
## to refinePhi0Lik, and stand the family M-step down for them.  What changes
## is only ownership --
##
##   * refinePhi0Lik runs for a declared-distribution fit (it was gated on
##     distribution == 4 || nonMuThetaRegress, and Bauer's models are prop());
##   * it no longer holds the declared thetas out of its own free list;
##   * the GLS holds them out unconditionally, since refinePhi0Lik owns them;
##   * the family MLE loop is skipped.
##
## The COPULA M-step still runs.  The correlation is a property of the latent
## block rather than of the mean function, and no observation-likelihood term
## identifies it.
##
## ---------------------------------------------------------------------------
## 3b.  How the step is actually taken: through the etas, off ONE solve
## ---------------------------------------------------------------------------
##
## A declared theta reaches the model ONLY through its own random effect, so
##
##   d(state)/d(theta_j) = sum_k d(state)/d(eta_k) * d(eta_k)/d(theta_j)
##
## The first factor is already in `ind->solve` from the solve that produced the
## prediction; the second is pure algebra on Q(phiU(z); args(theta)) with the
## latent FIXED.  So the theta derivative is a READ AND A MULTIPLY against the
## solved buffer, and the sensitivity system scales with the number of declared
## ETAS rather than with the number of parameters the declarations carry.
##
## Measured, 2-state ODE model, 2 declared distributions, 6 estimated non-mu
## thetas: 10 state-sensitivity ODEs the theta way, 4 the eta way, same 6 theta
## columns out of both.  On a linCmt() model there are no ODE states and both
## routes fall through to linCmtB's own parameter sensitivities -- there the
## columns come out BYTE IDENTICAL, which is the check that this is a
## reparameterization of one derivative and not a different one.
##
## Two things that are easy to get wrong here:
##
##   * the differentiation variable is the declared eta's VALUE, not its
##     latent.  The latent is held fixed during the step, so d/d(latent) is
##     identically zero and the gradient silently vanishes.
##   * saem does NOT integrate eta sensitivities by default -- its MCMC is
##     derivative-free -- and what this needs is not focei's full eta
##     sensitivity set either, only the declared (non-normal) etas.  So it is a
##     distinct, smaller model, not a reuse of something already being solved.
##
## ---------------------------------------------------------------------------
## 3.  The general objective (the heuristic route, still the default)
## ---------------------------------------------------------------------------
##
##   maximize over theta:
##     sum over (subject i, record j with evid == 0) of
##        log p_family( eta_ij ; args_ij(theta, covariates_ij) )
##
##   with   eta_ij = Q( phiU(w_i) ; args_ij(theta_old) )
##
## Gated on evid == 0, the way every other likelihood accumulation in saem is.
##
## Three properties:
##
##   * It SUBSUMES the current behaviour.  With no covariate, args_ij does not
##     depend on i or j, every record of a subject contributes the same term,
##     and the maximizer is the pooled fit of (b)+(c) -- reached directly,
##     without the second numerical solve of (d).
##
##   * It handles a TIME-VARYING covariate with no special case, exactly as the
##     normal case does.  There, a time-varying covariate makes each
##     observation's MEAN slightly different and the observation-based
##     likelihood is what gets optimized.  Here it makes each observation's
##     DISTRIBUTION slightly different and the same thing happens.  w_i is the
##     subject's percentile, held fixed, of a distribution whose parameters
##     move -- the analogue of exp(eta_i) being a fixed multiplier on a
##     time-varying typical value.  Nothing requires the covariate to be
##     constant within a subject.
##
##   * The gradient is available.  The argument expressions are symbolic, so
##     d(args)/d(theta) is known and the optimization can be gradient-based
##     rather than another derivative-free search.  d(log p)/d(args) is already
##     supplied per family by rxode2ll's exact derivatives (rxEtaDistGradD).
##
## ---------------------------------------------------------------------------
## 4.  ONE peer model, one solve, swapped in with odeSwap
## ---------------------------------------------------------------------------
##
## NOT by enumerating which symbols in a declaration are not thetas, resolving
## each to an lhs index, and finite-differencing the expressions.  That was the
## first attempt and it is the wrong shape twice over: it makes the R side
## describe the expressions to the C++ side, and it reimplements rxode2's own
## symbolic differentiation as a numeric approximation.
##
## Build ONE rxode2 model covering every declared distribution in the fit, and
## swap it in with odeSwap exactly as odeSlotThetaSens, odeSlotHess2 and
## odeSlotPred already are.  A single solve then produces everything the M-step
## needs, per OBSERVATION RECORD, in one pass.
##
## The model can go all the way to the objective itself.  rxode2 already has
## the families as model functions (rxode2ll's llikGamma() and friends) and
## already differentiates through phiU() and the inverse CDF exactly -- that is
## what the FOCEi path for these models relies on.  So the peer's lhs can be
##
##     rx_edll_              sum over declared families of log p( eta ; args )
##     rx_edll_dtheta_t      its derivative wrt each theta, symbolically
##
## and the M-step reduces to a weighted sum of lhs values over records.  The
## whole chain rule -- d(log p)/d(args) AND d(args)/d(theta) -- happens inside
## rxode2's symbolic engine rather than being assembled from two sources in
## C++.
##
## What that buys:
##
##   * A covariate needs NO handling at all.  The peer is evaluated per record
##     against the data, so a covariate -- fixed or time-varying -- is simply a
##     value the model already has.  No symbol list, no lhs resolution, no
##     branch separating the covariate case from the plain one.
##
##   * The same thetas are produced.  The peer is parameterized by the thetas
##     the fit already carries, so it composes with the existing phi
##     bookkeeping instead of introducing a parallel one.
##
##   * It stays ODE-free.  The peer has no states -- arithmetic on parameters
##     and covariates -- so swap/solve/read costs no integration, exactly like
##     the pred-only model.
##
##   * One solve per iteration for ALL declared distributions, rather than one
##     per family, and the per-record values are saved once and reused for
##     every candidate theta the optimizer tries.
##
##   * It obsoletes most of what the C++ side currently carries for this: the
##     RPN expression evaluator (etaDistExpr.h), rxEtaDistLoglikObj/Grad's
##     arithmetic, and the args->thetas inversion all exist to reconstruct in
##     C++ what rxode2 can emit directly.  Whatever survives should be the
##     accumulation, not the algebra.
##
##   * It reuses machinery that exists and is already threaded: the odeSwap
##     peer pool, odeSwapLhsIndex() for the offsets, and the per-record
##     calc_lhs() walk phi1PredAt() and the theta-sensitivity loop already do.
##
## SPLIT BY FAMILY.  One solve, but the peer emits a separate set of lhs for
## each declared eta/family in the model,
##
##     rx_edll_<k>              log p( eta_k ; args_k )
##     rx_edll_<k>_dtheta_t     its derivative wrt each of THAT family's thetas
##
## and each family is optimized on its own.  Not one summed objective, for
## three reasons:
##
##   * Each declared eta has its own family and its own thetas.  A joint
##     optimization would couple parameters that are not coupled -- the only
##     coupling between declared random effects is the copula, and that is
##     estimated separately and in closed form (section 5).
##
##   * A family can then decline independently.  The spread guard is per
##     family already; one eta whose chain has not settled must not hold back
##     another whose has, and must not drag a shared objective around.
##
##   * Smaller problems are better conditioned.  Each family carries one or two
##     thetas, so n1qn1 works on a 1-2 dimensional problem rather than the
##     concatenation of all of them.
##
## Cost, and what is deliberately NOT being optimized yet.
##
## Two schedules, at different levels, easily confused:
##
##   etaDistEvery   how many SAEM ITERATIONS between M-step firings.  Measured;
##                  20 is best on the models tried.
##   candidates     how many objective evaluations n1qn1 makes WITHIN one
##                  firing.  Each needs the peer at that candidate theta, so
##                  each is a peer solve -- population-wide, though ODE-free.
##
## They multiply: ~5 firings x ~30 candidates x per family.  Individually cheap,
## collectively not obviously so.
##
## Implement re-solve-per-candidate FIRST.  It is the general form -- it makes
## no assumption about the estimation method around it, so the same objective
## serves saem, imp and the focei family, which is the point of putting it in
## C++ at all.  It is also the version whose correctness can be tested in
## isolation, because the optimizer is doing ordinary work on an ordinary
## objective.
##
## The alternative -- one damped Newton step per firing, no inner search, the
## SAEM gain supplying the damping across iterations instead (NONMEM's shape
## for its non-mu thetas, eqs. 1.47-1.52) -- makes the per-candidate cost
## vanish, but it is SAEM-SPECIFIC: it borrows the gain sequence, which imp and
## focei do not have.  So it is a saem-only refinement to be MEASURED against
## the general form, not the thing to build first.  Measure it in saem and see
## whether it loses anything.
##
## Emitting enough to re-form the objective in C++ without re-solving is a
## third option and a false economy: it brings the algebra back to C++, which
## is exactly what this design removes.
##
## 5.  The copula
## ---------------------------------------------------------------------------
##
## Unchanged in shape.  rho is estimated from the correlation of the paired
## COMBINED latents, which is a closed form rather than a search, and is keyed
## off its own flag (etaDistCorMstep) because a model may declare distributions
## without correlating them.  It stays outside the objective above: the copula
## couples the latents, not the marginal families, and the two are separable by
## construction.
##
## ---------------------------------------------------------------------------
## 6.  Guards that must survive the rewrite
## ---------------------------------------------------------------------------
##
##   * The latent spread guard.  The latent is standard normal BY CONSTRUCTION,
##     so a pooled spread far from 1 means the chain has not mixed, not that
##     the family is wrong.  Fitting an unmixed chain is a runaway: it collapses
##     the distribution toward a point mass.  MUST be kept -- it is the only
##     thing standing between this step and that failure.
##
##   * The schedule (etaDistEvery).  The step needs roughly one mixing time
##     between updates; run every iteration it compounds its own output.
##
##   * Audibility.  A step that cannot run must SAY so.  An M-step that
##     silently does nothing is indistinguishable from a converged fit.
##
## ---------------------------------------------------------------------------
## 7.  Open questions -- to settle BEFORE trusting a test
## ---------------------------------------------------------------------------
##
##   Q1  eta_ij is computed at theta_old while the objective varies theta.
##       That is the standard EM two-argument form, but here the SAME
##       expressions appear on both sides, so the fixed point deserves a proof
##       rather than an assumption.  A toy shows it is stable given well-mixed
##       draws and divergent given conditional means or an unmixed chain --
##       that is evidence, not a proof.
##
##   Q2  Weighting.  A subject with 20 observations contributes 20 terms for
##       ONE random effect draw.  With no covariate that is a harmless constant
##       factor; with a time-varying one it is not, and it silently weights
##       subjects by their observation count.  Per-subject averaging, or
##       per-record with an explicit 1/n_i, has to be chosen deliberately.
##
##   Q3  Which records.  evid == 0 excludes doses, but a subject with no
##       observations still has a random effect.  Decide whether it contributes.
##
##   Q4  Identifiability.  A covariate on a distribution parameter is estimated
##       here from the random effects, while the same covariate may also enter
##       the structural model.  Whether both are identified is a modelling
##       question the step cannot answer, but it should not diverge silently
##       when they are not.
##
## ---------------------------------------------------------------------------
## 8.  Test plan
## ---------------------------------------------------------------------------
##
##   T1  No covariate: must reproduce the current inversion's answer to within
##       optimizer tolerance, on all four of Bauer's gamma datasets (relative
##       variance 0.09 / 0.5 / 1.0 / 2.0).
##   T2  No correlation: a single declared distribution, no copula -- the case
##       none of Bauer's datasets exercises, so simulated with known truth.
##   T3  Fixed covariate: recover a known covariate effect on a distribution
##       parameter from simulated data.
##   T4  Time-varying covariate: same, with the covariate varying within
##       subject; must not error and must not be biased by record count (Q2).
##   T5  Degenerate: the covariate effect is truly zero -- must return zero, not
##       drift.
##
## ===========================================================================

#' nlmixr2est's own Nelder-Mead, with two guards its interface needs
#'
#' `nmsimplex()` returns a LIST (`$par`, `$value`), and it derives the simplex
#' step as `-0.2 * start` -- so a coordinate whose starting value is 0 gets a
#' step of 0 and can never move.  Offset such coordinates before handing them
#' over.  Used instead of `stats::optim()` so every optimization in this package
#' goes through the same C simplex (neldermead_wrap -> nelder_fn).
#' @noRd
.etaDistNm <- function(start, fn, maxeval = 2000L, reltol = 1e-10) {
  .s <- as.numeric(start)
  .s[!is.finite(.s) | .s == 0] <- 1e-3
  .r <- tryCatch(nmsimplex(.s, fn, control = list(maxeval = maxeval,
                                                  reltol = reltol)),
                 error = function(e) NULL)
  if (is.null(.r) || is.null(.r$par) || anyNA(.r$par) ||
        !all(is.finite(.r$par))) return(NULL)
  list(par = as.numeric(.r$par), value = as.numeric(.r$value))
}

#' Maximum-likelihood M-step for a declared family, from sampled etas
#'
#' Uses the family's OWN density -- `etaDist` records a `d*()` call, so the
#' log-likelihood is that call with `log=TRUE`.  Family-agnostic, and a theta
#' may enter through an arbitrary expression such as `1/exp(lclrv)`.
#'
#' @param etaVals sampled values of the declared random effect
#' @param distCall the `d*()` call recorded on the iniDf
#' @param thetaNames free thetas to estimate
#' @param start their current values
#' @return named vector of updated thetas, or `NULL` if the fit failed
#' @noRd
.etaDistMstep <- function(etaVals, distCall, thetaNames, start) {
  .cl <- if (is.character(distCall)) str2lang(distCall) else distCall
  .dn <- as.character(.cl[[1]])
  if (!exists(.dn, mode = "function")) return(NULL)
  .e <- etaVals[is.finite(etaVals)]
  if (length(.e) < 2L) return(NULL)
  .nll <- function(p) {
    .tv <- stats::setNames(as.list(p), thetaNames)
    .args <- lapply(as.list(.cl)[-1], function(.a) {
      tryCatch(eval(.a, envir = .tv), error = function(e) NA_real_)
    })
    if (any(!vapply(.args, function(.a) is.numeric(.a) && length(.a) == 1L &&
                      is.finite(.a), logical(1)))) return(1e10)
    .ll <- tryCatch(do.call(.dn, c(list(.e), .args, list(log = TRUE))),
                    error = function(e) NULL)
    if (is.null(.ll) || any(!is.finite(.ll))) return(1e10)
    -sum(.ll)
  }
  ## nlmixr2est's own C simplex, matching _saemOpt()/refinePhi0Lik().  This is
  ## the R-side reference implementation used to validate the M-step; the
  ## in-loop version belongs in saem.cpp calling nelder_fn directly, with the
  ## family density coming from Rmath rather than an R callback on the hot path.
  .x <- .etaDistNm(start, .nll)
  if (is.null(.x) || .x$value >= 1e10) return(NULL)
  stats::setNames(.x$par, thetaNames)
}

#' Closed-form M-step for a Gaussian copula's correlation
#'
#' For a copula block the latent pair is bivariate normal with UNIT variances,
#' so the constrained MLE of the correlation is `sum(w1*w2)/n` -- a closed form,
#' not a search.  This is the same quantity the plan specifies as a simulation
#' CHECK (`cor(qnorm(pgamma(cl)), qnorm(pgamma(v1)))`); it works equally well as
#' an estimator.
#' @noRd
.etaDistCorMstep <- function(w1, w2) {
  .ok <- is.finite(w1) & is.finite(w2)
  if (sum(.ok) < 2L) return(NULL)
  .r <- sum(w1[.ok] * w2[.ok]) / sum(.ok)
  if (!is.finite(.r)) return(NULL)
  max(min(.r, 0.999), -0.999)
}

#' Family code for the C++ distribution M-step
#'
#' The code IS the row number in `lotri::lotriEtaDists()`, so the C++ dispatch
#' (`RXETADIST_*`, src/saem.cpp) and the catalog cannot drift apart -- adding a
#' family to the catalog shifts nothing already assigned.
#'
#' Returns `0L` for anything the C++ dispatch does not implement, which selects
#' the general R fallback (`.etaDistMstep()`).  That fallback evaluates the
#' declaration's own `d*()` call, so a family is never LESS supported than it
#' was; it just pays an R round trip per objective evaluation.
#' @noRd
.etaDistFamilyCode <- function(distCall) {
  .cl <- if (is.character(distCall)) str2lang(distCall) else distCall
  .nm <- as.character(.cl[[1]])
  .tab <- try(lotri::lotriEtaDists(), silent = TRUE)
  if (inherits(.tab, "try-error") || is.null(.tab$name)) return(0L)
  .w <- which(.tab$name == .nm)
  if (length(.w) != 1L) return(0L)
  as.integer(.w)
}

#' Support of a declared family ("real", "positive", "nonneg", "unit")
#'
#' Used to pick the surrogate in [etaDistInit()]: a positive-support family
#' gets a log-normal surrogate, a real-support one an ordinary normal.
#' @noRd
.etaDistSupport <- function(distCall) {
  .cl <- if (is.character(distCall)) str2lang(distCall) else distCall
  .nm <- as.character(.cl[[1]])
  .tab <- try(lotri::lotriEtaDists(), silent = TRUE)
  if (inherits(.tab, "try-error") || is.null(.tab$support)) return(NA_character_)
  .w <- which(.tab$name == .nm)
  if (length(.w) != 1L) return(NA_character_)
  .tab$support[.w]
}

#' Argument ROLES of a declared family, in the family's own argument order
#'
#' The role names what an argument DOES (`shape` vs `rate`), which is what a
#' covariate has to be attached to: a covariate on gamma's shape and one on its
#' rate are different models, and the SAME covariate may legitimately enter both
#' with different shapes.  Roles are unique within a family, so a role is a
#' usable group key -- `dbeta` has `shape1`/`shape2` rather than two `shape`s.
#'
#' Roles stay in R.  The C++ objective is role-agnostic (it maximizes over the
#' family's native parameters whatever they mean), and a role dispatch there
#' would be a second catalog to keep in sync with lotri's.
#'
#' `character(0)` when the family is unknown, or when lotri is too old to carry
#' the column -- callers must treat that as "no role information", not as "no
#' roles", and decline rather than group everything together.
#' @noRd
.etaDistRoles <- function(distCall) {
  .cl <- if (is.character(distCall)) str2lang(distCall) else distCall
  .nm <- as.character(.cl[[1]])
  .tab <- try(lotri::lotriEtaDists(), silent = TRUE)
  if (inherits(.tab, "try-error") || is.null(.tab$roles)) return(character(0))
  .w <- which(.tab$name == .nm)
  if (length(.w) != 1L) return(character(0))
  .r <- strsplit(.tab$roles[.w], ",", fixed = TRUE)[[1]]
  .r[nzchar(.r)]
}

#' Roles that refuse a covariate
#'
#' A SUPPORT endpoint (`dunif`'s bounds, the pareto minimum).  Letting it vary
#' by subject makes the density discontinuous in the parameter -- the likelihood
#' jumps the moment an eta crosses the moving endpoint -- so there is nothing
#' for a derivative-based M-step to follow.
#' @noRd
.etaDistRoleNoCovariate <- c("lower", "upper")

#' Map a declared family's ARGUMENT POSITIONS to their role groups
#'
#' Returns a named list: role -> integer positions in the family's argument
#' order.  With unique roles every group has exactly one member today; the list
#' shape is what lets a family with a genuinely shared role (should one ever be
#' added) group its arguments without changing any caller.
#'
#' `NULL` when the roles are unknown -- see [.etaDistRoles()].
#' @noRd
.etaDistRoleGroups <- function(distCall) {
  .r <- .etaDistRoles(distCall)
  if (length(.r) == 0L) return(NULL)
  split(seq_along(.r), factor(.r, levels = unique(.r)))
}

#' Solve a declared family's thetas so its ARGUMENTS take given values
#'
#' The C++ M-step estimates the family's NATIVE parameters (shape, rate, ...);
#' this turns those back into the user's thetas, which may enter through
#' arbitrary expressions (`shape = 1/exp(lclrv)`).  Called ONCE per iteration,
#' not per objective evaluation -- the expensive part (the likelihood over
#' every sampled eta) stays in C++.
#'
#' Exact when the map is invertible; when it is not, this returns the
#' least-squares closest thetas, which is the best that can be done without
#' constraining what `dist()` accepts.
#' @noRd
.etaDistArgsToThetas <- function(distCall, thetaNames, start, targetArgs) {
  .cl <- if (is.character(distCall)) str2lang(distCall) else distCall
  .ex <- as.list(.cl)[-1]
  if (length(.ex) != length(targetArgs)) return(NULL)
  .obj <- function(p) {
    .tv <- stats::setNames(as.list(p), thetaNames)
    .a <- vapply(.ex, function(.e) {
      .v <- tryCatch(eval(.e, envir = .tv), error = function(e) NA_real_)
      if (!is.numeric(.v) || length(.v) != 1L) NA_real_ else as.numeric(.v)
    }, numeric(1))
    if (anyNA(.a) || any(!is.finite(.a))) return(1e10)
    ## relative, so arguments on very different scales weigh comparably
    sum((log(pmax(abs(.a), 1e-300)) - log(pmax(abs(targetArgs), 1e-300)))^2) +
      sum((sign(.a) != sign(targetArgs)) * 1e3)
  }
  .x <- .etaDistNm(start, .obj)
  if (is.null(.x) || .x$value > 1e-6) return(NULL)
  stats::setNames(.x$par, thetaNames)
}

#' Metadata the C++ distribution M-step needs
#'
#' Returns the phi column of each declared eta's own latent normal, its family
#' code, which declared eta it is copula-correlated with, its current NATIVE
#' parameters and correlation -- plus the bookkeeping R needs to map the
#' updated parameters back onto thetas.
#' @noRd
.etaDistMstepInfo <- function(ui, etaTrans, etaNames, paramsToEstimate = NULL) {
  .c <- .etaDistMstepCore(ui)
  if (is.null(.c)) return(NULL)
  ## the expansion renames the declared eta's latent `rxz.<eta>`; saem indexes
  ## its phi columns through etaTrans rather than by eta number
  .direct <- .etaDistIsDirect(ui)
  .pfx <- .etaDistEtaPrefix(ui)
  .lat <- vapply(paste0(.pfx, .c$etas), function(.z) {
    .w <- which(etaNames == .z)
    if (length(.w) == 1L) as.integer(etaTrans[.w]) - 1L else -1L
  }, integer(1))
  if (any(.lat < 0)) return(NULL)                    # -> R fallback
  ## Each declared theta's PHI index.  saem's phi columns are the model
  ## parameters in `saemParamsToEstimate` order; a declared-distribution theta
  ## has no eta, so its phi column carries no random effect and lands in phi0,
  ## where .configsaem() finishes the mapping (phi index -> phi0 column).
  ## Match against the PHI parameters, which `saemParamsToEstimate` is NOT.
  ##
  ## That accessor INTERLEAVES each theta with its mu-referenced covariate
  ## coefficients (rxUiGet.saemParamsToEstimate builds [theta, cov1, ...] rows
  ## and flattens them), because it indexes MCOV.  A coefficient is carried by
  ## the COV/MCOV machinery and is NOT a phi parameter: with one covariate the
  ## list has 8 entries while covstruct stays 7x7 and nphi is 7.
  ##
  ## The indices built here are consumed as PHI indices (saem_fit.R matches them
  ## against i0, the phi columns with no random effect), so matching against the
  ## longer list puts every parameter at or after the coefficient one too high.
  ## Measured on a covariate model with deliberately distinct starting values:
  ## declaration 1's thetas (lclrv, lclm, bWT) mapped to phi0 columns (3, 0, 1),
  ## which are (lv1rv, lclm, lv1m) -- lclrv was never written, lv1m was
  ## clobbered, and the reported table showed lclrv holding lv1m's start value
  ## while only the residual error moved.  A range check cannot catch it: the
  ## indices stay within nphi, they just mean something else.
  ##
  ## Dropping the coefficients reproduces exactly the list a no-covariate model
  ## has.  A declaration whose theta IS a coefficient then fails to match and
  ## the step declines with a warning, which is right -- a coefficient has no
  ## phi0 column to be written back through, and saem's covariate machinery
  ## owns it.
  if (!is.null(paramsToEstimate)) {
    .muCov <- tryCatch(rxode2::rxUiDecompress(ui)$saemMuRefCovariateDataFrame,
                       error = function(e) NULL)
    if (!is.null(.muCov) && length(.muCov$covariateParameter) > 0L) {
      paramsToEstimate <- paramsToEstimate[
        !(paramsToEstimate %in% .muCov$covariateParameter)]
    }
  }
  .tp <- vector("list", .c$n)
  if (!is.null(paramsToEstimate)) {
    for (.i in seq_len(.c$n)) {
      .m <- match(.c$thetas[[.i]], paramsToEstimate)
      if (anyNA(.m)) return(NULL)        # not a plain phi theta -> fallback
      .tp[[.i]] <- as.integer(.m)
    }
  }
  ## the copula correlation is itself a theta and gets the same treatment,
  ## matched by the name the core resolved rather than by grep()ing rxCor.*
  .rn <- character(0); .rp <- integer(0)
  ## ...and PER FAMILY, so a model with more than one correlated pair is
  ## representable.  `.rp` is indexed by `.ck` (the families that have a
  ## partner); `.rpf` is indexed by family, NA where there is no partner, which
  ## is what lets each correlation keep its OWN phi0 column downstream instead
  ## of every pair sharing one.
  .rpf <- rep(NA_integer_, .c$n)
  .ck <- which(.c$corWith >= 0L)
  ## On the DIRECT route there is no `rxCor.*` theta to find: the expansion
  ## leaves the correlation in the omega, and `.etaDistMstepCore()` reads it
  ## from there.  Looking for the theta anyway returns NULL from this whole
  ## function, which downstream is indistinguishable from "no declarations" --
  ## so a correlated declared PAIR on the direct route was refused outright by
  ## the guard in saem.R, reporting that the distributions could not be
  ## resolved when in fact only their correlation's OWNER could not.
  ##
  ## Leaving the map empty is the right answer rather than a workaround: the
  ## copula M-step writes its update back through a phi0 column, and with no
  ## theta there is no column to write to.  The correlation is therefore used
  ## (the prior's copula term reads etaDistRho) but not ESTIMATED on this route
  ## -- .etaDistWarnCorFrozen() says so, because a rho sitting at its ini()
  ## value is otherwise indistinguishable from a converged one.
  if (.direct && length(.ck) > 0L) {
    .etaDistWarnCorFrozen(.c$etas[.ck])
  } else if (!is.null(paramsToEstimate) && length(.ck) > 0L) {
    .rc <- .c$corTheta[.ck]
    .m <- match(.rc, paramsToEstimate)
    if (anyNA(.m)) return(NULL)
    .rn <- .rc; .rp <- as.integer(.m)
    .rpf[.ck] <- as.integer(.m)
  }
  ## The ARGUMENT expressions and their theta names, so the C++ M-step can do
  ## the native-parameters -> thetas map itself instead of calling back into R.
  ## Deparsed here because this is where the dist() call is already parsed; the
  ## C++ side parses them once into RPN and declines anything outside its
  ## grammar, in which case the `map` closure below is still used.
  .exprs <- lapply(seq_len(.c$n), function(.k) {
    .cl <- .c$dist[[.k]]
    if (is.character(.cl)) .cl <- str2lang(.cl)
    vapply(as.list(.cl)[-1], function(.e) paste(deparse(.e), collapse = ""),
           character(1))
  })
  ## Follow saem's mu2 covariate rewrite into these expressions.
  ##
  ## mu2 hoists a covariate term out of the model into a generated column and
  ## DROPS the original covariate from the dataset: `bWT * log(WT/70)` becomes
  ## `nlmixrMuDerCov1 * bWT`, and data$data then holds nlmixrMuDerCov1 (which
  ## carries log(WT/70), verified to 15 digits) with no WT column at all.  The
  ## declaration's argument text is captured from the DECLARATION, so without
  ## this it still says log(WT/70) -- the expression the M-step parses and the
  ## data saem holds no longer share a symbol, the parse declines on the unknown
  ## name, and the covariate is silently invisible to the step.
  ##
  ## Rewriting here rather than teaching the M-step about mu2: the substitution
  ## is already recorded, the generated column already carries the transformed
  ## value, and doing it in one place keeps the two descriptions of the same
  ## model from drifting.
  .mu2 <- tryCatch(rxode2::rxUiDecompress(ui)$mu2RefCovariateReplaceDataFrame,
                   error = function(e) NULL)
  ## The table is EMPTY on the ui this receives -- the etaDist expansion
  ## rebuilds the ui and it does not survive (measured: 1 row while mu2 hoists,
  ## 0 rows here).  `.etaDistMu2Record()` puts the same substitution on
  ## `ui$meta`, which does survive, so fall back to that.  Without it the
  ## rewrite below silently never fires and the declaration keeps naming a
  ## covariate the fit data no longer has.
  if (is.null(.mu2) || NROW(.mu2) == 0L) {
    .st2 <- .etaDistMu2Get(ui)
    if (!is.null(.st2) && NROW(.st2) > 0L) {
      .mu2 <- data.frame(covariateParameter = .st2$covariateParameter,
                         modelExpression = .st2$modelExpression,
                         derived = .st2$derived,
                         stringsAsFactors = FALSE)
    }
  }
  .covK <- .c$cov
  if (!is.null(.mu2) && NROW(.mu2) > 0L) {
    .sub <- function(.txt) {
      for (.i in seq_len(NROW(.mu2))) {
        .o <- .mu2$modelExpression[.i]
        .n <- paste0(if (!is.null(.mu2$derived)) .mu2$derived[.i]
                     else paste0("nlmixrMuDerCov", .i),
                     " * ", .mu2$covariateParameter[.i])
        .txt <- gsub(.o, .n, .txt, fixed = TRUE)
        ## deparse() spacing need not match the recorded text exactly, so fall
        ## back to a whitespace-insensitive match rather than silently missing it
        if (!grepl(.n, .txt, fixed = TRUE)) {
          .txt <- gsub(gsub("[[:space:]]+", "", .o), .n,
                       .txt, fixed = TRUE)
        }
      }
      .txt
    }
    .exprs <- lapply(.exprs, function(.e) vapply(.e, .sub, character(1),
                                                 USE.NAMES = FALSE))
    ## Symbols come from the SUBSTITUTED text -- nlmixrMuDerCov1 is a real data
    ## column where WT no longer is.
    .thAll <- .c$iniDf$name[!is.na(.c$iniDf$ntheta)]
    .covK <- lapply(.exprs, function(.e) {
      setdiff(unique(unlist(lapply(.e, function(.z) all.vars(str2lang(.z))))),
              .thAll)
    })
  }
  ## Which of each declaration's thetas the PRIOR identifies rather than the
  ## data -- NoLimits' Q1/Q2 partition, see .etaDistThetaSplit().  Per
  ## declaration and per theta slot, so the C++ side can hold exactly those
  ## columns out of the observation-likelihood steps and give them to Q2.
  .split <- .etaDistThetaSplit(ui, unique(unlist(.c$thetas)))
  .q2 <- lapply(.c$thetas, function(.nm) as.integer(.nm %in% .split$q2))
  list(latent = .lat, fam = .c$fam, corWith = .c$corWith,
       direct = as.integer(.etaDistIsDirect(ui)),
       q2 = .q2, q2Names = .split$q2,
       usable = as.integer(.c$usable), cov = .covK,
       exprs = .exprs, exprThetas = .c$thetas,
       args = .c$args, rho = .c$rho,
       dist = .c$dist, thetas = .c$thetas, thetaPhi = .tp,
       corName = .rn, corPhi = .rp, corPhiByFam = .rpf, etas = .c$etas)
}

#' Read the declaration stash `.preProcessEtaDist()` left on the ui
#'
#' `NULL` when there is none -- an unexpanded ui, or a model with no
#' declaration.  Tolerant of a ui that is not an environment, so a caller that
#' hands in something unexpected gets the fallback rather than an error.
#'
#' @param ui decompressed rxode2 ui
#' @return the stash, or `NULL`
#' @noRd
.etaDistDeclGet <- function(ui) {
  tryCatch({
    .m <- rxode2::rxUiDecompress(ui)$meta
    if (!is.environment(.m)) return(NULL)
    if (!exists(".etaDistDecl", envir = .m, inherits = FALSE)) return(NULL)
    get(".etaDistDecl", envir = .m, inherits = FALSE)
  }, error = function(e) NULL)
}

#' Split the declared thetas into the ones the data identify and the ones the
#' PRIOR identifies
#'
#' NoLimits.jl's Q1/Q2 partition (`_partition_q1_q2_names`,
#' `src/estimation/common.jl:4557`), applied to our model text:
#'
#'   q2_candidates = setdiff(re_fe_syms, obs_fe)
#'
#' A parameter belongs to Q2 when it appears in a random-effect distribution
#' expression and in NO observation-side block.  The complete-data likelihood
#' `E[log p(y|eta,theta)] + E[log p(eta|theta)]` is then separable in it: it
#' enters only the second term, so its M-step needs no ODE and no observation
#' likelihood at all.
#'
#' Here the RE distribution expressions are the `rxEdA.<eta>.<role>` anchors the
#' expansion emits, so the test becomes: is this anchor READ by anything?
#'
#'   cdf     the decoder reads it (`gammapInv(rxEdA..., phiU(rxz...))/rxEdA...`)
#'           -> the theta is in the observation path        -> Q1
#'   direct  nothing reads it, the eta IS the random effect
#'           -> the theta parameterizes the prior and nothing else -> Q2
#'
#' One rule, no route-specific branching, and it is CORRECT rather than merely
#' convenient on both: on the cdf route the latent is a fixed N(0,1), so
#' `log p(z)` is theta-free and there is no prior term to maximize -- Q1 is the
#' whole story.  Measured, using the eta-prior objective on a cdf model instead
#' costs MARE 2.34% -> 17.57%, because the sample it scores was produced by
#' decoding with the current theta and the objective has a fixed point there.
#'
#' A theta appearing BOTH in an anchor and in the observation model is Q1.  That
#' is the non-separable case, and NoLimits does the same -- it empties the Q2 set
#' entirely when `extra_objective` couples the two (`saem.jl:3204-3211`).
#'
#' @param ui rxode2 ui, already expanded
#' @param thetas character vector of the declared thetas to classify; when
#'   missing every declared theta the stash names is classified
#' @return a list with `q1` and `q2`, character vectors
#' @noRd
#' @author Matthew L. Fidler
.etaDistThetaSplit <- function(ui, thetas = NULL) {
  .empty <- list(q1 = character(0), q2 = character(0))
  .ui <- tryCatch(rxode2::rxUiDecompress(ui), error = function(e) NULL)
  if (is.null(.ui)) return(.empty)
  .expr <- tryCatch(.ui$lstExpr, error = function(e) NULL)
  if (!is.list(.expr) || length(.expr) == 0L) return(.empty)
  .lhsOf <- function(.e) {
    if (is.call(.e) && length(.e) >= 3L && is.name(.e[[2]]) &&
          (identical(.e[[1]], quote(`<-`)) || identical(.e[[1]], quote(`=`)))) {
      as.character(.e[[2]])
    } else {
      NA_character_
    }
  }
  .rhsVars <- function(.e) {
    if (is.call(.e) && length(.e) >= 3L &&
          (identical(.e[[1]], quote(`<-`)) || identical(.e[[1]], quote(`=`)))) {
      all.vars(.e[[3]])
    } else {
      all.vars(.e)
    }
  }
  .lhs <- vapply(.expr, .lhsOf, character(1))
  .rhs <- lapply(.expr, .rhsVars)
  .isAnchor <- !is.na(.lhs) & grepl("^rxEdA[.]", .lhs)
  if (!any(.isAnchor)) return(.empty)
  ## every symbol read anywhere on a NON-anchor line is observation-side
  .obsSide <- unique(unlist(.rhs[!.isAnchor]))
  ## an anchor that nothing reads contributes no observation-side dependence
  .anchorRead <- vapply(.lhs[.isAnchor], function(.a) {
    any(vapply(.rhs, function(.v) .a %in% .v, logical(1)))
  }, logical(1))
  ## thetas feeding a READ anchor are in the observation path through it
  .viaRead <- unique(unlist(.rhs[.isAnchor][.anchorRead]))
  .viaDead <- unique(unlist(.rhs[.isAnchor][!.anchorRead]))
  if (is.null(thetas)) {
    .st <- .etaDistDeclGet(.ui)
    .all <- if (is.null(.st)) character(0) else {
      .tn <- .ui$iniDf$name[!is.na(.ui$iniDf$ntheta)]
      unique(unlist(lapply(.st$etaDist, function(.d)
        intersect(all.vars(str2lang(.d)), .tn))))
    }
    thetas <- .all
  }
  if (length(thetas) == 0L) return(.empty)
  ## Q2: reaches the model ONLY through an anchor nothing reads, and appears
  ## nowhere observation-side in its own right.
  .q2 <- thetas[thetas %in% .viaDead & !(thetas %in% .viaRead) &
                  !(thetas %in% .obsSide)]
  list(q1 = setdiff(thetas, .q2), q2 = .q2)
}

#' Say that a declared copula correlation is used but not estimated
#'
#' Only on the direct route, and only once per fit.  The rho reaches the prior
#' (the copula term in `rxEtaDistPairLogD`) but nothing updates it, because the
#' update is written back through a phi0 column and this route has no `rxCor.*`
#' theta to own one.  A frozen rho reported without comment is indistinguishable
#' from a converged one.
#'
#' @param etas the declared random effects whose correlation is affected
#' @return nothing, called for the message
#' @noRd
#' @author Matthew L. Fidler
.etaDistWarnCorFrozen <- function(etas) {
  message("the declared copula correlation for '", paste(etas, collapse="', '"),
          "' is used by the prior but held at its ini() value\n",
          "  etaDistParam=\"direct\" keeps the correlation in the omega rather ",
          "than in an rxCor.* theta, and the copula update writes back through ",
          "a theta\n",
          "  use etaDistParam=\"cdf\" to estimate it")
  invisible()
}

#' Which parameterization the expansion used
#'
#' `rxEtaDistExpand()` records it on `etaDistInfo$param`.  It is absent on a ui
#' expanded before the option existed, and on an UNEXPANDED one -- both mean the
#' CDF construction, which is the default and was the only route.
#'
#' This is read rather than the control, deliberately: the control says what was
#' ASKED for and the ui says what was BUILT, and an estimator that samples a
#' direct prior against a cdf-expanded model (or the reverse) is fitting the
#' wrong model silently.  Only the second question has a right answer here.
#'
#' @param ui rxode2 ui
#' @return `TRUE` when the model was expanded on the direct route
#' @noRd
#' @author Matthew L. Fidler
.etaDistIsDirect <- function(ui) {
  ## The STASH first.  `etaDistInfo` is where the expansion records it, but that
  ## is a ui environment variable and it does not survive to the estimator --
  ## saem.R read FALSE off a ui whose eta was already renamed `rxd.eta.cl`.  The
  ## stash lives in `ui$meta` and does survive.
  .st <- .etaDistDeclGet(ui)
  if (!is.null(.st) && !is.null(.st$param)) return(identical(.st$param, "direct"))
  tryCatch({
    .u <- rxode2::rxUiDecompress(ui)
    if (!exists("etaDistInfo", envir = .u, inherits = FALSE)) return(FALSE)
    .i <- get("etaDistInfo", envir = .u, inherits = FALSE)
    identical(.i$param, "direct")
  }, error = function(e) FALSE)
}

#' The eta-name prefix the expansion gave the declared random effect
#'
#' `rxz.` on the cdf route (the latent normal), `rxd.` on the direct route (the
#' declared eta itself).  Both are non-mu-referenced and both are renamed, so
#' the only thing an index map needs is which prefix to look for.
#'
#' @param ui rxode2 ui
#' @return the prefix, with its trailing dot
#' @noRd
#' @author Matthew L. Fidler
.etaDistEtaPrefix <- function(ui) {
  if (.etaDistIsDirect(ui)) "rxd." else "rxz."
}

#' Write the declaration stash where it will survive
#'
#' @param ui decompressed rxode2 ui
#' @param value the stash
#' @return `TRUE` if it was written
#' @noRd
.etaDistDeclSet <- function(ui, value) {
  tryCatch({
    .m <- ui$meta
    if (!is.environment(.m)) return(FALSE)
    assign(".etaDistDecl", value, envir = .m)
    TRUE
  }, error = function(e) FALSE)
}

#' Estimator-independent half of the declared-distribution M-step metadata
#'
#' Everything the M-step needs that does not depend on how a particular
#' estimator numbers its random effects and thetas: the family code, the current
#' native parameters, the copula pairing and each declaration's own thetas.
#' `.etaDistMstepInfo()` (saem) and `.etaDistMstepInfoFocei()` (imp/impmap) add
#' their own index maps on top of this.
#'
#' @param ui rxode2 ui, already expanded by `rxEtaDistExpand()`
#' @return a list, or `NULL` when the model carries no usable declaration
#' @noRd
.etaDistMstepCore <- function(ui) {
  .ui <- rxode2::rxUiDecompress(ui)
  .ini <- .ui$iniDf
  ## After rxEtaDistExpand() the declarations are gone from the iniDf, so the
  ## stash .preProcessEtaDist() took before expanding is the only record.  Fall
  ## back to reading the ui directly for an UNEXPANDED one (which is what the
  ## tests and etaDistInit() hand in).
  .st <- .etaDistDeclGet(.ui)
  if (is.null(.st)) {
    .d <- rxode2::rxUiEtaDists(.ui)
    if (nrow(.d) == 0L) return(NULL)
    .st <- .etaDistDeclStash(.ui, .d,
                             param = if (.etaDistIsDirect(.ui)) "direct" else "cdf")
    if (is.null(.st)) return(NULL)
  }
  .n <- length(.st$name)
  if (.n == 0L) return(NULL)
  .thNames <- .ini$name[!is.na(.ini$ntheta)]
  .thVals <- stats::setNames(as.list(.ini$est[!is.na(.ini$ntheta)]), .thNames)
  .fam <- integer(.n); .maxA <- 0L
  .args <- vector("list", .n); .tn <- vector("list", .n)
  .hasCov <- logical(.n)
  ## The covariate SYMBOLS each declaration reads, in the order all.vars() finds
  ## them.  Kept rather than just counted: the C++ argument parser resolves
  ## symbols against a flat name list and rxEtaDistLoglikObj() lays out
  ## vals[0..nth) thetas then vals[nth..nth+nSym) this record's symbols, so
  ## `c(thetas, cov)` is exactly the vector it needs.
  .cov <- vector("list", .n)
  for (.i in seq_len(.n)) {
    .fam[.i] <- .etaDistFamilyCode(.st$etaDist[.i])
    .cl <- str2lang(.st$etaDist[.i])
    ## A name in a declaration argument that is not a theta is DATA -- a
    ## covariate column.  That is the only thing separating "has no single
    ## population value to fit" from "evaluates to nonsense", and the two want
    ## opposite handling: the first is a supported model whose family MLE has
    ## to stand down for THIS declaration, the second is a broken declaration.
    .cov[[.i]] <- setdiff(all.vars(.cl), .thNames)
    .hasCov[.i] <- length(.cov[[.i]]) > 0L
    .a <- vapply(as.list(.cl)[-1], function(.x) {
      .v <- tryCatch(eval(.x, envir = .thVals), error = function(e) NA_real_)
      if (!is.numeric(.v) || length(.v) != 1L) NA_real_ else as.numeric(.v)
    }, numeric(1))
    .args[[.i]] <- .a
    .maxA <- max(.maxA, length(.a))
    .tn[[.i]] <- intersect(all.vars(.cl), .thNames)
  }
  ## Per DECLARATION, not per model.  This used to be `if (any(...)) return(NULL)`
  ## three times over, so one covariate-carrying or unimplemented declaration
  ## disabled the M-step for every OTHER declaration in the model too -- and
  ## silently, since a NULL here just means "no metadata" downstream.
  .usable <- .fam > 0L & !.hasCov & !vapply(.args, anyNA, logical(1))
  if (!any(.usable)) return(NULL)
  ## Current copula correlation.  Read from the rxCor.* theta the expansion
  ## created, not from the declaration's ini() value: that theta is what the
  ## model actually uses, and it moves during the fit.  It carries atanh(rho)
  ## for a pair (the k=2 case of the row-normalized Cholesky), so tanh() of it
  ## is the correlation.
  .cw <- as.integer(.st$corWith)
  .rho <- rep(0, .n)
  .corTheta <- rep(NA_character_, .n)
  ## On the DIRECT route there is no `rxCor.*` theta to read: the expansion
  ## leaves the correlation in the omega where the user wrote it, because
  ## nothing in the model text needs it (on the cdf route the theta exists
  ## precisely because the decoder BUILDS the latent from it).  Reading the
  ## missing theta returns NULL from this function, which reads downstream as
  ## "no declarations" -- so a correlated direct model would silently lose both
  ## the copula AND the family M-step, with no error.
  ##
  ## Both declared variances are fixed at the placeholder 1, so the off-diagonal
  ## IS the correlation.
  .direct <- .etaDistIsDirect(.ui)
  if (.direct) {
    .off <- .ini[!is.na(.ini$neta1) & .ini$neta1 != .ini$neta2, , drop = FALSE]
    .etaRow <- .ini[!is.na(.ini$neta1) & .ini$neta1 == .ini$neta2, , drop = FALSE]
    .etaNum <- stats::setNames(.etaRow$neta1, .etaRow$name)
    for (.i in seq_len(.n)) {
      if (.cw[.i] < 0L) next
      .n1 <- .etaNum[[paste0("rxd.", .st$name[.i])]]
      .n2 <- .etaNum[[paste0("rxd.", .st$name[.cw[.i] + 1L])]]
      if (is.null(.n1) || is.null(.n2)) return(NULL)
      .w <- which((.off$neta1 == .n1 & .off$neta2 == .n2) |
                    (.off$neta1 == .n2 & .off$neta2 == .n1))
      ## no off-diagonal row means the block was declared correlated but the
      ## covariance is zero -- an independent pair, not a failure
      .rho[.i] <- if (length(.w) == 1L) .off$est[.w] else 0
    }
  }
  for (.i in seq_len(.n)) {
    if (.direct) break
    if (.cw[.i] < 0L) next
    .nm <- .st$corTheta[.i]
    if (is.na(.nm)) return(NULL)
    ## tolerate either name order rather than assuming the expansion's
    .w <- which(.thNames == .nm)
    if (length(.w) != 1L) {
      .alt <- paste0("rxCor.", .st$name[.cw[.i] + 1L], ".", .st$name[.i])
      .w <- which(.thNames == .alt)
      if (length(.w) != 1L) return(NULL)
      .nm <- .alt
    }
    .corTheta[.i] <- .nm
    .rho[.i] <- tanh(.thVals[[.nm]])
  }
  .am <- matrix(0, nrow = .n, ncol = max(1L, .maxA))
  for (.i in seq_len(.n)) .am[.i, seq_along(.args[[.i]])] <- .args[[.i]]
  list(n = .n, fam = .fam, corWith = .cw, args = .am, rho = .rho,
       corTheta = .corTheta, dist = .st$etaDist, thetas = .tn,
       hasCov = .hasCov, usable = .usable, cov = .cov,
       etas = .st$name, iniDf = .ini)
}

#' Declared-distribution M-step metadata for the FOCEi-family estimators
#'
#' The imp/impmap flavour of [.etaDistMstepInfo()].  imp numbers its random
#' effects by `neta1` and its thetas by `ntheta`, so the index maps are direct
#' -- there is no phi/phi0 split to go through.
#'
#' The returned `map` closure is what turns a fitted set of NATIVE family
#' parameters back into the user's thetas; it remembers its last answer and
#' starts the solve there, so each call is warm and local.
#'
#' @inheritParams .etaDistMstepCore
#' @return a list for `impEtaDistMstep()` (src/imp.cpp), or `NULL`
#' @noRd
.etaDistMstepInfoFocei <- function(ui) {
  .c <- .etaDistMstepCore(ui)
  if (is.null(.c)) return(NULL)
  .ini <- .c$iniDf
  .etaRows <- .ini[!is.na(.ini$neta1) & .ini$neta1 == .ini$neta2, , drop = FALSE]
  .etaNames <- .etaRows[order(.etaRows$neta1), "name"]
  .th <- .ini[!is.na(.ini$ntheta), , drop = FALSE]
  .thNames <- .th[order(.th$ntheta), "name"]
  ## the expansion renames the declared eta's latent `rxz.<eta>`
  .lat <- as.integer(match(paste0(.etaDistEtaPrefix(ui), .c$etas), .etaNames) - 1L)
  if (anyNA(.lat)) return(NULL)
  .ti <- lapply(seq_len(.c$n), function(.k) {
    as.integer(match(.c$thetas[[.k]], .thNames) - 1L)
  })
  if (any(vapply(.ti, anyNA, logical(1)))) return(NULL)
  ## the copula correlation is itself a theta and gets the same treatment; one
  ## per correlated PAIR, in the order the driver walks them.  Matched by the
  ## name the core resolved, not by grep()ing rxCor.* and trusting the order --
  ## an unrelated rxCor.* (a second, undeclared block) would silently shift it.
  .ck <- which(.c$corWith >= 0L)
  .cti <- integer(0)
  if (length(.ck) > 0L) {
    .cti <- as.integer(match(.c$corTheta[.ck], .thNames) - 1L)
    if (anyNA(.cti)) return(NULL)
  }
  ## warm-start state for the native-parameters -> thetas solve
  .env <- new.env(parent = emptyenv())
  .env$cur <- lapply(seq_len(.c$n), function(.k) {
    .nm <- .c$thetas[[.k]]
    stats::setNames(vapply(.nm, function(.t) {
      .w <- which(.ini$name == .t)
      if (length(.w) == 1L) as.numeric(.ini$est[.w]) else NA_real_
    }, numeric(1)), .nm)
  })
  .map <- function(k, args) {
    .k <- as.integer(k)
    ## Decline LOUDLY for a declaration this M-step does not own.  There is no
    ## single population `a` to invert when an argument varies by subject, so a
    ## quiet best-effort answer here would be a wrong number rather than a
    ## missing one.
    if (!isTRUE(.c$usable[.k])) return(NULL)
    .st <- .env$cur[[.k]]
    if (is.null(.st) || anyNA(.st)) return(NULL)
    .s <- .etaDistArgsToThetas(.c$dist[.k], .c$thetas[[.k]], .st, as.numeric(args))
    if (is.null(.s)) return(NULL)
    .env$cur[[.k]] <- .s
    as.numeric(.s)
  }
  ## The ARGUMENT expressions and their theta names, so imp's C++ M-step can do
  ## the native-parameters -> thetas map itself rather than calling `map` back
  ## into R once per iteration per declared eta.  Same emission saem's
  ## .etaDistMstepInfo() makes, for the same reason and consumed the same way:
  ## C++ parses them once into RPN and DECLINES anything outside its grammar,
  ## in which case `map` is still there to fall back on.
  .exprs <- lapply(seq_len(.c$n), function(.k) {
    .cl <- .c$dist[[.k]]
    if (is.character(.cl)) .cl <- str2lang(.cl)
    vapply(as.list(.cl)[-1], function(.e) paste(deparse(.e), collapse = ""),
           character(1))
  })
  ## The thetas this M-step OWNS, so the caller can drop them from the Newton
  ## step's sensitivity list / the outer free-parameter vector.
  ##
  ## Only the usable declarations' thetas.  A theta held out here with nothing
  ## to update it stays at its ini() value for the whole fit and reports as an
  ## estimate -- so a theta that ANY unusable declaration also reads has to
  ## stay in the outer optimizer, even if a usable declaration reads it too.
  ## Sharing is why this is a setdiff and not just `.c$thetas[.c$usable]`.
  .own <- unique(unlist(.c$thetas[.c$usable]))
  .foreign <- unique(unlist(.c$thetas[!.c$usable]))
  list(latent = .lat, fam = as.integer(.c$fam), corWith = as.integer(.c$corWith),
       direct = as.integer(.etaDistIsDirect(ui)),
       args = .c$args, rho = as.numeric(.c$rho),
       usable = as.integer(.c$usable),
       exprs = .exprs, exprThetas = .c$thetas,
       thetaIdx = .ti, corThetaIdx = .cti, map = .map,
       thetaNames = unique(c(setdiff(.own, .foreign), .c$corTheta[.ck])))
}

#' Warn when a covariate-carrying declaration is left unidentified
#'
#' MEASURED, and it is silent otherwise.  With the M-step on, the COVARIATE
#' COEFFICIENT on a declaration argument is not estimated -- the declaration's
#' other parameters are fine (lclm +1.660 and lclrv -2.095 against truth +1.63
#' and -2.40 in the same run).  The C++ argument
#' parser resolves symbols against the declaration's THETA names only, so an
#' argument reading a data column fails to parse
#' (`etaDistExprParse`, src/etaDistExpr.h), the general family objective at
#' src/saem.cpp never runs, and the coefficient is left in the outer problem
#' with an objective that no longer identifies it.  On simulated data with a
#' known allometric effect (truth +0.75, started +0.20) focei returned -0.027
#' with the M-step ON and +0.598 with it OFF, against +0.574 from a log-normal
#' reference fit of the same data -- so the DEFAULT route is right and the
#' M-step route quietly returns a wrong number rather than declining.
#'
#' The fix is to pass the covariate names to the parser and supply their
#' per-record values (`nSym`/`rec`, `rxEtaDistLoglikObj()` in src/etaDistFam.cpp,
#' which already takes them); until that is wired, say so.
#'
#' @param edi metadata from `.etaDistMstepInfoFocei()`
#' @param ui the ui, for the declaration names
#' @return nothing; called for the warning
#' @noRd
.etaDistWarnCovMstep <- function(edi, ui) {
  .u <- as.integer(edi$usable)
  if (length(.u) == 0L || all(.u == 1L)) return(invisible())
  .c <- .etaDistMstepCore(ui)
  .nm <- if (!is.null(.c) && length(.c$etas) == length(.u)) .c$etas[.u == 0L] else "a declaration"
  warning("the declared-distribution M-step cannot estimate the covariate ",
          "coefficient on ", paste0("dist(", .nm, ")", collapse = ", "),
          ".  The family fit stands down for that declaration and the ",
          "coefficient is left in the outer problem, where it is not ",
          "identified: measured on a known allometric effect it returned ",
          "-0.027 for a truth of +0.75.  The other parameters are unaffected. ",
          "Use etaDistMstep=FALSE (the default), where the covariate is ",
          "evaluated per record by the ordinary solve and the same effect ",
          "recovers as +0.598.", call. = FALSE)
  invisible()
}

#' Wire the declared-distribution M-step into a FOCEi-family control
#'
#' Builds the metadata `foceiEtaDistMstep()` (src/inner.cpp) reads and the
#' per-theta hold-out mask that keeps those thetas out of the outer optimizer's
#' free-parameter vector, and assigns both onto the ui's control.
#'
#' Both halves are assigned together or neither is: a mask without metadata
#' would hold thetas out of the optimizer with nothing to update them, leaving
#' them silently at their `ini()` values.  When the model carries no usable
#' declaration this assigns `NULL`/`integer(0)`, and the fit estimates those
#' thetas in the outer problem exactly as before.
#'
#' @param ui rxode2 ui carrying the FOCEi control
#' @return `ui`, invisibly; called for the control assignment
#' @noRd
.foceiEtaDistSetup <- function(ui) {
  .info <- NULL
  .skip <- integer(0)
  if (isTRUE(rxode2::rxGetControl(ui, "etaDistMstep", FALSE))) {
    .edi <- .etaDistMstepInfoFocei(ui)
    if (is.null(.edi)) {
      .etaDistMstepWarnInert("focei")
    } else {
      .etaDistWarnCovMstep(.edi, ui)
      .ini <- rxode2::rxUiDecompress(ui)$iniDf
      .th <- .ini[!is.na(.ini$ntheta), , drop = FALSE]
      .thNames <- .th[order(.th$ntheta), "name"]
      .m <- match(.edi$thetaNames, .thNames)
      .m <- .m[!is.na(.m)]
      ## a theta the user fixed is already out of the free set; leaving it in
      ## the mask would double-count it in foceiSetupTheta_'s fixedn
      .fx <- .th[order(.th$ntheta), "fix"]
      .skip <- as.integer(seq_along(.thNames) %in% .m & !(!is.na(.fx) & .fx))
      .info <- .edi
    }
  }
  rxode2::rxAssignControlValue(ui, "foceiEtaDistInfo", .info)
  rxode2::rxAssignControlValue(ui, "foceiEtaDistThetaSkip", .skip)
  invisible(ui)
}

#' Say so when the declared-distribution M-step was asked for but cannot run
#'
#' `etaDistMstep=TRUE` is opt-in, so a user who sets it has a reason to think it
#' is running.  Every reason it can decline is a silent one -- no declaration in
#' the model, a family the C++ dispatch does not implement, a copula block wider
#' than a pair, a declaration whose arguments are not plain thetas -- and the
#' fit then proceeds by the ordinary route with estimates that look perfectly
#' reasonable.  That is exactly the failure that is hardest to notice: the
#' option appears to work because the fit converges.
#'
#' `warning()` is the established route onto the fit's `$runInfo` (collected in
#' nlmixr2Est.R and printed under "Information about run").
#'
#' @param what which estimator is reporting, for the message
#' @return `NULL`, called for the warning
#' @noRd
.etaDistMstepWarnInert <- function(what) {
  warning(paste0(what, ": etaDistMstep=TRUE was requested but the declared ",
                 "distributions could not be mapped, so those parameters are ",
                 "estimated the ordinary way (no declaration in the model, an ",
                 "unimplemented family, or a copula block wider than a pair)"),
          call. = FALSE)
  NULL
}

#' Record and read mu2's covariate substitution for the declared M-step
#'
#' mu2 rewrites a covariate term into a generated `nlmixrMuDerCov#` column and
#' the raw covariate then never reaches the fit data.  A `dist()` declaration
#' that reads the same covariate has to be rewritten the same way, or the
#' expression the M-step parses names a column that is not there.
#'
#' `ui$mu2RefCovariateReplaceDataFrame` cannot serve that downstream: the
#' etaDist expansion rebuilds the ui, and the table that has a row while the
#' hoist runs has none by the time the M-step metadata is built.  `ui$meta` is
#' an environment that survives the rebuild -- it already carries the
#' declaration stash -- so the substitution is recorded there instead.
#'
#' @param ui rxode2 ui
#' @param modelExpression the text mu2 replaced
#' @param covariateParameter the coefficient it replaced it with
#' @param derived the generated column name
#' @return nothing, called for the side effect
#' @noRd
.etaDistMu2Record <- function(ui, modelExpression, covariateParameter, derived) {
  tryCatch({
    .m <- rxode2::rxUiDecompress(ui)$meta
    if (!is.environment(.m)) return(invisible())
    .cur <- if (exists(".etaDistMu2Map", envir = .m, inherits = FALSE)) {
      get(".etaDistMu2Map", envir = .m, inherits = FALSE)
    } else {
      data.frame(modelExpression = character(0), covariateParameter = character(0),
                 derived = character(0), stringsAsFactors = FALSE)
    }
    if (derived %in% .cur$derived) return(invisible())
    assign(".etaDistMu2Map",
           rbind(.cur, data.frame(modelExpression = as.character(modelExpression),
                                  covariateParameter = as.character(covariateParameter),
                                  derived = as.character(derived),
                                  stringsAsFactors = FALSE)),
           envir = .m)
  }, error = function(e) NULL)
  invisible()
}

#' Read back what `.etaDistMu2Record()` stored
#'
#' @param ui rxode2 ui
#' @return the recorded substitution table, or `NULL` when there is none
#' @noRd
.etaDistMu2Get <- function(ui) {
  tryCatch({
    .m <- rxode2::rxUiDecompress(ui)$meta
    if (!is.environment(.m)) return(NULL)
    if (!exists(".etaDistMu2Map", envir = .m, inherits = FALSE)) return(NULL)
    get(".etaDistMu2Map", envir = .m, inherits = FALSE)
  }, error = function(e) NULL)
}
