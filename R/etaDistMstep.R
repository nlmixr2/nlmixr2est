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
## ---------------------------------------------------------------------------
## 3.  The general objective
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
  .lat <- vapply(paste0("rxz.", .c$etas), function(.z) {
    .w <- which(etaNames == .z)
    if (length(.w) == 1L) as.integer(etaTrans[.w]) - 1L else -1L
  }, integer(1))
  if (any(.lat < 0)) return(NULL)                    # -> R fallback
  ## Each declared theta's PHI index.  saem's phi columns are the model
  ## parameters in `saemParamsToEstimate` order; a declared-distribution theta
  ## has no eta, so its phi column carries no random effect and lands in phi0,
  ## where .configsaem() finishes the mapping (phi index -> phi0 column).
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
  .ck <- which(.c$corWith >= 0L)
  if (!is.null(paramsToEstimate) && length(.ck) > 0L) {
    .rc <- .c$corTheta[.ck]
    .m <- match(.rc, paramsToEstimate)
    if (anyNA(.m)) return(NULL)
    .rn <- .rc; .rp <- as.integer(.m)
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
  list(latent = .lat, fam = .c$fam, corWith = .c$corWith,
       exprs = .exprs, exprThetas = .c$thetas,
       args = .c$args, rho = .c$rho,
       dist = .c$dist, thetas = .c$thetas, thetaPhi = .tp,
       corName = .rn, corPhi = .rp, etas = .c$etas)
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
    .st <- .etaDistDeclStash(.ui, .d)
    if (is.null(.st)) return(NULL)
  }
  .n <- length(.st$name)
  if (.n == 0L) return(NULL)
  .thNames <- .ini$name[!is.na(.ini$ntheta)]
  .thVals <- stats::setNames(as.list(.ini$est[!is.na(.ini$ntheta)]), .thNames)
  .fam <- integer(.n); .maxA <- 0L
  .args <- vector("list", .n); .tn <- vector("list", .n)
  for (.i in seq_len(.n)) {
    .fam[.i] <- .etaDistFamilyCode(.st$etaDist[.i])
    .cl <- str2lang(.st$etaDist[.i])
    .a <- vapply(as.list(.cl)[-1], function(.x) {
      .v <- tryCatch(eval(.x, envir = .thVals), error = function(e) NA_real_)
      if (!is.numeric(.v) || length(.v) != 1L) NA_real_ else as.numeric(.v)
    }, numeric(1))
    .args[[.i]] <- .a
    .maxA <- max(.maxA, length(.a))
    .tn[[.i]] <- intersect(all.vars(.cl), .thNames)
  }
  if (any(.fam <= 0)) return(NULL)               # -> R fallback
  if (any(vapply(.args, anyNA, logical(1)))) return(NULL)
  ## Current copula correlation.  Read from the rxCor.* theta the expansion
  ## created, not from the declaration's ini() value: that theta is what the
  ## model actually uses, and it moves during the fit.  It carries atanh(rho)
  ## for a pair (the k=2 case of the row-normalized Cholesky), so tanh() of it
  ## is the correlation.
  .cw <- as.integer(.st$corWith)
  .rho <- rep(0, .n)
  .corTheta <- rep(NA_character_, .n)
  for (.i in seq_len(.n)) {
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
  .lat <- as.integer(match(paste0("rxz.", .c$etas), .etaNames) - 1L)
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
  list(latent = .lat, fam = as.integer(.c$fam), corWith = as.integer(.c$corWith),
       args = .c$args, rho = as.numeric(.c$rho),
       exprs = .exprs, exprThetas = .c$thetas,
       thetaIdx = .ti, corThetaIdx = .cti, map = .map,
       ## the thetas this M-step owns, so the caller can drop them from the
       ## Newton step's sensitivity list / the outer free-parameter vector
       thetaNames = unique(c(unlist(.c$thetas), .c$corTheta[.ck])))
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
