# est="impmap" -- importance-sampling EM.
#
# impmap centers a multivariate-normal proposal at each subject's MAP mode
# (produced by the mu-referenced FOCEI inner problem, muModel="lin"), draws
# importance samples, and updates the population parameters by an EM step.
# The numerical kernel lives in C++ (src/imp.cpp); this file is orchestration
# only: control construction, dispatch, and post-fit assembly.

# Importance-sampling / EM control names -- stripped when down-converting to a
# plain foceiControl for the MAP inner problem / output.
.impmapIsControlNames <- c("isample", "nIter", "mapIter", "gamma",
                           "gammaMethod", "gammaMethodUser", "gammaRule",
                           "df", "auto", "autoNonNormal",
                           "autoNonmemSparse", "autoDfPatience",
                           "iscaleMin", "iscaleMax", "iaccept",
                           "ctol", "nConvWindow", "impSeed", "impCov",
                           "qr", "qrShift", "qrRefresh", "sir", "sirSample",
                           # internal M-step index maps added in .impmapFamilyFit;
                           # not foceiControl() arguments, so they must be dropped
                           # when down-converting (e.g. .setOfvFo's do.call(foceiControl))
                           "impMuThetaIdx", "impMuEtaIdx", "impThetaSensIdx",
                           "impOmegaFixedEta")

#' Control options for the impmap (importance-sampling EM) estimation method
#'
#' A NONMEM-style Monte Carlo importance-sampling EM built on the mu-referenced
#' FOCEI MAP.  The proposal density for each subject is centered at the MAP mode
#' (`muModel="lin"`); mu-referenced population parameters are updated by the EM
#' gradient, while non-mu parameters (structural and residual error) are updated
#' by a symbolic-sensitivity Newton step -- the importance-sampling-weighted
#' score and Gauss-Newton information built from the analytic `d(f)/d(theta)` and
#' `d(V)/d(theta)` (exact censored partials for BLQ/M2/M3/M4 points).
#'
#' @inheritParams foceiControl
#' @param ... Parameters used in the default `foceiControl()`
#' @param isample Number of importance samples drawn per subject per iteration
#'   (NONMEM ISAMPLE).  Either a single count used for every subject, or a
#'   vector of length `nsub` giving a count **per subject**.
#'
#'   Per-subject counts are the NM7 Technical Guide's own remedy for poor
#'   coverage (its derivation is Gaussian throughout and never mentions a t
#'   proposal): a subject whose weights are badly behaved can be given more
#'   samples without charging every other subject for them.  Note this treats
#'   the symptom rather than the cause -- more draws from a proposal whose
#'   tails are too light still gives weights with infinite variance, which
#'   `fit$env$impPsisK` will show.  See `df` for the shape-based remedy.
#' @param nIter Maximum number of importance-sampling EM iterations.
#' @param mapIter Number of MAP re-centering iterations per EM step; `> 0`
#'   re-centers the proposal at the MAP mode each iteration.
#' @param gamma Initial proposal-variance inflation factor (NONMEM ISCALE); the
#'   proposal covariance is `gamma` times the inverse of the inner information
#'   matrix at the mode.
#' @param df Degrees of freedom of the importance-sampling proposal (NONMEM
#'   `DF`).  `0` (default) uses a multivariate **normal** proposal; any value
#'   `> 0` uses a multivariate **t** with that many degrees of freedom.
#'
#'   This changes the proposal's SHAPE rather than its width, and that is the
#'   distinction that matters.  `gamma` can only make a Gaussian proposal
#'   wider; it cannot give it heavier tails.  When the target posterior has
#'   heavier tails than the proposal, the importance weights have infinite
#'   variance -- and neither `xi` nor the Kish effective sample size can detect
#'   it, because the offending mass lies where the proposal rarely lands.  The
#'   Pareto k-hat diagnostic (`fit$env$impPsisK`) does detect it: `k > 0.7`
#'   means that subject's weights are unreliable.  A t proposal has polynomial
#'   tails that dominate a Gaussian target's, which bounds the weights.
#'
#'   NONMEM's guidance (Bauer, *NONMEM Tutorial Part II*) is to set a nonzero
#'   `DF` when there are fewer data points than etas, or for categorical data.
#'   Small values (3-8) are heavy; large values approach the Gaussian.
#' @param auto NONMEM `AUTO=1` equivalent: adapt the proposal degrees of
#'   freedom, the sample count and the acceptance target **per subject** rather
#'   than applying one global setting to everybody.
#'
#'   * **`df`** -- any subject gets a heavy-tailed t proposal when the model is
#'     not transformably normal, and any subject whose Pareto k-hat reports tail
#'     failure gets one on that evidence.  The tutorial's other trigger, "fewer
#'     data points than there are ETAs", is **not** applied on its own -- see
#'     `autoNonmemSparse`.  An escalation that fails to improve k-hat over two
#'     iterations is withdrawn, since a heavy tail the data creates is not
#'     repairable by proposal shape.
#'   * **`isample`** -- the total sample budget (`isample * nsub`) is
#'     reallocated toward subjects whose effective-sample fraction is lowest,
#'     the tutorial's "many ETAs or ... large stochastic fluctuations".  It is
#'     load-balancing, not a cost increase.  Note sample count is deliberately
#'     *not* driven by Pareto k-hat: a heavy tail is not repairable by more
#'     draws (see `df`).
#'   * **`iaccept`** -- lowered to 0.2 for the same sparse/categorical
#'     subjects, per the tutorial.
#'
#'   The concrete values (the `df` ladder, the k-hat thresholds, the budget
#'   reallocation rule) are **nlmixr2's choices**.  NONMEM does not publish what
#'   `AUTO=1` picks internally; only the triggers and `IACCEPT ~ 0.2` are
#'   documented.  Unlike NONMEM's AUTO, which its own tutorial warns "may result
#'   in lack of stochastic reproducibility", this remains seeded and
#'   thread-count independent.
#'
#'   Per-subject values are reported in `fit$env$impDfInd`,
#'   `fit$env$impNsampleInd` and `fit$env$impIacceptInd`.
#'
#'   **Measured trade-off.**  `auto = TRUE` (the default) improves the *tail*
#'   behaviour of the importance weights and the accuracy of `Omega`, at some
#'   cost in Monte-Carlo noise on the objective.  On theophylline with three
#'   ETAs -- a one-ETA model has no tail failure to fix, so it shows nothing
#'   either way -- against a reference computed at `isample = 8000`, over 8
#'   seeds:
#'
#'   \itemize{
#'     \item `auto = FALSE`: max Pareto k-hat 0.941, 2.38 subjects above 0.7;
#'       objective RMSE 0.0113, `Omega` RMSE 0.00406
#'     \item `auto = TRUE`: max Pareto k-hat 0.571, 0.25 subjects above 0.7;
#'       objective RMSE 0.0165, `Omega` RMSE 0.00301
#'   }
#'
#'   So `auto` removes the tail failure and estimates `Omega` about 26% more
#'   accurately, for about 46% more Monte-Carlo noise on the objective.  It is
#'   on by default because weights with infinite variance are a correctness
#'   problem -- their error is unbounded in the worst case -- whereas the extra
#'   noise is a bounded, measurable cost.  Setting `df` globally gives a better
#'   objective still on a model that is uniformly heavy-tailed, but costs 75%
#'   more objective RMSE on one that is not, which is why it is not the default.
#'
#'   Set `auto = FALSE` to recover the previous behaviour exactly (that path
#'   remains bit-identical to earlier versions), which is worth doing when you
#'   want the tightest possible objective on a model whose `fit$env$impPsisK`
#'   values are already comfortably below 0.7.
#' @param autoNonmemSparse Apply the NONMEM tutorial's `nobs < neta` trigger
#'   unconditionally, so a subject with fewer observations than random effects
#'   always gets a t proposal.  `FALSE` (default) leaves such subjects to the
#'   Pareto k-hat evidence like any other, and withdraws an escalation that does
#'   not improve k-hat within two iterations.
#'
#'   The default diverges from the tutorial deliberately, on measurement.  On a
#'   fixture built to the tutorial's own definition of sparse (2 observations, 3
#'   etas), 8 seeds against an `isample = 8000` reference, applying the rule made
#'   every number worse -- objective RMSE 0.1517 -> 0.2059, `Omega` RMSE 0.02379
#'   -> 0.03163, and *more* failing subjects (3.88 -> 4.88).  With fewer
#'   observations than random effects the individual posterior is not identified,
#'   so the heavy tail is structural and no proposal shape repairs it; the
#'   reference itself still reads max k-hat 0.794.
#'
#'   Set `TRUE` to get the documented NONMEM rule anyway.  NONMEM's own testing
#'   is not published and was not done on these models, so someone who measures
#'   the opposite on their own problem should be able to have it.
#' @param autoDfPatience Number of consecutive iterations an escalated `df` may
#'   fail to improve Pareto k-hat before that escalation is withdrawn and the
#'   subject returned to the proposal it had without it.  `0` never withdraws.
#'   Withdrawal is final, so escalation and withdrawal cannot oscillate, and it
#'   never goes below the `df` the model itself requires -- a non-normal endpoint
#'   keeps its t proposal.
#'
#'   A rung is judged against what the subject manages WITHOUT any escalation,
#'   measured while it sits there, so deterioration that happens before any
#'   escalation is tracked and an escalation is never credited for it.
#' @param gammaMethod How the proposal scale `gamma` is adapted during the EM.
#'
#'   `"auto"` (default) picks per model: `"individual"` when the model is not
#'   transformably normal -- a general log-likelihood (`ll()`) endpoint, or a
#'   count/categorical/time-to-event distribution -- and `"global"` otherwise.
#'   That is the split the hypothesis actually rests on: `gamma = 1` is already
#'   the efficient proposal when the individual posterior is close to Gaussian,
#'   so a normal model gains nothing from per-subject adaptation and would only
#'   pay for it in effective sample size, while a general-likelihood model is
#'   exactly where the posteriors go non-Gaussian and per-subject coverage
#'   starts to vary.  The test is `all(ui$predDf$distribution == "norm")`, the
#'   same line [rxode2::assertRxUiTransformNormal()] draws.  The resolved value
#'   is reported in the fit's `$runInfo`.
#'
#'   `"global"` keeps one scale shared by every subject, inflated (never
#'   relaxed) only when the *mean* Kish effective-sample fraction falls below
#'   `iaccept`.  It leaves `gamma` at its efficient starting value while
#'   coverage is healthy, which is the right behaviour when the individual
#'   posteriors are close to Gaussian.
#'
#'   `"individual"` gives every subject its own `gamma_i` and adapts it
#'   two-sided toward a target on that subject's `xi_i`, clamped to
#'   `[iscaleMin, iscaleMax]`.  This follows NONMEM's *objective* -- the NM7
#'   Technical Guide states that `gamma` is per-subject (eq. 1.90) and is
#'   "continually adjusted so that xi_i approximates IACCEPT" (note after
#'   eq. 1.76), bounded by `ISCALE_MIN`/`ISCALE_MAX` -- but the guide publishes
#'   no update formula, so the functional form used here is nlmixr2's own.
#'   Prefer this when the
#'   individual posteriors are heavy-tailed or the design is heterogeneous:
#'   a global scale is driven by the mean, so a handful of badly-covered
#'   subjects never trip it and their likelihood and `Omega` contributions end
#'   up carried by a few samples.
#'
#'   The two modes report *different* efficiency statistics -- `"individual"`
#'   targets `xi` (the mean normalized importance weight, NONMEM's `IACCEPT`
#'   quantity) while `"global"` targets the Kish effective-sample fraction.
#'   These are not comparable; the fit's `$runInfo` says which is in force.
#' @param gammaRule How the SHARED (`gammaMethod="global"`) proposal scale is
#'   adapted.  Ignored for `gammaMethod="individual"`, which always follows
#'   NONMEM's two-sided per-subject rule.
#'
#'   `"floor"` treats `iaccept` as a one-sided FLOOR on the mean Kish
#'   effective-sample fraction: `gamma` stays at its efficient starting value
#'   while coverage is healthy and is inflated only when coverage drops below the
#'   floor.  It never comes back down.
#'
#'   `"target"` (default) follows the NM7 Technical Guide, which says `gamma` is
#'   "continually adjusted so that xi_i approximates IACCEPT" -- adapted BOTH
#'   ways, on `xi` rather than the Kish fraction, using the same analytic
#'   inversion as the per-subject controller (`gamma * (xi/iaccept)^(2/neta)`,
#'   capped 1.25x each way and clamped to `[iscaleMin, iscaleMax]`).
#'
#'   The two rules settle at different operating points, so they are not
#'   interchangeable, and constants tuned against one do not carry over to the
#'   other.
#'
#'   **Measured trade-off** (theophylline, 6 seeds, `isample = 300`, RMSE against
#'   an `isample = 6000` reference):
#'
#'   \itemize{
#'     \item 3 ETAs -- `"floor"`: theta RMSE 0.00185, `Omega` RMSE 0.00298, max
#'       Pareto k-hat 0.604, 0.33 subjects above 0.7, converged 100% of seeds.
#'       `"target"`: theta RMSE 0.00257, `Omega` RMSE 0.00378, max k-hat -0.491,
#'       0 subjects above 0.7, converged 0%.
#'     \item 8 ETAs -- `"floor"` never adapts at all (`gamma` pinned at 1.0 with
#'       `xi` 1.35, i.e. a proposal far too narrow), because it watches the mean
#'       Kish fraction, which stays healthy while `xi` says the proposal is
#'       wrong.  `"target"` moves `gamma` to 1.25 and puts `xi` on 0.41.
#'   }
#'
#'   So `"target"` fixes the tail and is the only rule that reacts at high ETA
#'   dimension, but it pays for it: putting `xi` on `iaccept` deliberately widens
#'   the proposal (`gamma` 1.79 on the 3-ETA fixture), which costs effective
#'   sample size (0.96 to 0.70) and adds Monte-Carlo noise to the objective --
#'   enough that the fit often does not meet `ctol` within `nIter`.
#'
#'   `"target"` is the default: tuned against tuned, it wins every column at 3
#'   ETAs and is the only rule that adapts at all at 8 ETAs.  Its cost is on the
#'   1-ETA fixture, where it roughly doubles theta RMSE (0.00113 to 0.00224) and
#'   takes about twice as many iterations.  That trade follows the same weighting
#'   `auto` uses: weights with infinite variance are a correctness problem whose
#'   error is unbounded, while the extra Monte-Carlo noise is bounded and
#'   measurable.
#'
#'   Choose `"floor"` for the previous behaviour -- a proposal left at its
#'   efficient starting value while coverage is healthy.  It is also what the
#'   tail-machinery tests pin, because `"target"` repairs the tail itself and
#'   leaves the `df` ladder nothing to fix.
#'
#'   **The tuned constants travel with the rule.**  Selecting a rule also selects
#'   its tuned defaults (`nConvWindow` 20 for `"target"`, 10 for `"floor"`), so
#'   switching to the NONMEM method does not silently run NONMEM's law on the
#'   other rule's tuning.  An explicitly supplied value always wins.
#' @param iscaleMin,iscaleMax Lower/upper bounds for the adapted `gamma`
#'   (NONMEM ISCALE_MIN / ISCALE_MAX).  Both bounds are reachable under
#'   `gammaRule="target"`; under `"floor"` `gamma` only ever moves up, so only
#'   `iscaleMax` can bind.
#' @param iaccept Minimum importance-sampling effective-sample fraction
#'   (NONMEM IACCEPT).  The proposal scale `gamma` is kept at its efficient
#'   starting value while the achieved fraction stays at or above `iaccept`, and
#'   is inflated (toward `iscaleMax`) only when it drops below this floor.
#' @param ctol Convergence tolerance on the windowed objective-function change;
#'   `NULL` derives it from `sigdig`.
#' @param nConvWindow Length of the trailing iteration window used to average
#'   the objective-function change for convergence (NONMEM-style CTYPE).
#' @param muModel Mu-referencing variant for the MAP inner problem; for
#'   `impmapControl()` this is always `"lin"` and cannot be changed.
#' @param impSeed Base seed for the per-subject thread-safe (threefry) RNG
#'   streams; results are reproducible and independent of the thread count.
#' @param covMethod Covariance method.  `"imp"` (default) computes the
#'   Monte-Carlo importance-sampling observed-information covariance for the
#'   estimated thetas and Omega parameters (a finite-difference Hessian of the
#'   importance-sampling objective over fixed common-random-number samples),
#'   stashed as `$impCov` / `$impSe` and installed as the fit covariance; the
#'   theta standard errors match the Hessian-based FOCEI covariance, though the
#'   variance of a tightly-determined random effect (an Omega diagonal) can be
#'   over-estimated because the fixed samples barely span its prior variation.
#'   `"analytic"`, `"r,s"`, `"r"`, `"s"` instead compute the FOCEI covariance
#'   post-fit at the converged estimates (see [foceiControl()]); `""` skips the
#'   covariance step.
#' @param qr When `TRUE`, draw quasi-random (Sobol low-discrepancy) importance
#'   samples instead of pseudo-random Gaussian samples (QRPEM, Leary &
#'   Dunlavey PAGE 2012); the E-step integrals converge at O(1/N) instead of
#'   O(1/sqrt(N)).
#' @param qrShift Only used with `qr=TRUE`.  When `TRUE` each (iteration,
#'   subject) applies a random Cranley-Patterson shift to the Sobol points
#'   (seeded, thread-count independent); `FALSE` reuses one fixed Sobol point
#'   set everywhere (fully deterministic E-step, no RNG in the draw).
#' @param qrRefresh Only used with `qr=TRUE` and `qrShift=TRUE`.  When `TRUE`
#'   the shift is redrawn each iteration so residual quasi-random error
#'   averages out over the EM; `FALSE` draws one shift per subject at the fit
#'   start, making each EM iteration a deterministic map (smoothest objective
#'   trace).
#' @param sir When `TRUE`, accelerate the non-mu / residual-error M-step by
#'   SIR (sampling-importance-resampling): the theta-sensitivity Newton step
#'   uses `sirSample` equal-weight resampled points per subject instead of all
#'   `isample` weighted samples.
#' @param sirSample Number of SIR resampled points per subject; `NULL` uses
#'   `max(25, ceiling(isample/10))`.  Must be at most `isample`.
#' @return impmapControl object
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' impmapControl()
impmapControl <- function(sigdig=3,
                          ...,
                          isample=300L,
                          nIter=100L,
                          mapIter=1L,
                          gamma=1.0,
                          gammaMethod=c("auto", "global", "individual"),
                          gammaRule=c("target", "floor"),
                          df=0,
                          auto=TRUE,
                          autoNonmemSparse=FALSE,
                          autoDfPatience=2L,
                          iscaleMin=0.1,
                          iscaleMax=10.0,
                          iaccept=0.4,
                          ctol=NULL,
                          nConvWindow=10L,
                          impSeed=42L,
                          covMethod=c("imp", "analytic", "r,s", "r", "s", "sa", ""),
                          qr=FALSE,
                          qrShift=TRUE,
                          qrRefresh=TRUE,
                          sir=FALSE,
                          sirSample=NULL,
                          muModel=c("lin", "none")) {
  muModel <- match.arg(muModel)
  gammaMethod <- match.arg(gammaMethod)
  gammaRule <- match.arg(gammaRule)
  # RULE-DEPENDENT DEFAULTS.  The tuned constants belong to the RULE, so switching
  # to the NONMEM rule must bring its constants with it -- otherwise the user gets
  # NONMEM's law running on the floor rule's tuning, which is the unfair pairing
  # this whole comparison was redone to avoid.  Only defaults are substituted: an
  # explicitly supplied value always wins (and a control round-tripped through
  # do.call(impmapControl, ctl) supplies everything, so it is idempotent).
  #
  # nConvWindow: 10 for "floor", 20 for "target".  The target rule tracks a
  # Monte-Carlo statistic, so it needs a longer window to average the noise out.
  # Measured (3 ETAs, 6 seeds, isample=300, RMSE vs an isample=6000 reference):
  #
  #   target w=10   thetaRMSE 0.00196  omegaRMSE 0.00266  maxK -0.356  iter 16
  #   target w=20   thetaRMSE 0.00167  omegaRMSE 0.00196  maxK -0.429  iter 27
  #   target w=30   thetaRMSE 0.00199  omegaRMSE 0.00291  maxK -0.451  iter 34
  #
  # w=20 is the accuracy optimum; w=30 costs iterations for nothing.
  if (identical(gammaRule, "target") && missing(nConvWindow)) {
    nConvWindow <- 20L
  }
  checkmate::assertLogical(qr, any.missing=FALSE, len=1, .var.name="qr")
  checkmate::assertLogical(qrShift, any.missing=FALSE, len=1, .var.name="qrShift")
  checkmate::assertLogical(qrRefresh, any.missing=FALSE, len=1, .var.name="qrRefresh")
  checkmate::assertLogical(sir, any.missing=FALSE, len=1, .var.name="sir")
  # isample may be a single count or one count PER SUBJECT (NONMEM's per-subject
  # ISAMPLE): a badly covered subject can buy more samples without charging
  # every other subject for them.
  checkmate::assertIntegerish(isample, any.missing=FALSE, min.len=1, lower=1,
                              .var.name="isample")
  checkmate::assertIntegerish(impSeed, any.missing=FALSE, len=1, .var.name="impSeed")
  .isampleAll <- as.integer(isample)
  # A per-subject vector must be exactly nsub long.  nsub is not known here, so
  # the length is checked at fit time in the kernel (it used to fall back to a
  # uniform max() and ignore a wrong-length vector without a word).
  # the scalar used for sizing / sirSample defaults / reporting is the largest
  .isample <- max(.isampleAll)
  if (is.null(sirSample)) {
    .sirSample <- max(25L, as.integer(ceiling(.isample / 10)))
  } else {
    checkmate::assertIntegerish(sirSample, any.missing=FALSE, len=1, lower=1,
                                .var.name="sirSample")
    .sirSample <- as.integer(sirSample)
  }
  if (.sirSample > .isample) {
    stop("'sirSample' (", .sirSample, ") cannot exceed 'isample' (", .isample, ")",
         call.=FALSE)
  }
  # covMethod="imp" drives the Monte-Carlo importance-sampling covariance in the
  # C++ kernel (op_focei.impCov); every other token is the post-fit FOCEI
  # covariance, so hand foceiControl a valid token (the estimation pass forces
  # covMethod=0L regardless -- see .impmapFamilyFit).
  .dots <- list(...)
  .impCov <- isTRUE(.dots$impCov)   # may already be set on a round-tripped control
  .dots$impCov <- NULL              # internal field; do not forward to foceiControl
  # gammaMethodUser is stamped on the RUNTIME control by .impmapFamilyFit (it
  # records what the user asked for before "auto" was resolved).  A control that
  # has been round-tripped therefore carries it; keep it, but do not forward it
  # to foceiControl, which has no such argument.
  .gammaMethodUser <- .dots$gammaMethodUser
  .dots$gammaMethodUser <- NULL
  # autoNonNormal is likewise stamped on the RUNTIME control (it records whether
  # the resolved gammaMethod came out "individual"), and is not an argument here
  # either.  Without this, do.call(impmapControl, <a fit's own control>) died
  # with "unused argument: 'autoNonNormal'" -- so re-validating a completed fit's
  # control never worked.
  .autoNonNormal <- .dots$autoNonNormal
  .dots$autoNonNormal <- NULL
  if (is.character(covMethod)) {
    if (length(covMethod) == 1L && !nzchar(covMethod)) {
      covMethod <- ""
    } else {
      covMethod <- match.arg(covMethod)
    }
    .impCov <- identical(covMethod, "imp")
    .foceiCovMethod <- if (.impCov) "analytic" else covMethod
  } else {
    # round-trip: covMethod is already a resolved foceiControl integer slot;
    # keep the impCov flag from the incoming control (read above)
    .foceiCovMethod <- covMethod
  }
  .control <- do.call(foceiControl,
                      c(list(sigdig=sigdig), .dots,
                        list(covMethod=.foceiCovMethod, muModel="lin")))
  .control$impCov <- .impCov
  if (!is.null(.autoNonNormal)) .control$autoNonNormal <- .autoNonNormal
  .control$isample <- .isampleAll
  .control$nIter <- as.integer(nIter)
  .control$mapIter <- as.integer(mapIter)
  .control$gamma <- as.double(gamma)
  .control$gammaMethod <- gammaMethod
  .control$gammaRule <- gammaRule
  .control$df <- as.double(df)
  checkmate::assertLogical(auto, any.missing=FALSE, len=1, .var.name="auto")
  checkmate::assertLogical(autoNonmemSparse, any.missing=FALSE, len=1,
                           .var.name="autoNonmemSparse")
  checkmate::assertIntegerish(autoDfPatience, lower=0, len=1, any.missing=FALSE,
                              .var.name="autoDfPatience")
  .control$auto <- auto
  .control$autoNonmemSparse <- autoNonmemSparse
  .control$autoDfPatience <- as.integer(autoDfPatience)
  if (!is.null(.gammaMethodUser)) .control$gammaMethodUser <- .gammaMethodUser
  .control$iscaleMin <- as.double(iscaleMin)
  .control$iscaleMax <- as.double(iscaleMax)
  .control$iaccept <- as.double(iaccept)
  .control$ctol <- if (is.null(ctol)) NULL else as.double(ctol)
  .control$nConvWindow <- as.integer(nConvWindow)
  .control$impSeed <- as.integer(impSeed)
  .control$qr <- qr
  .control$qrShift <- qrShift
  .control$qrRefresh <- qrRefresh
  .control$sir <- sir
  .control$sirSample <- .sirSample
  class(.control) <- "impmapControl"
  .control
}

#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.impmapControl <- function(control, env) {
  assign("impmapControl", control, envir=env)
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.impmap <- function(control) {
  .ctl <- control[[1]]
  .cls <- class(control)[1]
  if (is.null(.ctl)) .ctl <- impmapControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list"))
    .ctl <- do.call("impmapControl", .ctl)
  if (inherits(.ctl, "foceiControl")) {
    .minfo(paste0("converting ", class(.ctl)[1], " to impmapControl"))
    class(.ctl) <- NULL
    .ctl <- do.call(impmapControl, .ctl)
  } else if (!inherits(.ctl, "impmapControl")) {
    .minfo(paste0("invalid control for `est=\"", .cls, "\"`, using default"))
    .ctl <- impmapControl()
  } else {
    .ctl <- do.call(impmapControl, .ctl)
  }
  .ctl
}

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.impmap <- function(x, ...) {
  .env <- x[[1]]
  if (exists("impmapControl", .env)) {
    .control <- get("impmapControl", .env)
    if (inherits(.control, "impmapControl")) return(.control)
  }
  if (exists("control", .env)) {
    .control <- get("control", .env)
    if (inherits(.control, "impmapControl")) return(.control)
  }
  stop("cannot find impmap related control object", call.=FALSE)
}

.impmapControlToFoceiControl <- function(env, assign=TRUE) {
  .impmapControl <- env$impmapControl
  .n <- setdiff(names(.impmapControl), .impmapIsControlNames)
  # np* internals (npBoxLower/npPoints/npResidFreeze ...) and the npag/npb user
  # knobs (points/cycles/gammaOptimize/... -- but NOT seed, a real foceiControl arg)
  # are nonparametric-engine control fields that foceiControl does not accept; drop
  # them so a downstream do.call(foceiControl, .) (e.g. .setOfvFo, general-likelihood
  # tables) does not error with "unused argument".
  .npKnobs <- c("points", "cycles", "gammaOptimize", "residOptimize", "muExpand",
                "gridWidth", "gridBounds", "dfScan",
                "alpha", "burnin", "nsamp", "nchains", "propSd", "est")
  .n <- .n[!grepl("^np[A-Z]", .n) & !(.n %in% .npKnobs)]
  .foceiControl <- setNames(lapply(.n, function(n) .impmapControl[[n]]), .n)
  class(.foceiControl) <- "foceiControl"
  if (assign) env$control <- .foceiControl
  .foceiControl
}

#' @rdname nmObjGetFoceiControl
#' @export
nmObjGetFoceiControl.impmap <- function(x, ...) {
  .env <- x[[1]]
  .impmapControlToFoceiControl(.env, assign=FALSE)
}

#' Resolve gammaMethod="auto" against the model
#'
#' `"individual"` when the model is NOT transformably normal -- a general
#' log-likelihood (`ll()`) endpoint, or a count / categorical / time-to-event
#' distribution -- and `"global"` otherwise.
#'
#' The rationale is the hypothesis the per-subject controller rests on: the
#' proposal is the Laplace approximation, so `gamma = 1` is already efficient
#' when the individual posterior is near-Gaussian.  A normal model therefore
#' gains nothing from per-subject adaptation and pays for it in effective
#' sample size (measured on theophylline: ESS 0.95 -> 0.70 for identical
#' estimates), whereas a general-likelihood model is where the posteriors are
#' genuinely non-Gaussian and per-subject coverage varies.
#'
#' The test is `all(predDf$distribution == "norm")`, the same line
#' `rxode2::assertRxUiTransformNormal()` draws between transformably-normal and
#' general likelihoods.
#'
#' @param gammaMethod Requested setting: "auto", "global" or "individual".
#' @param ui rxode2 ui object
#' @return "global" or "individual"
#' @noRd
.impmapResolveGammaMethod <- function(gammaMethod, ui) {
  if (!identical(gammaMethod, "auto")) return(gammaMethod)
  .dist <- tryCatch(ui$predDf$distribution, error=function(e) NULL)
  # No usable predDf (should not happen for a fittable model): fall back to the
  # conservative choice, which is the historical behaviour.  Same for an NA
  # distribution -- `all(NA == "x")` is NA, which would error an `if`.
  if (is.null(.dist) || length(.dist) == 0L) return("global")
  .dist <- as.character(.dist)
  if (anyNA(.dist)) return("global")
  # Canonicalize before comparing: rxode2 spells the Gaussian family both
  # "norm" and "dnorm" and treats them as identical (rxPreferredDistributionName
  # maps both to "dnorm"), so a literal == "norm" test would send an otherwise
  # ordinary normal model down the "individual" path.  That matters because
  # `linCmt() ~ add(add.sd) + dnorm()` -- the exact-likelihood (Laplace) form of
  # a plain normal endpoint -- reports "dnorm" while having a posterior every
  # bit as Gaussian as the add() form, so it should not pay for per-subject
  # adaptation.
  #
  # Note lognormal / boxCox / yeoJohnson residuals are carried in
  # predDf$transform with distribution still "norm", so they canonicalize to
  # "dnorm" here and correctly count as Gaussian.
  .canon <- tryCatch(rxode2::rxPreferredDistributionName(.dist),
                     error=function(e) .dist)
  if (all(.canon == "dnorm")) "global" else "individual"
}

#' Build the $runInfo note naming the efficiency statistic in force
#'
#' Both statistics are always computed and stashed on the fit -- NONMEM-style
#' `xi` in `$impXiTrace`/`$impXi` and the Kish effective-sample fraction in
#' `$impNeffFrac`/`$impNeff` -- but only one of them drives the proposal-scale
#' adaptation, and they are NOT the same quantity even though both are governed
#' by `iaccept`.  Measured on theophylline at gamma = 1, mean xi = 1.009 while
#' the mean Kish fraction = 0.997; at gamma = 4 they read 0.503 and 0.667.
#' Reading one as if it were the other, or comparing a "global" fit's number
#' against an "individual" fit's, is a live mistake, so every fit says which is
#' which.
#'
#' @param resolved Resolved gammaMethod: "global" or "individual".
#' @param user What was requested ("auto", "global", "individual").
#' @param iaccept Target/floor value.
#' @return A single-line character message.
#' @noRd
.impmapGammaRunInfo <- function(resolved, user, iaccept) {
  .why <- if (identical(user, "auto")) {
    if (identical(resolved, "individual")) {
      " (auto: not every endpoint is Gaussian)"
    } else {
      " (auto: all endpoints are Gaussian)"
    }
  } else {
    ""
  }
  if (identical(resolved, "individual")) {
    paste0("gammaMethod=\"individual\"", .why,
           ": the proposal scale is adapted per subject toward xi=", iaccept,
           " (NONMEM IACCEPT).  Sampling efficiency for this fit is xi",
           " ($impXi/$impXiTrace); the Kish effective-sample fraction",
           " ($impNeffFrac) is a DIFFERENT statistic and the two are not",
           " comparable -- nor is xi comparable across gammaMethod settings.")
  } else {
    paste0("gammaMethod=\"global\"", .why,
           ": one shared proposal scale, adapted only when the mean Kish",
           " effective-sample fraction falls below ", iaccept,
           ".  Sampling efficiency for this fit is that fraction",
           " ($impNeff/$impNeffFrac); NONMEM-style xi ($impXiTrace) is also",
           " reported but is a DIFFERENT statistic and the two are not",
           " comparable.")
  }
}

#' Install the impmap control into the ui
#'
#' @param env Environment with ui in it
#' @param ... Other arguments
#' @return Nothing, called for side effects
#' @author Matthew L. Fidler
#' @noRd
.impmapFamilyControl <- function(env, ...) {
  .ui <- env$ui
  .control <- env$control
  if (is.null(.control)) {
    .control <- impmapControl()
  }
  if (!inherits(.control, "impmapControl")) {
    .control <- do.call(nlmixr2est::impmapControl, .control)
  }
  assign("control", .control, envir=.ui)
}

#' Fit the impmap family of models
#'
#' @param env Environment from nlmixr2Est
#' @param ui rxode2 ui object
#' @param ... Other arguments
#' @return fit environment
#' @author Matthew L. Fidler
#' @noRd
.impmapFamilyFit <- function(env, ui, ...) {
  # With est="impmap", foceiFitCpp_ runs impOuter() (src/imp.cpp) in place of
  # foceiOuter().  The FOCEI outer optimizer is turned off (maxOuterIterations=0)
  # because impOuter drives its own EM iteration; the in-fit covariance is off
  # (it would bail on muModel="lin" anyway) -- the covariance is computed
  # post-fit on the full model at the converged estimates (.foceiRecomputeMuCov).
  .control <- ui$control
  .covMethodUser <- .control$covMethod  # restored on the fit env control below
  .control$maxOuterIterations <- 0L
  .control$covMethod <- 0L
  # Resolve gammaMethod="auto" here, where the ui (and therefore predDf) is in
  # scope; the C++ kernel only ever sees a concrete "global"/"individual".
  #
  # Resolve from the USER's original choice, not from whatever this control
  # currently holds: resolution overwrites gammaMethod in place, so a control
  # taken off a finished fit (fit$env$impmapControl) and reused on a different
  # model would otherwise carry the previous model's resolved value and never
  # re-resolve -- an "auto" control from a log-likelihood fit would silently
  # pin a normal model to "individual".  Keying off gammaMethodUser makes
  # resolution idempotent and re-runnable.
  .gmUser <- .control$gammaMethodUser
  if (is.null(.gmUser)) .gmUser <- .control$gammaMethod
  .control$gammaMethodUser <- .gmUser                # kept for $runInfo
  .control$gammaMethod <- .impmapResolveGammaMethod(.gmUser, ui)
  # AUTO's "or data are categorical" trigger: reuse the same transformably-normal
  # test the gammaMethod resolution uses, so there is one notion of model class.
  .control$autoNonNormal <- identical(.impmapResolveGammaMethod("auto", ui), "individual")
  # Say which efficiency statistic this fit's number actually is.  Both are
  # always stashed but they are not the same quantity, and only one drives the
  # adaptation -- warning() here is the established route onto $runInfo
  # (collected in nlmixr2Est.R and printed under "Information about run").
  warning(.impmapGammaRunInfo(.control$gammaMethod, .gmUser,
                              .control$iaccept),
          call.=FALSE)
  # 0-based index maps for the SIMPLE mu-referenced intercepts (theta = population
  # mean of an eta, no covariates): impOuter's M-step shifts each such theta by
  # the mean conditional eta.  Covariate mu-groups are excluded here because they
  # are handled by the regression update (updateMuGroups) instead.
  .env <- ui$foceiOptEnv  # builds foceiMuGroupTheta (the covariate-group thetas)
  .iniDf <- ui$iniDf
  .th <- .iniDf[!is.na(.iniDf$ntheta), ]
  .thNames <- .th[order(.th$ntheta), "name"]
  .etaRows <- .iniDf[!is.na(.iniDf$neta1) & .iniDf$neta1 == .iniDf$neta2, ]
  .etaNames <- .etaRows[order(.etaRows$neta1), "name"]
  .mr <- ui$muRefDataFrame
  .muThetaIdx <- as.integer(match(.mr$theta, .thNames) - 1L)
  .muEtaIdx <- as.integer(match(.mr$eta, .etaNames) - 1L)
  .covGroupTheta <- rxode2::rxGetControl(ui, "foceiMuGroupTheta", integer(0))
  .keep <- !is.na(.muThetaIdx) & !is.na(.muEtaIdx) &
    !(.muThetaIdx %in% .covGroupTheta)
  .control$impMuThetaIdx <- .muThetaIdx[.keep]
  .control$impMuEtaIdx <- .muEtaIdx[.keep]
  # 0-based theta indices of the estimated non-mu thetas (structural + residual
  # error) with sensitivity outputs in the sensitivity model; the M-step Newton
  # update maps its output columns back to these thetas.
  .control$impThetaSensIdx <- as.integer(.impmapEstTheta(ui)$all - 1L)
  # 0-based eta indices whose Omega diagonal is FIXED; the EM Omega update
  # restores their rows/columns to the starting value so fix()ed variances hold.
  .etaOrd <- .etaRows[order(.etaRows$neta1), ]
  .control$impOmegaFixedEta <- as.integer(which(isTRUE(.etaOrd$fix) | .etaOrd$fix) - 1L)
  assign("control", .control, envir=ui)
  # Seed the importance-sampling RNG from the control (impSeed) right before the
  # fit, mirroring saem's set.seed(seed).  The E-step draws through rxode2's
  # threefry engine (getRxSeed1), so a fixed rxseed makes the fit reproducible
  # and independent of the thread count / ambient RNG state; rxWithSeed restores
  # the prior seed state afterward so it does not leak globally.
  .impSeed <- if (is.null(.control$impSeed)) 42L else as.integer(.control$impSeed)
  # est is "impmap" or "imp" (the no-MAP-search variant); pass it through so the
  # C++ kernel (impOuter) selects the proposal accordingly.
  .est <- if (exists("est", envir=env)) get("est", envir=env) else "impmap"
  .fit <- rxode2::rxWithSeed(.impSeed, rxseed=.impSeed,
                             code=.foceiFamilyReturn(env, ui, ..., est=.est))
  # The MC covariance (impCov=TRUE) is published with theta row/column names but
  # the Omega parameters come out unnamed on this path; fill them in (defensively,
  # only when the counts line up) so vcov()/$cov and the correlation are labelled.
  .impmapNameCov(.fit, ui)
  .impRestoreCovMethod(.fit, .covMethodUser)
  # Capture THIS fit's pooled-solve layout before anything else runs.  The odeSwap
  # registry is process-global and describes the most recent registration, so the
  # objective recompute below -- which runs a nested focei fit -- re-registers the
  # slots and the global view stops describing this fit.  Stash it here so the
  # diagnostic travels with the fit and cannot be overwritten by a later one.
  tryCatch({
    .fenv0 <- .fit$env
    if (is.environment(.fenv0)) assign("odeSwapInfo", .odeSwapInfo(), envir=.fenv0)
  }, error=function(e) NULL)
  .impmapRecomputeObjf(.fit)
  # Tail-sensitive companion to xi / Kish ESS: computed post-fit from the
  # stashed final-iteration weights so it costs nothing during the EM.
  tryCatch({
    .fenv <- .fit$env
    if (is.environment(.fenv)) {
      assign("impPsisK", .impPsisKAll(.fenv), envir=.fenv)
    }
  }, error=function(e) NULL)
  .fit
}

#' Publish the objective as a FOCEi evaluation at the converged estimates.
#'
#' The in-C++ finalize (`impMapPass` -> `foceiOuterFinal`) computes the individual
#' objective from an eta-Hessian that never receives its data term on this path:
#' `fInd->a` (d(pred)/d(eta)) is left at allocation residue, so
#' `sum(cHff * a * a)` vanishes and `H` collapses to `Omega^-1` alone.  Measured on
#' theophylline with one random effect: `log|H|` reads 2.17048 for EVERY subject
#' instead of the correct 3.63-3.96, and the published objective comes out 19.96
#' too LOW (173.63 against 193.60).  Models with 2+ random effects are unaffected.
#'
#' Rather than repair that C++ state, take the objective the same way `setOfv()`
#' does for SAEM (`.setOfvFo`, `R/ofv.R`): re-evaluate through the ordinary
#' `nlmixr2()` FOCEi path at the converged estimates.  `maxOuterIterations = 0`
#' pins the thetas and Omega; the INNER problem is deliberately left to optimize,
#' because the FOCEi objective is defined with the etas at their conditional mode
#' -- freezing them at the importance-sampling conditional means is a different
#' (and worse) quantity, measured at 193.615 against 193.601 for the mode.
#'
#' The EM's own estimates are untouched: only the published objective is replaced.
#' `$impObj` (the importance-sampling estimate) is left alone, so the two remain
#' separately readable.
#'
#' @param fit completed impmap-family fit
#' @return invisibly TRUE when the objective was replaced
#' @noRd
.impmapRecomputeObjf <- function(fit) {
  .env <- tryCatch(fit$env, error=function(e) NULL)
  if (!is.environment(.env)) return(invisible(FALSE))
  # deep-copy the UI: the nested re-fit must not mutate THIS fit's UI
  .ui <- tryCatch(rxode2::rxUiDecompress(unserialize(serialize(fit$ui, NULL))),
                  error=function(e) NULL)
  if (is.null(.ui)) return(invisible(FALSE))
  .sigdig <- tryCatch(fit$foceiControl$sigdig, error=function(e) NULL)
  # A nested nlmixr2() calls .nlmixr2globalReset(), which clears nlmixr2global --
  # including the timing environment (dropping the outer fit's "other" timing row)
  # and nlmixr2objectName (which would report ".ui" instead of the user's symbol).
  # Snapshot and restore the whole thing, plus the mu-referencing global.
  .savedMuRef <- .muRefTrans$cur
  on.exit(.muRefTrans$cur <- .savedMuRef, add=TRUE)
  .savedGlobal <- as.list(nlmixr2global, all.names=TRUE)
  on.exit({
    rm(list=ls(nlmixr2global, all.names=TRUE), envir=nlmixr2global)
    for (.gn in names(.savedGlobal)) assign(.gn, .savedGlobal[[.gn]], envir=nlmixr2global)
  }, add=TRUE)
  .ctl <- try(foceiControl(print=0L, covMethod="", maxOuterIterations=0L,
                           calcTables=FALSE, compress=FALSE,
                           sigdig=if (is.null(.sigdig)) 4 else .sigdig),
              silent=TRUE)
  if (inherits(.ctl, "try-error")) return(invisible(FALSE))
  .f2 <- try(suppressMessages(suppressWarnings(
    nlmixr2(.ui, data=nlme::getData(fit), est="focei", control=.ctl))),
    silent=TRUE)
  if (inherits(.f2, "try-error")) return(invisible(FALSE))
  .e2 <- tryCatch(.f2$env, error=function(e) NULL)
  if (!is.environment(.e2)) return(invisible(FALSE))
  # carry the objective AND the per-subject quantities derived from the same
  # (correct) Hessian, so $etaObf/$phiH do not disagree with the published number
  for (.n in c("objective", "OBJF", "objf", "logLik", "AIC", "BIC",
               "etaObf", "etaObfFull", "phiH", "phiC", "phiR", "phiSE", "phiRSE")) {
    if (exists(.n, envir=.e2, inherits=FALSE)) {
      assign(.n, get(.n, envir=.e2), envir=.env)
    }
  }
  # objDf is MERGED, not replaced: the re-fit runs covMethod="" so its objDf has no
  # Condition#(Cov)/Condition#(Cor), and replacing wholesale would drop the columns
  # the imp covariance had already filled in.
  if (exists("objDf", envir=.e2, inherits=FALSE)) {
    .newObjDf <- get("objDf", envir=.e2)
    .oldObjDf <- tryCatch(get("objDf", envir=.env, inherits=FALSE), error=function(e) NULL)
    if (is.data.frame(.oldObjDf) && is.data.frame(.newObjDf) &&
          nrow(.oldObjDf) == nrow(.newObjDf)) {
      for (.c in intersect(names(.newObjDf), names(.oldObjDf))) {
        .oldObjDf[[.c]] <- .newObjDf[[.c]]
      }
      assign("objDf", .oldObjDf, envir=.env)
    } else {
      assign("objDf", .newObjDf, envir=.env)
    }
  }
  invisible(TRUE)
}

#' Restore the requested covMethod on the fit env's stored control
#'
#' The estimation pass forces covMethod=0L (the in-fit C++ step would bail on
#' muModel="lin"), and that runtime control is what gets stored on the fit env;
#' the post-fit recompute (.foceiRecomputeMuCov) reads the covMethod from there,
#' so put the requested choice back.
#' @noRd
.impRestoreCovMethod <- function(fit, covMethod) {
  .fenv <- tryCatch(fit$env, error = function(e) NULL)
  if (is.environment(.fenv) &&
        exists("impmapControl", envir = .fenv, inherits = FALSE)) {
    .ic <- get("impmapControl", envir = .fenv)
    .ic$covMethod <- covMethod
    assign("impmapControl", .ic, envir = .fenv)
  }
  invisible(fit)
}

#' Fill the Omega row/column names on the impmap covariance
#' @param fit impmap fit
#' @param ui rxode2 ui
#' @return Nothing, called for side effects
#' @noRd
.impmapNameCov <- function(fit, ui) {
  .fenv <- tryCatch(fit$env, error=function(e) NULL)
  if (is.null(.fenv) || is.null(.fenv$cov) || !is.matrix(.fenv$cov)) return(invisible())
  tryCatch({
    .dn <- dimnames(.fenv$cov)[[1]]
    if (is.null(.dn)) return(invisible())
    .empty <- which(is.na(.dn) | .dn == "")
    if (length(.empty) == 0L) return(invisible())
    .etaN <- .foceiEtaThetaMap(ui)$etaNames
    .op <- .foceiOmegaPairs(.fenv$omega, ui$iniDf)
    .omN <- .foceiOmegaCovNames(.op, .etaN)
    if (length(.omN) == length(.empty)) {
      .dn[.empty] <- .omN
      dimnames(.fenv$cov) <- list(.dn, .dn)
      if (!is.null(.fenv$fullCor) && is.matrix(.fenv$fullCor)) {
        dimnames(.fenv$fullCor) <- list(.dn, .dn)
      }
    }
  }, error=function(e) NULL)
  invisible()
}

#' @rdname nlmixr2Est
#' @export
nlmixr2Est.impmap <- function(env, ...) {
  .ui <- env$ui
  # General (dnorm/ll) likelihoods flow through the shared FOCEI inner problem
  # (impEvalJointLik = likInner0), so only require transformable normality when
  # the rxode2 build has no llik support -- mirrors nlmixr2Est.focei.
  if (!rxode2hasLlik()) {
    rxode2::assertRxUiTransformNormal(.ui, " for the estimation routine 'impmap'", .var.name=.ui$modelName)
  }
  rxode2::assertRxUiIovNoCor(.ui, " for the estimation routine 'impmap'",
                             .var.name=.ui$modelName)
  .control <- env$control
  .foceiFamilyControl(env, ..., type="impmapControl")
  on.exit({
    if (is.environment(.ui) && exists("control", envir=.ui, inherits=FALSE)) {
      rm("control", envir=.ui)
    }
  }, add=TRUE)
  env$impmapControl <- .control
  env$est <- "impmap"
  .ui <- env$ui
  .impmapFamilyFit(env, .ui, ...)
}
attr(nlmixr2Est.impmap, "covPresent") <- TRUE
attr(nlmixr2Est.impmap, "unbounded") <- .foUnbounded
attr(nlmixr2Est.impmap, "iov") <- TRUE
# Activates the mu2/mu3/mu4 covariate-rewriting hook (.uiApplyMu2hook, R/mu2.R),
# gated on muModel/muRefCovAlg, exactly as the mfocei family does.
attr(nlmixr2Est.impmap, "mu") <- function(control) {
  isTRUE(!identical(control$muModel, "none")) && isTRUE(control$muRefCovAlg)
}
