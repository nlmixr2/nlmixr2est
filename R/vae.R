#' Control for vae (variational autoencoder) estimation method in nlmixr2
#'
#' Variational-autoencoder NLME estimation (Rohleff et al., CPT:PSP 2025): an
#' LSTM encoder learns the individual posterior q(eta|y) and an rxode2 decoder
#' reconstructs the observations, trained on an ELBO / BICc-ELBO objective for
#' simultaneous population-parameter estimation and covariate selection.
#'
#' @inheritParams saemControl
#' @inheritParams foceiControl
#'
#' @param seed Random seed for the VAE training (encoder init, Adam,
#'   reparameterization sampling); default 42. Training is stochastic, so a fixed
#'   seed makes every fit reproducible.
#' @param itersBurnIn Number of burn-in iterations (encoder-only, tiny KL
#'   weight) before the main EM phase.
#' @param klWarmup Number of KL-annealing iterations over which the KL weight is
#'   ramped from a small value to 1 (prevents posterior collapse).
#' @param gammaIter Number of main iterations before the EMA-smoothing phase of
#'   the population-parameter update begins.
#' @param iters Total number of main-loop iterations (after burn-in).
#' @param nGradStep Number of Adam gradient steps per EM outer iteration
#'   (the reference `L_iter`).
#' @param hiddenDim LSTM hidden dimension (the reference `h_dim`).
#' @param learningRate Adam learning rate used in the main training phase.
#' @param burnInLearningRate Adam learning rate used during burn-in.
#' @param sigma0 Encoder prior standard deviation(s) at initialization (a small
#'   value giving a sharp initial posterior). `NULL` uses a small default per
#'   individual parameter. This is distinct from the `ini()` omega.
#' @param likelihood Inner likelihood used for the objective, EBEs, and
#'   gradients, all run through the same FOCEi inner interface: `"focei"`
#'   (default, with eta-epsilon interaction), `"foce"` (no interaction, NONMEM
#'   FOCE with R frozen at the population prediction), `"focep"` (FOCE+, no
#'   interaction but R evaluated at the live conditional eta), or `"laplace"`.
#' @param covariateSelection When `TRUE` (default) perform automated BICc-ELBO
#'   covariate selection during training; when `FALSE` fit only the covariate
#'   structure written in the model.  In the `FALSE` case the model-declared
#'   covariate coefficients (both linear `beta*WT` effects and transformed ones
#'   such as `beta*log(WT/70)`) are estimated in place by the regress M-step
#'   regardless of `nonMuTheta`; a `ini(... ~ fix())` coefficient stays fixed.
#' @param pinCovariates When `TRUE` (default) and the model already declares
#'   covariate effects, restrict the automatic covariate selection to only the
#'   covariate/parameter pairs written in the model -- the branch-and-bound
#'   search may still drop a declared covariate, but can never add one the model
#'   did not specify.  A declared covariate that is not a valid search candidate
#'   (time-varying, or a raw-linear form that does not match the `log`/centered
#'   encoding) is estimated in place by the regress M-step instead, with a note
#'   in `$runInfo`.  When the model declares no covariates there is nothing to
#'   pin and the full search runs.  Has no effect when `covariateSelection` is
#'   `FALSE`.
#' @param muRefCovAlg When `TRUE` (default) an algebraic/centered covariate
#'   effect written in the model (e.g. `wt.cl*(WT/70)` or `wt.cl*log(WT/70)`) is
#'   handled as a mu2/mu3 reference: the covariate expression -- including its
#'   centering -- is evaluated into an internal `nlmixrMuDerCov#` data column and
#'   the model uses the linear `wt.cl*nlmixrMuDerCov#` form during fitting, so the
#'   VAE covariate search never re-centers it.  The original expression is
#'   restored in the reported model.
#' @param nonMuTheta How to treat a structural population `theta` that has no
#'   random effect (is not mu-referenced) so it can still be estimated by the VAE
#'   (which only estimates parameters that occupy the latent space).  For the
#'   eta-injection modes a small eta is injected so the parameter enters the
#'   latent space, and the reported fixed effect is `theta + mean(eta)` with the
#'   temporary eta dropped from the output model.
#'
#'   * `"regress"` (default, matching `saemControl(nonMuTheta=)`): no eta is
#'     injected; instead each such theta is estimated directly, re-optimized every
#'     M-step by a bounded `bobyqa` regression against the full FOCEi outer
#'     objective (bounds from the `ini()` lower/upper), blended with the M-step gain.
#'     `mStepObjective` selects which objective that regression targets.  This
#'     recovers a no-random-effect population parameter without adding a spurious
#'     random effect.  `nonMuEtaOmega` is unused in this mode.
#'   * `"grad"`: same target as `"regress"` but stepped with the EXACT analytic
#'     outer gradient (Almquist sensitivity equations, the machinery behind
#'     `foceiControl(fast=TRUE)`) instead of a derivative-free search: one
#'     augmented sensitivity solve per M-step replaces the bobyqa sweep.  Both
#'     modes optimize the same full outer objective (with every mu-referenced
#'     theta held at its current M-step value), so this changes the optimizer,
#'     not the target.  It is also the more natural fit for the method: the
#'     gradient is handed to the SAME Adam machinery that moves the encoder
#'     weights, so the parameter is learned alongside the rest of the model on a
#'     shared schedule (same gain, same KL warmup gate), whereas `"regress"`
#'     pauses each M-step to run a separate derivative-free optimizer to
#'     convergence and adopts its answer.  This is NOT a speed option -- it is
#'     measurably SLOWER than
#'     `"regress"` (on `theo_sd`, 1.47x with one non-mu theta and 1.13x with
#'     three; the gap narrows as the number grows, since bobyqa's cost scales in
#'     it and a single solve does not, but it does not close).  Choose it for
#'     accuracy: the exact gradient lands closer to the maximum-likelihood value
#'     than the derivative-free search (`theo_sd` non-mu `tv`: 3.4294 vs 3.4324,
#'     against a FOCEi MLE of 3.4293).  Applies to a conditionally Gaussian model
#'     and to a single non-Gaussian (`ll()`/generalized) endpoint, which
#'     differentiates the log-density directly.  Falls back to `"regress"` when
#'     the model is out of analytic scope (`linCmt()`, IOV, `fo`, a
#'     multi-endpoint or censored `ll()` model, ...); `nonMuEtaOmega` is unused.
#'   * `"eta"`: inject the eta with an ESTIMATED omega (starting at
#'     `nonMuEtaOmega`); the typical value is estimated and appears in the
#'     iteration table.
#'   * `"fix"`: inject the eta with omega held FIXED at `nonMuEtaOmega` AND hold
#'     the typical-value theta fixed at its `ini()` value.  Nothing about the
#'     parameter is estimated, so it is not shown in the iteration table (it is
#'     reported at its `ini()` value, marked fixed, with the injected eta dropped).
#'   * `"none"`: leave non-mu-referenced thetas frozen at their `ini()` value (the
#'     historic behavior).
#' @param mStepObjective Objective the non-mu-referenced theta M-step
#'   (`nonMuTheta = "regress"` or `"grad"`) is optimized against.  It has no
#'   effect when there is no non-mu-referenced structural theta, and it never
#'   changes the encoder/ELBO training step or the covariate branch-and-bound
#'   criterion, both of which always follow the reference.
#'
#'   * `"outer"` (default): the full FOCEi outer objective -- the frozen-eta
#'     joint likelihood PLUS the Laplace determinant, `0.5*log|Omega^-1|` and the
#'     DV-transform Jacobian.  This is a deliberate deviation from Rohleff et al.
#'     (2025): it keeps the quantity being optimized equal to the objective the
#'     fit reports, and it is the functional the analytic outer gradient
#'     differentiates, so `nonMuTheta = "grad"` optimizes one target rather than
#'     stepping one and scoring another.
#'   * `"elbo"`: the reference behavior -- the plain variational bound
#'     (frozen-eta joint likelihood, no Laplace term), matching the M-step in
#'     Rohleff et al. (2025).  Use it to reproduce the reference implementation.
#'     The analytic outer gradient does not apply to this objective, so
#'     `nonMuTheta = "grad"` is downgraded to `"regress"` with a note in
#'     `$runInfo`.
#'
#'   The two objectives differ by terms that depend on the non-mu thetas through
#'   the eta Hessian, so they can land on different estimates, and -- because
#'   those estimates feed the latent means the covariate search regresses on --
#'   on different covariate sets.
#' @param nonMuEtaOmega Variance of the eta injected for a non-mu-referenced theta
#'   (starting value for `nonMuTheta="eta"`, fixed value for `nonMuTheta="fix"`;
#'   unused for `"regress"`).
#' @param covSelectAlpha Starting multiplier for the covariate-selection L0
#'   penalty, ramped linearly from `covSelectAlpha` down to `1` over the
#'   `klWarmup` warmup iterations and held at `1` afterward (matching the
#'   reference implementation's `linspace(alpha, 1, kl_iter)`).  Values `> 1`
#'   penalize covariate entry more heavily early in training; `1` disables the
#'   ramp.
#' @param covSelectSmooth When `TRUE` (default) the covariate selection regresses
#'   the SAEM sufficient statistic -- an exponential moving average of the
#'   posterior means, updated with the same gain as the M-step -- rather than the
#'   current posterior means.  This matches the reference implementation
#'   (Rohleff et al. 2025), which is the reason for the default.  In practice it
#'   changes little: `gamma` is exactly 1 until `gammaIter`, so the statistic
#'   equals the posterior mean for most of a run and is averaged only over the
#'   closing tail.  `FALSE` regresses the current posterior means.
#' @param gammaSeries Decaying step-size series used once the smoothing phase
#'   starts (after `gammaIter`); the gain is 1 throughout the EM phase either way.
#'
#'   * `"reference"` (default): `1/(iter - gammaIter)`, the textbook
#'     Kuhn-Lavielle series the reference implementation uses.  The first
#'     smoothing step is still a full replacement, and the decay follows.
#'   * `"saem"`: `1/(1 + iter - gammaIter)`, the CONTINUATION form
#'     \code{\link{saemControl}()} uses -- nlmixr2est's SAEM builds its series so
#'     it continues rather than repeating a gain of 1, so the decay begins at
#'     `1/2`.  Select this to match the step-size convention of the other
#'     nlmixr2 estimation methods rather than the reference.
#' @param sigma0Interp How `sigma0` is turned into the encoder's initial posterior
#'   spread.  The encoder head emits `logSigma` and forms `diag(L) = exp(logSigma)`,
#'   so `diag(L)` is the posterior standard deviation.
#'
#'   * `"sd"` (default): the bias is `log(sigma0)`, so the initial posterior SD is
#'     `sigma0` -- what the argument says it is.
#'   * `"reference"`: the bias is `log(sigma0^2)`, matching the reference
#'     implementation, whose initial posterior SD is therefore `sigma0` SQUARED
#'     (`1e-6` rather than `1e-3` for the first neonatal dimension).  The
#'     reference documents `sigma0` as a standard deviation, so this appears to be
#'     unintended there; it is offered only to reproduce its published behavior.
#' @param residRhoend Final trust-region radius (`rhoend`) of the bounded
#'   `bobyqa` that estimates the residual parameters -- its convergence
#'   tolerance.  `NULL` (default) derives it from `sigdig` (`10^(-sigdig)`), the
#'   same way every other optimizer tolerance in the package is derived, so
#'   `sigdig` stays the single knob that moves them together.  Set it explicitly
#'   when the residual step should converge tighter than the rest: it runs with
#'   the ODE frozen, so tightening it is far cheaper than tightening `rhoend`,
#'   which also tightens the structural regression that re-solves per candidate.
#' @param residOptimize How the residual-error parameters are estimated.
#'
#'   Residual forms the optimizer estimates: `add`, `prop`, `add + prop`, `pow`,
#'   `lnorm`, and a `boxCox` or `yeoJohnson` lambda (bounded to `(-2, 2)`).  For
#'   a transform-both-sides model the objective transforms `dv` only and carries
#'   the log-Jacobian, since `f` leaves the solve already on the transformed
#'   scale.
#'
#'   `nonMuTheta = "grad"` bypasses this entirely: the analytic outer gradient
#'   already carries a residual sigma and a transform lambda as its own
#'   directions, so those parameters are stepped by the gradient through Adam and
#'   the two-stage path never runs.  Which converges better is model-dependent.
#'
#'   * `"moment"`: the closed-form moment estimator.  For a model with
#'     a single additive error this is exactly the optimum (`sqrt(SSE/n)`); for
#'     any other error model it is either a different estimator or, for the forms
#'     with no closed form (`pow`, Box-Cox, Yeo-Johnson), no estimator at all --
#'     the parameter stays at its `ini()` value.  There is no moment estimator for
#'     a log-likelihood (`ll()`) parameter either, so those also stay at `ini()`;
#'     use `"twoStage"` for such a model.
#'   * `"twoStage"` (default): block coordinate descent, as `npag`'s
#'     `residOptimize = "alternate"` does.  Stage
#'     one optimizes the non-mu-referenced structural thetas with the residual
#'     parameters held, so it is driven by `(dv - f)`; stage two then holds those
#'     and optimizes the residual parameters alone against the extended
#'     least-squares objective `sum[(y-f)^2/r + log r]` over the CACHED `(y, f)`
#'     pairs.  Because `f` is fixed by stage one, stage two needs no ODE re-solve
#'     -- the same structure SAEM uses.  On `theo_sd` this beats the moment
#'     estimator on both a pure-additive model (objective 131.79 vs 131.81) and a
#'     combined one (121.03 vs 122.47).
#'
#'     Which parameters stage two owns is decided per parameter: an error
#'     parameter, or one that no `d/dt()` right-hand side, initial condition or
#'     dosing modifier can reach.  The second case is what a log-likelihood
#'     (`ll()`) or generalized endpoint needs -- its residual-like parameters are
#'     plain thetas with no error row, and on the error-only rule stage two was
#'     empty for such a model, silently making `"twoStage"` behave like
#'     `"optimize"`.  When no regressed theta qualifies (every one feeds the
#'     solve) stage two has nothing to do and `residOptimize` has no effect.
#'   * `"optimize"` (EXPERIMENTAL, diagnostic): a single JOINT solve over the
#'     structural and residual parameters together, against the full outer
#'     objective.  Fine with one free residual parameter, but with `add` and
#'     `prop` both free it diverges -- they are near-collinear, and routing the
#'     residual through the full outer objective lets the Laplace terms move with
#'     it at frozen etas (objective 320.7 against the moment estimator's 122.5).
#'     Retained for comparison; prefer `"twoStage"`.
#' @param omegaUpdate How the population variances are updated in the covariate
#'   M-step.  `"suffStat"` (default) follows the reference: `omega` is formed from
#'   the EMA sufficient statistics and ASSIGNED outright.  `"blend"` is the
#'   historic behavior, blending the freshly computed `omega` with the previous
#'   value at the M-step gain (so it is smoothed twice).  Applies to `omega` only.
#'   The residual error estimate is still EMA-smoothed on the standard-deviation
#'   scale, where the reference smooths the residual sum of squares and takes the
#'   root afterwards -- a known remaining difference.  Matching it would need
#'   per-endpoint sufficient statistics plus an optimizer branch for the error
#'   models with no closed form (`add + prop`, `add + pow`, Box-Cox /
#'   Yeo-Johnson), as \code{\link{saemControl}()} does.
#' @param inputScale Which observations the encoder-input centering and scaling
#'   are computed over.  `"reference"` (default) matches the reference
#'   implementation, which takes the mean and SD across the whole padded
#'   observation matrix, so the zero padding of subjects with fewer observations
#'   enters both statistics.  On a ragged dataset that is a materially different
#'   scale from `"observed"`, which uses only the observed values (on the neonatal
#'   case study the SD is 1582 against 506).  Affects only the encoder's inputs,
#'   never the likelihood.
#' @param covSelectMethod How the covariate M-step searches subsets.  `"bnb"` is
#'   the exact branch-and-bound; it becomes impractical past a few dozen candidate
#'   covariates.  `"l0learn"` has the suggested `L0Learn` package propose supports,
#'   which are then scored and polished with the same exact objective -- so the
#'   search is approximate but the scoring is not.  `"auto"` (default) uses
#'   `"l0learn"` for a latent dimension with at least `covSelectMaxExact` candidate
#'   covariates and `"bnb"` otherwise.  If `L0Learn` is not installed when the
#'   exact search would be impractical (`"auto"` at or over the threshold, or an
#'   explicit `"l0learn"`) this errors rather than run the slow exact search
#'   silently; install `L0Learn`, or set `covSelectMaxExact = Inf` to force the
#'   exact search everywhere.
#' @param covSelectMaxExact Candidate-covariate count at or above which
#'   `covSelectMethod = "auto"` switches a latent dimension to `L0Learn` (default
#'   `17`, the measured wall-clock crossover).  Counted after `pinCovariates`
#'   trimming, so it is the size of the search actually run.  `Inf` forces the
#'   exact branch-and-bound for every dimension.
#' @param bnbStrategy Frontier discipline for the exact branch-and-bound covariate
#'   selection: `"lifo"` (default, last-in-first-out depth-first search),
#'   `"fifo"` (first-in-first-out) or `"lc"` (least cost / best-first).  The
#'   solver is exact, so the selected covariates are identical for every strategy;
#'   only the search order (and thus efficiency) differs.
#' @param parEncoderBackward Parallelize the encoder backward (gradient) pass over
#'   subjects.  Defaults to `TRUE` unless `options(nlmixr2.identical = TRUE)` is
#'   set (which flips the default to `FALSE`); an explicit value here always wins.
#'   The encoder forward pass and the covariate branch-and-bound already run
#'   multi-threaded and are bit-identical to the serial run.  The backward gradient
#'   is a continuous cross-subject sum, so parallelizing it (per-thread partials
#'   reduced in thread order) makes the result deterministic for a fixed number of
#'   `cores` but no longer bit-identical to the serial path: the per-step gradient
#'   differs at ~1e-12, which compounds through the iterative SGD/EM training to a
#'   small final difference (well below any estimation tolerance), and results may
#'   differ across different `cores`.  When it is active (and `cores > 1`) a note is
#'   added to the fit's `$runInfo`.  Set this to `FALSE` -- or globally
#'   `options(nlmixr2.identical = TRUE)` -- for bit-identical, fully reproducible
#'   results.
#' @param objf Which objective-function value is active for AIC/BIC/BICc. Both
#'   the linearization and importance-sampling -2LL are always computed and
#'   stored; this selects the default active one.
#' @param covMethod Method for calculating the covariance at the VAE estimates,
#'   run through the FOCEi covariance step; the same choices as
#'   \code{\link{foceiControl}()}: \code{"analytic"} (default), \code{"r,s"},
#'   \code{"r"}, \code{"s"}, or \code{""} to skip.
#' @param nIsSample Number of importance-sampling draws for the IS -2LL.
#' @param rhoend Final trust-region radius (`rhoend`) of the inner bounded
#'   `bobyqa` used by the non-mu / covariate regress M-step.  `NULL` (default)
#'   derives it from `sigdig` (`10^(-sigdig)`, matching the optimizer convergence
#'   tolerance), or `1e-4` when `sigdig` is `NULL`.
#' @param returnVae When `TRUE` return the raw VAE training object instead of the
#'   nlmixr2 fit.
#'
#' @details
#'
#' Covariate selection -- MIQP vs. branch-and-bound.  Per latent parameter the
#' selection step minimizes the same L0/BIC objective
#' `RSS_S/omega + log(N)*|S|` over subsets `S` of the candidate covariates (`RSS_S`
#' is the residual sum of squares of the ordinary-least-squares fit on the
#' intercept plus `S`).  The reference implementation (Rohleff et al.) writes this
#' as a Mixed-Integer Quadratic Program (MIQP) -- binary include/exclude
#' indicators with big-M constraints -- and solves it with the commercial Gurobi
#' solver through `cvxpy`.  No MIQP-capable solver is freely available in R: Gurobi
#' is commercial/licensed, and the open QP solvers on CRAN (e.g. `osqp`) are
#' continuous-only and cannot represent the binary selection.  A continuous convex
#' relaxation (L1 / lasso) would be solvable but only approximates best subset.
#'
#' This package instead solves the identical L0/BIC objective EXACTLY with a
#' self-contained branch-and-bound: each candidate support's coefficients are the
#' closed-form OLS fit and branches are pruned by a valid lower bound (the RSS of
#' the OLS fit using all still-free covariates).  It therefore returns the same
#' optimum the MIQP would -- no commercial dependency and no relaxation/accuracy
#' loss -- and scales to a few dozen covariates.  The search is worst-case
#' exponential in the number of covariates, but the pruning makes the practical
#' (sparse) case fast (e.g. 32 candidate covariates in a fraction of a second).
#'
#' @return vae control structure (class `vaeControl`)
#' @export
#' @author Matthew L. Fidler
vaeControl <- function(seed = 42L,
                       itersBurnIn = 100L,
                       klWarmup = 50L,
                       gammaIter = 250L,
                       iters = 300L,
                       nGradStep = 5L,
                       hiddenDim = 25L,
                       learningRate = 5e-3,
                       burnInLearningRate = 8e-3,
                       sigma0 = NULL,
                       covariateSelection = TRUE,
                       pinCovariates = TRUE,
                       muRefCovAlg = TRUE,
                       covSelectAlpha = 2,
                       covSelectSmooth = TRUE,
                       gammaSeries = c("reference", "saem"),
                       sigma0Interp = c("sd", "reference"),
                       residOptimize = c("twoStage", "moment", "optimize"),
                       residRhoend = NULL,
                       omegaUpdate = c("suffStat", "blend"),
                       inputScale = c("reference", "observed"),
                       covSelectMethod = c("auto", "bnb", "l0learn"),
                       covSelectMaxExact = 17L,
                       bnbStrategy = c("lifo", "fifo", "lc"),
                       parEncoderBackward = !isTRUE(getOption("nlmixr2.identical", FALSE)),
                       nonMuTheta = c("regress", "grad", "eta", "fix", "none"),
                       nonMuEtaOmega = 0.01,
                       mStepObjective = c("outer", "elbo"),
                       likelihood = c("focei", "foce", "focep", "laplace"),
                       objf = c("importanceSampling", "linear"),
                       nIsSample = 3000L,
                       returnVae = FALSE,

                       print = 1L,
                       useColor = NULL,
                       printNcol = NULL,

                       covMethod = c("r,s", "analytic", "r", "s", "sa", "imp", ""),
                       optExpression = TRUE,
                       sumProd = FALSE,
                       literalFix = TRUE,
                       literalFixRes = TRUE,
                       addProp = c("combined2", "combined1"),
                       calcTables = TRUE,
                       compress = FALSE,
                       adjObf = TRUE,
                       ci = 0.95,
                       sigdig = 4,
                       sigdigTable = NULL,
                       rhoend = NULL,

                       stickyRecalcN = 4,
                       maxOdeRecalc = 5,
                       odeRecalcFactor = 10^(0.5),
                       indTolRelax = TRUE,
                       eventSens = c("jump", "fd"),
                       rxControl = NULL,
                       ...) {

  checkmate::assertIntegerish(seed, any.missing = FALSE, len = 1)
  checkmate::assertIntegerish(itersBurnIn, lower = 0, any.missing = FALSE, len = 1)
  checkmate::assertIntegerish(klWarmup, lower = 0, any.missing = FALSE, len = 1)
  checkmate::assertIntegerish(gammaIter, lower = 0, any.missing = FALSE, len = 1)
  checkmate::assertIntegerish(iters, lower = 1, any.missing = FALSE, len = 1)
  checkmate::assertIntegerish(nGradStep, lower = 1, any.missing = FALSE, len = 1)
  checkmate::assertIntegerish(hiddenDim, lower = 1, any.missing = FALSE, len = 1)
  checkmate::assertNumeric(learningRate, lower = 0, finite = TRUE, any.missing = FALSE, len = 1)
  checkmate::assertNumeric(burnInLearningRate, lower = 0, finite = TRUE, any.missing = FALSE, len = 1)
  if (!is.null(sigma0)) {
    checkmate::assertNumeric(sigma0, lower = 0, finite = TRUE, any.missing = FALSE, min.len = 1)
  }
  checkmate::assertLogical(covariateSelection, len = 1, any.missing = FALSE)
  checkmate::assertLogical(pinCovariates, len = 1, any.missing = FALSE)
  checkmate::assertLogical(muRefCovAlg, len = 1, any.missing = FALSE)
  checkmate::assertNumeric(covSelectAlpha, lower = 1, finite = TRUE, any.missing = FALSE, len = 1)
  checkmate::assertLogical(covSelectSmooth, len = 1, any.missing = FALSE)
  gammaSeries <- match.arg(gammaSeries)
  sigma0Interp <- match.arg(sigma0Interp)
  residOptimize <- match.arg(residOptimize)
  omegaUpdate <- match.arg(omegaUpdate)
  inputScale <- match.arg(inputScale)
  covSelectMethod <- match.arg(covSelectMethod)
  ## Inf is allowed: it forces the exact branch-and-bound everywhere (the
  ## threshold is never reached), so keep it numeric rather than coercing (which
  ## would make it NA)
  checkmate::assertNumeric(covSelectMaxExact, lower = 1, len = 1, any.missing = FALSE)
  if (is.finite(covSelectMaxExact)) covSelectMaxExact <- as.integer(covSelectMaxExact)
  bnbStrategy <- match.arg(bnbStrategy)
  checkmate::assertLogical(parEncoderBackward, len = 1, any.missing = FALSE)
  nonMuTheta <- match.arg(nonMuTheta)
  mStepObjective <- match.arg(mStepObjective)
  checkmate::assertNumeric(nonMuEtaOmega, lower = 0, finite = TRUE, any.missing = FALSE, len = 1)
  checkmate::assertIntegerish(nIsSample, lower = 1, any.missing = FALSE, len = 1)
  checkmate::assertLogical(returnVae, len = 1, any.missing = FALSE)
  checkmate::assertLogical(optExpression, len = 1, any.missing = FALSE)
  checkmate::assertLogical(sumProd, len = 1, any.missing = FALSE)
  checkmate::assertLogical(literalFix, len = 1, any.missing = FALSE)
  checkmate::assertLogical(literalFixRes, len = 1, any.missing = FALSE)
  checkmate::assertLogical(calcTables, len = 1, any.missing = FALSE)
  checkmate::assertLogical(compress, len = 1, any.missing = TRUE)
  checkmate::assertLogical(adjObf, len = 1, any.missing = TRUE)
  checkmate::assertIntegerish(stickyRecalcN, lower = 0, any.missing = FALSE, len = 1)
  checkmate::assertIntegerish(maxOdeRecalc, any.missing = FALSE, len = 1)
  checkmate::assertNumeric(odeRecalcFactor, lower = 1, len = 1, any.missing = FALSE)
  checkmate::assertLogical(indTolRelax, len = 1, any.missing = FALSE)
  likelihood <- match.arg(likelihood)
  objf <- match.arg(objf)
  # match.arg cannot match ""; treat it (skip covariance) like foceiControl does
  if (length(covMethod) == 1L && covMethod == "") {
    covMethod <- ""
  } else {
    covMethod <- match.arg(covMethod)
  }
  addProp <- match.arg(addProp)
  eventSens <- match.arg(eventSens)

  .xtra <- list(...)
  .bad <- names(.xtra)
  .bad <- .bad[!(.bad %in% c("genRxControl", "iterPrintControl"))]
  if (length(.bad) > 0) {
    stop("unused argument: ", paste(paste0("'", .bad, "'"), collapse = ", "),
         call. = FALSE)
  }

  .genRxControl <- FALSE
  if (!is.null(.xtra$genRxControl)) {
    .genRxControl <- .xtra$genRxControl
  }
  if (is.null(rxControl)) {
    if (!is.null(sigdig)) {
      rxControl <- .rxControlScaleSigdig(rxode2::rxControl(sigdig = sigdig), sigdig)
    } else {
      rxControl <- rxode2::rxControl(atol = 1e-4, rtol = 1e-4)
    }
    .genRxControl <- TRUE
  } else if (inherits(rxControl, "rxControl")) {
  } else if (is.list(rxControl)) {
    rxControl <- .rxControlScaleSigdig(do.call(rxode2::rxControl, rxControl), sigdig, skip = names(rxControl))
  } else {
    stop("solving options 'rxControl' needs to be generated from 'rxode2::rxControl'",
         call. = FALSE)
  }
  if (!is.null(sigdig)) {
    checkmate::assertNumeric(sigdig, lower = 1, finite = TRUE, any.missing = TRUE, len = 1)
    if (is.null(sigdigTable)) {
      sigdigTable <- round(sigdig)
    }
  }
  if (is.null(sigdigTable)) {
    sigdigTable <- 3
  }
  checkmate::assertIntegerish(sigdigTable, lower = 1, len = 1, any.missing = FALSE)

  .iterPrintControl <- .absorbIterPrintControl(print = print,
                                               printNcol = printNcol,
                                               useColor = useColor,
                                               iterPrintControl = .xtra$iterPrintControl)

  # inner bounded-bobyqa final trust-region radius for the non-mu/covariate
  # regress M-step; FOCEi mechanism from sigdig, else the sigdig=4 value
  if (is.null(rhoend)) rhoend <- if (!is.null(sigdig)) .sigdigOptTol(sigdig) else 1e-4
  checkmate::assertNumeric(rhoend, len=1, lower=0, finite=TRUE, any.missing=FALSE)
  ## Convergence tolerance of the RESIDUAL optimizer.  Derived from `sigdig` the
  ## same way every other optimizer tolerance in the package is (10^-sigdig), not
  ## inherited from an explicitly-set `rhoend` -- so `sigdig` remains the single
  ## knob that moves all of them together.
  if (is.null(residRhoend)) {
    residRhoend <- if (!is.null(sigdig)) .sigdigOptTol(sigdig) else 1e-4
  }
  checkmate::assertNumeric(residRhoend, len=1, lower=0, finite=TRUE, any.missing=FALSE)
  .ret <- list(seed = as.integer(seed),
               rhoend = as.numeric(rhoend),
               residRhoend = as.numeric(residRhoend),
               itersBurnIn = as.integer(itersBurnIn),
               klWarmup = as.integer(klWarmup),
               gammaIter = as.integer(gammaIter),
               iters = as.integer(iters),
               nGradStep = as.integer(nGradStep),
               hiddenDim = as.integer(hiddenDim),
               learningRate = learningRate,
               burnInLearningRate = burnInLearningRate,
               sigma0 = sigma0,
               covariateSelection = covariateSelection,
               pinCovariates = pinCovariates,
               muRefCovAlg = muRefCovAlg,
               covSelectAlpha = covSelectAlpha,
               covSelectSmooth = covSelectSmooth,
               gammaSeries = gammaSeries,
               sigma0Interp = sigma0Interp,
               residOptimize = residOptimize,
               omegaUpdate = omegaUpdate,
               inputScale = inputScale,
               covSelectMethod = covSelectMethod,
               covSelectMaxExact = covSelectMaxExact,
               bnbStrategy = bnbStrategy,
               parEncoderBackward = parEncoderBackward,
               nonMuTheta = nonMuTheta,
               nonMuEtaOmega = nonMuEtaOmega,
               mStepObjective = mStepObjective,
               likelihood = likelihood,
               objf = objf,
               nIsSample = as.integer(nIsSample),
               returnVae = returnVae,
               covMethod = covMethod,
               optExpression = optExpression,
               sumProd = sumProd,
               literalFix = literalFix,
               literalFixRes = literalFixRes,
               addProp = addProp,
               calcTables = calcTables,
               compress = compress,
               adjObf = adjObf,
               ci = ci,
               sigdig = sigdig,
               sigdigTable = sigdigTable,
               stickyRecalcN = as.integer(stickyRecalcN),
               maxOdeRecalc = as.integer(maxOdeRecalc),
               odeRecalcFactor = odeRecalcFactor,
               indTolRelax = indTolRelax,
               eventSens = eventSens,
               iterPrintControl = .iterPrintControl,
               rxControl = rxControl,
               genRxControl = .genRxControl)
  class(.ret) <- "vaeControl"
  .ret
}

#' @export
rxUiDeparse.vaeControl <- function(object, var) {
  .default <- vaeControl()
  .w <- .deparseDifferent(.default, object, "genRxControl")
  .deparseFinal(.default, object, .w, var)
}

#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.vaeControl <- function(control, env) {
  assign("vaeControl", control, envir = env)
}

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.vae <- function(x, ...) {
  .env <- x[[1]]
  if (exists("vaeControl", .env)) {
    .control <- get("vaeControl", .env)
    if (inherits(.control, "vaeControl")) return(.control)
  }
  if (exists("control", .env)) {
    .control <- get("control", .env)
    if (inherits(.control, "vaeControl")) return(.control)
  }
  stop("cannot find vae related control object", call. = FALSE)
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.vae <- function(control) {
  .ctl <- control[[1]]
  if (is.null(.ctl)) .ctl <- vaeControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list")) .ctl <- do.call("vaeControl", .ctl)
  if (!inherits(.ctl, "vaeControl")) {
    .minfo("invalid control for `est=\"vae\"`, using default")
    .ctl <- vaeControl()
  } else {
    .ctl <- do.call(vaeControl, .ctl)
  }
  .ctl
}

#' @rdname nlmixr2Est
#' @export
nlmixr2Est.vae <- function(env, ...) {
  .ui <- env$ui
  rxode2::assertRxUiRandomOnIdOnly(.ui, " for the estimation routine 'vae'", .var.name = .ui$modelName)
  ## mu-referencing is NOT required: a non-mu-referenced eta is modeled as
  ## theta+eta with theta forced to 0 (see .vaeDataPrep isFree handling)
  ## absorb the validated control (set by getValidNlmixrControl before dispatch)
  if (exists("control", envir = env) && inherits(env$control, "vaeControl")) {
    assign("vaeControl", env$control, envir = env)
  } else {
    assign("vaeControl", vaeControl(), envir = env)
  }
  ## nonMuTheta="grad" needs the analytic outer gradient; out of analytic scope it
  ## must not silently become a no-op, so downgrade to the bobyqa regression once,
  ## up front, and say so in $runInfo.
  if (identical(env$vaeControl$nonMuTheta, "grad")) {
    ## the analytic gradient differentiates the OUTER objective; under the
    ## reference ELBO M-step it would step one functional and score another
    .elbo <- identical(env$vaeControl$mStepObjective, "elbo")
    if (.elbo || !.vaeGradInScope(.ui)) {
      .ctl <- env$vaeControl
      .ctl$nonMuTheta <- "regress"
      assign("vaeControl", .ctl, envir = env)
      assign("control", .ctl, envir = env)
      assign("control", .ctl, envir = .ui)   # the ui copy .analyticGradCaller reads
      warning(if (.elbo) "mStepObjective=\"elbo\": used nonMuTheta=\"regress\""
              else "analytic gradient out of scope; used nonMuTheta=\"regress\"",
              call. = FALSE)
    }
  }
  ## Seed the ENTIRE estimation ONCE here (encoder init, Adam, reparam sampling,
  ## and any random draws in the model / residual-table simulation) and restore
  ## the caller's global RNG state afterward -- a fit never perturbs it, and a
  ## given seed makes the whole process reproducible.
  rxode2::rxWithSeed(env$vaeControl$seed, {
    .fit <- .vaeFitModel(env)
    if (isTRUE(env$vaeControl$returnVae)) .fit else .vaeToFit(env, .fit)
  })
}
attr(nlmixr2Est.vae, "covPresent") <- TRUE
attr(nlmixr2Est.vae, "unbounded") <- FALSE
## enable the IOV preprocessing hook (.uiApplyIov): occasion-level random effects
## are materialized into the model so theta+eta+iov reads as a mu-referenced
## theta+eta expression (the ID-level eta the encoder learns) plus per-occasion
## deviations handled by the inner problem
attr(nlmixr2Est.vae, "iov") <- TRUE
## enable the mu2/mu3/mu4 covariate-rewriting hook (.uiApplyMu2hook, R/mu2.R) so a
## centered/algebraic covariate (e.g. wt.cl*(WT/70) or wt.cl*log(WT/70)) is turned
## into a linear nlmixrMuDerCov# data column -- the centering is carried by the
## mu2/mu3 data, not re-applied by the VAE covariate search -- gated on muRefCovAlg
attr(nlmixr2Est.vae, "mu") <- function(control) isTRUE(control$muRefCovAlg)
