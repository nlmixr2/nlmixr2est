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
#'   reparameterization sampling).
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
#'   covariate selection during training; when `FALSE` fit the given fixed
#'   covariate structure only (faster population-only mode).
#' @param objf Which objective-function value is active for AIC/BIC/BICc. Both
#'   the linearization and importance-sampling -2LL are always computed and
#'   stored; this selects the default active one.
#' @param nIsSample Number of importance-sampling draws for the IS -2LL.
#' @param returnVae When `TRUE` return the raw VAE training object instead of the
#'   nlmixr2 fit.
#'
#' @return vae control structure (class `vaeControl`)
#' @export
#' @author Matthew L. Fidler
vaeControl <- function(seed = 1L,
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
                       likelihood = c("focei", "foce", "focep", "laplace"),
                       objf = c("importanceSampling", "linear"),
                       nIsSample = 3000L,
                       returnVae = FALSE,

                       print = 1L,
                       useColor = NULL,
                       printNcol = NULL,

                       covMethod = c("linear", ""),
                       optExpression = TRUE,
                       sumProd = FALSE,
                       literalFix = TRUE,
                       literalFixRes = TRUE,
                       addProp = c("combined2", "combined1"),
                       calcTables = TRUE,
                       compress = FALSE,
                       adjObf = TRUE,
                       ci = 0.95,
                       sigdig = NULL,
                       sigdigTable = NULL,

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
  covMethod <- match.arg(covMethod)
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
      rxControl <- rxode2::rxControl(sigdig = sigdig)
    } else {
      rxControl <- rxode2::rxControl(atol = 1e-4, rtol = 1e-4)
    }
    .genRxControl <- TRUE
  } else if (inherits(rxControl, "rxControl")) {
  } else if (is.list(rxControl)) {
    rxControl <- do.call(rxode2::rxControl, rxControl)
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

  .ret <- list(seed = as.integer(seed),
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
  .fit <- .vaeFitModel(env)
  if (isTRUE(env$vaeControl$returnVae)) return(.fit)
  .vaeToFit(env, .fit)
}
attr(nlmixr2Est.vae, "covPresent") <- TRUE
attr(nlmixr2Est.vae, "unbounded") <- FALSE
