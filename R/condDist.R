#' Add each subject's conditional distribution of the etas to a fit
#'
#' Samples each subject's conditional (posterior) distribution of the random
#' effects, `p(eta | y, theta, Omega, sigma)`, at the fit's final population
#' estimates, and stores its mean and SD as `fit$condDist`.  It applies to a
#' fit from any mixed-effects estimation method.
#'
#' `method = "imp"` draws importance samples from a proposal centered at each
#' subject's mode, holding every population parameter at the fit's estimate
#' (the frozen `est = "imp"` E-step), and weights them by the posterior.
#' Mixture, IOV and non-parametric fits are not supported.
#'
#' @param fit nlmixr2 fit with random effects
#' @param method `"auto"` (default) or `"imp"`: the sampler
#' @param nsamp number of samples per subject (default 300)
#' @param keepSamples Boolean (default `FALSE`); keep `nsamp` equally weighted
#'   draws per subject as the `"samples"` attribute of `fit$condDist`, an
#'   `ID x eta x draw` array
#' @param seed random seed of the sampler
#' @return the fit, with `fit$condDist`: a data frame `ID`, `ETA`, `MEAN`,
#'   `SD`.  Its `"method"` attribute names the sampler and its `"ess"`
#'   attribute holds each subject's effective sample size.
#' @export
#' @author Matthew L. Fidler
#' @examples
#' \donttest{
#' one.cmt <- function() {
#'   ini({
#'     tka <- 0.45
#'     tcl <- 1
#'     tv <- 3.45
#'     eta.ka ~ 0.6
#'     eta.cl ~ 0.3
#'     eta.v ~ 0.1
#'     add.sd <- 0.7
#'   })
#'   model({
#'     ka <- exp(tka + eta.ka)
#'     cl <- exp(tcl + eta.cl)
#'     v <- exp(tv + eta.v)
#'     linCmt() ~ add(add.sd)
#'   })
#' }
#' fit <- nlmixr2(one.cmt, nlmixr2data::theo_sd, est = "focei", control = list(print = 0))
#' fit <- addCondDist(fit)
#' head(fit$condDist)
#' }
addCondDist <- function(fit, method = c("auto", "imp"), nsamp = 300L, keepSamples = FALSE, seed = 42L) {
  method <- match.arg(method)
  checkmate::assertIntegerish(nsamp, lower = 2, len = 1, any.missing = FALSE)
  checkmate::assertLogical(keepSamples, len = 1, any.missing = FALSE)
  checkmate::assertIntegerish(seed, len = 1, any.missing = FALSE)
  .env <- tryCatch(fit$env, error = function(e) NULL)
  if (!is.environment(.env)) {
    stop("'fit' needs to be an nlmixr2 fit", call. = FALSE)
  }
  .eta <- tryCatch(fit$eta, error = function(e) NULL)
  if (!is.data.frame(.eta) || ncol(.eta) < 2L) {
    stop("'fit' has no random effects", call. = FALSE)
  }
  if (length(tryCatch(fit$ui$mixProbs, error = function(e) NULL)) > 0L) {
    stop("conditional distribution is not available for mixture models", call. = FALSE)
  }
  if (isTRUE(tryCatch(.uiHasIov(fit$ui), error = function(e) FALSE))) {
    stop("conditional distribution is not available for IOV models", call. = FALSE)
  }
  # their eta distribution is the fitted discrete support, not Omega
  if (tryCatch(fit$est, error = function(e) "") %in% c("npag", "mnpag", "inpag", "npb", "mnpb", "inpb")) {
    stop("conditional distribution is not available for non-parametric fits", call. = FALSE)
  }
  .r <- .condDistImp(fit, as.integer(round(nsamp)), as.integer(seed))
  .etaNames <- names(.eta)[-1]
  if (!identical(dim(.r$mean), c(nrow(.eta), length(.etaNames)))) {
    stop("conditional distribution has an unexpected shape", call. = FALSE)
  }
  .ret <- .condDistSummary(.r$mean, .r$sd, .eta$ID, .etaNames, if (keepSamples) .r$draws)
  attr(.ret, "method") <- "imp"
  attr(.ret, "ess") <- .r$ess
  assign("condDist", .ret, envir = .env)
  fit
}

#' Importance-sampled conditional distribution at a fit's estimates
#'
#' @param fit nlmixr2 fit
#' @param nsamp samples per subject
#' @param seed sampler seed
#' @return list: `mean`, `sd` (`N x neta`), `ess` (length `N`) and `draws`,
#'   an `N x neta x nsamp` array of equally weighted draws
#' @noRd
.condDistImp <- function(fit, nsamp, seed) {
  .fc <- .condDistLikControl(fit)
  # mapIter = 1 centers each E-step's proposal at the subject's mode with the
  # Laplace covariance (mapIter = 0 keeps the prior as the proposal); the
  # E-steps before the last adapt its scale
  .ctl <- impmapControl(
    nIter = 5L,
    mapIter = 1L,
    isample = nsamp,
    impSeed = seed,
    sirSample = min(25L, nsamp),
    covMethod = "",
    calcTables = FALSE,
    print = 0L,
    addProp = .fc$addProp,
    rxControl = .fc$rxControl,
    impFrozen = TRUE
  )
  .fit2 <- .nativePinnedRefit(fit, "imp", .ctl, useEtaMat = TRUE)
  if (inherits(.fit2, "try-error")) {
    stop("importance sampling at the fit's estimates failed", call. = FALSE)
  }
  .e <- .fit2$env
  .u <- rxode2::rxWithSeed(seed, rxode2::rxunif(n = length(.e$impSamples)))
  .weightedDraws(.e$impSamples, .e$impWeights, nsamp, .u)
}

#' The likelihood settings of a fit that a refit at its estimates must keep
#' @param fit nlmixr2 fit
#' @return list: `addProp` and `rxControl` (`NULL` for the engine default)
#' @noRd
.condDistLikControl <- function(fit) {
  .fc <- tryCatch(fit$foceiControl, error = function(e) NULL)
  .addProp <- if (is.list(.fc)) .fc$addProp
  list(
    addProp = if (is.character(.addProp) && length(.addProp) == 1L) .addProp else "combined2",
    rxControl = if (is.list(.fc) && inherits(.fc$rxControl, "rxControl")) .fc$rxControl
  )
}

#' Moments and equally weighted draws of per-subject weighted samples
#'
#' Resampling is systematic: one uniform offset per subject, shared by its
#' `ndraw` equally spaced strata.
#' @param samples list, per subject, of `nsample x neta` matrices
#' @param weights list, per subject, of the samples' weights
#' @param ndraw equally weighted draws per subject
#' @param u per-subject offsets in `[0, 1)`
#' @return list: `mean`, `sd`, `ess` and `draws` (`N x neta x ndraw`)
#' @noRd
.weightedDraws <- function(samples, weights, ndraw, u = rep(0.5, length(samples))) {
  .n <- length(samples)
  .neta <- ncol(samples[[1]])
  .mean <- .sd <- matrix(NA_real_, .n, .neta)
  .ess <- numeric(.n)
  .draws <- array(NA_real_, c(.n, .neta, ndraw))
  for (.i in seq_len(.n)) {
    .u <- (seq_len(ndraw) - 1 + u[.i]) / ndraw
    .s <- samples[[.i]]
    .w <- as.numeric(weights[[.i]])
    .w <- .w / sum(.w)
    .m <- colSums(.s * .w)
    .mean[.i, ] <- .m
    .sd[.i, ] <- sqrt(colSums(.w * sweep(.s, 2L, .m)^2))
    .ess[.i] <- 1 / sum(.w^2)
    .idx <- pmin(findInterval(.u, cumsum(.w)) + 1L, nrow(.s))
    .draws[.i, , ] <- t(.s[.idx, , drop = FALSE])
  }
  list(mean = .mean, sd = .sd, ess = .ess, draws = .draws)
}

#' Conditional moments as the long `fit$condDist` data frame
#'
#' @param mean,sd `N x neta` matrices
#' @param id subject IDs, in fit order
#' @param etaNames eta names
#' @param draws `NULL`, or the `N x neta x ndraw` draws to keep as the
#'   `"samples"` attribute
#' @return data frame `ID`, `ETA`, `MEAN`, `SD`, by ID then eta
#' @noRd
.condDistSummary <- function(mean, sd, id, etaNames, draws = NULL) {
  .n <- length(etaNames)
  .ret <- data.frame(
    ID = rep(id, each = .n),
    ETA = factor(rep(etaNames, length(id)), levels = etaNames),
    MEAN = as.vector(t(mean)),
    SD = as.vector(t(sd))
  )
  if (!is.null(draws)) {
    dimnames(draws) <- list(as.character(id), etaNames, NULL)
    attr(.ret, "samples") <- draws
  }
  .ret
}

#' @rdname nmObjGet
#' @export
nmObjGet.condDist <- function(x, ...) {
  .env <- x[[1]]$env
  if (exists("condDist", envir = .env, inherits = FALSE)) {
    return(get("condDist", envir = .env, inherits = FALSE))
  }
  NULL
}
attr(nmObjGet.condDist, "desc") <- "per-subject conditional eta mean/SD (addCondDist())"
