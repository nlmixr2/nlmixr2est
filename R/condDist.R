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
#' Mixture, IOV, non-parametric and repeated (`same()`) Omega block fits are
#' not supported.
#'
#' `method = "saem"` uses the SAEM MCMC kernels.  For a SAEM fit, the chains of
#' its own last iterations are the draws when they are enough (`nsamp` draws
#' from `nmc` chains, within the last half of the second-phase iterations,
#' where the parameters have nearly settled), so nothing is simulated;
#' otherwise the kernels run with every population parameter held
#' at the fit's estimates, continuing a SAEM fit's chains or, for any other
#' fit, starting cold, and `burn` iterations are discarded; it needs a model
#' SAEM can fit (every eta on a population parameter of its own).  An eta SAEM does
#' not sample for a subject (uninformative) has `SD` `NA` and `MEAN` 0, its
#' prior, as in `fit$eta`; that is exact when the eta is uncorrelated with the
#' others.  A model SAEM fits with temporary etas (a likelihood parameter with
#' no eta) is sampled with `"imp"`, since those chains also sampled parameters
#' the reported model holds fixed.
#'
#' @param fit nlmixr2 fit with random effects
#' @param method the sampler: `"auto"` (default) is `"saem"` for a SAEM fit
#'   and `"imp"` otherwise
#' @param nsamp number of samples per subject (default 300)
#' @param keepSamples Boolean (default `FALSE`); keep `nsamp` equally weighted
#'   draws per subject as the `"samples"` attribute of `fit$condDist`, an
#'   `ID x eta x draw` array
#' @param burn MCMC iterations discarded by a held `"saem"` run (default 50)
#' @param seed random seed of the sampler
#' @return the fit, with `fit$condDistShrinkage`, the eta shrinkage (%) from
#'   the conditional variances, `100 * mean(SD^2) / omega` (an unsampled
#'   subject's `SD^2` taken as `omega`), and
#'   `fit$condDist`: a data frame `ID`, `ETA`, `MEAN`,
#'   `SD`.  Its `"method"` attribute names the sampler (`"imp"`, `"saem"`
#'   for a SAEM fit's own chains or `"mcmc"` for a held SAEM run) and, for
#'   `"imp"`, its `"ess"` attribute holds each subject's effective sample size.
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
addCondDist <- function(
  fit,
  method = c("auto", "imp", "saem"),
  nsamp = 300L,
  keepSamples = FALSE,
  burn = 50L,
  seed = 42L
) {
  method <- match.arg(method)
  checkmate::assertIntegerish(nsamp, lower = 2, upper = 1e7, len = 1, any.missing = FALSE)
  checkmate::assertIntegerish(burn, lower = 0, upper = 1e7, len = 1, any.missing = FALSE)
  checkmate::assertLogical(keepSamples, len = 1, any.missing = FALSE)
  checkmate::assertIntegerish(seed, len = 1, any.missing = FALSE)
  .eta <- .condDistCheckFit(fit)
  .env <- fit$env
  # a SAEM fit with temporary etas sampled parameters the reported model holds
  # fixed, so its chains are not this model's conditional distribution
  .pseudo <- length(tryCatch(fit$env$saemPseudoEtas, error = function(e) NULL)) > 0L
  if (method == "auto") {
    .isSaem <- identical(tryCatch(fit$est, error = function(e) NULL), "saem")
    method <- if (.isSaem && !.pseudo) "saem" else "imp"
  }
  if (method == "saem" && .pseudo) {
    stop("SAEM sampled temporary etas for this fit; use method = \"imp\"", call. = FALSE)
  }
  nsamp <- as.integer(round(nsamp))
  if (method == "imp") {
    .r <- .condDistImp(fit, nsamp, if (keepSamples) nsamp else 0L, as.integer(seed))
  } else {
    .r <- .condDistSaem(fit, nsamp, as.integer(round(burn)), as.integer(seed))
    .r$mean <- apply(.r$draws, c(1L, 2L), mean)
    .r$sd <- apply(.r$draws, c(1L, 2L), stats::sd)
    # not sampled: the eta stays at its prior mean, as in fit$eta
    .r$mean[is.na(.r$mean)] <- 0
    if (!keepSamples) {
      .r$draws <- NULL
    }
  }
  .etaNames <- names(.eta)[-1]
  if (!identical(dim(.r$mean), c(nrow(.eta), length(.etaNames)))) {
    stop("conditional distribution has an unexpected shape", call. = FALSE)
  }
  .ret <- .condDistSummary(.r$mean, .r$sd, .eta$ID, .etaNames, .r$draws)
  attr(.ret, "method") <- if (is.null(.r$method)) method else .r$method
  attr(.ret, "ess") <- .r$ess
  assign("condDist", .ret, envir = .env)
  assign("condDistShrinkage", .condDistShrinkage(.ret, fit$omega), envir = .env)
  fit
}

#' Refuse fits whose conditional eta distribution cannot be sampled
#' @param fit nlmixr2 fit
#' @return the fit's `$eta` data frame
#' @noRd
.condDistCheckFit <- function(fit) {
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
  if (isTRUE(tryCatch(fit$est, error = function(e) "") %in% c("npag", "mnpag", "inpag", "npb", "mnpb", "inpb"))) {
    stop("conditional distribution is not available for non-parametric fits", call. = FALSE)
  }
  .om <- tryCatch(as.matrix(fit$omega), error = function(e) NULL)
  if (!is.matrix(.om) || !all(is.finite(.om)) || inherits(try(chol(.om), silent = TRUE), "try-error")) {
    stop("conditional distribution needs a positive-definite Omega", call. = FALSE)
  }
  .eta
}

#' Distribution-based eta shrinkage from the conditional SDs
#'
#' `100 * mean_i(Var(eta_i | y_i)) / omega`, which by the law of total variance
#' equals the variance shrinkage `1 - Var(E(eta | y)) / omega` without reading
#' the spread of shrunk point estimates.  A subject whose eta is not sampled
#' (`SD` `NA`) keeps its prior variance, so counts as fully shrunk (exact
#' when the eta is uncorrelated with the others); an eta with no positive
#' variance is `NA`.
#' @param condDist `fit$condDist`
#' @param omega the fit's Omega
#' @return named numeric, percent per eta
#' @noRd
.condDistShrinkage <- function(condDist, omega) {
  .eta <- levels(condDist$ETA)
  vapply(
    stats::setNames(.eta, .eta),
    function(e) {
      .om <- if (is.matrix(omega) && e %in% rownames(omega)) omega[e, e] else NA_real_
      if (!is.finite(.om) || .om <= 0) {
        return(NA_real_)
      }
      .v <- condDist$SD[condDist$ETA == e]^2
      .v[is.na(.v)] <- .om
      100 * mean(.v) / .om
    },
    numeric(1)
  )
}

#' Importance-sampled conditional distribution at a fit's estimates
#'
#' @param fit nlmixr2 fit
#' @param nsamp samples per subject
#' @param ndraw equally weighted draws per subject to return (0 for none)
#' @param seed sampler seed
#' @return list: `mean`, `sd` (`N x neta`), `ess` (length `N`) and `draws`,
#'   an `N x neta x ndraw` array of equally weighted draws or `NULL`
#' @noRd
.condDistImp <- function(fit, nsamp, ndraw, seed) {
  # imp raises an Omega variance below 1e-6 to 1e-6, which would change the prior
  .om <- tryCatch(diag(as.matrix(fit$omega)), error = function(e) NULL)
  if (any(.om < 1e-6)) {
    stop("importance sampling needs every Omega variance >= 1e-6; try method = \"saem\"", call. = FALSE)
  }
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
    stop("importance sampling at the fit's estimates failed: ", .tryErrorMessage(.fit2), call. = FALSE)
  }
  .e <- .fit2$env
  .u <- rxode2::rxWithSeed(seed, rxode2::rxunif(n = length(.e$impSamples)))
  .weightedDraws(.e$impSamples, .e$impWeights, ndraw, .u)
}

#' The message of a `try-error`
#' @param x a `try-error`
#' @return character
#' @noRd
.tryErrorMessage <- function(x) {
  .c <- attr(x, "condition")
  trimws(if (inherits(.c, "condition")) conditionMessage(.c) else as.character(x))
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
#' @param ndraw equally weighted draws per subject (0 for none)
#' @param u per-subject offsets in `[0, 1)`
#' @return list: `mean`, `sd`, `ess` and `draws` (`N x neta x ndraw`, or
#'   `NULL` when `ndraw` is 0)
#' @noRd
.weightedDraws <- function(samples, weights, ndraw, u = rep(0.5, length(samples))) {
  .n <- length(samples)
  .neta <- ncol(samples[[1]])
  .mean <- .sd <- matrix(NA_real_, .n, .neta)
  .ess <- numeric(.n)
  .draws <- if (ndraw > 0L) array(NA_real_, c(.n, .neta, ndraw))
  for (.i in seq_len(.n)) {
    .s <- samples[[.i]]
    .w <- as.numeric(weights[[.i]])
    .w <- .w / sum(.w)
    .m <- colSums(.s * .w)
    .mean[.i, ] <- .m
    .sd[.i, ] <- sqrt(colSums(.w * sweep(.s, 2L, .m)^2))
    .ess[.i] <- 1 / sum(.w^2)
    if (ndraw > 0L) {
      .u <- (seq_len(ndraw) - 1 + u[.i]) / ndraw
      .idx <- pmin(findInterval(.u, cumsum(.w)) + 1L, nrow(.s))
      .draws[.i, , ] <- t(.s[.idx, , drop = FALSE])
    }
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

#' @rdname nmObjGet
#' @export
nmObjGet.condDistShrinkage <- function(x, ...) {
  .env <- x[[1]]$env
  if (exists("condDistShrinkage", envir = .env, inherits = FALSE)) {
    return(get("condDistShrinkage", envir = .env, inherits = FALSE))
  }
  NULL
}
attr(nmObjGet.condDistShrinkage, "desc") <- "eta shrinkage (%) from the conditional variances (addCondDist())"
