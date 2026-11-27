#' Draw each subject's etas from their conditional distribution
#'
#' Draws `n` samples per subject from the conditional (posterior)
#' distribution of the random effects, `p(eta | y, theta, Omega, sigma)`, at
#' the fit's final population estimates.  Unlike the EBEs, the draws carry
#' each subject's uncertainty, so they can be used to simulate that subject
#' rather than a new one.
#'
#' `method = "laplace"` draws from a normal centered at each subject's
#' conditional mode with the inverse of the inner Hessian (`fit$phiC`) as its
#' covariance.  A fit whose etas are not conditional modes (for example SAEM
#' or emvi) gets both from a posthoc FOCEi step at its estimates.
#' `method = "imp"` resamples the importance samples of [addCondDist()].
#' `method = "auto"` is `"laplace"`.
#'
#' Mixture, IOV and non-parametric fits are
#' not supported, nor is an Omega that is not positive definite.
#'
#' @param fit nlmixr2 fit with random effects
#' @param n number of draws per subject
#' @param method sampler: `"auto"`, `"laplace"` or `"imp"`
#' @param seed random seed; the same seed gives the same draws whatever the
#'   number of cores
#' @return a data frame with columns `ID`, `draw` and one column per eta,
#'   ordered by `draw` then `ID`.  Passed as `params` to
#'   `rxode2::rxSolve(fit, params = draws, nStud = n, simVariability = FALSE)`
#'   it simulates each subject once per draw (`sim.id` is the draw).  Its
#'   `"method"` attribute names the sampler.
#' @seealso [addCondDist()]
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
#' smp <- etaSamples(fit, n = 100)
#' # each subject's conditional SD of eta.cl
#' aggregate(eta.cl ~ ID, smp, sd)
#' # simulate each subject with its own uncertainty
#' s <- rxode2::rxSolve(fit, params = smp, nStud = 100, simVariability = FALSE)
#' }
etaSamples <- function(fit, n = 500L, method = c("auto", "laplace", "imp"), seed = 42L) {
  method <- match.arg(method)
  checkmate::assertIntegerish(n, lower = 1, len = 1, any.missing = FALSE)
  checkmate::assertIntegerish(seed, len = 1, any.missing = FALSE)
  .eta <- .condDistCheckFit(fit)
  n <- as.integer(round(n))
  seed <- as.integer(seed)
  if (method == "auto") {
    method <- "laplace"
  }
  .draws <- switch(
    method,
    laplace = .etaSamplesLaplace(fit, n, seed),
    imp = .condDistImp(fit, max(n, 300L), n, seed)$draws
  )
  .etaNames <- names(.eta)[-1]
  if (!identical(dim(.draws), c(nrow(.eta), length(.etaNames), n))) {
    stop("eta samples have an unexpected shape", call. = FALSE)
  }
  .ret <- .etaSamplesDf(.draws, .eta$ID, .etaNames)
  attr(.ret, "method") <- method
  .ret
}

#' Long data frame of per-subject eta draws
#' @param draws `N x neta x n` array
#' @param id subject IDs, in fit order
#' @param etaNames eta names
#' @return data frame `ID`, `draw`, one column per eta, by draw then ID
#' @noRd
.etaSamplesDf <- function(draws, id, etaNames) {
  .d <- dim(draws)
  # draw-major rows: subject varies fastest, as rxSolve() reads nStud params
  .m <- matrix(aperm(draws, c(1L, 3L, 2L)), .d[1] * .d[3], .d[2])
  colnames(.m) <- etaNames
  cbind(
    data.frame(ID = rep(id, .d[3]), draw = rep(seq_len(.d[3]), each = .d[1])),
    as.data.frame(.m)
  )
}

#' Estimation methods whose etas and `phiC` are the conditional modes and the
#' inverse inner Hessian there
#' @noRd
.etaSamplesModeEst <- c("focei", "foce", "laplace", "agq", "posthoc")

#' Each subject's conditional mode and Laplace covariance
#' @param fit nlmixr2 fit
#' @return list: `eta` (`N x neta` matrix) and `cov` (list of `neta x neta`
#'   matrices, in fit order)
#' @noRd
.etaSamplesLaplaceCov <- function(fit) {
  .etaNames <- names(fit$eta)[-1]
  .ok <- function(f) {
    .c <- tryCatch(f$phiC, error = function(e) NULL)
    is.list(.c) &&
      length(.c) == nrow(f$eta) &&
      all(vapply(.c, function(m) is.matrix(m) && identical(colnames(m), .etaNames) && all(is.finite(m)), logical(1)))
  }
  .src <- fit
  if (!(tryCatch(fit$est, error = function(e) "") %in% .etaSamplesModeEst) || !.ok(fit)) {
    .fc <- .condDistLikControl(fit)
    .ctl <- foceiControl(
      print = 0L,
      covMethod = "",
      calcTables = FALSE,
      addProp = .fc$addProp,
      rxControl = .fc$rxControl
    )
    .src <- .nativePinnedRefit(fit, "posthoc", .ctl, useEtaMat = TRUE)
    if (inherits(.src, "try-error")) {
      stop("posthoc step at the fit's estimates failed: ", .tryErrorMessage(.src), call. = FALSE)
    }
    if (!.ok(.src)) {
      stop("posthoc step gave no usable inner Hessian", call. = FALSE)
    }
  }
  list(eta = as.matrix(.src$eta[, .etaNames, drop = FALSE]), cov = .src$phiC)
}

#' Laplace draws: normal at each subject's mode with the inverse inner Hessian
#' @param fit nlmixr2 fit
#' @param n draws per subject
#' @param seed seed
#' @return `N x neta x n` array
#' @noRd
.etaSamplesLaplace <- function(fit, n, seed) {
  .lc <- .etaSamplesLaplaceCov(fit)
  .nid <- nrow(.lc$eta)
  .neta <- ncol(.lc$eta)
  .z <- rxode2::rxWithSeed(seed, rxode2::rxnorm(n = .nid * .neta * n))
  dim(.z) <- c(.neta, n, .nid)
  .ret <- array(NA_real_, c(.nid, .neta, n))
  for (.i in seq_len(.nid)) {
    .r <- .etaSamplesChol(.lc$cov[[.i]])
    .ret[.i, , ] <- .lc$eta[.i, ] + crossprod(.r, matrix(.z[,, .i], .neta, n))
  }
  .ret
}

#' Upper Cholesky factor of a covariance, allowing zero-variance etas
#' @param cov covariance matrix
#' @return `R` with `t(R) %*% R == cov`
#' @noRd
.etaSamplesChol <- function(cov) {
  .r <- suppressWarnings(chol(cov, pivot = TRUE))
  if (attr(.r, "rank") < nrow(cov)) {
    .r[-seq_len(attr(.r, "rank")), ] <- 0
  }
  .r[, order(attr(.r, "pivot")), drop = FALSE]
}
