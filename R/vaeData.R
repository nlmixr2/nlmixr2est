# vaeData.R -- VAE data preparation. Builds, from the rxode2 ui + data:
#  - encoder inputs: standardized (time, DV) sequences padded to Tmax + lengths
#  - per-subject decoder inputs: event table, obs times, observed DV
#  - eta<->theta (z_pop) mapping and initial z_pop / omega / a from ini()
# Encoder-input standardization follows the reference (time/max, (DV-mean)/sd).

# Data columns never treated as covariate candidates by the VAE search
.vaeReservedCols <- c("ID", "TIME", "DV", "EVID", "AMT", "CMT", "MDV", "SS", "II",
                      "ADDL", "RATE", "DUR", "DVID", "CENS", "LIMIT", "OCC")

#' Discover + encode subject-level covariates for the VAE search
#' @param d normalized data.frame (upper-case column names)
#' @param ids unique subject ids in estimation order
#' @return list(covNames, covMat, covType, covPop, tvExcl)
#' @noRd
.vaeCovariateSearch <- function(d, ids) {
  N <- length(ids)
  ## auto-discover subject-level covariate candidates (constant within ID),
  ## excluding reserved data columns. Paper encoding: continuous -> log(v/mean),
  ## categorical (<=2 levels) -> centered.
  .cand <- setdiff(names(d), .vaeReservedCols)
  .isSubjConst <- vapply(.cand, function(nm) {
    all(vapply(ids, function(id) {
      v <- d[[nm]][d$ID == id]; length(unique(v)) == 1L
    }, logical(1))) && is.numeric(d[[nm]])
  }, logical(1))
  # The VAE covariate search absorbs covariates as subject-level (time-invariant)
  # effects, so a covariate that varies within a subject (time-varying) cannot be
  # searched.  Unlike saem/mu-focei (whose covariates are declared in the model,
  # detectable via .nlmixrTimeVaryingCovariates), the VAE search scans every
  # numeric data column, so time-varying ones are those that are not
  # subject-constant.  They are reported back (tvExcl) so callers can warn.
  .covNames <- .cand[.isSubjConst]
  .numCand <- .cand[vapply(.cand, function(nm) is.numeric(d[[nm]]), logical(1))]
  .tvExcl <- setdiff(.numCand, .covNames)
  .covVal <- vapply(.covNames, function(nm) {
    vapply(ids, function(id) d[[nm]][d$ID == id][1], numeric(1))
  }, numeric(N))
  if (length(.covNames) == 1L) .covVal <- matrix(.covVal, N, 1L, dimnames = list(NULL, .covNames))
  .covType <- character(length(.covNames)); .covPop <- numeric(length(.covNames))
  .covMat <- matrix(0, N, length(.covNames), dimnames = list(NULL, .covNames))
  for (j in seq_along(.covNames)) {
    v <- .covVal[, j]
    if (length(unique(v)) > 2L && all(v > 0)) {
      .covType[j] <- "continuous"; .covPop[j] <- mean(v); .covMat[, j] <- log(v / .covPop[j])
    } else {
      .covType[j] <- "categorical"; .covPop[j] <- mean(v); .covMat[, j] <- v - .covPop[j]
    }
  }
  list(covNames = .covNames, covMat = .covMat, covType = .covType,
       covPop = .covPop, tvExcl = .tvExcl)
}

#' Covariates explored by the VAE covariate search
#'
#' Returns the subject-level covariates that `nlmixr2(..., est = "vae")` would
#' explore during automated covariate selection, using the same discovery rules
#' as the fit: every non-reserved numeric data column that is constant within
#' each subject is a candidate; a candidate with more than two unique values
#' (all positive) is treated as continuous (encoded `log(value/mean)`),
#' anything else as categorical (mean-centered).  Time-varying numeric columns
#' cannot be searched and are excluded with a warning.
#'
#' @param data estimation dataset containing at least an `ID` column; column
#'   names are matched case-insensitively, as in the VAE fit
#' @param warn when `TRUE` (default) warn about time-varying numeric columns
#'   excluded from the search; when `FALSE` exclude them silently
#' @return a data frame with one row per explored covariate and columns
#'   `covariate` (upper-cased column name), `type` (`"continuous"` or
#'   `"categorical"`) and `center` (the population value the covariate is
#'   centered at); zero rows when no covariates qualify
#' @export
#' @author Matthew L. Fidler
#' @examples
#' d <- data.frame(id = rep(1:3, each = 2), time = rep(0:1, 3), dv = rnorm(6),
#'                 wt = rep(c(70, 80, 60), each = 2),
#'                 sex = rep(c(0, 1, 0), each = 2))
#' vaeCovariates(d)
vaeCovariates <- function(data, warn = TRUE) {
  checkmate::assertLogical(warn, len = 1, any.missing = FALSE)
  d <- as.data.frame(data)
  names(d) <- toupper(names(d))
  if (is.null(d$ID)) {
    stop("'data' must contain an ID column", call. = FALSE)
  }
  .cov <- .vaeCovariateSearch(d, unique(d$ID))
  if (warn && length(.cov$tvExcl) > 0L) {
    warning("time-varying covariate(s) were excluded from automatic covariate search: ",
            paste(.cov$tvExcl, collapse = ", "), call. = FALSE)
  }
  data.frame(covariate = .cov$covNames, type = .cov$covType,
             center = .cov$covPop, row.names = NULL)
}

#' Prepare VAE inputs from a ui + data
#' @param ui rxode2 ui object
#' @param data estimation data (ID/TIME/DV/EVID/... columns)
#' @return list of prepared VAE inputs
#' @noRd
.vaeDataPrep <- function(ui, data) {
  .idf <- ui$iniDf
  .map <- .foceiEtaThetaMap(ui)
  .etaNames <- .map$etaNames
  .neta <- length(.etaNames)
  if (.neta == 0L) stop("est=\"vae\" requires at least one random effect", call. = FALSE)

  ## full theta vector (THETA_i_ in ntheta order), from ini estimates
  .thRows <- .idf[!is.na(.idf$ntheta), , drop = FALSE]
  .thRows <- .thRows[order(.thRows$ntheta), , drop = FALSE]
  .th <- setNames(as.numeric(.thRows$est), paste0("THETA_", seq_len(nrow(.thRows)), "_"))
  ## structural theta index (in the full theta vector) paired with each eta.
  ## A random effect that is not mu-referenced to a single theta -- a mixture eta
  ## (mix(exp(lke1+eta.ke),p,exp(lke2+eta.ke))), an eta on a fixed (literalFix-ed)
  ## theta, or a genuinely free eta -- is modeled as theta+eta with theta forced
  ## to 0: it centers at 0, is held there by the M-step, and is excluded from
  ## covariate selection; the rest of the model (literal / component thetas /
  ## covariate expression) carries its structure.
  .zPopThetaIdx <- match(.map$thetaForEta, .thRows$name)
  .isFree <- is.na(.zPopThetaIdx)
  .zPop <- numeric(.neta)                                      # structural population means (transformed)
  .zPop[!.isFree] <- as.numeric(.th[.zPopThetaIdx[!.isFree]])

  ## omega init (diagonal) for the etas + which variances are FIXED (held by the
  ## M-step, not estimated)
  .omega <- vapply(.etaNames, function(nm) {
    .r <- .idf[!is.na(.idf$neta1) & .idf$neta1 == .idf$neta2 & .idf$name == nm, , drop = FALSE]
    as.numeric(.r$est[1])
  }, numeric(1))
  .omegaFix <- vapply(.etaNames, function(nm) {
    .r <- .idf[!is.na(.idf$neta1) & .idf$neta1 == .idf$neta2 & .idf$name == nm, , drop = FALSE]
    isTRUE(as.logical(.r$fix[1]))
  }, logical(1))
  ## structural-theta bounds per eta (Inf/-Inf when unbounded or free): the M-step
  ## clamps the population estimate to [lower, upper], giving the constrained
  ## estimate (at the bound when the unconstrained optimum is outside).
  .zPopLower <- rep(-Inf, .neta); .zPopUpper <- rep(Inf, .neta)
  .zPopLower[!.isFree] <- as.numeric(.thRows$lower[.zPopThetaIdx[!.isFree]])
  .zPopUpper[!.isFree] <- as.numeric(.thRows$upper[.zPopThetaIdx[!.isFree]])

  ## residual error params (all of them, in theta order): value, theta index,
  ## type (add/prop/...), and bounds. Combined models have >1 row; log-likelihood
  ## models may have none. `a` is the (named) error-param vector.
  .errRow <- .idf[!is.na(.idf$err) & !is.na(.idf$ntheta), , drop = FALSE]
  .errRow <- .errRow[order(.errRow$ntheta), , drop = FALSE]
  .a <- if (nrow(.errRow) > 0) setNames(as.numeric(.errRow$est), .errRow$name) else numeric(0)
  .errThetaIdx <- as.integer(.errRow$ntheta)
  .errType <- as.character(.errRow$err)
  .errLower <- as.numeric(.errRow$lower); .errUpper <- as.numeric(.errRow$upper)

  ## normalize data columns
  d <- as.data.frame(data)
  names(d) <- toupper(names(d))
  if (is.null(d$EVID)) d$EVID <- ifelse(is.na(d$AMT) | d$AMT == 0, 0L, 1L)
  .ids <- unique(d$ID)
  N <- length(.ids)

  ## per-subject decoder inputs + gather all obs for standardization
  subj <- vector("list", N)
  .allTime <- numeric(0); .allDv <- numeric(0)
  for (i in seq_len(N)) {
    .di <- d[d$ID == .ids[i], , drop = FALSE]
    .obs <- .di[.di$EVID == 0, , drop = FALSE]
    .times <- .obs$TIME
    .y <- .obs$DV
    ## M2/M3/M4 censoring columns (0 / NA when absent)
    .cens <- if (is.null(.obs$CENS)) integer(length(.y)) else as.integer(.obs$CENS)
    .limit <- if (is.null(.obs$LIMIT)) rep(NA_real_, length(.y)) else as.numeric(.obs$LIMIT)
    subj[[i]] <- list(ev = .di, times = .times, y = .y, n = length(.times),
                      cens = .cens, limit = .limit)
    .allTime <- c(.allTime, .times); .allDv <- c(.allDv, .y)
  }
  .tMax <- max(.allTime); .dvMean <- mean(.allDv); .dvSd <- stats::sd(.allDv)

  ## encoder inputs: [N, Tmax, 2] standardized (time, DV), padded; lengths
  Tmax <- max(vapply(subj, function(s) s$n, integer(1)))
  dataIn <- array(0, c(N, Tmax, 2L))
  lengths <- integer(N)
  for (i in seq_len(N)) {
    s <- subj[[i]]; ni <- s$n; lengths[i] <- ni
    dataIn[i, seq_len(ni), 1L] <- s$times / .tMax
    dataIn[i, seq_len(ni), 2L] <- (s$y - .dvMean) / .dvSd
  }
  covIn <- matrix(0, N, 0L)                     # encoder-head covariates (unused for now)

  ## subject-level covariate discovery + encoding (shared with vaeCovariates())
  .cov <- .vaeCovariateSearch(d, .ids)
  if (length(.cov$tvExcl) > 0L) {
    warning("time-varying covariate(s) were excluded from automatic covariate search: ",
            paste(.cov$tvExcl, collapse = ", "), call. = FALSE)
  }

  list(N = N, neta = .neta, zDim = .neta, etaNames = .etaNames,
       th = .th, zPopThetaIdx = .zPopThetaIdx, isFree = .isFree, omegaFix = .omegaFix,
       zPopLower = .zPopLower, zPopUpper = .zPopUpper,
       errThetaIdx = .errThetaIdx, errType = .errType,
       errLower = .errLower, errUpper = .errUpper,
       zPop = .zPop, omega = .omega, a = .a,
       subj = subj, dataIn = dataIn, lengths = lengths, covIn = covIn,
       covNames = .cov$covNames, covMat = .cov$covMat, covType = .cov$covType,
       covPop = .cov$covPop,
       tMax = .tMax, dvMean = .dvMean, dvSd = .dvSd, Nobs = length(.allDv))
}
