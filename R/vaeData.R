# vaeData.R -- VAE data preparation. Builds, from the rxode2 ui + data:
#  - encoder inputs: standardized (time, DV) sequences padded to Tmax + lengths
#  - per-subject decoder inputs: event table, obs times, observed DV
#  - eta<->theta (z_pop) mapping and initial z_pop / omega / a from ini()
# Encoder-input standardization follows the reference (time/max, (DV-mean)/sd).

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

  ## auto-discover subject-level covariate candidates (constant within ID),
  ## excluding reserved data columns. Paper encoding: continuous -> log(v/mean),
  ## categorical (<=2 levels) -> centered.
  .reserved <- c("ID", "TIME", "DV", "EVID", "AMT", "CMT", "MDV", "SS", "II",
                 "ADDL", "RATE", "DUR", "DVID", "CENS", "LIMIT", "OCC")
  .cand <- setdiff(names(d), .reserved)
  .isSubjConst <- vapply(.cand, function(nm) {
    all(vapply(.ids, function(id) {
      v <- d[[nm]][d$ID == id]; length(unique(v)) == 1L
    }, logical(1))) && is.numeric(d[[nm]])
  }, logical(1))
  .covNames <- .cand[.isSubjConst]
  .covVal <- vapply(.covNames, function(nm) {
    vapply(.ids, function(id) d[[nm]][d$ID == id][1], numeric(1))
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

  list(N = N, neta = .neta, zDim = .neta, etaNames = .etaNames,
       th = .th, zPopThetaIdx = .zPopThetaIdx, isFree = .isFree, omegaFix = .omegaFix,
       zPopLower = .zPopLower, zPopUpper = .zPopUpper,
       errThetaIdx = .errThetaIdx, errType = .errType,
       errLower = .errLower, errUpper = .errUpper,
       zPop = .zPop, omega = .omega, a = .a,
       subj = subj, dataIn = dataIn, lengths = lengths, covIn = covIn,
       covNames = .covNames, covMat = .covMat, covType = .covType, covPop = .covPop,
       tMax = .tMax, dvMean = .dvMean, dvSd = .dvSd, Nobs = length(.allDv))
}
