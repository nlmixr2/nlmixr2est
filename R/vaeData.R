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
  ## A mixture random effect (e.g. eta.ke in mix(exp(lke1+eta.ke),p,exp(lke2+
  ## eta.ke))) has no single structural theta -- its prior is the component-
  ## independent N(0,omega), so it centers at 0 and the fixed component thetas in
  ## the full theta vector carry the structure (the inner problem selects them).
  .zPopThetaIdx <- match(.map$thetaForEta, .thRows$name)
  .isMix <- is.na(.zPopThetaIdx)
  .nMix <- tryCatch(as.integer(ui$saemNMix), error = function(e) 1L)
  if (any(.isMix) && (is.na(.nMix) || .nMix < 2L))
    stop("est=\"vae\" needs every random effect mu-referenced to a theta", call. = FALSE)
  .zPop <- numeric(.neta)                                      # structural population means (transformed)
  .zPop[!.isMix] <- as.numeric(.th[.zPopThetaIdx[!.isMix]])

  ## omega init (diagonal) for the etas, and residual a init (error theta est)
  .omega <- vapply(.etaNames, function(nm) {
    .r <- .idf[!is.na(.idf$neta1) & .idf$neta1 == .idf$neta2 & .idf$name == nm, , drop = FALSE]
    as.numeric(.r$est[1])
  }, numeric(1))
  .errRow <- .idf[!is.na(.idf$err) & !is.na(.idf$ntheta), , drop = FALSE]
  .a <- if (nrow(.errRow) > 0) as.numeric(.errRow$est[1]) else 1.0
  .aThetaIdx <- if (nrow(.errRow) > 0) as.integer(.errRow$ntheta[1]) else NA_integer_

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
       th = .th, zPopThetaIdx = .zPopThetaIdx, aThetaIdx = .aThetaIdx, isMix = .isMix,
       zPop = .zPop, omega = .omega, a = .a,
       subj = subj, dataIn = dataIn, lengths = lengths, covIn = covIn,
       covNames = .covNames, covMat = .covMat, covType = .covType, covPop = .covPop,
       tMax = .tMax, dvMean = .dvMean, dvSd = .dvSd, Nobs = length(.allDv))
}
