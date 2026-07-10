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
  ## structural theta index (in the full theta vector) paired with each eta
  .zPopThetaIdx <- match(.map$thetaForEta, .thRows$name)
  if (anyNA(.zPopThetaIdx)) stop("est=\"vae\" needs every random effect mu-referenced to a theta", call. = FALSE)
  .zPop <- as.numeric(.th[.zPopThetaIdx])                      # structural population means (transformed)

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
    subj[[i]] <- list(ev = .di, times = .times, y = .y, n = length(.times))
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
  covIn <- matrix(0, N, 0L)                     # Milestone A: no encoder covariates

  list(N = N, neta = .neta, zDim = .neta, etaNames = .etaNames,
       th = .th, zPopThetaIdx = .zPopThetaIdx, aThetaIdx = .aThetaIdx,
       zPop = .zPop, omega = .omega, a = .a,
       subj = subj, dataIn = dataIn, lengths = lengths, covIn = covIn,
       tMax = .tMax, dvMean = .dvMean, dvSd = .dvSd, Nobs = length(.allDv))
}
