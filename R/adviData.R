# adviData.R -- ADVI data preparation.  Builds, from the rxode2 ui + data:
#  - the eta<->theta (population) map and initial theta / omega from ini()
#  - per-subject observation/event inputs (times, DV, censoring)
#  - the latent-parameter classification that drives the gradient path:
#      * mu-referenced structural thetas    -> gradient via d/dtheta = d/deta
#      * non-mu structural thetas           -> theta-sensitivity ODE (impThetaScore)
#      * residual-error (sigma) thetas      -> algebraic d(V)/d(theta) (impThetaScore)
# ADVI dispatches with attr `unbounded=TRUE`, so bounded population parameters are
# already on the unconstrained real scale by the time this runs (the variational
# family lives in that real coordinate space, per Kucukelbir 2017 Sec 2.3).

#' Classify the population parameters for the ADVI gradient path.
#' @param ui rxode2 ui object
#' @return list(muRefThetaIdx (per eta, 1-based ntheta or NA), struct, sigma, all)
#'   where `struct`/`sigma`/`all` are the non-mu theta-sensitivity indices from
#'   `.impmapEstTheta` and `muRefThetaIdx[k]` is the theta mu-referenced by eta k.
#' @noRd
.adviClassifyPars <- function(ui) {
  .map <- .foceiEtaThetaMap(ui)
  .thRows <- ui$iniDf[!is.na(ui$iniDf$ntheta), , drop = FALSE]
  .thRows <- .thRows[order(.thRows$ntheta), , drop = FALSE]
  .muRefThetaIdx <- match(.map$thetaForEta, .thRows$name) # 1-based ntheta, NA if free eta
  .sens <- .impmapEstTheta(ui)                            # $struct, $sigma, $all (1-based ntheta)
  list(etaNames = .map$etaNames, muRefThetaIdx = .muRefThetaIdx,
       struct = .sens$struct, sigma = .sens$sigma, thetaSensIdx = .sens$all)
}

#' Prepare ADVI inputs from a ui + data.
#' @param ui rxode2 ui object
#' @param data estimation data (ID/TIME/DV/EVID/... columns)
#' @return list of prepared ADVI inputs (see field comments below)
#' @noRd
.adviDataPrep <- function(ui, data) {
  .map <- .foceiEtaThetaMap(ui)
  .etaNames <- .map$etaNames
  .neta <- length(.etaNames)
  if (.neta == 0L) stop("est=\"advi\" requires at least one random effect", call. = FALSE)
  .idf <- ui$iniDf

  ## full theta vector (THETA_i_ in ntheta order), from ini estimates
  .thRows <- .idf[!is.na(.idf$ntheta), , drop = FALSE]
  .thRows <- .thRows[order(.thRows$ntheta), , drop = FALSE]
  .th <- setNames(as.numeric(.thRows$est), paste0("THETA_", seq_len(nrow(.thRows)), "_"))

  ## population parameter classification for the gradient path
  .cls <- .adviClassifyPars(ui)

  ## omega init (diagonal) + which variances are FIXED
  .omega <- vapply(.etaNames, function(nm) {
    .r <- .idf[!is.na(.idf$neta1) & .idf$neta1 == .idf$neta2 & .idf$name == nm, , drop = FALSE]
    as.numeric(.r$est[1])
  }, numeric(1))
  .omegaFix <- vapply(.etaNames, function(nm) {
    .r <- .idf[!is.na(.idf$neta1) & .idf$neta1 == .idf$neta2 & .idf$name == nm, , drop = FALSE]
    isTRUE(as.logical(.r$fix[1]))
  }, logical(1))

  ## which thetas are FIXED (held at ini, never gradient-updated)
  .thetaFix <- as.logical(.thRows$fix); .thetaFix[is.na(.thetaFix)] <- FALSE

  ## residual error params (in theta order): value, theta index, type
  .errRow <- .idf[!is.na(.idf$err) & !is.na(.idf$ntheta), , drop = FALSE]
  .errRow <- .errRow[order(.errRow$ntheta), , drop = FALSE]
  .a <- if (nrow(.errRow) > 0) setNames(as.numeric(.errRow$est), .errRow$name) else numeric(0)
  .errThetaIdx <- as.integer(.errRow$ntheta)
  .errType <- as.character(.errRow$err)

  ## normalize data columns + per-subject observation/event inputs
  d <- as.data.frame(data)
  names(d) <- toupper(names(d))
  if (is.null(d$EVID)) d$EVID <- ifelse(is.na(d$AMT) | d$AMT == 0, 0L, 1L)
  .ids <- unique(d$ID)
  N <- length(.ids)
  subj <- vector("list", N); .nobs <- 0L
  for (i in seq_len(N)) {
    .di <- d[d$ID == .ids[i], , drop = FALSE]
    .obs <- .di[.di$EVID == 0, , drop = FALSE]
    .cens <- if (is.null(.obs$CENS)) integer(nrow(.obs)) else as.integer(.obs$CENS)
    .limit <- if (is.null(.obs$LIMIT)) rep(NA_real_, nrow(.obs)) else as.numeric(.obs$LIMIT)
    subj[[i]] <- list(ev = .di, times = .obs$TIME, y = .obs$DV, n = nrow(.obs),
                      cens = .cens, limit = .limit)
    .nobs <- .nobs + nrow(.obs)
  }

  list(N = N, neta = .neta, zDim = .neta, etaNames = .etaNames, ids = .ids,
       th = .th, theta = as.numeric(.th), ntheta = length(.th), thetaFix = .thetaFix,
       thetaRealNames = .thRows$name,   # actual parameter names in ntheta order
       muRefThetaIdx = .cls$muRefThetaIdx, isFree = is.na(.cls$muRefThetaIdx),
       structIdx = .cls$struct, sigmaIdx = .cls$sigma, thetaSensIdx = .cls$thetaSensIdx,
       omega = .omega, omegaFix = .omegaFix,
       errThetaIdx = .errThetaIdx, errType = .errType, a = .a,
       subj = subj, Nobs = .nobs)
}
