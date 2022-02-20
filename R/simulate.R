.getSimModel <- function(obj, hideIpred=FALSE) {
  .lines <- rxode2::rxCombineErrorLines(obj$ui)
  .f <- function(x) {
    if (is.atomic(x) || is.name(x) || is.pairlist(x)) {
      return(x)
    } else if (is.call(x)) {
      if (identical(x[[1]], quote(`<-`)) ||
            identical(x[[1]], quote(`=`))) {
        if (identical(x[[2]], quote(`ipredSim`))) {
          x[[2]] <- quote(`ipred`)
          if (hideIpred) {
            x[[1]] <- quote(`~`)
          } else {
            x[[1]] <- quote(`<-`)
          }
        } else if (identical(x[[2]], quote(`sim`))) {
          x[[2]] <- quote(`sim`)
          x[[1]] <- quote(`<-`)
        } else if (length(x[[2]]) == 1L) {
          x[[1]] <- quote(`~`)
        } else {
          if (identical(x[[2]][[1]], quote(`/`))) {
            x[[1]] <- quote(`~`)
          } else {
            x[[1]] <- quote(`<-`)
          }
        }
      }
      return(as.call(lapply(x, .f)))
    }
  }
  .f(.lines)
}

.simInfo <- function(object) {
  .mod <- .getSimModel(object, hideIpred=FALSE)
  .omega <- object$omega
  .etaN <- dimnames(.omega)[[1]]
  .params <- nlme::fixed.effects(object)
  .params <- .params
  .dfObs <- object$nobs
  .nlmixr2Data <- object$origData
  .dfSub <- object$nsub
  .env <- object$env
  if (exists("cov", .env)) {
    .thetaMat <- nlme::getVarCov(object)
  } else {
    .thetaMat <- NULL
  }
  if (all(is.na(object$ui$ini$neta1))) {
    .omega <- NULL
    .dfSub <- 0
  }
  .sigma <- object$ui$simulationSigma
  return(list(
    rx = .mod, params = .params, events = .nlmixr2Data, thetaMat = .thetaMat,
    omega = .omega, sigma = .sigma, dfObs = .dfObs, dfSub = .dfSub
  ))
}

#' @importFrom rxode2 rxSolve
#' @export
rxode2::rxSolve
