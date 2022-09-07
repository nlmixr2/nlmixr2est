#' Expands the simulation model to add tad data item
#'
#' 
#' @param obj Object to expand
#' @return quoted model
#' @author Matthew L. Fidler
#' @noRd
.expandSimModelAddTad <- function(obj) {
  .ret <- obj
  .tmp <- .ret[[2]]
  .idx <- NULL
  .idxDvid <- NULL
  .tmp <- lapply(seq(2, length(.tmp)), function(i) {
    if (identical(.tmp[[i]][[1]], quote(`cmt`))) {
      .idx <<- i - 1
    }
    if (identical(.tmp[[i]][[1]], quote(`dvid`))) {
      .idxDvid <<- i - 1
    }
    .tmp[[i]]
  })
  if (is.null(.idx) && !is.null(.idxDvid)) {
    # use dvid() instead of cmt()
    .idx <- .idxDvid
  } else if (is.null(.idx) && is.null(.idxDvid)) {
    # simply append to the end.
    .ret[[2]] <- as.call(c(list(quote(`{`)),
                           .tmp,
                           list(str2lang("tad <- tad()"))))
    return(.ret)

  }
  .ret[[2]] <- as.call(lapply(seq(1, length(.tmp)+2), function(i) {
    if (i == 1) {
      quote(`{`)
    } else if (i-1 == .idx) {
      str2lang("tad <- tad()")
    } else if (i-1 < .idx) {
      .tmp[[i-1]]
    } else {
      .tmp[[i-2]]
    }
  }))
  .ret
}
#' Get the simulation model for VPC and NPDE
#'
#' 
#' @param obj nlmixr fit object
#' @param hideIpred Hide the ipred (by default FALSE)
#' @param tad Include `tad` calculation (by default FALSE)
#' @return quoted simulation model (simply need to evaluate it)
#' @author Matthew L. Fidler
#' @noRd
.getSimModel <- function(obj, hideIpred=FALSE, tad=TRUE) {
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
  .ret <- .f(.lines)
  if (tad) {
    .ret <- .expandSimModelAddTad(.ret)
  }
  .ret
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
