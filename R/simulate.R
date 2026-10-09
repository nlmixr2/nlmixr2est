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
  .tmp <- lapply(seq(2, length(.tmp)), function(i) .tmp[[i]])
  for (i in seq_along(.tmp)) {
    if (identical(.tmp[[i]][[1]], quote(`cmt`))) {
      .idx <- i
    }
    if (identical(.tmp[[i]][[1]], quote(`dvid`))) {
      .idxDvid <- i
    }
  }
  if (is.null(.idx) && !is.null(.idxDvid)) {
    # use dvid() instead of cmt()
    .idx <- .idxDvid
  } else if (is.null(.idx) && is.null(.idxDvid)) {
    # simply append to the end.
    .ret[[2]] <- as.call(c(list(quote(`{`)), .tmp, list(str2lang("tad <- tad()"))))
    return(.ret)
  }
  .ret[[2]] <- as.call(lapply(seq_len(length(.tmp) + 2), function(i) {
    if (i == 1) {
      quote(`{`)
    } else if (i - 1 == .idx) {
      str2lang("tad <- tad()")
    } else if (i - 1 < .idx) {
      .tmp[[i - 1]]
    } else {
      .tmp[[i - 2]]
    }
  }))
  .ret
}
#' Variables used inside a history function of a model
#'
#' @param x quoted model
#' @return character vector of the variables referenced inside `lag()`,
#'   `lead()`, `first()`, `last()` or `diff()`
#' @author Matthew L. Fidler
#' @noRd
.simModelLaggedVars <- function(x) {
  if (!is.call(x)) {
    return(character(0))
  }
  .ret <- unlist(lapply(as.list(x)[-1], .simModelLaggedVars))
  if (
    is.name(x[[1]]) &&
      as.character(x[[1]]) %in% c("lag", "lead", "first", "last", "diff") &&
      length(x) >= 2L &&
      is.name(x[[2]])
  ) {
    .ret <- c(as.character(x[[2]]), .ret)
  }
  unique(as.character(.ret))
}
#' Get the simulation model for VPC and NPDE
#'
#'
#' @param obj nlmixr fit object
#' @param hideIpred Hide the ipred (by default FALSE)
#' @param tad Include `tad` calculation (by default FALSE)
#' @return quoted simulation model (simply need to evaluate it); its
#'   `"lagged"` attribute names the extra outputs kept for `lag()`
#' @author Matthew L. Fidler
#' @noRd
.getSimModel <- function(obj, hideIpred = FALSE, tad = TRUE) {
  .lines <- rxode2::getBaseSimModel(obj)
  # rxode2 only allows a history function of a real lhs, so these stay `<-`
  .lagged <- .simModelLaggedVars(.lines)
  .acc <- new.env(parent = emptyenv())
  .acc$keptLhs <- character(0)
  .f <- function(x) {
    if (is.atomic(x) || is.name(x) || is.pairlist(x)) {
      return(x)
    } else if (is.call(x)) {
      if (
        identical(x[[1]], quote(`<-`)) ||
          identical(x[[1]], quote(`=`))
      ) {
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
          if (as.character(x[[2]]) %in% .lagged) {
            x[[1]] <- quote(`<-`)
            .acc$keptLhs <- c(.acc$keptLhs, as.character(x[[2]]))
          } else {
            x[[1]] <- quote(`~`)
          }
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
  # outputs only so `lag()` compiles; callers drop them from the solve
  attr(.ret, "lagged") <- unique(.acc$keptLhs)
  .ret
}

.simInfo <- function(object) {
  .env <- new.env(parent = emptyenv())
  .env$ui <- object$ui
  .env$data <- object$origData
  suppressMessages(.preProcessHooksRun(.env, "rxSolve"))
  .mod <- .getSimModel(.env$ui, hideIpred = FALSE)
  .omega <- .env$ui$omega
  .etaN <- dimnames(.omega)[[1]]
  .params <- nlme::fixed.effects(.env$ui)
  .params <- .params
  .dfObs <- object$nobs
  .nlmixr2Data <- .env$data
  .dfSub <- object$nsub
  .env <- object$env
  if (exists("cov", .env)) {
    .thetaMat <- nlme::getVarCov(object)
  } else {
    .thetaMat <- NULL
  }
  if (all(is.na(.env$ui$ini$neta1))) {
    .omega <- NULL
    .dfSub <- 0
  }
  .sigma <- .env$ui$simulationSigma
  return(list(
    rx = .mod,
    params = .params,
    events = .nlmixr2Data,
    thetaMat = .thetaMat,
    omega = .omega,
    sigma = .sigma,
    dfObs = .dfObs,
    dfSub = .dfSub
  ))
}

#' @importFrom rxode2 rxSolve
#' @export
rxode2::rxSolve
