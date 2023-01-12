#' Method for getting simulation rxode2 classic models based on fits
#'
#' @param x list where first element is the fit.  The class represents the estimation method.
#' @return  model for fit$simulationModel
#' @author Matthew L. Fidler
#' @export
#' @keywords internal
getBaseSimModelFit <- function(x) {
  UseMethod("getBaseSimModelFit")
}

#attr(rxUiGet.foceiModelDigest, "desc") <- "Get the md5 digest for the focei model"
#' @export
rxUiGet.foceiSimModelCache <- function(x, ...) {
  file.path(rxode2::rxTempDir(),
            paste0("focei-sim-", rxUiGet.foceiModelDigest(x, ...), ".qs"))
}


#' @rdname getBaseSimModelFit
#' @export
getBaseSimModelFit.focei <- function(x) {
  obj <- x[[1]]
  .cacheFile <- obj$ui$foceiSimModelCache
  if (file.exists(.cacheFile)) {
    .ret <- qs::qread(.cacheFile)
    return(.ret)
  }
  .base <- rxode2::rxCombineErrorLines(obj$ui)
  .base[[1]] <- quote(`rxModelVars`)
  .b2 <- .base[[2]]
  .params <- .b2[[2]]
  .b2 <- .b2[-2]
  .end <- NULL
  .cur <- .b2[[length(.b2)]][[1]]
  while (identical(.cur, quote(`dvid`)) ||
           identical(.cur, quote(`cmt`))) {
             .end <- c(.b2[[length(.b2)]], .end)
             .b2 <- .b2[-length(.b2)]
             .cur <- .b2[[length(.b2)]][[1]]
           }
  .base[[2]] <- .b2
  .base <- eval(.base)
  # focei models are always pruned
  .s <- .loadSymengine(rxode2::rxPrune(.base))
  .msuccess("done")
  .prd <- get("rx_pred_", envir = .s)
  .prd <- paste0("rx_pred_ ~", rxode2::rxFromSE(.prd))
  .r <- get("rx_r_", envir = .s)
  .r <- paste0("rx_r_~", rxode2::rxFromSE(.r))
  .yj <- paste(get("rx_yj_", envir = .s))
  .yj <- paste0("rx_yj_~", rxode2::rxFromSE(.yj))
  .lambda <- paste(get("rx_lambda_", envir = .s))
  .lambda <- paste0("rx_lambda_~", rxode2::rxFromSE(.lambda))
  .hi <- paste(get("rx_hi_", envir = .s))
  .hi <- paste0("rx_hi_~", rxode2::rxFromSE(.hi))
  .low <- paste(get("rx_low_", envir = .s))
  .low <- paste0("rx_low_~", rxode2::rxFromSE(.low))
  .lhs0 <- .s$..lhs0
  .ipredSim <- get("ipredSim", envir = .s)
  .ipredSim <- paste0("ipredSim =", rxode2::rxFromSE(.ipredSim))
  .sim <- get("sim", envir = .s)
  .sim <- paste0("sim =", rxode2::rxFromSE(.sim))
  if (is.null(.lhs0)) .lhs0 <- ""
  .lhs <- .s$..lhs
  if (is.null(.lhs)) .lhs <- ""
  .ddt <- .s$..ddt
  if (is.null(.ddt)) .ddt <- ""
  .base <- paste(c(
    .s$..stateInfo["state"],
    .lhs0,
    .ddt,
    .yj,
    .lambda,
    .hi,
    .low,
    .prd,
    .r,
    .ipredSim,
    .sim,
    .lhs,
    .s$..stateInfo["statef"],
    .s$..stateInfo["dvid"]
  ), collapse = "\n")

  .ctl <- obj$foceiControl
  if (.ctl$sumProd) {
    .malert("stabilizing round off errors in simulation problem...")
    .base <- rxode2::rxSumProdModel(.base)
    .msuccess("done")
  }
  if (.ctl$optExpression) {
    .base <- rxode2::rxOptExpr(.base, ifelse(.getRxPredLlikOption(),
                                             "simulation llik model",
                                             "simulation model"))
    .msuccess("done")
  }
  .ret <- eval(parse(text=paste0("quote(rxode2({\n", .base, "\n}))")))
  qs::qsave(.ret, .cacheFile)
  .ret
}

#' @rdname getBaseSimModelFit
#' @export
getBaseSimModelFit.foce <- getBaseSimModelFit.focei

#' @rdname getBaseSimModelFit
#' @export
getBaseSimModelFit.fo <- getBaseSimModelFit.focei
#' @rdname getBaseSimModelFit
#' @export
getBaseSimModelFit.foi <- getBaseSimModelFit.focei

#' @rdname getBaseSimModelFit
#' @export
getBaseSimModelFit.posthoc <- getBaseSimModelFit.focei

#' @rdname getBaseSimModelFit
#' @export
getBaseSimModelFit.default <- function(x) {
  .obj <- x[[1]]
  .ui <- .obj$ui
  getBaseSimModel(.ui)
}

getBaseSimModel.nlmixr2FitCoreSilent <- function(obj) {
  .est <- obj$est
  .ret <- list(obj)
  class(.ret) <- c(.est, "getBaseSimModelFit")
  return(getBaseSimModelFit(.ret))
}

getBaseSimModel.nlmixr2FitData <- getBaseSimModel.nlmixr2FitCoreSilent
getBaseSimModel.nlmixr2FitCore <- getBaseSimModel.nlmixr2FitCoreSilent
