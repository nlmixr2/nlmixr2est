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

.isAssignExpr <- function(x) {
  length(x) == 3L &&
    (identical(x[[1]],quote(`=`)) ||
       identical(x[[1]],quote(`<-`)))
}

.replaceThetaEtaWithNamed <- function(x, iniDf) {
  if (is.call(x)) {
    if (length(x) == 3L &&
          is.numeric(x[[3]]) &&
          identical(x[[1]], quote(`[`))) {
      if (identical(x[[2]], quote(`THETA`))) {
        return(str2lang(iniDf[which(iniDf$ntheta== x[[3]]), "name"]))
      }
      if (identical(x[[2]], quote(`ETA`))) {
        return(str2lang(iniDf[which(iniDf$neta1 == x[[3]] & iniDf$neta2 == x[[3]]),
                              "name"]))
      }
    }
    return(as.call(lapply(x, .replaceThetaEtaWithNamed, iniDf)))
  }
  x
}

#' @rdname getBaseSimModelFit
#' @export
getBaseSimModelFit.focei <- function(x) {
  obj <- x[[1]]
    getBaseSimModelFit.default(x) # fall back to basic with new method
  if (all(obj$ui$predDf$distribution == "norm")) {
    .expr <- eval(parse(text=paste0("quote(rxode2({",
                                    rxode2::rxNorm(obj$foceiModel$predOnly),
                                    "}))")))
    .e2 <- .expr[[2]]
    .e2 <- lapply(seq_along(.e2), function(i) {
      .cur <- .e2[[i]]
      .replaceThetaEtaWithNamed(.cur, obj$ui$iniDf)
    })
    .w <- vapply(seq_along(.e2), function(i) {
      .cur <- .e2[[i]]
      if (.isAssignExpr(.cur) &&
            identical(.cur[[2]], .cur[[3]])) {
        return(FALSE)
      }
      return(TRUE)
    }, logical(1), USE.NAMES=FALSE)
    .e2 <- .e2[.w]
    .w <- which(vapply(seq_along(.e2), function(i) {
      .cur <- .e2[[i]]
      if (.isAssignExpr(.cur) &&
            identical(.cur[[2]], quote(`rx_r_`))) return(TRUE)
      FALSE
    }, logical(1), USE.NAMES=TRUE))
    .w <- seq(1, .w)
    .e21 <- .e2[.w]
    .e22 <- .e2[-.w]
    .e2 <- as.call(c(.e21,
      list(quote(ipredSim <- rxTBSi(rx_pred_, rx_lambda_, rx_yj_, rx_low_, rx_hi_)),
           eval(parse(text=paste0("quote(sim <- rxTBSi(rx_pred_ + sqrt(rx_r_) *(",
                                  paste(paste0("error.", obj$predDf$var, "*(CMT==", obj$predDf$cmt, ")"),
                                        collapse = "+"),
                                  "), rx_lambda_, rx_yj_, rx_low_, rx_hi_))")))),
      .e22))
    .expr[[2]] <- .e2
  }
  getBaseSimModelFit.default(x)
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
  rxode2::getBaseSimModel(.ui)
}

getBaseSimModel.nlmixr2FitCoreSilent <- function(obj) {
  .est <- obj$est
  .ret <- list(obj)
  class(.ret) <- c(.est, "getBaseSimModelFit")
  return(getBaseSimModelFit(.ret))
}

getBaseSimModel.nlmixr2FitData <- getBaseSimModel.nlmixr2FitCoreSilent
getBaseSimModel.nlmixr2FitCore <- getBaseSimModel.nlmixr2FitCoreSilent
