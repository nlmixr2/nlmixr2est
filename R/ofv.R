.setOfvFo <- function(fit, type = c("focei", "foce", "fo")) {
  .type <- match.arg(type)
  nlmixrWithTiming(
    paste0(.type, "Lik"),
    {
      .foceiControl <- fit$foceiControl
      .foceiControl$etaMat <- fit$etaMat
      .foceiControl$fo <- FALSE
      .foceiControl$maxOuterIterations <- 0L
      .foceiControl$maxInnerIterations <- 0L
      .foceiControl$calcTables <- FALSE
      .foceiControl$covMethod <- 0L
      .foceiControl$compress <- FALSE
      .foceiControl$nAGQ <- 0L # focei/foce/fo objectives, not the fit's quadrature
      if (.type == "focei") {
        .foceiControl$interaction <- TRUE
        .rn <- "FOCEi"
      } else if (.type == "foce") {
        .foceiControl$interaction <- FALSE
        .rn <- "FOCE"
      } else {
        .foceiControl$interaction <- FALSE
        .foceiControl$fo <- TRUE
        .foceiControl$etaMat <- NULL
        .rn <- "FO"
      }
      .inObjDf <- fit$objDf
      if (any(rownames(.inObjDf) == .rn)) {
        return(fit)
      }
      .foceiControl <- do.call(foceiControl, .foceiControl)
      .newFit <- .nlmixr2PriorGateBypass(
        nlmixr2(fit, nlme::getData(fit), "focei", control = .foceiControl)
      )
      .env <- fit$env
      .addFoceiInfoToFit(.env, .newFit)
      .ob1 <- .newFit$objDf
      .etaObf <- .newFit$etaObf
      nlmixrAddObjectiveFunctionDataFrame(fit, .ob1, .rn, .etaObf)
      invisible(fit)
    },
    envir = fit
  )
}

#' Add an importance-sampling objective by an E-step-only run at the fit's estimates
#'
#' @param fit nlmixr2 fit
#' @param type "imp" or "impmap"
#' @return fit, invisibly
#' @noRd
.setOfvImp <- function(fit, type = c("imp", "impmap")) {
  .type <- match.arg(type)
  .rn <- toupper(.type)
  if (any(rownames(fit$objDf) == .rn)) {
    return(invisible(fit))
  }
  nlmixrWithTiming(
    paste0(.type, "Lik"),
    {
      .sigdig <- fit$foceiControl$sigdig
      .ctl <- impmapControl(
        print = 0L,
        nIter = 0L,
        covMethod = "",
        calcTables = FALSE,
        compress = FALSE,
        sigdig = if (is.null(.sigdig)) 3 else .sigdig
      )
      .ctl$etaMat <- fit$etaMat
      .newFit <- .nlmixr2PriorGateBypass(
        nlmixr2(fit, nlme::getData(fit), .type, control = .ctl)
      )
      # impObj omits the normal constant, like an adjusted objective
      .objf <- .newFit$env$impObj
      .setOfvRow(fit, .rn, .objf + fit$env$nobs * log(2 * pi), .objf)
    },
    envir = fit
  )
}

#' Add an objective function row computed from a -2 log-likelihood
#'
#' AIC and BIC use the fit's `df` and `nobs`.  OBJF omits the normal constant
#' unless `adjObf` (on the fit environment, else its control) is `FALSE`.
#' @param fit nlmixr2 fit
#' @param type objective function type (the row name)
#' @param m2ll -2 log-likelihood, with the normal constant
#' @param adjObjf the objective without the normal constant
#' @return fit, invisibly
#' @noRd
.setOfvRow <- function(fit, type, m2ll, adjObjf = m2ll - fit$env$nobs * log(2 * pi)) {
  .env <- fit$env
  .df <- attr(get("logLik", .env), "df")
  .adj <- .env$adjObf
  if (is.null(.adj)) {
    .adj <- fit$control$adjObf
  }
  .tmp <- data.frame(
    OBJF = if (isFALSE(.adj)) m2ll else adjObjf,
    AIC = m2ll + 2 * .df,
    BIC = m2ll + log(.env$nobs) * .df,
    "Log-likelihood" = -m2ll / 2,
    check.names = FALSE
  )
  nlmixrAddObjectiveFunctionDataFrame(fit, .tmp, type)
  invisible(fit)
}

##' Set/get Objective function type for a nlmixr2 object
##'
##' @param x nlmixr2 fit object
##' @param type Type of objective function to use for AIC, BIC, and
##'     $objective.  `"imp"` and `"impmap"` add an importance-sampling
##'     objective from an E-step-only run (`nIter=0`) at the fit's estimates.
##' @return Nothing
##' @author Matthew L. Fidler
##' @export
setOfv <- function(x, type) {
  assertNlmixrFit(x)
  .objDf <- x$objDf
  .w <- which(tolower(row.names(.objDf)) == tolower(type))
  if (length(.w) != 1) {
    return(.setOfvAdd(x, type))
  }
  .env <- x$env
  .objf <- .objDf[.w, "OBJF"]
  .lik <- .objDf[.w, "Log-likelihood"]
  attr(.lik, "df") <- attr(get("logLik", .env), "df")
  attr(.lik, "nobs") <- attr(get("logLik", .env), "nobs")
  class(.lik) <- "logLik"
  .bic <- .objDf[.w, "BIC"]
  .aic <- .objDf[.w, "AIC"]
  assign("OBJF", .objf, .env)
  assign("objf", .objf, .env)
  assign("objective", .objf, .env)
  assign("logLik", .lik, .env)
  assign("AIC", .aic, .env)
  assign("BIC", .bic, .env)
  if (!is.null(x$saem)) {
    .setSaemExtra(.env, type)
  }
  .env$ofvType <- type
  invisible(x)
}

#' Compute and add an objective function type the fit does not have yet
#' @param x nlmixr2 fit
#' @param type objective function type
#' @return fit, invisibly
#' @noRd
.setOfvAdd <- function(x, type) {
  .type <- tolower(type)
  if (any(.type == c("focei", "foce", "fo"))) {
    return(.setOfvFo(x, .type))
  }
  if (any(.type == c("imp", "impmap"))) {
    return(.setOfvImp(x, .type))
  }
  if (is.null(x$saem)) {
    stop("cannot switch objective function to '", type, "' type", call. = FALSE)
  }
  .setOfvSaemQuad(x, type)
}

#' Add a saem laplace/gauss quadrature objective function
#' @param x saem fit
#' @param type "laplace<nsd>" or "gauss<nnodes>_<nsd>"
#' @return fit, invisibly
#' @noRd
.setOfvSaemQuad <- function(x, type) {
  nlmixrWithTiming(
    paste0(type, "Lik"),
    {
      .q <- .saemParseLikName(type)
      if (is.null(.q)) {
        stop("cannot switch objective function to '", type, "' type", call. = FALSE)
      }
      .setOfvRow(x, type, calc.2LL(x$saem, nnodes.gq = .q[1], nsd.gq = .q[2], x$phiM))
    },
    envir = x
  )
}

##' @rdname setOfv
##' @export
getOfvType <- function(x) {
  return(x$ofvType)
}
