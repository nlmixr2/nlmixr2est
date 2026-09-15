.setOfvFo <- function(fit, type = c("focei", "foce", "fo")) {
  .type <- match.arg(type)
  nlmixrWithTiming(paste0(.type, "Lik"), {
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
      nlmixr2(fit, nlme::getData(fit), "focei", control = .foceiControl))
    .env <- fit$env
    .addFoceiInfoToFit(.env, .newFit)
    .ob1 <- .newFit$objDf
    .etaObf <- .newFit$etaObf
    nlmixrAddObjectiveFunctionDataFrame(fit, .ob1, .rn, .etaObf)
    invisible(fit)
  }, envir=fit)
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
  nlmixrWithTiming(paste0(.type, "Lik"), {
    .sigdig <- fit$foceiControl$sigdig
    .ctl <- impmapControl(print = 0L, nIter = 0L, covMethod = "",
                          calcTables = FALSE, compress = FALSE,
                          sigdig = if (is.null(.sigdig)) 3 else .sigdig)
    .ctl$etaMat <- fit$etaMat
    .newFit <- .nlmixr2PriorGateBypass(
      nlmixr2(fit, nlme::getData(fit), .type, control = .ctl))
    .newEnv <- .newFit$env
    .env <- fit$env
    .df <- attr(get("logLik", .env), "df")
    .nobs <- .env$nobs
    # impObj omits the normal constant, like an adjusted objective
    .objf <- .newEnv$impObj
    .m2ll <- .objf + .nobs * log(2 * pi)
    .adj <- .env$adjObf  # saem stores it on the env, the focei family on the control
    if (is.null(.adj)) .adj <- fit$foceiControl$adjObf
    if (isFALSE(.adj)) .objf <- .m2ll
    .tmp <- data.frame(
      OBJF = .objf, AIC = .m2ll + 2 * .df, BIC = .m2ll + log(.nobs) * .df,
      "Log-likelihood" = -.m2ll / 2, check.names = FALSE
    )
    nlmixrAddObjectiveFunctionDataFrame(fit, .tmp, .rn)
    invisible(fit)
  }, envir=fit)
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
    stop("cannot switch objective function to '", type, "' type",
         call. = FALSE)
  }
  .setOfvSaemQuad(x, type)
}

#' Add a saem laplace/gauss quadrature objective function
#' @param x saem fit
#' @param type "laplace<nsd>" or "gauss<nnodes>_<nsd>"
#' @return fit, invisibly
#' @noRd
.setOfvSaemQuad <- function(x, type) {
  nlmixrWithTiming(paste0(type, "Lik"), {
    .env <- x$env
    .reg <- rex::rex(start, "laplace", capture(.regNum), end)
    .regG <- rex::rex(start, "gauss", capture(.regNum), "_", capture(.regNum), end)
    if (regexpr(.reg, type, perl = TRUE) != -1) {
      .nnode <- 1
      .nsd <- as.numeric(sub(.reg, "\\1", type, perl = TRUE))
    } else if (regexpr(.regG, type, perl = TRUE) != -1) {
      .nnode <- as.numeric(sub(.regG, "\\1", type, perl = TRUE))
      .nsd <- as.numeric(sub(.regG, "\\2", type, perl = TRUE))
    } else {
      stop("cannot switch objective function to '", type, "' type", call. = FALSE)
    }
    .saemObf <- calc.2LL(x$saem, nnodes.gq = .nnode, nsd.gq = .nsd, x$phiM)
    .llik <- -.saemObf / 2
    .nobs <- .env$nobs
    attr(.llik, "df") <- attr(get("logLik", .env), "df")
    .objf <- ifelse(.env$adjObf, .saemObf - .nobs * log(2 * pi), .saemObf)
    .tmp <- data.frame(
      OBJF = .objf, AIC = .saemObf + 2 * attr(get("logLik", .env), "df"),
      BIC = .saemObf + log(.env$nobs) * attr(get("logLik", .env), "df"),
      "Log-likelihood" = as.numeric(.llik), check.names = FALSE
    )
    nlmixrAddObjectiveFunctionDataFrame(x, .tmp, type)
    return(invisible(x))
  }, envir=x)
}

##' @rdname setOfv
##' @export
getOfvType <- function(x) {
  return(x$ofvType)
}
