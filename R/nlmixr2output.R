.getBackTransformationFunction <- function(par, ui) {
  # This has a specified back-transformation
  .w <- which(ui$iniDf$name == par)
  if (length(.w) == 1L) {
    .b <- ui$iniDf$backTransform
    if (!is.na(.b)) {
      return(.b)
    }
  }
  # If this has extra information in the mu-ref
  # (ie exp(tv + eta.v + 2)), then this shouldn't
  # have a default back-transformation
  .w <- which(ui$muRefExtra$parameter == par)
  if (length(.w) == 1L) {
    return("")
  }
  # Covariates should be reported without back-transformation
  .w <- which(ui$muRefCovariateDataFrame$covariateParameter == par)
  if (length(.w) == 1L) {
    return("")
  }
}
#' Get the parameter label and apply to parameter dataset
#'
#' @param ret Final fit environment
#' @param ui UI fit information
#' @return Nothing, called for side effects
#' @author Matthew L. Fidler
#' @noRd
.updateParFixedAddParameterLabel <- function(ret, ui) {
  .lab <- paste(ui$iniDf$label[!is.na(ui$iniDf$ntheta)])
  .lab[.lab == "NA"] <- ""
  .lab <- gsub(" *$", "", gsub("^ *", "", .lab))
  if (!all(.lab == "")) {
    ret$popDf <- data.frame(Parameter = .lab, ret$popDf, check.names = FALSE)
    ret$popDfSig <- data.frame(Parameter = .lab, ret$popDfSig, check.names = FALSE)
  }
}
#'  This applies the manually specified back-transformations
#'
#' @param .ret focei environment
#' @return Nothing, called for side effecs
#' @author Matthew L. Fidler
#' @examples
.updateParFixedApplyManualBacktransformations <- function(.ret, .ui) {
  .qn <- qnorm(1.0-(1-0.95)/2)
  .btName <- names(.ret$popDfSig)[4]
  .sigdig <- rxode2::rxGetControl(.ui, "sigdig", 3L)
  .fmt <- paste0("%", .sigdig, "g (%", .sigdig, "g, %", .sigdig, "g)")
  .fmt2 <- paste0("%", .sigdig, "g")
  lapply(seq_along(.ret$popDf$Estimate), function(i) {
    theta <- row.names(.ret$popDf)[i]
    .w <- which(.ui$iniDf$name == theta)
    if (length(.w) == 1L) {
      .b <- .ui$iniDf$backTransform[.w]
      if (!is.na(.b)) {
        .est <- .ret$popDf$Estimate[i]
        .se <- .ret$popDf$SE[i]
        .bt <- .ret$popDf[["Back-transformed"]]
        .bt[i] <- get(.b, envir=globalenv())(.est)
        .ret$popDf[["Back-transformed"]] <- .bt
        if (is.na(.se)) {
          .i1 <- .ret$popDf[["CI Lower"]]
          .low <- get(.b, envir=globalenv())(.est - .se * .qn)
          .i1[i] <- .low
          .ret$popDf[["CI Lower"]] <- .i1
          .i1 <- .ret$popDf[["CI Upper"]]
          .hi <- get(.b, envir=globalenv())(.est + .se * .qn)
          .i1[i] <- .hi
          .ret$popDf[["CI Upper"]] <- .i1
          .bt2 <- .ret$popDfSig[[.btName]]
          .bt2[i] <- sprintf(.fmt, .est, .low, .hi)
          .ret$popDfSig[[.btName]] <- .bt2
        } else {
          .bt2 <- .ret$popDfSig[[.btName]]
          .bt2[i] <- sprintf(.fmt, .est, .low, .hi)
          .ret$popDfSig[[.btName]] <- .bt2
          .bt2 <- .ret$popDfSig[[.btName]]
          .bt2[i] <- sprintf(.fmt2, .est)
          .ret$popDfSig[[.btName]] <- .bt2
        }
      }
    }
  })
}

#' This gest the CV/SD for a single ETA
#'
#'
#' @param .eta Eta Name
#' @param .env Environment where the indicators of `.sdOnly`, `.cvOnly` are stored so the column name can be changed to match the data
#' @param .ome Omega fixed vector
#' @param .muRefCurEval The current mu ref evaluation.  This determines if the ETA is logit normal and %CV should be calculated.
#' @param .sigdig is the number of significant digits used in the evaulation
#' @return Data frame row with ch= the charaacter representation and v is the vector representation of the CV or sd
#' @author Matthew L. Fidler
#' @noRd
.updateParFixedGetEtaRow <- function(.eta, .env, .ome, .omegaFix, .muRefCurEval, .sigdig) {
  .v <- .ome[.eta, .eta]
  .w <- which(.muRefCurEval$parameter == .eta)
  if (.muRefCurEval$curEval[.w] == "exp") {
    assign(".sdOnly", FALSE, envir=.env)
    return(data.frame(
      ch = paste0(
        ifelse(.omegaFix[.eta], "fix(", ""),
        formatC(signif(sqrt(exp(.v) - 1) * 100, digits = .sigdig),
                digits = .sigdig, format = "fg", flag = "#"),
        ifelse(.omegaFix[.eta], ")", "")
      ),
      v = sqrt(exp(.v) - 1) * 100))
  } else {
    assign(".cvOnly", FALSE, envir=.env)
    return(data.frame(
      ch = paste0(
        ifelse(.omegaFix[.eta], "fix(", ""),
        formatC(signif(sqrt(.w), digits = .sigdig),
                digits = .sigdig, format = "fg", flag = "#"),
        ifelse(.omegaFix[.eta], ")", "")),
      v = .v))
  }
}

#'  This will add the between subject varaibility to the mu-referenced theta.  It also expands the table to include non-mu referenced ETAs
#'
#'
#' @param .ret The focei return environment
#' @param .ui The rxode2 ui environment
#' @return Nothing called for side effects on popDf and popDfSig in the .ret environment
#' @author Matthew L. Fidler
#' @noRd
.updateParFixedAddBsv <- function(.ret, .ui) {
  .ome <- .ret$omega
  .omegaFix <- as.data.frame(.ui$ini)
  .omegaFix <- .omegaFix[is.na(.omegaFix$ntheta), ]
  .omegaFix <- setNames(.omegaFix$fix, paste(.omegaFix$name))
  .sigdig <- rxode2::rxGetControl(.ui, "sigdig", 3L)

  .muRefDataFrame <- .ui$muRefDataFrame
  .muRefCurEval   <- .ui$muRefCurEval
  .env <- new.env(parent=emptyenv())
  .env$.cvOnly <- TRUE
  .env$.sdOnly <- TRUE
  .env$.muRefVars <- NULL
  .cvp <- lapply(row.names(.ret$popDfSig), function(x) {
    .w <- which(.muRefDataFrame$theta == x)
    if (length(.w) != 1) {
      return(data.frame(ch = " ", v = NA_real_))
    }
    .eta <- .muRefDataFrame$eta[.w]
    assign(".muRefVars", c(.env$.muRefVars, .eta), envir=.env)
    .updateParFixedGetEtaRow(.eta, .env, .ome, .omegaFix, .muRefCurEval, .sigdig)
  })
  .cvp <- do.call("rbind", .cvp)
  .nonMuRef <- setdiff(dimnames(.ome)[[1]], .env$.muRefVars)
  if (length(.nonMuRef) > 0) {
    .ret$popDfSig2 <- as.data.frame(lapply(names(.ret$popDfSig), function(x) { rep("", length(.nonMuRef))}))
    names(.ret$popDfSig2) <- names(.ret$popDfSig)
    .ret$popDf2 <- as.data.frame(lapply(names(.ret$popDf), function(x) { rep(NA_real_, length(.nonMuRef))}))
    names(.ret$popDf2) <- names(.ret$popDf)
    row.names(.ret$popDfSig2) <- .nonMuRef
    row.names(.ret$popDf2) <- .nonMuRef
  }
  .ret$popDfSig <- data.frame(.ret$popDfSig, "BSD" = .cvp$ch, check.names = FALSE)
  .ret$popDf <- data.frame(.ret$popDf, "BSD" = .cvp$v, check.names = FALSE)
  if (length(.nonMuRef) > 0) {
    .cvp <- lapply(row.names(.ret$popDfSig2), function(x) {
      .updateParFixedGetEtaRow(x, .env, .ome, .omegaFix, .muRefCurEval, .sigdig)
    })
    .cvp <- do.call("rbind", .cvp)
    .ret$popDfSig2 <- data.frame(.ret$popDfSig2, "BSD" = .cvp$ch, check.names = FALSE)
    .ret$popDf2 <- data.frame(.ret$popDf2, "BSD" = .cvp$v, check.names = FALSE)
    .ret$popDfSig <- rbind(.ret$popDfSig, .ret$popDfSig2)
    .ret$popDf <- rbind(.ret$popDf, .ret$popDf2)
    rm(list=c("popDf2", "popDfSig2"), envir=.ret)
  }
  .w <- which(names(.ret$popDfSig) == "BSD")
  if (length(.w) == 1) {
    names(.ret$popDfSig)[.w] <- ifelse(.env$.sdOnly, "BSV(SD)", ifelse(.env$.cvOnly, "BSV(CV%)", "BSV(CV% or SD)"))
  }
  .w <- which(names(.ret$popDf) == "BSD")
  if (length(.w) == 1) {
    names(.ret$popDf)[.w] <- ifelse(.env$.sdOnly, "BSV(SD)", ifelse(.env$.cvOnly, "BSV(CV%)", "BSV(CV% or SD)"))
  }
}

.updateParFixedAddShrinkage <- function(.ret, .ui) {
  .shrink <- .ret$shrink
  .errs <- as.data.frame(.ui$ini)
  .errs <- paste(.errs[which(!is.na(.errs$err)), "name"])
  .muRefDataFrame <- .ui$muRefDataFrame
  .sigdig <- rxode2::rxGetControl(.ui, "sigdig", 3L)
  .sh <- lapply(row.names(.ret$popDfSig), function(x) {
    .w <- which(.muRefDataFrame$theta == x)
    if (length(.w) != 1) {
      .w <- which(names(.shrink) == x)
      if (length(.w) != 1) {
        if (any(x == .errs)) {
          x <- "IWRES"
        } else {
          return(data.frame(ch = " ", v = NA_real_))
        }
      }
      .eta <- x
    } else {
      .eta <- .muRefDataFrame$eta[.w]
    }
    .v <- .shrink[7, .eta]
    if (length(.v) != 1) {
      return(data.frame(ch = " ", v = NA_real_))
    }
    if (is.na(.v)) {
      return(data.frame(ch = " ", v = NA_real_))
    }
    .t <- ">"
    if (.v < 0) {
    } else if (.v < 20) {
      .t <- "<"
    } else if (.v < 30) {
      .t <- "="
    }
    return(data.frame(
      ch = sprintf("%s%%%s", formatC(signif(.v, digits = .sigdig),
                                     digits = .sigdig, format = "fg", flag = "#"
                                     ), .t),
      v = .v
    ))
  })
  .sh <- do.call("rbind", .sh)
  .ret$popDfSig <- data.frame(.ret$popDfSig, "Shrink(SD)%" = .sh$ch, check.names = FALSE)
  .ret$popDf <- data.frame(.ret$popDf, "Shrink(SD)%" = .sh$v, check.names = FALSE)
}
#' Create the parFixed dataset
#'
#' @param .ret focei style environment
#' @return Nothing, called for side effects
#' @author Matthew L. Fidler
#' @noRd
.updateParFixed <- function(.ret) {
  .ui <- .ret$ui
  .updateParFixedApplyManualBacktransformations(.ret, .ui)
  .updateParFixedAddParameterLabel(.ret, .ui)
  .updateParFixedAddBsv(.ret, .ui)
  .updateParFixedAddShrinkage(.ret, .ui)
  .ret$parFixed <- .ret$popDfSig
  .ret$parFixedDf <- .ret$popDf
  rm(list=c("popDfSig", "popDf"), envir=.ret)
  class(.ret$parFixed) <- c("nlmixr2ParFixed", "data.frame")
}

.nmObjEnsureObjective <- function(obj) {
  .env <- obj
  if (!is.null(obj$saem)) {
    .tmp <- obj$saem
    .curObj <- get("objective", .env)
    if (is.na(.curObj)) {
      .nnodes <- 3
      if (exists("nnodes.gq", .env)) {
        .nnodes <- .env$nnodes.gq
      }
      .nsd <- 1.6
      if (exists("nsd.gq", .env)) {
        .nsd <- .env$nsd.gq
      }
      if (.nnodes == 1) {
        .tmp <- try(setOfv(obj, paste0("laplace", .nsd)), silent = TRUE)
      } else {
        .tmp <- try(setOfv(obj, paste0("gauss", .nnodes, "_", .nsd)), silent = TRUE)
      }
      if (inherits(.tmp, "try-error")) {
        message("gaussian quadrature failed, changed to focei")
        setOfv(obj, "focei")
      }
    }
  }
}

#' Get an item from a nlmixr core object
#'
#' @param x A specialized list with:
#' - First argument is a nlmixrFitCore environment
#' - Second argument is if the exact argument is requested
#' - The class would be the requested argument name followed by the class "nmObjGet"
#' @param ... Other arguments
#' @return Value of argument or NULL
#' @author Matthew L. Fidler
#' @keywords internal
nmObjGet <- function(x, ...) {
  if (!inherits(x, "nmObjGet")) {
    stop("'x' is wrong type for 'nmObjGet'", call.=FALSE)
  }
  .arg <- class(x)[1]
  if (any(.arg == c(
    "logLik", "value", "obf", "ofv",
    "objf", "OBJF", "objective", "AIC",
    "BIC"))) {
    .nmObjEnsureObjective(x[[1]])
  }
  UseMethod("nmObjGet")
}

#' @rdname nmObjGet
#' @export
nmObjGet.default <- function(x, ...) {
  .arg <- class(x)[1]
  .env <- x[[1]]
  if (exists(.arg, envir = .env)) {
    .ret <- get(.arg, envir = .env)
    if (inherits(.ret, "raw")) .ret <- qs::qdeserialize(.ret)
    return(.ret)
  }
  .lst <- list(get("ui", envir=.env), x[[2]])
  class(.lst) <- c(.arg, "rxUiGet")
  .ret <- rxUiGet(.lst)
  .ret
}

#' @rdname nmObjGet
#' @export
nmObjGet.cor <- function(x, ...) {
  .obj <- x[[1]]
  .cov <- .obj$cov
  .sd2 <- sqrt(diag(.cov))
  .cor <- stats::cov2cor(.cov)
  dimnames(.cor) <- dimnames(.cov)
  diag(.cor) <- .sd2
  .cor
}
attr(nmObjGet.cor, "desc") <- "correlation matrix of theta, calculated from covariance of theta"

#' @rdname nmObjGet
#' @export
nmObjGet.omegaR <- function(x, ...) {
  .obj <- x[[1]]
  .cov <- .obj$omega
  .sd2 <- sqrt(diag(.cov))
  .cor <- stats::cov2cor(.cov)
  dimnames(.cor) <- dimnames(.cov)
  diag(.cor) <- .sd2
  .cor
}
attr(nmObjGet.omegaR, "desc") <- "correlation matrix of omega"

#' @rdname nmObjGet
#' @export
nmObjGet.dataSav <- function(x, ...){
  .obj <- x[[1]]
  .objEnv <- .obj$env
  if (exists("dataSav", .objEnv)) return(get("dataSav", envir=.objEnv))
  .data <- .obj$origData
  .env <- new.env(emptyenv())
  .foceiPreProcessData(.data, .env, .obj$ui)
  .env$dataSav
}
#attr(nmObjGet.dataSav, "desc") <- "data that focei sees for optimization"

#' @rdname nmObjGet
#' @export
nmObjGet.saemTransformedData <- function(x, ...) {
  .dataSav <- nmObjGet.dataSav(x, ...)
  .obj <- x[[1]]
  .ui <- .obj$ui
  .saemGetDataForFit(.dataSav, .ui)
}
#attr(nmObjGet.saemTransformedData, "desc") <- "data that saem sees for optimization"


#' @rdname nmObjGet
#' @export
nmObjGet.parHistStacked <- function(x, ...) {
  .obj <- x[[1]]
  .env <- .obj$env
  if (exists("parHist", envir=.env)) {
    .parHist <- .env$parHist
    .iter <- .parHist$iter
    .ret <- data.frame(iter=.iter, stack(.parHist[, -1]))
    names(.ret) <- sub("values", "val",
                       sub("ind", "par", names(.ret)))
    .ret
  } else {
    NULL
  }
}
attr(nmObjGet.parHistStacked, "desc") <- "stacked parameter history"



#' @rdname nmObjGet
#' @export
nmObjGet.md5 <- function(x, ...) {
  .nlmixr2Md5(x[[1]])
}

#' @rdname nmObjGet
#' @export
nmObjGet.notes <- function(x, ...) {
  .notesFit(x[[1]])
}

#' @rdname nmObjGet
#' @export
nmObjGet.sigma <- function(x, ...) {
  .sigma(x[[1]])
}

#' @rdname nmObjGet
#' @export
nmObjGet.coefficients <- function(x, ...) {
  list(fixed = fixef(x[[1]]),
       random = ranef(x[[1]]))
}

#' @rdname nmObjGet
#' @export
nmObjGet.env <- function(x, ...) {
  x[[1]]
}

#' @rdname nmObjGet
#' @export
nmObjGet.condition <- function(x, ...) {
  .env <- x[[1]]
  .objDf <- .env$objDf
  if (any(names(.objDf) == "Condition Number")) {
    .cn <- .objDf[, "Condition Number"]
    .cn <- .cn[!is.na(.cn)]
    return(.cn)
  }
  return(NULL)
}

#' @rdname nmObjGet
#' @export
nmObjGet.simInfo <- function(x, ...) {
  .simInfo(x[[1]])
}

#' @rdname nmObjGet
#' @export
nmObjGet.seed <- function(x, ...) {
  .env <- x[[1]]
  if (exists("saem", .env)) {
    attr(get("saem", .env), "saem.cfg")$seed
  }
  NULL
}
#' @rdname nmObjGet
#' @export
nmObjGet.saemCfg <- function(x, ...) {
  .env <- x[[1]]
  if (exists("saem", .env)) {
    return(attr(get("saem", .env), "saem.cfg"))
  }
}

.nmObjBackward <- c(
  "value"="objf",
  "obf"="objf",
  "ofv"="objf",
  "par.hist"="parHist",
  "par.hist.stacked"="parHistStacked",
  "omega.R"="omegaR",
  "par.fixed"="parFixed",
  "eta"="ranef",
  "theta"="fixef",
  "varFix"="cov",
  "thetaMat"="cov",
  "modelName"="model.name",
  "dataName"="data.name",
  "saem.cfg"="saemCfg",
  "objf"="objective",
  "OBJF"="objective"
)

#' @export
`$.nlmixr2FitCore` <- function(obj, arg, exact = FALSE) {
  .env <- obj
  .arg <- .nmObjBackward[arg]
  if (is.na(.arg)) .arg <- arg
  .lst <- list(obj, exact)
  class(.lst) <- c(.arg, "nmObjGet")
  .ret <- nmObjGet(.lst)
  if (!is.null(.ret)) return(.ret)
}

 #' @export
`$.nlmixr2FitCoreSilent` <- `$.nlmixr2FitCore`

#' @export
`$.nlmixr2FitData` <- function(obj, arg, exact = FALSE) {
  .ret <- obj[[arg]]
  if (arg == "md5") {
    return(.nlmixr2Md5(obj))
  } else if (is.null(.ret)) {
    .cls <- class(obj)
    .env <- attr(.cls, ".foceiEnv")
    .ret <- `$.nlmixr2FitCore`(.env, arg, exact)
  }
  return(.ret)
}


#' @importFrom nlme VarCorr
#' @export
VarCorr.nlmixr2FitCore <- function(x, sigma = NULL, ...) {
  .ret <- x$nlme
  if (is.null(.ret)) {
    .var <- diag(x$omega)
    .ret <- data.frame(
      Variance = .var, StdDev = sqrt(.var),
      row.names = names(.var)
    )
    .ret <- .ret[!is.na(.ret[, 1]), ]
    return(.ret)
  } else {
    VarCorr(.ret, ...)
  }
}

#' @export
VarCorr.nlmixr2FitCoreSilent <- VarCorr.nlmixr2FitCore

.sigma <- function(x) {
  .ret <- x$nlme
  if (is.null(.ret)) {
    if (exists("uif", envir = x$env)) {
      .df <- as.data.frame(x$uif$ini)
      .errs <- paste(.df[which(!is.na(.df$err)), "name"])
      return(fixef(x)[.errs])
    }
  } else {
    return(.ret$sigma)
  }
}

#' @export
str.nlmixr2FitData <- function(object, ...) {
  NextMethod(object)
  .env <- object$env
  ## cat(" $ par.hist         : Parameter history (if available)\n")
  ## cat(" $ par.hist.stacked : Parameter history in stacked form for easy plotting (if available)\n")
  cat(" $ omega            : Omega matrix\n")
  cat(" $ omegaR           : Omega Correlation matrix\n")
  cat(" $ shrink           : Shrinkage table, includes skewness, kurtosis, and eta p-values\n")
  cat(" $ parFixed         : Fixed Effect Parameter Table\n")
  cat(" $ theta            : Fixed Parameter Estimates\n")
  cat(" $ eta              : Individual Parameter Estimates\n")
  cat(" $ seed             : Seed (if applicable)\n")
  cat(" $ coefficients     : Fixed and random coefficients\n")
  if (exists("uif", envir = object$env)) {
    cat(" $ meta             : Model meta information environment\n")
    cat(" $ modelName        : Model name (from R function)\n")
    cat(" $ dataName         : Name of R data input\n")
    cat(" $ simInfo          : rxode2 list for simulation\n")
    cat(" $ sigma            : List of sigma components and their values\n")
  }
}


#' Extract residuals from the FOCEI fit
#'
#' @param object focei.fit object
#' @param ... Additional arguments
#' @param type Residuals type fitted.
#' @return residuals
#' @author Matthew L. Fidler
#' @export
residuals.nlmixr2FitData <- function(object, ..., type = c("ires", "res", "iwres", "wres", "cwres", "cpred", "cres")) {
  return(object[, toupper(match.arg(type))])
}

#' Return the objective function
#'
#' @param x object to return objective function value
#' @param type Objective function type value to retrieve or add.
#'
#' \itemize{
#'
#' \item{focei} For most models you can specify "focei" and it will
#' add the focei objective function.
#'
#' \item{nlme} This switches/chooses the nlme objective function if
#'    applicable.  This objective function cannot be added if it
#'    isn't present.
#'
#' \item{fo} FO objective function value. Cannot be generated
#'
#' \item{foce} FOCE object function value. Cannot be generated
#'
#' \item{laplace#} This adds/retrieves  the Laplace objective function value.
#' The \code{#} represents the number of standard deviations
#' requested when expanding the Gaussian Quadrature.  This can
#' currently only be used with saem fits.
#'
#' \item{gauss#.#} This adds/retrieves the Gaussian Quadrature
#' approximation of the objective function.  The first number is the
#' number of nodes to use in the approximation. The second number is
#' the number of standard deviations to expand upon.
#'
#' }
#'
#' @param ... Other arguments sent to ofv for other methods.
#'
#' @return Objective function value
#'
#' @author Matthew Fidler
#'
#' @export
ofv <- function(x, type, ...) {
  UseMethod("ofv")
}

#' @export
ofv.nlmixr2FitData <- function(x, type, ...) {
  if (!missing(type)) setOfv(x, type)
  return(x$ofv)
}

#' @export
logLik.nlmixr2FitData <- function(object, ...) {
  .objName <- substitute(object)
  .lst <- list(...)
  if (!is.null(.lst$type)) {
    .new <- setOfv(object, .lst$type)
    .parent <- globalenv()
    .bound <- do.call("c", lapply(ls(.parent, all.names = TRUE), function(.cur) {
      if (.cur == .objName && identical(.parent[[.cur]]$env, object$env)) {
        return(.cur)
      }
      return(NULL)
    }))
    if (length(.bound) == 1) {
      if (exists(.bound, envir = .parent)) {
        assign(.bound, .new, envir = .parent)
      }
    }
    return(get("logLik", .new$env))
  } else {
    return(object$logLik)
  }
}

#' @export
logLik.nlmixr2FitCore <- function(object, ...) {
  object$logLik
}

#' @export
nobs.nlmixr2FitCore <- function(object, ...) {
  object$nobs
}

#' @export
vcov.nlmixr2FitCore <- function(object, ...) {
  object$cov
}
#' This gets the parsed data in the lower-level manner that nlmixr2 expects.
#'
#' @param object nlmixr2 Object
#'
#' @return Gets the parsed data
#'
#' @export
#'
#' @author Matthew L. Fidler
#' @keywords internal
.nmGetData <- function(object, keep=NULL) {
  if (is.null(keep)) keep <- character(0)
  .uif <- object$uif
  .tmp <- deparse(body(.uif$theta.pars))[-1]
  .tmp <- .tmp[-length(.tmp)]
  return(rxode2::etTrans(object$origData, paste(paste(.tmp, collapse = "\n"), "\n", .uif$rxode), TRUE, TRUE, TRUE, keep=keep))
}

#' @export
getData.nlmixr2FitCore <- function(object) {
  object$origData
}

#' @export
ranef.nlmixr2FitCore <- function(object, ...) {
  object$ranef
}

#' @export
fixef.nlmixr2FitCore <- function(object, ...) {
  object$fixef
}
#' @export
fixef.nlmixr2FitCoreSilent <- fixef.nlmixr2FitCore

#' @export
ranef.nlmixr2FitCoreSilent <- ranef.nlmixr2FitCore

#' @export
getData.nlmixr2FitCoreSilent <- getData.nlmixr2FitCore

#' @export
logLik.nlmixr2FitCoreSilent <- logLik.nlmixr2FitCore

#' @export
nobs.nlmixr2FitCoreSilent <- nobs.nlmixr2FitCore

#' @export
vcov.nlmixr2FitCoreSilent <- vcov.nlmixr2FitCore

#' Update model to have final parameter estimates for piping and save orig data
#'
#' @param x Data to fix
#' @return Nothing, called for side effects
#' @noRd
.nlmixr2FitUpdateParams <- function(x) {
  # Update initial estimates to match current initial estimates
  .ui <- x$ui
  if (exists("fullTheta", x)) {
    .thetas <- x$fullTheta
  } else {
    .thetas <- x$theta
  }
  for (.n in names(.thetas)) {
    .ui$iniDf$est[.ui$iniDf$name == .n] <- .thetas[.n]
  }
  .omega <- x$omega
  for (.i in seq_along(.ui$iniDf$neta1)) {
    if (!is.na(.ui$iniDf$neta1[.i])) {
      .ui$iniDf$est[.i] <- .omega[.ui$iniDf$neta1[.i], .ui$iniDf$neta2[.i]]
    }
  }
  assign("ui", .ui, envir=x)
}
