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
  NULL
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
#' Apply the manual translation for only one item
#'
#' @param i Which parameter should be updated
#' @param .ret return environment
#' @param .ui user interface function
#' @param .qn The qn for the
#' @param .btName Back-transform name
#' @param .fmt format for one estimate
#' @param .fmt2 format for estimate and ci
#' @return nothing, called for side effects
#' @keywords internal
#' @author Matthew L. Fidler
#' @noRd
.updateParFixedApplyManualBacktransformationsI <- function(i, .ret, .ui,
                                                           .qn, .btName,
                                                           .fmt, .fmt2) {
  theta <- row.names(.ret$popDf)[i]
  .w <- which(.ui$iniDf$name == theta)
  if (length(.w) == 1L) {
    .b <- .ui$iniDf$backTransform[.w]
    if (!is.na(.b)) {
      .bfun <- try(get(.b, envir=nlmixr2global$nlmixrEvalEnv$envir, mode="function"), silent=TRUE)
      if (inherits(.bfun, "try-error")) {
        warning("unknown function '", .b, "' for manual backtransform, revert to nlmixr2 back-transformation detection for, '", theta, "'",
                call.=FALSE)
        return(invisible())
      }
      .est <- .ret$popDf$Estimate[i]
      .se <- .ret$popDf$SE[i]
      .bt <- .ret$popDf[["Back-transformed"]]
      .bt[i] <- .bfun(.est)
      .ret$popDf[["Back-transformed"]] <- .bt
      if (!is.na(.se)) {
        .i1 <- .ret$popDf[["CI Lower"]]
        .low <- .bfun(.est - .se * .qn)
        .i1[i] <- .low
        .ret$popDf[["CI Lower"]] <- .i1
        .i1 <- .ret$popDf[["CI Upper"]]
        .hi <- .bfun(.est + .se * .qn)
        .i1[i] <- .hi
        .ret$popDf[["CI Upper"]] <- .i1
        .bt2 <- .ret$popDfSig[[.btName]]
        .bt2[i] <- sprintf(.fmt, .bfun(.est),
                           .low, .hi)
        .ret$popDfSig[[.btName]] <- .bt2
      } else {
        .bt2 <- .ret$popDfSig[[.btName]]
        .bt2[i] <- sprintf(.fmt2, .bfun(.est))
        .ret$popDfSig[[.btName]] <- .bt2
      }
    }
  }
}

#'  This applies the manually specified back-transformations
#'
#' @param .ret focei environment
#' @return Nothing, called for side effects
#' @author Matthew L. Fidler
#' @noRd
.updateParFixedApplyManualBacktransformations <- function(.ret, .ui) {
  .qn <- qnorm(1.0-(1-.ret$control$ci)/2)
  .n <- names(.ret$popDfSig)
  if (length(.n) >= 4) {
    .btName <- names(.ret$popDfSig)[4]
  } else if (length(.n) >= 2) {
    .btName <- names(.ret$popDfSig)[2]
  } else {
    return(NULL)
  }
  .sigdig <- rxode2::rxGetControl(.ui, "sigdig", 3L)
  .fmt <- paste0("%", .sigdig, "g (%", .sigdig, "g, %", .sigdig, "g)")
  .fmt2 <- paste0("%", .sigdig, "g")
  lapply(seq_along(.ret$popDf$Estimate), .updateParFixedApplyManualBacktransformationsI,
         .ret=.ret, .ui=.ui, .qn=.qn, .btName=.btName,
         .fmt=.fmt, .fmt2=.fmt2)
}

#' This gets the CV/SD for a single ETA
#'
#'
#' @param .eta Eta Name
#' @param .env Environment where the indicators of `.sdOnly`, `.cvOnly` are stored so the column name can be changed to match the data
#' @param .ome Omega fixed vector
#' @param .muRefCurEval The current mu ref evaluation.  This determines if the ETA is logit normal and %CV should be calculated.
#' @param .sigdig is the number of significant digits used in the evaluation
#' @return Data frame row with ch= the character representation and v is the vector representation of the CV or sd
#' @author Matthew L. Fidler and Bill Denney
#' @noRd
.updateParFixedGetEtaRow <- function(.eta, .env, .ome, .omegaFix, .muRefCurEval, .sigdig) {
  if (is.null(.ome)) {
    # This can happen if there are no BSV parameters in a model
    return("")
  } else if (!(.eta %in% rownames(.ome))) {
    # This can happen if .eta is a fixed BSV parameter
    return("")
  }
  .v <- .ome[.eta, .eta]
  .w <- which(.muRefCurEval$parameter == .eta)
  if (length(.w) == 1L && .muRefCurEval$curEval[.w] == "exp") {
    assign(".sdOnly", FALSE, envir=.env)
    .valNumber <- sqrt(exp(.v) - 1) * 100
    .valCharPrep <- .valNumber
  } else {
    assign(".cvOnly", FALSE, envir=.env)
    .valNumber <- .v
    .valCharPrep <- sqrt(.v)
  }
  if (.eta %in% names(.omegaFix) && .omegaFix[.eta]) {
    .charPrefix <- "fix("
    .charSuffix <- ")"
  } else {
    .charPrefix <- ""
    .charSuffix <- ""
  }
  .valChar <-
    formatC(
      signif(.valCharPrep, digits = .sigdig),
      digits = .sigdig, format = "fg", flag = "#"
    )
  data.frame(
    ch = paste0(.charPrefix, .valChar, .charSuffix),
    v = .valNumber
  )
}

#'  This will add the between subject variability to the mu-referenced theta.  It also expands the table to include non-mu referenced ETAs
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
    data.frame(
      ch = sprintf("%s%%%s", formatC(signif(.v, digits = .sigdig),
                                     digits = .sigdig, format = "fg", flag = "#"
                                     ), .t),
      v = .v
    )
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
  if (!is.null(nlmixr2global$nlmixr2EstEnv$uiUnfix)) {
    .ui <- nlmixr2global$nlmixr2EstEnv$uiUnfix
    .theta <- .ui$theta
    .tn <- names(.theta)
    .fmt <- paste0("%.", .ret$control$sigdig, "g")
    .row.names <- row.names(.ret$popDf)
    .popDf <-
      data.frame(
        `Estimate`=vapply(.tn,
                          function(n) {
                            .w <- which(.row.names == n)
                            if (length(.w) ==1L) return(setNames(.ret$popDf[.w, "Estimate"], NULL))
                            .theta[n]
                          }, double(1), USE.NAMES = FALSE),
        `SE`=vapply(.tn,
                    function(n) {
                      .ret <- .ret$popDf[n, "SE"]
                      setNames(.ret, NULL)
                    }, double(1), USE.NAMES = FALSE),
        `%RSE`=vapply(.tn,
                    function(n) {
                      .ret <- .ret$popDf[n, "%RSE"]
                      setNames(.ret, NULL)
                    }, double(1), USE.NAMES = FALSE),
        `Back-transformed`=vapply(.tn,
                                  function(n) {
                                    .w <- which(.row.names == n)
                                    if (length(.w) ==1L) return(setNames(.ret$popDf[.w, "Back-transformed"], NULL))
                                    .theta[n]
                                  }, double(1), USE.NAMES = FALSE),
        `CI Lower`=vapply(.tn,
                      function(n) {
                        .ret <- .ret$popDf[n, "CI Lower"]
                        setNames(.ret, NULL)
                      }, double(1), USE.NAMES = FALSE),
        `CI Upper`=vapply(.tn,
                          function(n) {
                            .ret <- .ret$popDf[n, "CI Upper"]
                            setNames(.ret, NULL)
                           }, double(1), USE.NAMES = FALSE),
        row.names = .tn,
        check.rows = FALSE, check.names = FALSE
      )
    if (any(names(.ret$popDfSig) == "SE")) {
      .popDfSig <-
        data.frame(
          `Est.`=vapply(.tn,
                        function(n) {
                          .w <- which(.row.names == n)
                          if (length(.w) ==1L) return(setNames(.ret$popDfSig[.w, "Est."], NULL))
                          sprintf(.fmt, .theta[n])
                        }, character(1), USE.NAMES = FALSE),
          `SE`=vapply(.tn,
                      function(n) {
                        .w <- which(.row.names == n)
                        if (length(.w) == 1L) return(setNames(.ret$popDfSig[.w, "SE"], NULL))
                        "FIXED"
                      }, character(1), USE.NAMES = FALSE),
          `%RSE`=vapply(.tn,
                        function(n) {
                          .w <- which(.row.names == n)
                          if (length(.w) == 1L) return(setNames(.ret$popDfSig[.w, "%RSE"], NULL))
                          "FIXED"
                        }, character(1), USE.NAMES = FALSE),
          `Back-transformed`=vapply(.tn,
                                    function(n) {
                                      .w <- which(.row.names == n)
                                      if (length(.w) == 1L) return(setNames(.ret$popDfSig[.w, 4], NULL))
                                      sprintf(.fmt, .theta[n])
                                    }, character(1), USE.NAMES = FALSE),
          row.names = .tn,
          check.rows = FALSE, check.names = FALSE
        )
      names(.popDfSig)[4] <- names(.ret$popDfSig)[4]
    } else {
      .popDfSig <-
        data.frame(
          `Est.`=vapply(.tn,
                        function(n) {
                          .w <- which(.row.names == n)
                          if (length(.w) ==1L) return(setNames(.ret$popDfSig[.w, "Est."], NULL))
                          sprintf(.fmt, .theta[n])
                        }, character(1), USE.NAMES = FALSE),
          `Back-transformed`=vapply(.tn,
                                    function(n) {
                                      .w <- which(.row.names == n)
                                      if (length(.w) == 1L) return(setNames(.ret$popDfSig[.w, 2], NULL))
                                      sprintf(.fmt, .theta[n])
                                    }, character(1), USE.NAMES = FALSE),
          row.names = .tn,
          check.rows = FALSE, check.names = FALSE
        )
      names(.popDfSig)[2] <- names(.ret$popDfSig)[2]
    }
    .ret$popDfSig <- .popDfSig
    .ret$popDf <- .popDf
  }
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
      if (exists("nnodesGq", .env)) {
        .nnodes <- .env$nnodesGq
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
  "model.name"="modelName",
  "dataName"="data.name",
  "saem.cfg"="saemCfg",
  "objf"="objective",
  "OBJF"="objective",
  "theta"="fixef",
  "etaR"="phiR",
  "etaH"="phiH",
  "etaC"="phiC",
  "etaSE"="phiSE",
  "etaRSE"="phiRSE",
  "uiIni"="iniUi"
)

#' @export
`$.nlmixr2FitCore` <- function(obj, arg, exact = FALSE) {
  rxode2::.udfEnvSet(parent.frame(1))
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
  rxode2::.udfEnvSet(parent.frame(1))
  .ret <- obj[[arg]]
  if (arg == "md5") {
    return(.nlmixr2Md5(obj))
  } else if (is.null(.ret)) {
    .lst <- list(obj, exact)
    class(.lst) <- c(arg, "nmObjGetData")
    .ret <- nmObjGetData(.lst)
    if (is.null(.ret)) {
      .cls <- class(obj)
      .env <- attr(.cls, ".foceiEnv")
      .ret <- `$.nlmixr2FitCore`(.env, arg, exact)
    }
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
    if (exists("ui", envir = x$env)) {
      .df <- as.data.frame(x$uif$ini)
      .errs <- paste(.df[which(!is.na(.df$err)), "name"])
      return(fixef(x)[.errs])
    }
  }
  .ret$sigma
}

#' @export
str.nlmixr2FitData <- function(object, ...) {
  NextMethod(object)
  .s <- .nmObjGetSupportedDollars()
  cat(paste(strtrim(paste(vapply(names(.s), function(x){
    .nchar <- nchar(x)
    if (.nchar >= 10) {
      paste0(" $ ", x, ": ")
    } else {
      paste0(" $ ",x, paste(rep(" ", 10 - .nchar), collapse=""), ": ")
    }
  }, character(1), USE.NAMES=FALSE), .s), 128), collapse="\n"))
  cat("\n")
  invisible()
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
  x$ofv
}

#' @export
logLik.nlmixr2FitData <- function(object, ...) {
  .objName <- substitute(object)
  .lst <- list(...)
  if (!is.null(.lst$type)) {
    .new <- setOfv(object, .lst$type)
    .parent <- globalenv()
    .bound <- do.call("c", lapply(ls(.parent, all.names = TRUE),
      function(.cur) {
       if (.cur == .objName && identical(.parent[[.cur]]$env, object$env)) {
         return(.cur)
       }
       NULL
      }))
    if (length(.bound) == 1) {
      if (exists(.bound, envir = .parent)) {
        assign(.bound, .new, envir = .parent)
      }
    }
    get("logLik", .new$env)
  } else {
    object$logLik
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
#' @keywords internal
#' @export
.nlmixr2FitUpdateParams <- function(x) {
  # Update initial estimates to match current initial estimates
  .ui <- x$ui
  .iniDf <- .ui$iniDf
  assign("iniDf0", nlmixr2global$nlmixr2EstEnv$iniDf0, envir=x)
  if (exists("fullTheta", x)) {
    .thetas <- x$fullTheta
  } else if (exists("fixef", x)) {
    .thetas <- get("fixef", x)
  } else {
    .thetas <- x$theta
  }
  if (is.null(names(.thetas))) {
    .thetaNames <- .iniDf$name[which(!is.na(.iniDf$ntheta))]
    if (length(.thetaNames) > length(.thetas)) stop("corrupted rxode2 ui", call.=FALSE)
    .thetas <- .thetas[seq_along(.thetaNames)]
    names(.thetas) <- .thetaNames
  }
  for (.n in names(.thetas)) {
    .iniDf$est[.iniDf$name == .n] <- .thetas[.n]
  }
  # In the case of nlme, it estimates the whole covariance matrix, even if you don't want it to.
  # Allow the omega to expand the initial estimates if needed.
  .omega <- x$omega
  if (is.null(.omega)) {
    .ui <- rxode2::rxUiDecompress(.ui)
    assign("iniDf", .iniDf, envir=.ui)
    .ui <- rxode2::rxUiCompress(.ui)
    assign("ui", .ui, envir=x)
  } else {
    .fixComps <- .iniDf[is.na(.iniDf$ntheta),]
    .fixComps <- setNames(.fixComps$fix, .fixComps$name)
    .lotri <- lotri::as.lotri(.iniDf)
    attr(.omega, "lotriEst") <- attr(.lotri, "lotriEst")
    class(.omega) <- class(.lotri)
    .iniDf1 <- .iniDf[is.na(.iniDf$neta1), ]
    .iniDf2 <- as.data.frame(.omega)
    .iniDf2 <- .iniDf2[!is.na(.iniDf2$neta1), ]
    .iniDf2$err <- NA_character_
    .names <- names(.iniDf)
    .iniDf <- rbind(.iniDf1, .iniDf2)
    .iniDf <- .iniDf[, .names]
    for (.n in names(.fixComps)) {
      .w  <- which(.iniDf$name == .n)
      if (length(.w) == 1L) .iniDf[.w, "fix"] <- .fixComps[.n]
    }
    .ui <- rxode2::rxUiDecompress(.ui)
    assign("iniDf", .iniDf, envir=.ui)
    .ui <- rxode2::rxUiCompress(.ui)
    assign("ui", .ui, envir=x)
  }
}
