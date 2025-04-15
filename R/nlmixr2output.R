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
#' @param popDf data.frame of parameter estimates
#' @param iniDf The ini data.frame from the rxUi object
#' @return `popDf` with a "Parameter" column added if labels are present
#' @author Matthew L. Fidler
#' @noRd
.updateParFixedAddParameterLabel <- function(popDf, iniDf) {
  .lab <- paste(iniDf$label[!is.na(iniDf$ntheta)])
  .lab[.lab == "NA"] <- ""
  .lab <- gsub(" *$", "", gsub("^ *", "", .lab))
  if (!all(.lab == "")) {
    popDf <- data.frame(Parameter = .lab, popDf, check.names = FALSE)
  }
  popDf
}
#' Apply the manual translation for only one item
#'
#' @param theta The current parameter name as a character string
#' @inheritParams .updateParFixedApplyManualBacktransformations
#' @returns popDf with back-transform and CIs updated for `theta`
#' @keywords internal
#' @author Matthew L. Fidler
#' @noRd
.updateParFixedApplyManualBacktransformationsI <- function(theta, popDf, btFun, ci, btEnv) {
  .qn <- qnorm(1.0-(1-ci)/2)
  .bfun <- NULL
  if (!is.na(btFun)) {
    .bfun <- try(get(btFun, envir = btEnv, mode="function"), silent=TRUE)
    if (inherits(.bfun, "try-error")) {
      warning("unknown function '", btFun, "' for manual backtransform, revert to nlmixr2 back-transformation detection for, '", theta, "'",
              call.=FALSE)
      .bfun <- NULL
    }
  }
  if (!is.null(.bfun)) {
    .est <- popDf[theta, "Estimate"]
    .se <- popDf[theta, "SE"]
    popDf[theta, "Back-transformed"] <- .bfun(.est)
    if (!is.na(.se)) {
      popDf[theta, "CI Lower"] <- .bfun(.est - .se * .qn)
      popDf[theta, "CI Upper"] <- .bfun(.est + .se * .qn)
    }
  }
  popDf
}

#'  This applies the manually specified back-transformations
#'
#' @param popDf data.frame of parameter estimates
#' @param iniDf The ini data.frame from the rxUi object
#' @param ci The confidence-interval level
#' @param env The environment to find back-transform functions
#' @returns popDf with back-transform and CIs updated
#' @author Matthew L. Fidler
#' @noRd
.updateParFixedApplyManualBacktransformations <- function(popDf, iniDf, ci, btEnv) {
  if (nrow(popDf) == 0) {
    return(popDf)
  }
  btFun =
    stats::setNames(
      .ret$iniDf$backTransform[!is.na(iniDf$ntheta)],
      .ret$iniDf$name[!is.na(iniDf$ntheta)]
    )
  for (currentTheta in rownames(popDf)) {
    if (currentTheta %in% names(btFun)) {
      popDf <-
        .updateParFixedApplyManualBacktransformationsI(
          theta = currentTheta,
          popDf = popDf,
          btFun = btFun[[currentTheta]],
          ci = ci,
          btEnv = btEnv
        )
    }
  }
  popDf
}

#' This gets the CV/SD for a single ETA
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

#'  This will add the between subject variability to the mu-referenced theta.
#'  It also expands the table to include non-mu referenced ETAs
#'
#' @param .ret The focei return environment
#' @param .ui The rxode2 ui environment
#' @returns Nothing called for side effects on popDf in the .ret environment
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
  .cvp <- lapply(row.names(.ret$popDf), function(x) {
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
    .ret$popDf2 <- as.data.frame(lapply(names(.ret$popDf), function(x) { rep(NA_real_, length(.nonMuRef))}))
    names(.ret$popDf2) <- names(.ret$popDf)
    row.names(.ret$popDf2) <- .nonMuRef
  }
  .ret$popDf <- data.frame(.ret$popDf, "BSD" = .cvp$v, check.names = FALSE)
  if (length(.nonMuRef) > 0) {
    .cvp <- lapply(row.names(.ret$popDf2), function(x) {
      .updateParFixedGetEtaRow(x, .env, .ome, .omegaFix, .muRefCurEval, .sigdig)
    })
    .cvp <- do.call("rbind", .cvp)
    .ret$popDf2 <- data.frame(.ret$popDf2, "BSD" = .cvp$v, check.names = FALSE)
    .ret$popDf <- rbind(.ret$popDf, .ret$popDf2)
    rm(list=c("popDf2", "popDfSig2"), envir=.ret) # popDfSig2 may no longer exist
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
  .sh <- lapply(row.names(.ret$popDf), function(x) {
    .w <- which(.muRefDataFrame$theta == x)
    if (length(.w) != 1) {
      .w <- which(names(.shrink) == x)
      if (length(.w) != 1) {
        if (any(x == .errs)) {
          x <- "IWRES"
        } else {
          return(data.frame(v = NA_real_))
        }
      }
      .eta <- x
    } else {
      .eta <- .muRefDataFrame$eta[.w]
    }
    .v <- .shrink[7, .eta]
    if (length(.v) != 1) {
      return(data.frame(v = NA_real_))
    }
    data.frame(v = .v)
  })
  .sh <- do.call("rbind", .sh)
  .ret$popDf <- data.frame(.ret$popDf, "Shrink(SD)%" = .sh$v, check.names = FALSE)
}

#' Create the $popDfSig data.frame with all formatting
#'
#' @param df The $popDf data.frame
#' @param digits The number of significant digits
#' @param ci The confidence interval (for the column name of the back-transformed column)
#' @param fixedNames Character vector of parameters that are fixed
#' @returns `df` with formatting applied
#' @noRd
.updateParFixedApplySig <- function(df, digits, ci, fixedNames) {
  if (is.null(digits)) {
    # The FO method does not have a `control` element (see to
    # https://github.com/nlmixr2/nlmixr2est/pull/509#issuecomment-2802688590)
    digits <- 3
  }
  if (is.null(ci)) {
    # The FO method does not have a `control` element (see to
    # https://github.com/nlmixr2/nlmixr2est/pull/509#issuecomment-2802688590)
    ci <- 0.95
  }
  ret <- df
  colNumEst <- which(names(ret) %in% "Estimate")
  names(ret)[colNumEst] <- "Est."
  for (nm in names(ret)) {
    if (!is.character(ret[[nm]])) {
      ret[[nm]] <- formatMinWidth(ret[[nm]], digits = digits, naValue = "")
    }
  }
  if ("CI Lower" %in% names(ret)) {
    colNumBt <- which(startsWith(names(ret), "Back"))
    ret[[colNumBt]] <-
      ifelse(
        ret$`CI Lower` == "",
        ret[[colNumBt]],
        sprintf("%s (%s, %s)", ret[[colNumBt]], ret$`CI Lower`, ret$`CI Upper`)
      )
    names(ret)[colNumBt] <- sprintf("Back-transformed(%g%%CI)", 100*ci)
    ret$`CI Lower` <- NULL
    ret$`CI Upper` <- NULL
  }
  # Add the suffix to the shrinkage
  if ("Shrink(SD)%" %in% names(df)) {
    # Only add the suffix if they are not NA
    shrinkMask <- !is.na(df$`Shrink(SD)%`)
    if (any(shrinkMask)) {
      shrinkSuffix <-
        as.character(cut(
          df$`Shrink(SD)%`,
          breaks = c(-Inf, 0, 20, 30, Inf),
          labels = c(">", "<", "=", ">")
        ))
      ret$`Shrink(SD)%`[shrinkMask] <- paste0(ret$`Shrink(SD)%`[shrinkMask], shrinkSuffix[shrinkMask])
    }
  }
  # Add SE and RSE FIXED
  if (length(fixedNames) > 0) {
    ret[fixedNames, "SE"] <- "FIXED"
    ret[fixedNames, "%RSE"] <- "FIXED"
  }
  ret
}

#' Create the parFixed dataset
#'
#' @param .ret focei style environment
#' @return Nothing, called for side effects
#' @author Matthew L. Fidler
#' @noRd
.updateParFixed <- function(.ret) {
  .ui <- .ret$ui
  .fixedNames <- character()
  if (!is.null(nlmixr2global$nlmixr2EstEnv$uiUnfix)) {
    .ui <- nlmixr2global$nlmixr2EstEnv$uiUnfix
    .theta <- .ui$theta
    .tn <- names(.theta)

    .popDfEst <- .ret$popDf
    .popDfEst$Estimate <- unname(.popDfEst$Estimate)
    .popDfEst$SE <- unname(.popDfEst$SE)
    .fixedNames <- setdiff(names(.theta), rownames(.popDfEst))
    .popDfFixed <-
      data.frame(
        Estimate = unname(.theta[.fixedNames]),
        SE = NA_real_,
        `%RSE` = NA_real_,
        `Back-transformed` = unname(.theta[.fixedNames]),
        `CI Lower` = NA_real_,
        `CI Upper` = NA_real_,
        row.names = .fixedNames,
        check.names = FALSE,
        check.rows = FALSE
      )
    # Combine estimated and fixed parameters, then order them by the theta names
    .popDf <- rbind(.popDfEst, .popDfFixed)[.tn, ]

    # Show the fixed values in the model
    .ret$popDf <- .popDf
  }
  popDf <- .ret$popDf
  popDf <-
    .updateParFixedApplyManualBacktransformations(
      popDf = popDf,
      iniDf = .ui$iniDf,
      ci = .ret$control$ci,
      btEnv = nlmixr2global$nlmixrEvalEnv$envir
    )
  popDf <- .updateParFixedAddParameterLabel(popDf, iniDf = .ui$iniDf)
  .ret$popDf <- popDf
  .updateParFixedAddBsv(.ret, .ui)
  .updateParFixedAddShrinkage(.ret, .ui)
  # Applying significant digits happens via .updateParFixedApplySig
  # (.ret$popDfSig is ignored)
  .ret$parFixed <-
    .updateParFixedApplySig(
      .ret$popDf,
      digits = .ret$control$sigdig,
      ci = .ret$control$ci,
      fixedNames = .fixedNames
    )
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
  nmObjGet(.lst)
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
  .ret
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
    .ret
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
  object[, toupper(match.arg(type))]
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
#' @noRd
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

#' Format numeric values to minimize the printing width
#'
#' Special values (`NaN`, `Inf`, `-Inf`, and `0`) are returned as their
#' character representation without additional modification. `NA` is returned as
#' the value of the `naValue` argument.
#'
#' @param x The numeric vector to convert
#' @param digits The number of significant digits to show
#' @param naValue The value to return if `is.na(x)`
#' @returns A character vector converting the numbers with minimum width
#' @examples
#' formatMinWidth(x = -123456*10^(-10:10))
#' @export
formatMinWidth <- function(x, digits = 3, naValue = "NA") {
  checkmate::assert_numeric(x)
  maskSpecialValue <- x %in% c(NA, NaN, Inf, -Inf, 0)
  signifX <- signif(x[!maskSpecialValue], digits = digits)
  # Generate scientific notation values
  formatSN <- paste0("%1.", digits - 1, "e")
  valueSN <- sprintf(formatSN, signifX)
  # Drop the + and leading zeros for positive exponent SN values
  valueSN <- gsub(x = valueSN, pattern = "e\\+0*", replacement = "e")
  # Drop the leading zeros for negative exponent SN values
  valueSN <- gsub(x = valueSN, pattern = "e\\-0*", replacement = "e-")

  firstDigit <- floor(log10(abs(signifX)))
  lastDigit <- firstDigit - digits + 1
  formatNormal <-
    paste0(
      "%",
      ifelse(
        firstDigit < 0,
        "0",
        firstDigit
      ),
      ifelse(
        lastDigit < 0,
        paste0(".", abs(lastDigit)),
        ".0"
      ),
      "f"
    )
  valueNormal <- sprintf(formatNormal, signifX)

  # Prepare the return value with special values
  ret <- rep(NA_character_, length(x))
  if (any(maskSpecialValue)) {
    ret[maskSpecialValue] <- as.character(x[maskSpecialValue])
    maskNA <- !is.nan(x) & is.na(x)
    if (any(maskNA)) {
      ret[maskNA] <- naValue
    }
  }

  # Place in the normal values
  ret[!maskSpecialValue] <-
    ifelse(
      nchar(valueSN) < nchar(valueNormal),
      valueSN,
      valueNormal
    )
  ret
}
