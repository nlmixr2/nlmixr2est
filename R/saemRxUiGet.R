#' @export
rxUiGet.saemMuRefCovariateDataFrame <- function(x, ...) {
  .ui <- x[[1]]
  if (!rxode2::rxGetControl(x[[1]], "muRefCov", getOption("nlmixr2.saemMuRefCov", TRUE))) {
    return(data.frame(theta=character(0), covariate=character(0), covariateParameter=character(0)))
  }
  if (exists("muRefFinal", .ui)) {
    .cov <- .ui$muRefFinal
  } else {
    .cov <- .ui$muRefCovariateDataFrame
  }
  .iniDf <- .ui$iniDf
  .rm <- NULL
  for (.i in seq_along(.cov$covariateParameter)) {
    .cp <- .cov$covariateParameter[.i]
    .cv <- .cov$covariate[.i]
    .w <- which(.iniDf$name == .cp)
    if (.iniDf$fix[.w]) {
      .rm <- c(.rm, -.i)
    }
  }
  if (!is.null(.rm)) {
    .cov <- .cov[.rm,]
  }
  .cov
}
attr(rxUiGet.saemMuRefCovariateDataFrame, "rstudio") <- NA


#' @export
rxUiGet.saemInParsAndMuRefCovariates <- function(x, ...) {
  .ui <- x[[1]]
  # mu ref final removes time varying covariates
  .muRefFinal <- rxUiGet.saemMuRefCovariateDataFrame(x, ...)
  .cov <- .ui$covariates
  .muCov <- unique(.muRefFinal$covariate)
  .cov <- .cov[!(.cov %in% .muCov)]
  if (length(.ui$predDf$cond) > 1) {
    .cov <- unique(c("CMT", .cov))
  }
  list(inPars=.cov, covars=.muCov)
}
#attr(rxUiGet.saemInParsAndMuRefCovariates, "desc") <- "Get inPars and covars for saem"

#' @export
rxUiGet.saemInPars <- function(x, ...) {
  .ret <- rxUiGet.saemInParsAndMuRefCovariates(x, ...)
  .ret$inPars
}
#attr(rxUiGet.saemInPars, "desc") <- "get inPars"
attr(rxUiGet.saemInPars, "rstudio") <- "char"

#' @export
rxUiGet.saemCovars <- function(x, ...) {
  .ret <- rxUiGet.saemInParsAndMuRefCovariates(x, ...)
  .ret$covars
}
#attr(rxUiGet.saemInPars, "desc") <- "get saemn mu-referenced non-time varying covariates"
attr(rxUiGet.saemCovars, "rstudio") <- "char"

#' @export
rxUiGet.saemFunctionModPredQuote <- function(x, ...) {
  .mod <- rxode2::rxode2(rxUiGet.saemModel(x, ...))
  .fnPred <- bquote(function(a, b, c) {
    rxode2::rxLoad(.(.mod))
    rxode2::rxLock(.(.mod))
    rxode2::rxAllowUnload(FALSE)
    on.exit({
      rxode2::rxUnlock(.(.mod))
      rxode2::rxAllowUnload(TRUE)
      rxode2::rxSolveFree()
    })
    .Call(`_nlmixr2est_saem_do_pred`, a, b, c)
  })
  list(.mod, .fnPred)
}


#' @export
rxUiGet.saemFunction <- function(x, ...) {
  # This function depends on the number of time varying covariates in the data
  .ui <- x[[1]]
  .mod <- rxUiGet.saemFunctionModPredQuote(x, ...)
  .fnPred <- .mod[[2]]
  .mod <- .mod[[1]]
  .fn <- bquote(function(a, b, c) {
    rxode2::rxLoad(.(.mod))
    rxode2::rxLock(.(.mod))
    on.exit({
      rxode2::rxUnlock(.(.mod))
      rxode2::rxAllowUnload(TRUE)
      rxode2::rxSolveFree()
    })
    if (missing(b) && missing(c)) {
      .ret <- .Call(`_nlmixr2est_saem_fit`, a, PACKAGE = "nlmixr2est")
      attr(.ret, "dopred") <- .(.fnPred)
      return(.ret)
    } else {
      .curFn <- .(.fnPred)
      return(.curFn(a, b, c))
    }
  })
  .inPars <- rxUiGet.saemInPars(x, ...)
  .param <- rxode2::rxParam(.mod)
  .estParam <- rxUiGet.saemParamsToEstimate(x, ...)
  .cov <- rxUiGet.saemMuRefCovariateDataFrame(x, ...)
  if (length(.cov$covariateParameter) > 0) {
    .param <- .param[!(.param %in% .cov$covariateParameter)]
    .estParam <- .estParam[!(.estParam %in% .cov$covariateParameter)]
  }
  .parmUpdate <- vapply(.param, function(x) {
    if (x %in% .estParam) {
      return(1L)
    } else {
      return(0L)
    }
  }, integer(1), USE.NAMES=FALSE)
  .nendpnt <- length(.ui$predDf$cond)
  .fn <- eval(.fn)
  attr(.fn, "form") <- "ode" ## Not sure this is necessary any more
  attr(.fn, "neq") <- length(rxode2::rxState(.mod))
  attr(.fn, "nlhs") <- length(rxode2::rxLhs(.mod))
  attr(.fn, "nrhs") <- sum(.parmUpdate)
  attr(.fn, "paramUpdate") <- .parmUpdate
  attr(.fn, "rx") <- .mod
  attr(.fn, "inPars") <- .inPars
  attr(.fn, "nendpnt") <- .nendpnt
  .fn
}
attr(rxUiGet.saemFunction, "rstudio") <- function(){}

#' @export
rxUiGet.saemFixed <- function(x, ...) {
  .ui <- x[[1]]
  .df <- .ui$iniDf
  .dft <- .df[!is.na(.df$ntheta), ]
  .fixError <- .dft[!is.na(.dft$err), ]
  .dft <- .dft[is.na(.dft$err), ]
  .dft <- setNames(.dft$fix, paste(.dft$name))
  .cov <- rxUiGet.saemMuRefCovariateDataFrame(x, ...)
  if (length(.cov$theta) > 0) {
    .theta <- .dft
    .theta <- .theta[!(names(.theta) %in% .cov$covariateParameter)]
    .allCovs <- rxUiGet.saemCovars(x, ...)
    .lc <- length(.allCovs)
    .m <- matrix(rep(NA_character_, .lc * length(.theta)), ncol = .lc)
    dimnames(.m) <- list(names(.theta), .allCovs)
    for (.c in seq_along(.cov$covariateParameter)) {
      .curTheta <- .cov[.c, "theta"]
      .curCov <- .cov[.c, "covariate"]
      .curPar <- .cov[.c, "covariateParameter"]
      .m[.curTheta, .curCov] <- .curPar
    }
    .m <- cbind(matrix(names(.theta), ncol=1), .m)
    .m <- as.vector(t(.m))
    .dft <- .dft[.m[!is.na(.m)]]
  }
  .extra <- .ui$nonMuEtas
  .extra <- setNames(rep(TRUE, length(.extra)), .ui$nonMuEtas)
  c(.dft, .extra)
}
#attr(rxUiGet.saemFixed, "desc") <- "Get the saem fixed parameters"
attr(rxUiGet.saemFixed, "rstudio") <- function(){}

.saemEtaTrans <- function(x, ..., nonMu=FALSE) {
  .ui <- x[[1]]
  .etas <- .ui$iniDf[!is.na(.ui$iniDf$neta1), ]
  .etas <- .etas$name[.etas$neta1 == .etas$neta2]
  .thetas <- rxUiGet.saemParamsToEstimateCov(x, ...)
  if (nonMu) {
    .thetas <- .thetas[!(.thetas %in% .ui$nonMuEtas)]
  }
  .muRefDataFrame <- .ui$muRefDataFrame
  vapply(.etas, function(eta) {
    .w <- which(eta == .muRefDataFrame$eta)
    if (length(.w) == 1L) {
      .muTheta <- .muRefDataFrame$theta[.w]
      .w <- which(.muTheta == .thetas)
      if (length(.w) == 1L) return(.w)
    }
    if (nonMu && eta %in% .ui$nonMuEtas) {
      .w <- which(eta == .etas)
      if (length(.w) == 1L) return(-.w)
    }
    .w <- which(eta == .thetas)
    if (length(.w) == 1L) return(.w)
    return(NA_integer_)
  }, integer(1), USE.NAMES=FALSE)
}

#' @export
rxUiGet.saemEtaTrans <- function(x, ...) {
  .saemEtaTrans(x, ...)
}
attr(rxUiGet.saemEtaTrans, "rstudio") <- c(1L, 3L)

#' @export
rxUiGet.saemEtaTransPred <- function(x, ...) {
  .saemEtaTrans(x, ..., nonMu=TRUE)
}
#attr(rxUiGet.saemEtaTrans, "desc") <- "Get the saem eta to theta translation"
attr(rxUiGet.saemEtaTransPred, "rstudio") <- c(1L, 3L)

#' @export
rxUiGet.saemOmegaTrans <- function(x, ...) {
  .etaTrans <- rxUiGet.saemEtaTrans(x, ...)
  .o <- order(.etaTrans)
  .etaTrans2 <- .etaTrans
  .c <- 1
  for (i in .o) {
    .etaTrans2[i] <- .c
    .c <- .c + 1
  }
  .etaTrans2
}
#attr(rxUiGet.saemOmegaTrans, "desc") <- "Get the saem omega to UI omega translation"
attr(rxUiGet.saemOmegaTrans, "rstudio") <- c(1L, 3L)


#' @export
rxUiGet.saemModelOmega <- function(x, ...) {
  .ui <- x[[1]]
  .thetas <- rxUiGet.saemParamsToEstimateCov(x, ...)
  .etaTrans <- rxUiGet.saemEtaTrans(x, ...)
  .dm <- length(.thetas)
  .mat <- matrix(rep(0, .dm * .dm), .dm)
  .iniDf <- .ui$iniDf
  .etd <- .iniDf[which(!is.na(.iniDf$neta1)), ]
  for (i in seq_along(.etd$neta1)) {
    .mat[.etaTrans[.etd$neta1[i]], .etaTrans[.etd$neta2[i]]] <-
      .mat[.etaTrans[.etd$neta2[i]], .etaTrans[.etd$neta1[i]]] <- 1
  }
  .mat
}
#attr(rxUiGet.saemModelOmega, "desc") <- "Get the saem model omega"
attr(rxUiGet.saemModelOmega, "rstudio") <- lotri::lotri(a+b~c(1, 0.1, 1))

#' @export
rxUiGet.saemModelOmegaFixed <- function(x, ...) {
  .ui <- x[[1]]
  .thetas <- rxUiGet.saemParamsToEstimateCov(x, ...)
  .etaTrans <- rxUiGet.saemEtaTrans(x, ...)
  .dm <- length(.thetas)
  .mat <- matrix(rep(0, .dm * .dm), .dm)
  .iniDf <- .ui$iniDf
  .etd <- .iniDf[which(!is.na(.iniDf$neta1)), ]
  for (i in seq_along(.etd$neta1)) {
    .mat[.etaTrans[.etd$neta1[i]], .etaTrans[.etd$neta2[i]]] <-
      .mat[.etaTrans[.etd$neta2[i]], .etaTrans[.etd$neta1[i]]] <- as.integer(.etd$fix[i])
  }
  .mat
}
#attr(rxUiGet.saemModelOmegaFixed, "desc") <- "Get the indicator for saem model omega fixed components"
attr(rxUiGet.saemModelOmegaFixed, "rstudio") <- lotri::lotri(a+b~c(1, 0.1, 1))

#' @export
rxUiGet.saemModelOmegaFixedValues <- function(x, ...) {
  .ui <- x[[1]]
  .thetas <- rxUiGet.saemParamsToEstimateCov(x, ...)
  .etaTrans <- rxUiGet.saemEtaTrans(x, ...)
  .dm <- length(.thetas)
  .mat <- matrix(rep(0, .dm * .dm), .dm)
  .iniDf <- .ui$iniDf
  .etd <- .iniDf[which(!is.na(.iniDf$neta1)), ]
  for (i in seq_along(.etd$neta1)) {
    .mat[.etaTrans[.etd$neta1[i]], .etaTrans[.etd$neta2[i]]] <-
      .mat[.etaTrans[.etd$neta2[i]], .etaTrans[.etd$neta1[i]]] <- .etd$est[i]
  }
  .mat
}
#attr(rxUiGet.saemModelOmegaFixedValues, "desc") <- "Get the omega values may be fixed"
attr(rxUiGet.saemModelOmegaFixedValues, "rstudio") <- lotri::lotri(a+b~c(1, 0.1, 1))

#' @export
rxUiGet.saemLow <- function(x, ...) {
  .ui <- x[[1]]
  .ui$predDf$trLow
}
#attr(rxUiGet.saemLow, "desc") <- "Get the saem error transformation lower boundary"
attr(rxUiGet.saemLow, "rstudio") <- -Inf

#' @export
rxUiGet.saemHi <- function(x, ...) {
  .ui <- x[[1]]
  .ui$predDf$trHi
}
#attr(rxUiGet.saemHi, "desc") <- "Get the saem error transformation higher boundary"
attr(rxUiGet.saemHi, "rstudio") <- Inf

#' @export
rxUiGet.saemPropT <- function(x, ...) {
  .ui <- x[[1]]
  as.integer((.ui$predDf$errTypeF=="transformed")*1L)
}
#attr(rxUiGet.saemPropT, "desc") <- "Get the saem transformation type for the function"
attr(rxUiGet.saemPropT, "rstudio") <- 1L

#' @export
rxUiGet.saemYj <- function(x, ...) {
  .ui <- x[[1]]
  as.integer(.ui$predDf$transform) - 1
}
#attr(rxUiGet.saemYj, "desc") <- "Get the saem transformation type"
attr(rxUiGet.saemYj, "rstudio") <- 1L

#' @export
rxUiGet.saemResMod <- function(x, ...) {
  .ui <- x[[1]]
  .predDf <- .ui$predDf
  vapply(seq_along(.predDf$errType),
         function(i) {
           .errType <- as.integer(.predDf$errType[i])
           .hasLambda <- .predDf$transform[i] %in% c("boxCox", "yeoJohnson",
                                                     "logit + yeoJohnson",
                                                     "probit + yeoJohnson",
                                                     "logit + boxCox",
                                                     "probit + boxCox")
           if (.hasLambda) {
             return(.errType + 5L)
           } else {
             return(.errType)
           }
         }, integer(1), USE.NAMES=FALSE)
}
#attr(rxUiGet.saemResMod, "desc") <- "saem res.mod component"
attr(rxUiGet.saemResMod, "rstudio") <- c(1L, 2L)

#' @export
rxUiGet.saemModNumEst <- function(x, ...) {
  .resMod <- rxUiGet.saemResMod(x, ...)
  vapply(.resMod, function(i) {
    switch(i,
           1L, # add = 1
           1L, # prop = 2
           2L, # pow = 3
           2L, # add + prop = 4
           3L, # add + pow = 5
           2L, # add + lambda = 6
           2L, # prop + lambda = 7
           3L, # pow + lambda = 8
           3L, # add + prop + lambda = 9
           4L # add + pow + lambda = 10
           )
  }, integer(1), USE.NAMES=TRUE)
}
#attr(rxUiGet.saemModNumEst, "desc") <- "saem number of parameters that can be estimated for each component"
attr(rxUiGet.saemModNumEst, "rstudio") <- c(1L, 2L)

#' @export
rxUiGet.saemModResOffset <- function(x, ...) { # res_offset
  cumsum(c(0, rxUiGet.saemModNumEst(x, ...)))
}
#attr(rxUiGet.saemModResOffset, "desc") <- "saem residual parameters offset"
attr(rxUiGet.saemModResOffset, "rstudio") <- c(1, 2)

#' @export
rxUiGet.saemModResTotalResiduals <- function(x, ...) { # res_offset
  sum(rxUiGet.saemModNumEst(x, ...))
}
#attr(rxUiGet.saemModResTotalResiduals, "desc") <- "saem total number of residuals"
attr(rxUiGet.saemModResTotalResiduals, "rstudio") <- c(1, 2)

#' Get the ares name based on the condition
#'
#' @param iniDf Ini data frame
#' @param cond Condition of the parameter
#' @param types The types of errors that support this distribution
#' @return Name of the ares for saem
#' @author Matthew L. Fidler
#' @noRd
.saemGetIniDfResNameFromType <- function(iniDf, cond, types, column="name") {
  .ini <- iniDf[iniDf$condition == cond, ]
  .ini[.ini$err %in% types, column]
}

.saemGetIniDfAResName <- function(iniDf, cond, column="name") {
  .saemGetIniDfResNameFromType(iniDf, cond, c("add", "lnorm", "probitNorm", "logitNorm"), column=column)
}

.saemGetIniDfBResName <- function(iniDf, cond, column="name") {
  .saemGetIniDfResNameFromType(iniDf, cond, c("prop", "propT", "propF", "pow", "powT", "powF"), column=column)
}

.saemGetIniDfCResName <- function(iniDf, cond, column="name") {
  .saemGetIniDfResNameFromType(iniDf, cond, c("pow2", "powF2", "powT2"), column=column)
}

.saemGetIniDfLResName <- function(iniDf, cond, column="name") {
  .saemGetIniDfResNameFromType(iniDf, cond, c("boxCox", "yeoJohnson"), column=column)
}
#'  Get saem residual item
#' @param ui rxode2 UI
#' @param column column name
#' @return Residual names
#' @author Matthew L. Fidler
#' @noRd
.saemGetResItem <- function(ui, column="name") {
  .predDf <- ui$predDf
  .iniDf <- ui$iniDf
  .numEst <- ui$saemModNumEst
  .resMod <- ui$saemResMod
  do.call("c", lapply(seq_along(.numEst),
         function(i) {
           .num <- .numEst[i]
           .cond <- .predDf$cond[i]
           .ret <- switch(.resMod[i],
                          .saemGetIniDfAResName(.iniDf, .cond, column), # add = 1
                          .saemGetIniDfBResName(.iniDf, .cond, column), # prop = 2
                          c(.saemGetIniDfBResName(.iniDf, .cond, column),
                            .saemGetIniDfCResName(.iniDf, .cond, column)), # pow = 3

                          c(.saemGetIniDfAResName(.iniDf, .cond, column),
                            .saemGetIniDfBResName(.iniDf, .cond, column)), # add + prop = 4

                          c(.saemGetIniDfAResName(.iniDf, .cond, column),
                            .saemGetIniDfBResName(.iniDf, .cond, column),
                            .saemGetIniDfCResName(.iniDf, .cond, column)), # add + pow = 5

                          c(.saemGetIniDfAResName(.iniDf, .cond, column),
                            .saemGetIniDfLResName(.iniDf, .cond, column)), # add + lambda = 6

                          c(.saemGetIniDfBResName(.iniDf, .cond, column),
                            .saemGetIniDfLResName(.iniDf, .cond, column)), # prop + lambda = 7

                          c(.saemGetIniDfBResName(.iniDf, .cond, column),
                            .saemGetIniDfCResName(.iniDf, .cond, column),
                            .saemGetIniDfLResName(.iniDf, .cond, column)), # pow + lambda = 8

                          c(.saemGetIniDfAResName(.iniDf, .cond, column),
                            .saemGetIniDfBResName(.iniDf, .cond, column),
                            .saemGetIniDfLResName(.iniDf, .cond, column)), # add + prop + lambda = 9

                          c(.saemGetIniDfAResName(.iniDf, .cond, column),
                            .saemGetIniDfBResName(.iniDf, .cond, column),
                            .saemGetIniDfCResName(.iniDf, .cond, column),
                            .saemGetIniDfLResName(.iniDf, .cond, column)) # add + pow + lambda = 10
                          )
           if (.num != length(.ret)) {
             stop("endpoint '", .cond, "' for saem cannot locate the residual error(s) correctly", call.=FALSE)
           }
           .ret
         }))
}

#' @export
rxUiGet.saemResNames <- function(x, ...) {
  .ui <- x[[1]]
  .saemGetResItem(.ui, column="name")
}
#attr(rxUiGet.saemResNames, "desc") <- "Get error names for SAEM"
attr(rxUiGet.saemResNames, "rstudio") <- "add.err"

#' @export
rxUiGet.saemResFixed <- function(x, ...) {
  .ui <- x[[1]]
  as.integer(.saemGetResItem(.ui, column="fix"))
}
#attr(rxUiGet.saemResFixed, "desc") <- "Integer vector of residual fixed components"
attr(rxUiGet.saemResFixed, "rstudio") <- c(1L, 2L)

#' @export
rxUiGet.saemParHistResNames <- function(x, ...) {
  .fix <- (rxUiGet.saemResFixed(x, ...) == 0L)
  rxUiGet.saemResNames(x, ...)[.fix]
}
#attr(rxUiGet.saemParHistResNames, "desc") <-"Get the SAEM parameter history residual names"
attr(rxUiGet.saemParHistResNames, "rstudio") <- "add.sd"

#' @export
rxUiGet.saemResValue <- function(x, ...) {
  .ui <- x[[1]]
  .saemGetResItem(.ui, column="est")
}
attr(rxUiGet.saemResValue, "rstudio") <- 0.7

#' @export
rxUiGet.saemEtaNames <- function(x, ...) {
  .ui <- x[[1]]
  .etaNames <- .ui$iniDf[!is.na(.ui$iniDf$neta1), ]
  .etaNames <- .etaNames[.etaNames$neta1 == .etaNames$neta2, "name"]
  ## .etaTrans <- rxUiGet.saemOmegaTrans(x, ...)
  .etaTrans <- rxUiGet.saemEtaTrans(x, ...)
  .names <- rxUiGet.saemParamsToEstimateCov(x, ...)
  .names <- rep("", length(.names))
  for (.i in seq_along(.etaTrans)) {
    .names[.etaTrans[.i]] <- .etaNames[.i]
  }
  .names <- .names[.names != ""]
  .names
}
#attr(rxUiGet.saemParHistEtaNames, "desc") <- "Get ETA names for SAEM based on theta order"
attr(rxUiGet.saemEtaNames,"rstudio") <- "eta.names"

#' @export
rxUiGet.saemParHistOmegaKeep <- function(x, ...) {
  .ui <- x[[1]]
  .etaNames <- .ui$iniDf[!is.na(.ui$iniDf$neta1), ]
  .etaNames <- .etaNames[.etaNames$neta1 == .etaNames$neta2,]
  .names <- rxUiGet.saemEtaNames(x, ...)
  vapply(.names, function(etaName) {
    .w <- which(.etaNames$name == etaName)
    if (length(.w) == 1) {
      return(1L - as.integer(.etaNames$fix[.w]))
    } else {
      stop("cannot figure out saemParHistOmegaKeep", call.=FALSE)
    }
  }, integer(1))
}
#attr(rxUiGet.saemOmegaKeep, "desc") <- "Get the etas that are kept for SAEM based on theta order"
attr(rxUiGet.saemParHistOmegaKeep, "rstudio") <- c("eta.ka"=1)

#' @export
rxUiGet.saemParHistEtaNames <- function(x, ...) {
  .ui <- x[[1]]
  .names <- rxUiGet.saemParHistOmegaKeep(x, ...)
  .names <- .names[.names == 1L]
  if (length(.names) == 0) return(NULL)
  paste0("V(", names(.names), ")")
}
#attr(rxUiGet.saemParHistEtaNames, "desc") <- "Get the parameter history eta names"
attr(rxUiGet.saemParHistEtaNames, "rstudio") <- "V(ka)"

#' @export
rxUiGet.saemParHistNames <- function(x, ...) {
  #join_cols(join_cols(Plambda, Gamma2_phi1.diag()), vcsig2).t();
  .plambda <- rxUiGet.saemParamsToEstimate(x, ...)
  .plambda <- .plambda[!rxUiGet.saemFixed(x, ...)]
  c(.plambda, rxUiGet.saemParHistEtaNames(x, ...), rxUiGet.saemParHistResNames(x, ...))
}
attr(rxUiGet.saemParHistNames, "rstudio") <- c("ka", "add.sd")

#' @export
rxUiGet.saemAres <- function(x, ...) {
  .ui <- x[[1]]
  .predDf <- .ui$predDf
  .ini <- .ui$iniDf
  .ini <- .ini[!is.na(.ini$err), ]
  return(vapply(.predDf$cond, function(x) {
    .tmp <- .ini[which(.ini$condition == x), ]
    .w <- which(vapply(.tmp$err, function(x) {
      x %in% c(
        "add", "norm", "dnorm", "dpois",
        "pois", "dbinom", "binom", "dbern", "bern",
        "lnorm", "dlnorm", "logn", "dlogn")
    }, logical(1), USE.NAMES=FALSE))
    if (length(.w) == 1) {
      return(.tmp$est[.w])
    } else {
      return(10)
    }
  }, numeric(1), USE.NAMES=FALSE))
}
#attr(rxUiGet.saemAres, "desc") <- "ares initial estimates for saem"
attr(rxUiGet.saemAres, "rstudio") <- 0.7

#' @export
rxUiGet.saemBres <- function(x, ...) {
  .ui <- x[[1]]
  .predDf <- .ui$predDf
  .ini <- .ui$iniDf
  .ini <- .ini[!is.na(.ini$err), ]
  return(vapply(.predDf$cond, function(x) {
    .tmp <- .ini[which(.ini$condition == x), ]
    .w <- which(vapply(.tmp$err,
                       function(x) (any(x == "prop") || any(x == "propT")),
                       logical(1), USE.NAMES=FALSE))
    if (length(.w) == 1) {
      return(.tmp$est[.w])
    } else {
      .w <- which(vapply(.tmp$err,
                         function(x) (any(x == "pow") || any(x == "powT")),
                         logical(1), USE.NAMES=FALSE))
      if (length(.w) == 1) {
        return(.tmp$est[.w])
      } else {
        return(1)
      }
    }
  }, numeric(1), USE.NAMES=FALSE))
}
#attr(rxUiGet.saemBres, "desc") <- "bres initial estimates for saem"
attr(rxUiGet.saemBres, "rstudio") <- 0.7

#' @export
rxUiGet.saemCres <- function(x, ...) {
  .ui <- x[[1]]
  .predDf <- .ui$predDf
  .ini <- .ui$iniDf
  .ini <- .ini[!is.na(.ini$err), ]
  return(vapply(.predDf$cond, function(x) {
    .tmp <- .ini[which(.ini$condition == x), ]
    .w <- which(vapply(.tmp$err, function(x) (any(x == "pow2") || any(x == "powT2")),
                       logical(1), USE.NAMES=FALSE))
    if (length(.w) == 1) {
      return(.tmp$est[.w])
    } else {
      return(1)
    }
  }, numeric(1), USE.NAMES=FALSE))
}
#attr(rxUiGet.saemCres, "desc") <- "cres initial estimates for saem"
attr(rxUiGet.saemCres, "rstudio") <- 0.7

#' @export
rxUiGet.saemLres <- function(x, ...) {
  .ui <- x[[1]]
 .predDf <- .ui$predDf
  .ini <- .ui$iniDf
  .ini <- .ini[!is.na(.ini$err), ]
  return(vapply(.predDf$cond, function(x) {
    .tmp <- .ini[which(.ini$condition == x), ]
    .boxCox <- which(.tmp$err == "boxCox")
    if (length(.boxCox) == 1L) {
      return(.tmp$est[.boxCox])
    }
    .yeoJohnson <- which(.tmp$err == "yeoJohnson")
    if (length(.yeoJohnson) == 1L) {
      return(.tmp$est[.yeoJohnson])
    }
    return(1.0)
  }, numeric(1), USE.NAMES=FALSE))
}
#attr(rxUiGet.saemLres, "desc") <- "lres (lambda) initial estimates for saem"
attr(rxUiGet.saemLres, "rstudio") <- 0.7

#' @export
rxUiGet.saemLogEta <- function(x, ...) {
  .ui <- x[[1]]
  .thetas <- rxUiGet.saemParamsToEstimate(x, ...)
  .ce <- .ui$muRefCurEval
  .cov <- rxUiGet.saemMuRefCovariateDataFrame(x, ...)
  .thetas <- .thetas[!(.thetas %in% .cov$covariateParameter)]
  vapply(.thetas, function(x) {
    .w <- which(.ce$parameter == x)
    if (length(.w) == 1L) return(.ce$curEval[.w] == "exp")
    FALSE
  }, logical(1))
}
#attr(rxUiGet.saemLogEta, "desc") <- "saem's log.eta for saem"
attr(rxUiGet.saemLogEta, "rstudio") <- c(tka=TRUE)

#' @export
rxUiGet.saemModelList <- function(x, ...) {
  .ui <- x[[1]]
  .mod <- list(saem_mod = .ui$saemFunction)
  .covars <- rxUiGet.saemCovars(x, ...)
  if (length(.covars) > 0) {
    .mod$covars <- .covars
  }
  .mod$res.mod <- rxUiGet.saemResMod(x, ...)
  .mod$log.eta <- rxUiGet.saemLogEta(x, ...)
  .mod$ares    <- rxUiGet.saemAres(x, ...)
  .mod$bres    <- rxUiGet.saemBres(x, ...)
  .mod$omega   <- rxUiGet.saemModelOmega(x, ...)
  .mod
}
#attr(rxUiGet.saemModelList "desc") <- "saem's log.eta for saem"

#' @export
rxUiGet.saemInitTheta <- function(x, ...) {
  .logEta <- rxUiGet.saemLogEta(x, ...)
  .names <- names(.logEta)
  .ui <- x[[1]]
  .iniDf <- .ui$iniDf
  .est <- setNames(.iniDf[!is.na(.iniDf$ntheta) & is.na(.iniDf$err), "est"],
                   .iniDf[!is.na(.iniDf$ntheta) & is.na(.iniDf$err), "name"])
  .cov <- rxUiGet.saemMuRefCovariateDataFrame(x, ...)
  .est <- .est[!(names(.est) %in% .cov$covariateParameter)]
  .etaNames <- .iniDf[is.na(.iniDf$ntheta), ]
  .etaNames <- .iniDf[.iniDf$neta1 == .iniDf$neta2, "name"]
  .fixed <- rxUiGet.saemFixed(x, ...)
  .theta <- .fixed
  .theta <- .theta[!(names(.theta) %in% .cov$covariateParameter)]
  .logEta <- .logEta[!(names(.logEta) %in% .cov$covariateParameter)]
  .n <- vapply(.theta, function(x) ifelse(x, "FIXED", ""),
               character(1), USE.NAMES=FALSE)
  .ret <- vapply(seq_along(.logEta),
                   function(i) {
                    .isEta <- any(.names[i] %in% .etaNames)
                    if (.logEta[i]) {
                      if (.isEta) {
                        return(1)
                      } else {
                        return(exp(.est[i]))
                      }
                    } else {
                      if (.isEta) {
                        return(0)
                      } else {
                        return(.est[i])
                      }
                    }
                  }, numeric(1), USE.NAMES=FALSE)
  if (length(.cov$theta) > 0) {
    .allCovs <- rxUiGet.saemCovars(x, ...)
    .lc <- length(.allCovs)
    .m <- matrix(rep(NA, .lc * length(.ret)), ncol = .lc)
    dimnames(.m) <- list(names(.theta), .allCovs)
    for (.c in seq_along(.cov$theta)) {
      .curTheta <- .cov[.c, "theta"]
      .curCov <- .cov[.c, "covariate"]
      .curPar <- .cov[.c, "covariateParameter"]
      .w <- which(.iniDf$name == .curPar)
      .est <- .iniDf$est[.w]
      .m[.curTheta, .curCov] <- .est
    }
    .ret <- c(.ret, as.vector(.m))
    .ret <- setNames(.ret, c(.n, rep("", length(.ret) - length(.n))))
    #.ret <- setNames(, c(.n, rep("", .lc + 1)))
  } else {
    .ret <- setNames(.ret, .n)
  }
  .ret
}
#attr(rxUiGet.saemInitTheta, "desc") <- "initialization for saem's theta"
attr(rxUiGet.saemInitTheta, "rstudio") <- c(" "=1, tcl=1)

#' @export
rxUiGet.saemInitOmega <- function(x, ...) {
  .ui <- x[[1]]
  .iniDf <- .ui$iniDf
  .eta <- .iniDf[is.na(.iniDf$ntheta), ]
  .eta <- .eta[.eta$neta1 == .eta$neta2, ]
  .eta <- setNames(.eta$est, .eta$name)
  .pars <- rxUiGet.saemParamsToEstimateCov(x, ...)
  .ret <- rep(1.0, length(.pars))
  .etaTrans <- rxUiGet.saemEtaTrans(x, ...)
  for (i in seq_along(.etaTrans)) {
    .ret[.etaTrans[i]] <- .eta[i]
  }
  .ret <- setNames(.ret, .pars)
  .cov <- rxUiGet.saemMuRefCovariateDataFrame(x, ...)
  if (length(.cov$covariateParameter) > 0) {
    .ret <- .ret[!(names(.ret) %in% .cov$covariateParameter)]
  }
  .ret
}
#attr(rxUiGet.saemInitOmega, "desc") <- "initialization for saem's omega"
attr(rxUiGet.saemInitOmega, "rstudio") <- c(tka=0.6)

#' @export
rxUiGet.saemInit <- function(x, ...) {
  list(theta=rxUiGet.saemInitTheta(x, ...),
       omega=rxUiGet.saemInitOmega(x, ...))
}
#attr(rxUiGet.saemInit, "desc") <- "initialization for saem's theta and omega"

#' @export
rxUiGet.saemThetaDataFrame <- function(x, ...) {
  .ui <- x[[1]]
  .theta <- .ui$theta
  .fixed <- .ui$iniDf[!is.na(.ui$iniDf$ntheta), "fix"]
  data.frame(lower= -Inf, theta=.theta, fixed=.fixed, upper=Inf, row.names=names(.theta))
}
#attr(rxUiGet.saemThetaDataFrame, "desc") <- "Get theta data frame"
attr(rxUiGet.saemThetaDataFrame, "rstudio") <- NA

#' @export
rxUiGet.saemParHistThetaKeep <- function(x, ...) {
  1L-as.integer(rxUiGet.saemFixed(x, ...))
}
#attr(rxUiGet.saemParHistThetaKeep, "desc") <- "The thetas that are kept in the parameter history"
attr(rxUiGet.saemParHistThetaKeep, "rstudio") <- c(1L, 2L)

#' @export
rxUiGet.saemAddProp <- function(x, ...) {
  .ui <- x[[1]]
  .parDf <- .ui$predDf
  .addProp <- as.integer(.parDf$addProp)
  .w <- which(.addProp == 3L)
  if (length(.w) > 0) {
    .default <- c(combined1=1, combined2=2)[rxode2::rxGetControl(.ui, "addProp", "combined2")]
    .addProp[.w] <- .default
  }
  .addProp
}
#attr(rxUiGet.saemParHistThetaKeep, "desc") <- "Get the saem addProp integer vector"
attr(rxUiGet.saemAddProp, "rstudio") <- 2
