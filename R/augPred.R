#' Augment Prediction for Ipred Model
#'
#' This function augments the prediction for an individual prediction
#' (Ipred) model.  It retrieves the simulation model from the fit
#' object and evaluates the model variables.
#'
#' @param fit The fitted model object from which to retrieve the simulation model.
#' @return The evaluated model variables for the Ipred model.
#' @details
#'
#' The function performs the following steps:
#'
#' - Retrieves the simulation model from the provided `fit` object using `.getSimModel` with `hideIpred` and `tad` set to `FALSE`.
#'
#' - Evaluates the model variables using `rxModelVars`.
#'
.augPredIpredModel <- function(fit) {
  .ipredModel <- .getSimModel(fit, hideIpred=FALSE,tad=FALSE)
  eval(as.call(list(quote(`rxModelVars`), .ipredModel[[-1]])))
}

#' This expands data based on predicted data for `augPred` style solves
#'
#'
#' This method will:
#'
#' - Interpolates covariate values based on the specified
#' interpolation method.
#'
#' - Combines the original data with the expanded data and ensures all
#' necessary columns are present.
#'
#' - Orders the resulting data frame by ID and TIME.
#'
#' @noRd
.augPredExpandData <- function(fit,
                               covsInterpolation = c("locf", "nocb", "linear", "midpoint"),
                               minimum = NULL, maximum = NULL, length.out = 51L) {
  .origData <- rxode2::etTrans(fit$dataSav,
                               .augPredIpredModel(fit),
                               addCmt=TRUE, keepDosingOnly=TRUE, allTimeVar=TRUE,
                               addlKeepsCov = fit$control$rxControl$addlKeepsCov,
                               addlDropSs = fit$control$rxControl$addlDropSs,
                               ssAtDoseTime = fit$control$rxControl$ssAtDoseTime)
  .predDf <- fit$ui$predDf
  .range <- range(.origData$TIME)
  .covs <- fit$ui$allCovs
  .obsData <- .origData[.origData$EVID %in% c(0, 2), ]
  .allCmt <- unique(.obsData$CMT)
  .idLvl <- attr(class(.origData), ".rxode2.lst")$idLvl
  if (is.null(minimum)) {
    minimum <- .range[1]
  }
  if (is.null(maximum)) {
    maximum <- .range[2]
  }
  .fs <- c(locf = 0, nocb = 1, midpoint = 0.5, linear = 0)
  .base <- expand.grid(TIME=seq(minimum, maximum, length.out=length.out),
                       EVID=2, AMT=NA_real_, II=NA_real_, DV=NA_real_, CMT=.allCmt)
  .covsi <- match.arg(covsInterpolation)
  .ret0 <- c(list(as.data.frame(.origData)),
                    lapply(seq_along(.idLvl), function(id) {
                      .cur <- .origData[.origData$ID == id, ]
                      if (length(.covs) > 0) {
                        cbind(data.frame(ID=id, .base),
                              setNames(data.frame(lapply(.covs, function(cov){
                                suppressWarnings({
                                  .fun <- stats::approxfun(.cur$TIME, .cur[[cov]],
                                                           method = ifelse(.covsi == "linear", "linear", "constant"),
                                                           rule = 2,
                                                           f = .fs[.covsi])
                                  .fun(.base$TIME)
                                })
                              })), .covs))
                      } else {
                        data.frame(ID=id, .base)
                      }
                    }))
  .u <- unique(unlist(lapply(.ret0, function(x){
    names(x)
  })))
  .ret0 <- lapply(.ret0, function(x) {
    .d <- setdiff(.u, names(x))
    if (any(.d == "CENS")) {
      x$CENS <- 0
    }
    if (any(.d == "LIMIT")) {
      x$LIMIT <- NA_real_
    }
    .d <- setdiff(.d, c("CENS", "LIMIT"))
    for (.c in .d) {
      x[[.c]] <- NA_real_
    }
    x
  })
  .ret <- do.call("rbind", .ret0)
  attr(.ret$ID, "levels") <- .idLvl
  class(.ret$ID) <- "factor"
  .ret <- .ret[order(.ret$ID, .ret$TIME), ]
  .ret
}
#' Get the ipred parameters for the fit
#'
#' @param fit the fit to extract the ipred parameters from
#' @return data.frame of parameters for ipred
#' @noRd
#' @author Matthew L. Fidler
.nlmixrGetIpredParams <- function(fit) {
  .si <- fit$simInfo
  .sigma <- .si$sigma
  .omega <- .si$omega
  if (is.null(.omega)) {
    .params <- data.frame(t(fit$theta),
                          t(setNames(rep(0, dim(.sigma)[1]), dimnames(.sigma)[[2]])))
    .params <- setNames(as.numeric(.params), names(.params))
  } else {
    .params <- data.frame(t(fit$theta),fit$eta[, -1, drop = FALSE],
                          t(setNames(rep(0, dim(.sigma)[1]), dimnames(.sigma)[[2]])))
  }
  .params
}

#' Augmented Prediction for nlmixr2 fit
#'
#' @param fit Nlmixr2 fit object
#' @inheritParams nlme::augPred
#' @inheritParams rxode2::rxSolve
#' @return Stacked data.frame with observations, individual/population predictions.
#' @author Matthew L. Fidler
#' @export
nlmixr2AugPredSolve <- function(fit, covsInterpolation = c("locf", "nocb", "linear", "midpoint"),
                                minimum = NULL, maximum = NULL, length.out = 51L, ...) {
  .si <- fit$simInfo
  .env <- new.env(parent=emptyenv())
  .env$ui <- fit$ui
  .env$data <- fit$origData
  suppressMessages(.preProcessHooksRun(.env, "rxSolve"))
  .rx <- .getSimModel(.env$ui, hideIpred=TRUE)
  .rx <- eval(.rx)
  .sigma <- .si$sigma
  .omega <- .si$omega
  .params <- .nlmixrGetIpredParams(fit)
  .events <- .augPredExpandData(fit, covsInterpolation = covsInterpolation,
                                minimum = minimum, maximum = maximum,
                                length.out = length.out)
  # ipred
  .sim <- rxode2::rxSolve(object=.rx, .params, .events,
                          keepInterpolation="na",
                          keep=c("DV", "CMT"), returnType="data.frame")
  # now do pred
  if (is.null(.omega)) {
    names(.sim) <- sub("sim", "pred", names(.sim))
    .stk <- stack(.sim[, c("pred", "DV")])
  } else {
    names(.sim) <- sub("sim", "ipred", names(.sim))
    .params <- c(t(fit$theta),t(rep(0, dim(.omega)[1])),
                 t(rep(0, dim(.sigma)[1])))
    .params <- setNames(.params, c(names(fit$theta),
                                   dimnames(.omega)[[2]],
                                   dimnames(.sigma)[[2]]))
    .sim2 <- rxode2::rxSolve(object=.rx, params=.params, events=.events,
                             returnType="data.frame")
    .sim$pred <- .sim2$sim
    .stk <- stack(.sim[, c("ipred", "pred", "DV")])
  }

  .stk$id <- .sim$id
  .stk$time <- .sim$time
  .stk$cmt <- as.integer(.sim$CMT)
  .ipredModel <- .augPredIpredModel(fit)
  .lvl <- c(.ipredModel$state, .ipredModel$stateExtra)
  if (length(.lvl) == 1L && .lvl == "rxLinCmt") {
    if (rxModelVars(fit)$flags["ka"] == c(ka=1L)) {
      .lvl <- c("depot", "central")
    } else {
      .lvl <- "central"
    }
  }
  levels(.stk$cmt) <- .lvl
  class(.stk$cmt) <- "factor"
  .stk <- .stk[!is.na(.stk$values), ]
  class(.stk) <- c("nlmixr2AugPred", "data.frame")
  if (is.null(.omega)) {
    levels(.stk$ind) <- sub("pred", "Population",
                            sub("DV", "Observed", levels(.stk$ind)))
  } else {
    levels(.stk$ind) <- sub("pred", "Population",
                            sub("ipred", "Individual",
                                sub("DV", "Observed", levels(.stk$ind))))
  }
  .stk$Endpoint <- factor(paste(.stk$cmt))
  .stk <- .stk[, names(.stk) != "cmt"]
  class(.stk) <- c("nlmixr2AugPred", "data.frame")
  .stk
}

#' @rdname nlmixr2AugPredSolve
#' @export
augPred.nlmixr2FitData <- function(object, primary = NULL, minimum = NULL, maximum = NULL,
                                                    length.out = 51, ...) {
  nlmixr2AugPredSolve(
    fit=object, minimum = minimum, maximum = maximum,
    length.out = length.out, ...
  )
}
