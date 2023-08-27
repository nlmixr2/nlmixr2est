.augPredIpredModel <- function(fit) {
  .ipredModel <- .getSimModel(fit, hideIpred=FALSE,tad=FALSE)
  eval(as.call(list(quote(`rxModelVars`), .ipredModel[[-1]])))
}

.augPredExpandData <- function(fit, covsInterpolation = c("locf", "nocb", "linear", "midpoint"),
                               minimum = NULL, maximum = NULL, length.out = 51L) {
  .origData <- rxode2::etTrans(fit$dataSav, .augPredIpredModel(fit), addCmt=TRUE, keepDosingOnly=TRUE, allTimeVar=TRUE,
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
  .rx <- .getSimModel(fit, hideIpred=TRUE)
  .rx <- eval(.rx)
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

  .events <- .augPredExpandData(fit, covsInterpolation = covsInterpolation,
                                minimum = minimum, maximum = maximum,
                                length.out = length.out)

  # ipred
  .sim <- rxode2::rxSolve(object=.rx, .params, .events,
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
