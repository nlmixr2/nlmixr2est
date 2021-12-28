.augPredExpandData <- function(fit, covsInterpolation = c("locf", "nocb", "linear", "midpoint"),
                               minimum = NULL, maximum = NULL, length.out = 51L) {

  .origData <- rxode2::etTrans(fit$origData, fit$model$pred.only, addCmt=TRUE, keepDosingOnly=TRUE, allTimeVar=TRUE)
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
  .ret <- do.call("rbind",
                  c(list(as.data.frame(.origData)),
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
                    })))
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
  .si$object <- fit
  .si$rx <- .getSimModel(fit, hideIpred=TRUE)
  .si$nsim <- 1
  .si$keep <- c("DV", "CMT")
  .si$events <- .augPredExpandData(fit, covsInterpolation = covsInterpolation,
                                   minimum = minimum, maximum = maximum,
                                   length.out = length.out)
  .si$modelName <- "augPred ipred"
  .si$dfObs <- 0
  .si$dfSub <- 0
  .si$thetaMat <- NA
  .si$returnType <- "data.frame"
  .params <- .si$params
  .omega <- .si$omega
  .sigma <- .si$sigma
  .si$omega <- NA
  .si$sigma <- NA
  .si$params <- data.frame(t(fit$theta),fit$eta,
                           t(setNames(rep(0, dim(.sigma)[1]), dimnames(.sigma)[[2]])))
  rxode2::.setWarnIdSort(FALSE)
  on.exit(rxode2::.setWarnIdSort(TRUE))
  .sim <- do.call("nlmixr2Sim", .si)
  names(.sim) <- sub("sim", "ipred", names(.sim))
  # now do pred
  .si$keep <- NULL
  .si$params <- c(
      .params, setNames(rep(0, dim(.omega)[1]), dimnames(.omega)[[2]]),
      setNames(rep(0, dim(.sigma)[1]), dimnames(.sigma)[[2]]))
  .si$modelName <- "augPred pred"
  .sim2 <- do.call("nlmixr2Sim", .si)
  .sim$pred <- .sim2$sim
  .stk <- stack(.sim[, c("ipred", "pred", "DV")])
  .stk$id <- .sim$id
  .stk$time <- .sim$time
  .stk$cmt <- as.integer(.sim$CMT)
  levels(.stk$cmt) <- c(fit$model$pred.only$state, fit$model$pred.only$stateExtra)
  class(.stk$cmt) <- "factor"
  .stk <- .stk[!is.na(.stk$values), ]
  class(.stk) <- c("nlmixr2AugPred", "data.frame")
  levels(.stk$ind) <- sub("pred", "Population",
                          sub("ipred", "Individual",
                              sub("DV", "Observed", levels(.stk$ind))))

  .stk$Endpoint <- factor(paste(.stk$cmt))
  .stk <- .stk[, names(.stk) != "cmt"]
  class(.stk) <- c("nlmixr2AugPred", "data.frame")
  .stk
}


.augPredEndpoint <- NULL

#' @rdname nlmixr2AugPredSolve
#' @export
augPred.nlmixr2FitData <- memoise::memoise(function(object, primary = NULL, minimum = NULL, maximum = NULL,
                                                    length.out = 51, ...) {
  nlmixr2AugPredSolve(
    fit=object, minimum = minimum, maximum = maximum,
    length.out = length.out, ...
  )
})

#' @export
plot.nlmixr2AugPred <- function(x, y, ...) {
  if (any(names(x) == "Endpoint")) {
    for (.tmp in levels(x$Endpoint)) {
      assignInMyNamespace(".augPredEndpoint", .tmp)
      .x <- x[x$Endpoint == .tmp, names(x) != "Endpoint"]
      plot.nlmixr2AugPred(.x)
    }
  } else {
    ids <- unique(x$id)
    time <- values <- ind <- id <- NULL # Rcheck fix
    for (i in seq(1, length(ids), by = 16)) {
      tmp <- ids[seq(i, i + 15)]
      tmp <- tmp[!is.na(tmp)]
      d1 <- x[x$id %in% tmp, ]
      dobs <- d1[d1$ind == "Observed", ]
      dpred <- d1[d1$ind != "Observed", ]
      p3 <- ggplot(d1, aes(time, values, col = ind)) +
        geom_line(data = dpred, size = 1.2) +
        geom_point(data = dobs) +
        facet_wrap(~id) + rxode2::rxTheme() + ggplot2::ggtitle(label=.augPredEndpoint)
      print(p3)
    }
  }
}
