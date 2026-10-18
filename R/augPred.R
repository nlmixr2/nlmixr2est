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
#' - Retrieves the simulation model from the provided `fit` object using `.getSimModel` with
#'   `hideIpred` and `tad` set to `FALSE`.
#'
#' - Evaluates the model variables using `rxModelVars`.
#'
.augPredIpredModel <- function(fit) {
  .ipredModel <- .getSimModel(fit, hideIpred = FALSE, tad = FALSE)
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
.augPredExpandData <- function(
  fit,
  covsInterpolation = c("locf", "nocb", "linear", "midpoint"),
  minimum = NULL,
  maximum = NULL,
  length.out = 51L
) {
  # a fit with IOV: its occasion variables, which the user's model does not
  # read, ride along as covariates, so each time takes the occasion the
  # fit's covariate interpolation and record order give it
  .occCols <- as.character(names(fit$env$iov))
  .origData <- rxode2::etTrans(
    .augPredDataSavOcc(fit, .occCols),
    .augPredIpredModel(fit),
    addCmt = TRUE,
    keepDosingOnly = TRUE,
    allTimeVar = TRUE,
    keep = .occCols,
    addlKeepsCov = fit$control$rxControl$addlKeepsCov,
    addlDropSs = fit$control$rxControl$addlDropSs,
    ssAtDoseTime = fit$control$rxControl$ssAtDoseTime
  )
  .idLvl <- attr(class(.origData), ".rxode2.lst")$idLvl
  if (length(.occCols) > 0L) {
    .origData <- cbind(as.data.frame(.origData), attr(class(.origData), ".rxode2.lst")$keepL[[1]])
  }
  .predDf <- fit$ui$predDf
  .range <- range(.origData$TIME)
  .covs <- c(fit$ui$allCovs, .occCols)
  .obsData <- .origData[.origData$EVID %in% c(0, 2), ]
  .allCmt <- unique(.obsData$CMT)
  if (is.null(minimum)) {
    minimum <- .range[1]
  }
  if (is.null(maximum)) {
    maximum <- .range[2]
  }
  .fs <- c(locf = 0, nocb = 1, midpoint = 0, linear = 0)
  .base <- expand.grid(
    TIME = seq(minimum, maximum, length.out = length.out),
    EVID = 2,
    AMT = NA_real_,
    II = NA_real_,
    DV = NA_real_,
    CMT = .allCmt
  )
  .covsi <- match.arg(covsInterpolation)
  .ret0 <- c(
    list(as.data.frame(.origData)),
    lapply(seq_along(.idLvl), function(id) {
      .cur <- .origData[.origData$ID == id, ]
      if (length(.covs) > 0) {
        cbind(
          data.frame(ID = id, .base),
          setNames(
            data.frame(lapply(.covs, function(cov) {
              suppressWarnings({
                .x <- .cur$TIME
                .y <- .cur[[cov]]
                if (.covsi == "midpoint" && length(.x) > 1L) {
                  # rxode2 steps to the next value at the midpoint, it does
                  # not average the two
                  .o <- order(.x)
                  .x <- .x[.o]
                  .y <- .y[.o]
                  .x <- c(.x[1], (.x[-1] + .x[-length(.x)]) / 2)
                }
                # tied records keep their record order (etTrans()'s), so an
                # occasion steps through them as rxode2 does, not to a mean
                .fun <- stats::approxfun(
                  .x,
                  .y,
                  method = ifelse(.covsi == "linear", "linear", "constant"),
                  rule = 2,
                  f = .fs[.covsi],
                  ties = if (cov %in% .occCols) "ordered" else mean
                )
                .fun(.base$TIME)
              })
            })),
            .covs
          )
        )
      } else {
        data.frame(ID = id, .base)
      }
    })
  )
  .u <- unique(unlist(lapply(.ret0, function(x) {
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
    .params <- data.frame(t(fit$theta), t(setNames(rep(0, dim(.sigma)[1]), dimnames(.sigma)[[2]])))
    .params <- setNames(as.numeric(.params), names(.params))
  } else {
    .params <- data.frame(
      t(fit$theta),
      fit$eta[, -1, drop = FALSE],
      t(setNames(rep(0, dim(.sigma)[1]), dimnames(.sigma)[[2]]))
    )
  }
  .params
}

#' The fit's processed data with its occasion columns
#'
#' `fit$dataSav` has only the user's model's columns; an occasion variable is
#' not one of them.  Its rows are rows of `fit$origData` (`nlmixrRowNums`).
#'
#' @param fit nlmixr2 fit
#' @param occCols the occasion variables (`names(fit$env$iov)`)
#' @return `fit$dataSav`, with each occasion variable as a column
#' @noRd
#' @author Matthew L. Fidler
.augPredDataSavOcc <- function(fit, occCols) {
  .ds <- fit$dataSav
  if (length(occCols) == 0L) {
    return(.ds)
  }
  .orig <- as.data.frame(fit$origData)
  .w <- match(tolower(occCols), tolower(names(.orig)))
  for (.i in seq_along(occCols)) {
    .ds[[occCols[.i]]] <- .orig[[.w[.i]]][.ds$nlmixrRowNums]
  }
  .ds
}

#' The fit's occasion etas for rows of data
#'
#' A fit with IOV estimated one value of each occasion eta per subject and
#' occasion (`fit$env$iov`, a table per occasion variable).  Each row takes
#' the value for its own subject and occasion; a pair the fit did not
#' estimate (a new subject or occasion) takes 0, the population value.
#'
#' @param fit nlmixr2 fit with IOV
#' @param id subject label of each row
#' @param occ list named by the occasion variables of `fit$env$iov`: the
#'   occasion of each row
#' @return a list: `etas`, a data frame with one column per occasion eta,
#'   and `unseen`, the number of (subject, occasion) pairs the fit did not
#'   estimate
#' @noRd
#' @author Matthew L. Fidler
.iovEtaRows <- function(fit, id, occ) {
  .iov <- fit$env$iov
  .etas <- list()
  .unseen <- character(0)
  for (.lvl in names(.iov)) {
    .tab <- .iov[[.lvl]]
    .key <- paste(id, occ[[.lvl]])
    .row <- match(.key, paste(as.character(.tab$ID), .tab[[.lvl]]))
    if (anyNA(.row)) {
      .unseen <- c(.unseen, paste(.lvl, .key[is.na(.row)]))
    }
    for (.eta in setdiff(names(.tab), c("ID", .lvl))) {
      .v <- .tab[[.eta]][.row]
      .v[is.na(.v)] <- 0
      .etas[[.eta]] <- .v
    }
  }
  list(etas = as.data.frame(.etas), unseen = length(unique(.unseen)))
}

#' augPred()'s events with the fit's occasion etas
#'
#' Each row takes its subject's estimate for its occasion, which
#' `.augPredExpandData()` gives every row of a fit with IOV, as columns the
#' model reads as those etas.  A grid point between two occasions under
#' linear covariate interpolation has no occasion of its own and takes 0.
#'
#' @param fit nlmixr2 fit
#' @param events augPred()'s events (`.augPredExpandData()`), whose `ID` is
#'   a factor of the fit's subject labels
#' @return `events`, with one column per occasion eta when the fit has IOV
#' @noRd
#' @author Matthew L. Fidler
.augPredIovEvents <- function(fit, events) {
  .occCols <- names(fit$env$iov)
  if (length(.occCols) == 0L) {
    return(events)
  }
  .r <- .iovEtaRows(fit, as.character(events$ID), as.list(events[.occCols]))
  events[names(.r$etas)] <- .r$etas
  events
}

#' Data with the fit's occasion etas, for individual predictions
#'
#' The model reads each occasion eta from its column, so none is simulated.
#' A (subject, occasion) pair the fit did not estimate gets 0, with a
#' warning.
#'
#' @param fit nlmixr2 fit with IOV
#' @param data data to predict, with the subject and occasion columns
#' @return `data` as a data frame, with one column per occasion eta
#' @noRd
#' @author Matthew L. Fidler
.nlmixrIovEtaData <- function(fit, data) {
  .d <- as.data.frame(data)
  .occCols <- names(fit$env$iov)
  .w <- match(tolower(c("ID", .occCols)), tolower(names(.d)))
  if (anyNA(.w)) {
    stop(
      "the data needs ",
      paste(c("ID", .occCols)[is.na(.w)], collapse = ", "),
      " for the fit's occasion etas",
      call. = FALSE
    )
  }
  .occ <- stats::setNames(as.list(.d[.w[-1]]), .occCols)
  .r <- .iovEtaRows(fit, as.character(.d[[.w[1]]]), .occ)
  if (.r$unseen > 0L) {
    warning(
      .r$unseen,
      " subject-occasion pair(s) not in the fit get occasion etas of 0",
      call. = FALSE
    )
  }
  .d[names(.r$etas)] <- .r$etas
  .d
}

#' Name of the covariate interpolation a fit was estimated with
#'
#' @param fit nlmixr2 fit
#' @return interpolation name; `"locf"` when the fit did not set one
#' @noRd
.augPredCovsInterpolation <- function(fit) {
  .codes <- c(linear = 0L, locf = 1L, nocb = 2L, midpoint = 3L)
  .ci <- .residCovsInterpolation(fit)$covsInterpolation
  if (is.character(.ci) && length(.ci) == 1L && .ci %in% names(.codes)) {
    return(.ci)
  }
  .name <- names(.codes)[match(as.integer(.ci), .codes)]
  if (length(.name) != 1L || is.na(.name)) {
    return("locf")
  }
  .name
}

#' Augmented Prediction for nlmixr2 fit
#'
#' @param fit Nlmixr2 fit object
#' @param covsInterpolation covariate interpolation; by default the one
#'   the fit was estimated with
#' @inheritParams nlme::augPred
#' @inheritParams rxode2::rxSolve
#' @details
#'
#' For a fit with inter-occasion variability, each time takes its subject's
#' estimate for the occasion the fit's covariate interpolation gives that
#' time, in the fit's record order where records tie.  Under `"midpoint"`
#' or `"linear"` interpolation the individual curve differs from the fit's
#' predictions after an occasion changes, as `predict()` does.
#'
#' For a `linCmt()` model with time-varying covariates or occasions the
#' curves are approximate: under the default `"locf"` interpolation,
#' rxode2's `linCmt()` takes a covariate's value at the end of each interval,
#' as `"nocb"` does, so the times `augPred()` adds move where a change takes
#' effect.
#'
#' @return Stacked data.frame with observations, individual/population predictions.
#' @author Matthew L. Fidler
#' @export
nlmixr2AugPredSolve <- function(
  fit,
  covsInterpolation = c("locf", "nocb", "linear", "midpoint"),
  minimum = NULL,
  maximum = NULL,
  length.out = 51L,
  ...
) {
  if (missing(covsInterpolation)) {
    covsInterpolation <- .augPredCovsInterpolation(fit)
  }
  covsInterpolation <- match.arg(covsInterpolation)
  .naInterpolation <- .residCovsInterpolation(fit)$naInterpolation
  .nonmem <- .residCovsInterpolation(fit)$nonmem
  .si <- fit$simInfo
  .env <- new.env(parent = emptyenv())
  .env$ui <- fit$ui
  .env$data <- fit$origData
  suppressMessages(.preProcessHooksRun(.env, "rxSolve"))
  .rx <- .getSimModel(.env$ui, hideIpred = TRUE)
  .drop <- setdiff(attr(.rx, "lagged"), c("DV", "CMT"))
  if (length(.drop) == 0L) {
    .drop <- NULL
  }
  .rx <- eval(.rx)
  .sigma <- .si$sigma
  .omega <- .si$omega
  .params <- .nlmixrGetIpredParams(fit)
  .events <- .augPredExpandData(
    fit,
    covsInterpolation = covsInterpolation,
    minimum = minimum,
    maximum = maximum,
    length.out = length.out
  )
  .tolFactor <- fit$env$tolFactor
  # ipred
  .sim <- .rxSolveNonmem(
    object = .rx,
    .params,
    .augPredIovEvents(fit, .events),
    keepInterpolation = "na",
    covsInterpolation = covsInterpolation,
    naInterpolation = .naInterpolation,
    nonmem = .nonmem,
    tolFactor = .tolFactor,
    keep = c("DV", "CMT"),
    drop = .drop,
    returnType = "data.frame"
  )
  # now do pred
  if (is.null(.omega)) {
    names(.sim) <- sub("^sim$", "pred", names(.sim))
    .stk <- stack(.sim[, c("pred", "DV")])
  } else {
    names(.sim) <- sub("^sim$", "ipred", names(.sim))
    # every eta, of every level of a nested omega, is zero
    .etaNames <- dimnames(lotri::lotriMat(.omega))[[2]]
    .params <- c(t(fit$theta), t(rep(0, length(.etaNames))), t(rep(0, dim(.sigma)[1])))
    .params <- setNames(.params, c(names(fit$theta), .etaNames, dimnames(.sigma)[[2]]))
    .sim2 <- .rxSolveNonmem(
      object = .rx,
      params = .params,
      events = .events,
      covsInterpolation = covsInterpolation,
      naInterpolation = .naInterpolation,
      nonmem = .nonmem,
      tolFactor = .tolFactor,
      drop = .drop,
      returnType = "data.frame"
    )
    .sim$pred <- .sim2$sim
    .stk <- stack(.sim[, c("ipred", "pred", "DV")])
  }

  .stk$id <- .sim$id
  # rxSolve returns integer-coded ids; restore the fit's original ID labels
  # (e.g. character/factor subject ids) so augPred identifies the actual subject
  .idLvl <- fit$idLvl
  if (!is.null(.idLvl) && is.factor(.stk$id)) {
    levels(.stk$id) <- .idLvl[as.integer(levels(.stk$id))]
  }
  .stk$time <- .sim$time
  .stk$cmt <- as.integer(.sim$CMT)
  .ipredModel <- .augPredIpredModel(fit)
  .lvl <- c(.ipredModel$state, .ipredModel$stateExtra)
  if (length(.lvl) == 1L && .lvl == "rxLinCmt") {
    if (rxModelVars(fit)$flags["ka"] == c(ka = 1L)) {
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
    levels(.stk$ind) <- sub("pred", "Population", sub("DV", "Observed", levels(.stk$ind)))
  } else {
    levels(.stk$ind) <- sub("pred", "Population", sub("ipred", "Individual", sub("DV", "Observed", levels(.stk$ind))))
  }
  .stk$Endpoint <- factor(paste(.stk$cmt))
  .stk <- .stk[, names(.stk) != "cmt"]
  class(.stk) <- c("nlmixr2AugPred", "data.frame")
  .stk
}

#' @rdname nlmixr2AugPredSolve
#' @export
augPred.nlmixr2FitData <- function(object, primary = NULL, minimum = NULL, maximum = NULL, length.out = 51, ...) {
  ## event bus: one solveComplete linked to the fit (see rxEvents.R)
  force(object)
  .nlmixr2EventEnter()
  .evCall <- sys.call()
  on.exit(.nlmixr2EventExitSolve(returnValue(), object, .evCall, "augPred"), add = TRUE)
  nlmixr2AugPredSolve(
    fit = object,
    minimum = minimum,
    maximum = maximum,
    length.out = length.out,
    ...
  )
}
