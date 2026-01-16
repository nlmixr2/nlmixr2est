#' Get Control Settings for nlmixr
#'
#' This function retrieves and sets control settings for the `nlmixr`
#' package from the given environment  It ensures that the control
#' settings are valid and, if necessary, uses default settings from
#' `rxode2::rxControl()`.
#'
#' @param env Environment from which to retrieve control settings.
#' @return A list of control settings for `rxSolve`.
#' @details
#' The function performs the following steps:
#'
#' - Retrieves the `ui` object from the provided environment.
#'
#' - Checks if a `control` object exists in the environment and
#' retrieves it.
#'
#' - Validates if the retrieved `control` object is of class
#' `rxControl`. If not, it attempts to retrieve the `rxControl`
#' element from the `control` object.
#'
#' - If the `rxControl` object is still not valid, it uses default
#' solving options from `rxode2::rxControl()`.
#'
#' - Determines if the model is a prediction model based on the
#' `omega` and `sigma` values.
#'
#' - If additional simulation information (`nlmixr2global$nlmixr2SimInfo`) is
#' available, it updates the `rxControl` object with population
#' uncertainty, number of observations, number of subjects, and
#' diagonal `sigma` based on the fitted model.
#'
#' - Checks if a `table` object exists in the environment.  If it
#' does, adjust the rxode2 solving control options by preferring
#' non-default values from table as well as combining `keep` and
#' `drop` from `tableControl()`.  If `cores` is non-NULL, use that
#' instead of the value from `rxControl()`.
#'
#' @noRd
.rxSolveGetControlForNlmixr <- function(env) {
  .ui <- get("ui", envir=env)
  if (exists("control", envir=env)) {
    .rxControl <- get("control", envir=env)
  }
  if (!inherits(.rxControl, "rxControl")) {
    .rxControl <- try(.rxControl$rxControl)
    if (!inherits(.rxControl, "rxControl")) {
      .minfo("using default solving options `rxode2::rxControl()`")
      .rxControl <- rxode2::rxControl()
    }
  }
  .isPred <- FALSE
  if (length(.rxControl$omega) == 1 && length(.rxControl$sigma) == 1) {
    .isPred <- is.na(.rxControl$omega) & (is.na(.rxControl$sigma))
  }
  if (!.isPred) {
    if (is.na(.rxControl$simVariability)) {
      if (.rxControl$nStud == 1) {
        .isPred <- TRUE
      }
    } else if (!.rxControl$simVariability) {
      .isPred <- TRUE
    }
  }
  if (!is.null(nlmixr2global$nlmixr2SimInfo)) {
    .thetaMat <- nlmixr2global$nlmixr2SimInfo$thetaMat
    if (is.null(.rxControl$thetaMat) & !.isPred) {
      .minfo("using population uncertainty from fitted model (`thetaMat`)")
      .rxControl$thetaMat <- .thetaMat
    }
    if (.rxControl$dfObs == 0L & !.isPred) {
      .minfo(paste0("using `dfObs=", nlmixr2global$nlmixr2SimInfo$dfObs,
             "` from the number of observations in fitted model"))
      .rxControl$dfObs <- nlmixr2global$nlmixr2SimInfo$dfObs
    }
    if (.rxControl$dfSub == 0L & !.isPred) {
      .minfo(paste0("using `dfSub=", nlmixr2global$nlmixr2SimInfo$dfSub,
             "` from the number of subjects in fitted model"))
      .rxControl$dfSub <- nlmixr2global$nlmixr2SimInfo$dfSub
    }

    if (is.null(.rxControl$sigma) & !.isPred) {
      .minfo("using diagonal `sigma` based on model")
      .rxControl$sigma <- nlmixr2global$nlmixr2SimInfo$sigma
    }
  }
  if (exists("table", envir=env) &&
        !is.null(env$table)) {
    .table <- env$table
    if (checkmate::testLogical(.table$covariates, any.missing=FALSE, len=1) &&
          !.table$covariates && .rxControl$addCov) {
      .rxControl$addCov <- FALSE
    }
    if (checkmate::testLogical(.table$addDosing, any.missing=FALSE, len=1) &&
        .table$addDosing && !.rxControl$addDosing) {
      .rxControl$addDosing <- TRUE
    }
    if (checkmate::testLogical(.table$subsetNonmem, any.missing=FALSE, len=1) &&
          !.table$subsetNonmem && .rxControl$subsetNonmem) {
      .rxControl$subsetNonmem <- FALSE
    }
    if (checkmate::testIntegerish(.table$cores, len=1, lower=1, any.missing=FALSE)) {
      .rxControl$cores <-.table$cores
    }
    if (checkmate::testCharacter(.table$keep, any.missing=FALSE)) {
      .keep <- unique(c(.table$keep, .rxControl$keep))
      .rxControl$keep <- .keep
    }
    if (checkmate::testCharacter(.table$drop, any.missing=FALSE)) {
      .drop <- unique(c(.table$drop, .rxControl$drop))
      .rxControl$drop <- .drop
    }
  }
  .rxControl
}

##' @rdname nmObjGet
##' @export
nmObjGet.rxControlWithVar <- function(x, ...) {
  .tmp <- x[[1]]
  nlmixr2global$nlmixr2SimInfo <- .tmp$simInfo
  .env <- .tmp$env
  if (exists("control", .env)) {
    .oldControl <- get("control", .env)
    on.exit({
      nlmixr2global$nlmixr2SimInfo <- NULL
      assign("control", .oldControl, envir=.env)})
    if (!inherits(.oldControl, "rxControl")) {
      .rxControl <- nmObjGet.rxControl(x, ...)
    } else {
      .rxControl <- .oldControl
    }
    assign("control", .rxControl, envir=.env)
  } else {
    .rxControl <- nmObjGet.rxControl(x, ...)
    assign("control", .rxControl, envir=.env)
    on.exit({
      nlmixr2global$nlmixr2SimInfo <- NULL
      if (exists("control", envir=.env)) {
        rm(list="control", envir=.env)
      }
    })
  }
  .rxSolveGetControlForNlmixr(.env)
}

#'@rdname nlmixr2Est
#'@export
nlmixr2Est.rxSolve <- function(env, ...) {
  .events <- get("data", envir=env)
  do.call(rxode2::rxSolve, c(list(object = get("ui", envir=env), params = NULL,
                                  events = .events, inits = NULL), .rxSolveGetControlForNlmixr(env),
                             list(theta = NULL, eta = NULL)))
}
attr(nlmixr2Est.rxSolve, "covPresent") <- TRUE

#'@rdname nlmixr2Est
#'@export
nlmixr2Est.simulate <- function(env, ...) {
  .rxControl <- .rxSolveGetControlForNlmixr(env)
  .events <- get("data", envir=env)
  do.call(rxode2::rxSolve, c(list(object = get("ui", envir=env), params = NULL,
                                  events = .events, inits = NULL), .rxSolveGetControlForNlmixr(env),
                             list(theta = NULL, eta = NULL)))
}
attr(nlmixr2Est.simulate, "covPresent") <- TRUE

#'@rdname nlmixr2Est
#'@export
nlmixr2Est.simulation <- function(env, ...) {
  .nlmixr2clearPipe()
  nlmixr2global$nlmixr2SimInfo <- NULL
  on.exit({
    .nlmixr2clearPipe()
    nlmixr2global$nlmixr2SimInfo <- NULL
  })
  .rxControl <- .rxSolveGetControlForNlmixr(env)
  env$control <- .rxControl
  .events <- get("data", envir=env)
  do.call(rxode2::rxSolve, c(list(object = get("ui", envir=env), params = NULL,
                                  events = .events, inits = NULL), .rxControl,
                             list(theta = NULL, eta = NULL)))
}
attr(nlmixr2Est.simulation, "covPresent") <- TRUE


#'@rdname nlmixr2Est
#'@export
nlmixr2Est.predict <- function(env, ...) {
  .nlmixr2clearPipe()
  nlmixr2global$nlmixr2SimInfo <- NULL
  on.exit({
    .nlmixr2clearPipe()
    nlmixr2global$nlmixr2SimInfo <- NULL
  })
  .rxControl <- .rxSolveGetControlForNlmixr(env)
  .rxControl$omega <- NA
  .rxControl$sigma <- NA
  .events <- get("data", envir=env)
  if (is.na(.rxControl$simVariability)) {
    .rxControl$simVariability <- FALSE
  }
  nlmixr2(object=get("ui", envir=env), data=.events,
          est="rxSolve", control=.rxControl)
}
attr(nlmixr2Est.predict, "covPresent") <- TRUE

#' Get new data
#'
#'
#' @param both the adjusted control
#' @return both adjusted for single data frame to be newdata
#' @noRd
#' @author Matthew L. Fidler
.getNewData <- function(both) {
  .both <- both
  if (!any(names(.both$rest) == "newdata")) {
    .w <- which(vapply(seq_along(.both$rest),
                       function(i) {
                         inherits(.both$rest[[i]], "data.frame")
                       }, logical(1), USE.NAMES=FALSE))
    if (length(.w) == 1L) {
      names(.both$rest)[1] <- "newdata"
    }
  }
  .both
}
#' Predict method for nlmixr2 fit core objects
#'
#' This function generates predictions from an `nlmixr2FitCore` object.
#' It allows for both population-level and individual-level predictions
#' based on the specified `level` parameter.
#'
#' @param object nlmixr2 fit core object to predict
#'
#' @param ... additional arguments passed to rxode2::rxSolve or
#'   nlmixr2; matching other `predict` methods, these can include
#'   `newdata` and `rxControl` settings
#'
#' @param level the prediction level; one of `"population"` (default) or
#'  `"individual"`; numeric values `0` and `1` are also accepted
#'
#' @return A data frame with predictions
#'
#' @export
#'
#' @examples
#'
#' \donttest{
#'
#' one.compartment <- function() {
#'  ini({
#'   tka <- log(1)
#'   tcl <- log(10)
#'   tv <- log(35)
#'   eta.ka ~ 0.1
#'   eta.cl ~ 0.1
#'   eta.v ~ 0.1
#'   add.sd <- 0.1
#'  })
#'  model({
#'   ka <- exp(tka + eta.ka)
#'   cl <- exp(tcl + eta.cl)
#'   v <- exp(tv + eta.v)
#'   d/dt(depot) = -ka * depot
#'   d/dt(center) = ka * depot - cl / v * center
#'   cp = center / v
#'   cp ~ add(add.sd)
#'  })
#' }
#'
#' # The fit is performed by the function nlmixr/nlmix2 specifying
#' # the model, data and estimate
#' fit <- nlmixr2(one.compartment, theo_sd, est = "focei",
#'                foceiControl(maxOuterIterations = 0L))
#'
#' # Population predictions
#' ppred <- predict(fit, theo_sd, level="population")
#'
#' # Individual predictions
#' ipred <- predict(fit, theo_sd, level="individual")
#'
#' }
#'
predict.nlmixr2FitCore <- function(object, ...,
                                   level = c("population", "individual")) {
  if (checkmate::testNumeric(level, len=1)) {
    level <- switch(as.character(level),
                    "0" = "population",
                    "1" = "individual",
                    "bad")
    if (identical(level, "bad")) {
      stop("level numeric must be 0 (population) or 1 (individual)",
           call.=FALSE)
    }
  }
  if (identical(level, "ipred")) {
    level <- "individual"
  } else if (identical(level, "pred") || identical(level, "ppred")) {
    level <- "population"
  } else {
    level <- match.arg(level)
  }
  .nlmixr2clearPipe()
  nlmixr2global$nlmixr2SimInfo <- NULL
  on.exit({
    .nlmixr2clearPipe()
    nlmixr2global$nlmixr2SimInfo <- NULL
  })
  .env <- nlmixr2global$nlmixrEvalEnv$envir
  if (!is.environment(.env)) {
    .env <- parent.frame(1)
  }
  .est <- if (identical(level, "population")) "predict" else "ipred"
  if (.est == "ipred") {
    .minfo("individual predictions requested (`level=\"individual\"`)")
  } else {
    .minfo("population predictions requested (`level=\"population\"`)")
  }

  .both <- .getNewData(.getControlFromDots(rxode2::rxControl(envir=.env), ...))
  if (.est != "ipred") {
    .both$ctl$omega <- NA
    .both$ctl$sigma <- NA
  }
  .env <- nlmixr2global$nlmixrEvalEnv$envir
  if (!is.environment(.env)) {
    .env <- parent.frame(1)
  }
  .rxControl <- do.call(rxode2::rxControl, .both$ctl)
  .rxControl$envir <- .env
  .data <- getData(object)

  if (inherits(.both$rest$newdata, "data.frame")) {
    .minfo("using new data for predictions")
    .data <- .both$rest$newdata
  }
  if (.est == "ipred") {
    .params <- .nlmixrGetIpredParams(object)
    do.call(rxode2::rxSolve,
            c(list(object, .params, .data), .rxControl))
  } else {
    nlmixr2(object=object, data=.data,
            est=.est, control=.rxControl)
  }
}

#' @export
simulate.nlmixr2FitCore <- function(object, ...) {
  .nlmixr2clearPipe()
  nlmixr2global$nlmixr2SimInfo <- NULL
  on.exit({
    .nlmixr2clearPipe()
    nlmixr2global$nlmixr2SimInfo <- NULL
  })
  .env <- nlmixr2global$nlmixrEvalEnv$envir
  if (!is.environment(.env)) {
    .env <- parent.frame(1)
  }
  .both <- .getNewData(.getControlFromDots(rxode2::rxControl(envir=.env), ...))
  .rxControl <- do.call(rxode2::rxControl, .both$ctl)
  .rxControl$envir <- .env
  if (inherits(.both$rest$newdata, "data.frame")) {
    nlmixr2(object=object, data=.both$rest$newdata,
            est="rxSolve", control=.rxControl)
  } else {
    nlmixr2(object=object, est="rxSolve", control=.rxControl)
  }
}
