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
  if (!is.null(.nlmixr2SimInfo)) {
    .thetaMat <- .nlmixr2SimInfo$thetaMat
    if (is.null(.rxControl$thetaMat) & !.isPred) {
      .minfo("using population uncertainty from fitted model (`thetaMat`)")
      .rxControl$thetaMat <- .thetaMat
    }
    if (.rxControl$dfObs == 0L & !.isPred) {
      .minfo(paste0("using `dfObs=", .nlmixr2SimInfo$dfObs,
             "` from the number of observations in fitted model"))
      .rxControl$dfObs <- .nlmixr2SimInfo$dfObs
    }
    if (.rxControl$dfSub == 0L & !.isPred) {
      .minfo(paste0("using `dfSub=", .nlmixr2SimInfo$dfSub,
             "` from the number of subjects in fitted model"))
      .rxControl$dfSub <- .nlmixr2SimInfo$dfSub
    }

    if (is.null(.rxControl$sigma) & !.isPred) {
      .minfo("using diagonal `sigma` based on model")
      .rxControl$sigma <- .nlmixr2SimInfo$sigma
    }
  }
  .rxControl
}

##' @rdname nmObjGet
##' @export
nmObjGet.rxControlWithVar <- function(x, ...) {
  .tmp <- x[[1]]
  assignInMyNamespace(".nlmixr2SimInfo", .tmp$simInfo)
  .env <- .tmp$env
  if (exists("control", .env)) {
    .oldControl <- get("control", .env)
    on.exit({
      assignInMyNamespace(".nlmixr2SimInfo", NULL)
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
      assignInMyNamespace(".nlmixr2SimInfo", NULL)
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

#'@rdname nlmixr2Est
#'@export
nlmixr2Est.simulate <- function(env, ...) {
  .rxControl <- .rxSolveGetControlForNlmixr(env)
  .events <- get("data", envir=env)
  do.call(rxode2::rxSolve, c(list(object = get("ui", envir=env), params = NULL,
                                  events = .events, inits = NULL), .rxSolveGetControlForNlmixr(env),
                             list(theta = NULL, eta = NULL)))
}

#'@rdname nlmixr2Est
#'@export
nlmixr2Est.simulation <- function(env, ...) {
  .rxControl <- .rxSolveGetControlForNlmixr(env)
  .events <- get("data", envir=env)
  do.call(rxode2::rxSolve, c(list(object = get("ui", envir=env), params = NULL,
                                  events = .events, inits = NULL), .rxSolveGetControlForNlmixr(env),
                             list(theta = NULL, eta = NULL)))
}


#'@rdname nlmixr2Est
#'@export
nlmixr2Est.predict <- function(env, ...) {
  .rxControl <- .rxSolveGetControlForNlmixr(env)
  .events <- get("data", envir=env)
  if (is.na(.rxControl$simVariability)) {
    .rxControl$simVariability <- FALSE
  }
  do.call(rxode2::rxSolve, c(list(object = get("ui", envir=env), params = NULL,
                                  events = .events, inits = NULL), .rxSolveGetControlForNlmixr(env),
                             list(theta = NULL, eta = NULL)))
}

#' @export
predict.nlmixr2FitCore <- function(object, ...) {
  .env <- .nlmixrEvalEnv$envir
  if (!is.environment(.env)) {
    .env <- parent.frame(1)
  }
  .both <- .getControlFromDots(rxode2::rxControl(envir=.env), ...)
  .both$ctl$omega <- NA
  .both$ctl$sigma <- NA
  .env <- .nlmixrEvalEnv$envir
  if (!is.environment(.env)) {
    .env <- parent.frame(1)
  }
  .rxControl <- do.call(rxode2::rxControl, .both$ctl)
  .rxControl$envir <- .env
  if (inherits(.both$rest$newdata, "data.frame")) {
    nlmixr2(object=object, data=.both$rest$newdata,
            est="rxSolve", control=.rxControl)
  } else {
    nlmixr2(object=object, est="rxSolve", control=.rxControl)
  }
}

#' @export
simulate.nlmixr2FitCore <- function(object, ...) {
  .env <- .nlmixrEvalEnv$envir
  if (!is.environment(.env)) {
    .env <- parent.frame(1)
  }
  .both <- .getControlFromDots(rxode2::rxControl(envir=.env), ...)
  .rxControl <- do.call(rxode2::rxControl, .both$ctl)
  .rxControl$envir <- .env
  if (inherits(.both$rest$newdata, "data.frame")) {
    nlmixr2(object=object, data=.both$rest$newdata,
            est="rxSolve", control=.rxControl)
  } else {
    nlmixr2(object=object, est="rxSolve", control=.rxControl)
  }
}
