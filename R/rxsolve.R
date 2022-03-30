.rxSolveGetControlForNlmixr <- function(env) {
  .ui <- get("ui", envir=env)
  .rxControl <- get("control", envir=env)
  if (!inherits(.rxControl, "rxControl")) {
    .rxControl <- rxode2::rxControl()
  }
  .events <- get("data", envir=env)
  if (!is.null(.nlmixr2SimInfo)) {
    .thetaMat <- .nlmixr2SimInfo$thetaMat
    if (is.null(.rxControl$thetaMat)) {
      .minfo("using population uncertainty from fitted model (`thetaMat`)")
      .rxControl$thetaMat <- .thetaMat
    }
    if (.rxControl$dfObs == 0L) {
      .minfo(paste0("using `dfObs=", .nlmixr2SimInfo$dfObs,
             "` from the number of observations in fitted model"))
      .rxControl$dfObs <- .nlmixr2SimInfo$dfObs
    }
    if (.rxControl$dfSub == 0L) {
      .minfo(paste0("using `dfSub=", .nlmixr2SimInfo$dfSub,
             "` from the number of subjects in fitted model"))
      .rxControl$dfSub <- .nlmixr2SimInfo$dfSub
    }

    if (is.null(.rxControl$sigma)) {
      .minfo("using diagonal `sigma` based on model")
      .rxControl$sigma <- .nlmixr2SimInfo$sigma
    }
  }
  .rxControl
}

##' @rdname nmObjGet
##' @export
nmObjGet.rxControlWithVar <- function(x, ...) {
  .rxSolveGetControlForNlmixr(x[[1]])
}

#'@rdname nlmixr2Est
#'@export
nlmixr2Est.rxSolve <- function(env, ...) {
  do.call(rxode2::rxSolve, c(list(object = get("ui", envir=env), params = NULL,
                                  events = .events, inits = NULL), .rxSolveGetControlForNlmixr(env),
                             list(theta = NULL, eta = NULL)))
}

#'@rdname nlmixr2Est
#'@export
nlmixr2Est.simulation <- function(env, ...) {
  .rxControl <- .rxSolveGetControlForNlmixr(env)
  if (is.na(.rxControl$simVariability)) {
    .rxControl$simVariability <- TRUE
  }
  do.call(rxode2::rxSolve, c(list(object = get("ui", envir=env), params = NULL,
                                  events = .events, inits = NULL), .rxSolveGetControlForNlmixr(env),
                             list(theta = NULL, eta = NULL)))
}

#'@rdname nlmixr2Est
#'@export
nlmixr2Est.predict <- function(env, ...) {
  .rxControl <- .rxSolveGetControlForNlmixr(env)
  if (is.na(.rxControl$simVariability)) {
    .rxControl$simVariability <- FALSE
  }
  do.call(rxode2::rxSolve, c(list(object = get("ui", envir=env), params = NULL,
                                  events = .events, inits = NULL), .rxSolveGetControlForNlmixr(env),
                             list(theta = NULL, eta = NULL)))
}

