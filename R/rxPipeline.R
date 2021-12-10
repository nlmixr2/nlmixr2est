##' @export
et.nlmixr2FitData <- function(x, ..., envir = parent.frame()) {
  .si <- x$simInfo
  rxode2::.clearPipe(
    rx = rxode2::rxode2(.si$rx),
    ## events=nlmixr2Data(getData(x)),
    params = .si$params,
    thetaMat = .si$thetaMat,
    dfObs = .si$dfObs,
    omega = .si$omega,
    dfSub = .si$dfSub,
    sigma = .si$sigma
  )
  do.call(rxode2::et, c(list(...), list(envir = envir)), envir = envir)
}

##' @export
rxParams.nlmixr2FitData <- function(obj, ...) {
  .si <- obj$simInfo
  rxode2::.clearPipe(
    rx = rxode2::rxode2(.si$rx),
    events = nlmixr2Data(getData(obj)),
    thetaMat = .si$thetaMat,
    dfObs = .si$dfObs,
    omega = .si$omega,
    dfSub = .si$dfSub,
    sigma = .si$sigma
  )
  do.call(rxode2::rxParams, list(...))
}
