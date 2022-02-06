.nlmixr2Time <- NULL
.currentTimingEnvironment <- NULL
.extraTimingTable <- NULL

.finalizeOverallTiming <- function() {
  on.exit({
    assignInMyNamespace(".nlmixr2Time", NULL)
    assignInMyNamespace(".currentTimingEnvironment", NULL)
    assignInMyNamespace(".extraTimingTable", NULL)
  })
  if (is.environment(.currentTimingEnvironment) &
        inherits(.nlmixr2Time, "proc_time")) {
    .time <- get("time", envir=.currentTimingEnvironment)
    if (inherits(.extraTimingTable, "data.frame")) {
      .time <- cbind(.time, .extraTimingTable)
    }
    # Keep non-zero times
    .keep <- vapply(seq_along(names(.time)),
                    function(i){
                      !(.time[[i]] < 5e-5)
                    }, logical(1), USE.NAMES=FALSE)
    .time <- .time[, .keep]
    .sum <- sum(vapply(seq_along(names(.time)),
                       function(i) {
                         .time[[i]]
                       }, numeric(1), USE.NAMES=TRUE))
    .other <- (proc.time() - .nlmixr2Time)["elapsed"] - .sum
    if (.other > 5e-5) {
      .time <- cbind(.time, data.frame(other=.other, row.names="elapsed"))
    }
    assign("time", .time, envir=.currentTimingEnvironment)
  }
}

.timingStack <- NULL

.nlmixrPushTimingStack <- function(name) {
  assignInMyNamespace(".timingStack", c(.timingStack, setNames(0, name)))
}

.nlmixrPopTimingStack <- function(preTiming) {
  .time <- (proc.time() - preTiming)["elapsed"]
  .lastTime <- setNames(.timingStack[length(.timingStack)], NULL)
  if (length(.timingStack) == 1L) {
    assign(".timingStack", NULL)
  } else {
    assign(".timingStack", .timingStack[-length(.timingStack)] + .time)
  }
  .time - .lastTime
}

.nlmixrFinalizeTiming <- function(name, preTiming, envir=NULL) {
  .time <- .nlmixrPopTimingStack(preTiming)
  .df <- list(.time)
  names(.df) <- name
  .df <- as.data.frame(.df, check.names=FALSE, row.names="elapsed")
  if (is.environment(.currentTimingEnvironment) & !is.environment(envir)) {
    envir <- .currentTimingEnvironment
  }
  if (is.environment(envir)) {
    .time <- get("time", envir=envir)
    assign("time",cbind(.time, .df), envir=envir)
  } else {
    if (inherits(.extraTimingTable, "data.frame")) {
      assignInMyNamespace(".extraTimingTable", cbind(.extraTimingTable, .df))
    } else {
      assignInMyNamespace(".extraTimingTable", .df)
    }
  }
}

.nlmixrWithTiming <- function(name, code, envir=NULL) {
  .pt <- proc.time()
  force(name)
  force(envir)
  .nlmixrPushTimingStack(name)
  on.exit(.nlmixrFinalizeTiming(name, .pt, envir), add=TRUE)
  force(code)
}
