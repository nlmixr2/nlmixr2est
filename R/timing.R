#' Push the nlmixr timing stack for a nested nlmixr call
#'
#' @return Nothing, called for side effects
#' @author Matthew L. Fidler
#' @noRd
.pushNlmixr2timing <- function() {
  nlmixr2global$timingStackNlmixr <-
    c(nlmixr2global$timingStackNlmixr,
      list(list(nlmixr2global$nlmixr2Time,
                nlmixr2global$currentTimingEnvironment,
                nlmixr2global$extraTimingTable,
                nlmixr2global$timingStack)))
  nlmixr2global$nlmixr2Time <- NULL
  nlmixr2global$currentTimingEnvironment <- NULL
  nlmixr2global$extraTimingTable <- NULL
  nlmixr2global$timingStack <- NULL
}
#' Pop the full nlmixr timing stack (if needed)
#'
#' @return Nothing, called for side effects
#' @author Matthew L. Fidler
#' @noRd
.popNlmixr2Timing <- function() {
  .l <- length(nlmixr2global$timingStackNlmixr)
  if (.l == 0) {
    nlmixr2global$nlmixr2Time <- NULL
    nlmixr2global$currentTimingEnvironment <- NULL
    nlmixr2global$extraTimingTable <- NULL
    nlmixr2global$timingStack <- NULL
  } else {
    .cur <- nlmixr2global$timingStackNlmixr[[.l]]
    if (.l == 1) {
      nlmixr2global$timingStackNlmixr <- NULL
    } else {
      nlmixr2global$timingStackNlmixr <- nlmixr2global$timingStackNlmixr[[-.l]]
    }
    nlmixr2global$nlmixr2Time <- .cur[[1]]
    nlmixr2global$currentTimingEnvironment <- .cur[[2]]
    nlmixr2global$extraTimingTable <- .cur[[3]]
    nlmixr2global$timingStack <- .cur[[4]]
  }
}

.finalizeOverallTiming <- function() {
  on.exit({
    .popNlmixr2Timing()
  })
  if (is.environment(nlmixr2global$currentTimingEnvironment) &
        inherits(nlmixr2global$nlmixr2Time, "proc_time")) {
    .time <- .nlmixrMergeTimeWithExtraTime(get("time", envir=nlmixr2global$currentTimingEnvironment))
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
    .other <- (proc.time() - nlmixr2global$nlmixr2Time)["elapsed"] - .sum
    if (.other > 5e-5) {
      .time <- cbind(.time, data.frame(other=.other, row.names="elapsed"))
    }
    assign("time", .time, envir=nlmixr2global$currentTimingEnvironment)
  }
}


.nlmixrPushTimingStack <- function(name) {
  nlmixr2global$timingStack <- c(nlmixr2global$timingStack, setNames(0, name))
}

.nlmixrPopTimingStack <- function(preTiming) {
  .lastTime <- setNames(nlmixr2global$timingStack[length(nlmixr2global$timingStack)], NULL)
  if (length(nlmixr2global$timingStack) == 1L) {
    nlmixr2global$timingStack <- NULL
    .time <- setNames((proc.time() - preTiming)["elapsed"], NULL)
  } else {
    .time <- setNames((proc.time() - preTiming)["elapsed"], NULL)
    nlmixr2global$timingStack <- nlmixr2global$timingStack[-length(nlmixr2global$timingStack)] + .time
  }
  .time - .lastTime
}
#' Construct a final table, add the timing at the last second
#'
#' @param time Initial time table
#' @param name Character name of timing
#' @param preTiming Pre-timing to calculate the total time
#' @return New timing table
#' @author Matthew L. Fidler
#' @noRd
.nlmixrFinalizeTimingConstructTable <- function(time, name, preTiming) {
  .w <- which(names(time) == name)
  .time <- time
  if (length(.w) == 1){
    .amt <- .nlmixrPopTimingStack(preTiming)
    if (length(.amt) == 1) {
      if (!is.na(.amt)) {
        .time[, .w] <- .time[, .w] + .amt
      }
    }
  } else {
    .df <- list(0)
    names(.df) <- name
    .df <- as.data.frame(.df, check.names=FALSE, row.names="elapsed")
    .amt <- .nlmixrPopTimingStack(preTiming)
    if (length(.amt) == 1) {
      if (!is.na(.amt)) {
        .df[, 1] <- .amt
        .time <- cbind(.time, .df)
      }
    }
  }
  .time
}
#' Merge the timing table with the internal extraTiming table
#'
#' @param time timing table
#' @return Merge table if needed
#' @author Matthew L. Fidler
#' @noRd
.nlmixrMergeTimeWithExtraTime <- function(time) {
  if (!inherits(nlmixr2global$extraTimingTable, "data.frame")) return(time)
  on.exit({nlmixr2global$extraTimingTable <- NULL})
  .time <- time
  .df <- nlmixr2global$extraTimingTable
  .dropNames <- NULL
  for (.n in names(.df)) {
    .w <- which(names(time) == .n)
    if (length(.w) == 1L) {
      .time[, .w] <- .time[, .w] + .df[, .n]
      .dropNames <- c(.dropNames, .n)
    }
  }
  .df <- .df[, !(names(.df) %fin% .dropNames), drop = FALSE]
  cbind(.time, .df)
}
#' Finalized the timer and integrate into the appropriate place
#'
#' @param name Name of timer
#' @param preTiming Pre-code evaluation time
#' @param envir Environment
#' @return Nothing called for side effects
#' @author Matthew L. Fidler
#' @noRd
.nlmixrFinalizeTiming <- function(name, preTiming, envir=NULL) {
  if (inherits(envir, "nlmixr2FitData")) {
    envir <- envir$env
  }
  if (is.environment(nlmixr2global$currentTimingEnvironment) & !is.environment(envir)) {
    envir <- nlmixr2global$currentTimingEnvironment
  }
  if (is.environment(envir)) {
    .time <- .nlmixrMergeTimeWithExtraTime(get("time", envir=envir))
    assign("time", .nlmixrFinalizeTimingConstructTable(.time, name, preTiming), envir=envir)
  } else {
    if (inherits(nlmixr2global$extraTimingTable, "data.frame")) {
      nlmixr2global$extraTimingTable <-
                          .nlmixrFinalizeTimingConstructTable(nlmixr2global$extraTimingTable, name, preTiming)
    } else {
      .df <- list(0)
      names(.df) <- name
      .df <- as.data.frame(.df, check.names=FALSE, row.names="elapsed")
      .amt <- .nlmixrPopTimingStack(preTiming)
      if (length(.amt) == 1) {
        if (!is.na(.amt)) {
          .df[, 1] <- .amt
          nlmixr2global$extraTimingTable <- .df
        }
      }
    }
  }
}

#' Time a part of a nlmixr operation and add to nlmixr object
#'
#' @param name Name of the timing to be integrated
#' @param code Code to be evaluated and timed
#' @param envir can be either the nlmixr2 fit data, the nlmixr2 fit
#'   environment or NULL, which implies it is going to be added to the
#'   nlmixr fit when it is finalized.  If the function is being called
#'   after a fit is created, please supply this environmental variable
#' @return Result of code
#' @author Matthew L. Fidler
#' @examples
#'
#' \donttest{
#'
#' one.cmt <- function() {
#'  ini({
#'    ## You may label each parameter with a comment
#'    tka <- 0.45 # Ka
#'    tcl <- log(c(0, 2.7, 100)) # Log Cl
#'    ## This works with interactive models
#'    ## You may also label the preceding line with label("label text")
#'    tv <- 3.45; label("log V")
#'    ## the label("Label name") works with all models
#'    eta.ka ~ 0.6
#'    eta.cl ~ 0.3
#'    eta.v ~ 0.1
#'    add.sd <- 0.7
#'  })
#'  model({
#'    ka <- exp(tka + eta.ka)
#'    cl <- exp(tcl + eta.cl)
#'    v <- exp(tv + eta.v)
#'    linCmt() ~ add(add.sd)
#'  })
#' }
#' fit <- nlmixr(one.cmt, theo_sd, est="saem")
#'
#' nlmixrWithTiming("time1", {
#'    Sys.sleep(1)
#'    # note this can be nested, time1 will exclude the timing from time2
#'    nlmixrWithTiming("time2", {
#'       Sys.sleep(1)
#'    }, envir=fit)
#' }, envir=fit)
#'
#' print(fit)
#'
#' }
#'
#' @export
nlmixrWithTiming <- function(name, code, envir=NULL) {
  .pt <- proc.time()
  checkmate::assertCharacter(name, len=1, all.missing=FALSE)
  force(name)
  force(envir)
  if (is.null(envir)){
  } else if (inherits(envir, "nlmixr2FitData")) {
  } else if (is.environment(envir)) {
  } else {
    stop("'envir' must be NULL, a nlmixr2 object or an environment",
         call.=FALSE)
  }
  .nlmixrPushTimingStack(name)
  on.exit(.nlmixrFinalizeTiming(name, .pt, envir), add=TRUE)
  force(code)
}
#' Manually add time to a nlmixr2 object
#'
#' @param object nlmixr2 object
#' @param name string of the timing name
#' @param time time (in seconds)
#' @return Nothing, called for side effects
#' @author Matthew L. Fidler
#' @examples
#'
#' \donttest{
#'
#' one.cmt <- function() {
#'  ini({
#'    ## You may label each parameter with a comment
#'    tka <- 0.45 # Ka
#'    tcl <- log(c(0, 2.7, 100)) # Log Cl
#'    ## This works with interactive models
#'    ## You may also label the preceding line with label("label text")
#'    tv <- 3.45; label("log V")
#'    ## the label("Label name") works with all models
#'    eta.ka ~ 0.6
#'    eta.cl ~ 0.3
#'    eta.v ~ 0.1
#'    add.sd <- 0.7
#'  })
#'  model({
#'    ka <- exp(tka + eta.ka)
#'    cl <- exp(tcl + eta.cl)
#'    v <- exp(tv + eta.v)
#'    linCmt() ~ add(add.sd)
#'  })
#' }
#'
#' fit <- nlmixr(one.cmt, theo_sd, est="saem")
#'
#' # will add to the current setup
#' nlmixrAddTiming(fit, "setup", 3)
#'
#' # Add a new item to the timing dataframe
#' nlmixrAddTiming(fit, "new", 3)
#'
#' }
#'
#' @export
nlmixrAddTiming <- function(object, name, time) {
  .env <- object
  if (inherits(object, "nlmixr2FitData")) {
    .env <- object$env
  }
  .time <- get("time", envir=.env)
  .w <- which(names(.time) == name)
  if (length(.w) == 1L) {
    .time[, .w] <- .time[, .w] + time
  } else {
    if (!is.na(time)) {
      .df <- list(time)
      names(.df) <- name
      .df <- as.data.frame(.df, check.names=FALSE, row.names="elapsed")
      .time <- cbind(.time, .df)
    }
  }
  assign("time", .time, envir=.env)
  invisible()
}
