#' VPC simulation
#'
#' @param object This is the nlmixr2 fit object
#' @param ... Other arguments sent to `rxSolve()`
#' @param keep Column names to keep in the output simulated dataset
#' @param n Number of simulations
#' @param pred Should predictions be added to the simulation
#' @param seed Seed to set for the VPC simulation
#' @param nretry Number of times to retry the simulation if there is
#'   NA values in the simulation
#' @param normRelated should the VPC style simulation be for normal
#'   related variables only
#' @param minN With retries, the minimum number of studies to
#'   restimulate (by default 10)
#' @return data frame of the VPC simulation
#' @author Matthew L. Fidler
#' @export
#' @examples
#'
#' \donttest{
#'
#' one.cmt <- function() {
#'  ini({
#'    ## You may label each parameter with a comment
#'    tka <- 0.45 # Log Ka
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
#' fit <- nlmixr(one.cmt, theo_sd, est="focei")
#'
#' head(vpcSim(fit, pred=TRUE))
#'
#' }
vpcSim <- function(object, ..., keep=NULL, n=300,
                   pred=FALSE, seed=1009, nretry=50, minN=10,
                   normRelated=TRUE) {
  checkmate::assertIntegerish(minN, len=1, any.missing=FALSE, lower=2)
  checkmate::assertLogical(pred, len=1, any.missing=FALSE)
  checkmate::assertIntegerish(nretry, len=1, any.missing=FALSE, lower=0)
  checkmate::assertLogical(normRelated, len=1, any.missing=FALSE)
  checkmate::assertCharacter(keep, null.ok=TRUE, pattern="^[.]*[a-zA-Z]+[a-zA-Z0-9._]*$")
  checkmate::assertIntegerish(seed)
  nlmixr2global$finalUiCompressed <- FALSE
  on.exit(nlmixr2global$finalUiCompressed <- TRUE)
  set.seed(seed)
  .si <- object$simInfo
  .env <- new.env(parent=emptyenv())
  .env$ui <- object$ui
  .env$data <- object$origData
  suppressMessages(.preProcessHooksRun(.env, "rxSolve"))
  .si$object <- eval(.getSimModel(.env$ui, hideIpred=FALSE))
  .w <- which(names(.si) == "rx")
  .si <- .si[-.w]
  .si$nsim <- n
  .si <- c(.si, list(...))
  .pt <- proc.time()
  .si$keep <- unique(c(keep, "nlmixrRowNums"))
  .data <- .si$events
  .data$nlmixrRowNums <- seq_along(.data[, 1])
  if (normRelated) {
    .ui <- .env$ui
    .predDf <-.ui$predDf
    if (all(.predDf$dist %fin% c("norm", "dnorm","t", "cauchy"))) {
    } else {
      if (is.null(.data$CMT)) {
        .ds <- object$dataSav
        .data$nlmixrRowNums <- seq_along(.data[,1])
        .ds <- .ds[, c("CMT", "nlmixrRowNums")]
        .data <- merge(.data, .ds, by ="nlmixrRowNums")
        .data <- .data[order(.data$nlmixrRowNums),]
      }
      .lst <- .Call(`_nlmixr2est_filterNormalLikeAndDoses`,
                    .data$CMT, .predDf$distribution, .predDf$cmt)
      .lst$nlmixrRowNums <- .data[.lst$filter, "nlmixrRowNums"]
      if (.lst$nnorm == 0L) {
        stop("need normal data for vpcSim (or use normRelated=FALSE)")
      }
      .data <- .data[.lst$filter, ]
    }
  }
  .si$events <- .data
  .si$thetaMat <- NULL
  .si$dfSub <- NULL
  .si$dfObs <- NULL
  .si$returnType <- "data.frame.TBS"
  .sim <- do.call(rxode2::rxSolve, .si)
  if (!("sim.id" %fin% names(.sim))) {
    .sim2$sim.id <- 1
  }
  # now look for how many have missing values
  .w <- which(is.na(.sim$ipred))
  .nretry <- 0
  while (length(.w) > 0 && .nretry < nretry) {
    .w <- which(is.na(.sim$ipred))
    .simIds <- unique(.sim$sim.id[.w])
    .sim <- .sim[!(.sim$sim.id %fin% .simIds),, drop = FALSE]
    if (length(.sim$sim.id) == 0) {
      warning("when filtering for simulations, could not find any though some were flagged, be cautious with results",
              call.=FALSE)
      break
    }
    .sim$sim.id <- as.integer(factor(.sim$sim.id))
    .mx <- max(.sim$sim.id)
    .si$nsim <- n - .mx
    .adjust <- FALSE
    if (.si$nsim < minN) {
      .si$nsim <- minN
      .adjust <- TRUE
    }
    .sim2 <- do.call(rxode2::rxSolve, .si)
    if (!("sim.id" %fin% names(.sim2))) {
      .sim2$sim.id <- .simIds[1]
    }
    if (.adjust) {
      # Select simulations without NA ipreds in them
      .w <- which(is.na(.sim2$ipred))
      .simIds <- unique(.sim2$sim.id[.w])
      .allSimIds <- unique(.sim2$sim.id)
      .simIds <- .allSimIds[!(.allSimIds %fin% .simIds)]
      .simIds <- .simIds[seq_len(min(n - .mx, length(.simIds)))]
      .sim2 <- .sim2[.sim2$sim.id %fin% .simIds, ]
    }
    .sim2$sim.id <- .sim2$sim.id + .mx
    .sim <- rbind(.sim, .sim2)
    .w <- which(is.na(.sim$ipred))
    .nretry <- .nretry + 1
  }
  if (.nretry != 0) {
    if (length(.w) == 0) {
      warning("'NA' values in vpc or npde simulation, re-simulated until all simulations were successful",
              call.=FALSE)
    } else {
      warning("'NA' values in vpc or npde simulation",
              call.=FALSE)
    }
  }
  if (pred) {
    .si$nsim <- n # restore for pred
    .si2 <- .si
    .si2$params <- c(
      .si$params, setNames(rep(0, dim(.si$omega)[1]), dimnames(.si$omega)[[2]]),
      setNames(rep(0, dim(.si$sigma)[1]), dimnames(.si$sigma)[[2]]))
    .si2$omega <- NULL
    .si2$sigma <- NULL
    .si2$returnType <- "data.frame"
    .si2$nStud <- 1
    .si2$nsim <- NULL
    nlmixr2global$lastPredSimulationInfo <- .si2
    .sim2 <- do.call(rxode2::rxSolve, .si2)
    .sim$pred <- .sim2$sim
  }
  .sim <- vpcNameDataCmts(object, .sim)
  .cls <- c("nlmixr2vpcSim", class(.sim))
  .fit <- object
  .cls0 <- c("rxHidden", class(.fit))
  attr(.cls0, ".foceiEnv") <- attr(class(.fit), ".foceiEnv")
  class(.fit) <- .cls0
  attr(.cls, "fit") <- .fit
  class(.sim) <- .cls
  return(.sim)
}
#' Name the data and compartments
#'
#' @param object nlmixr2 fit object
#' @param data dataset to name `dvid` and `cmt` columns to correspond with the model
#' @return Updated object/data
#' @author Matthew L. Fidler
#' @keywords internal
#' @export
vpcNameDataCmts <- function(object, data) {
  nlmixr2global$finalUiCompressed <- FALSE
  on.exit(nlmixr2global$finalUiCompressed <- TRUE)
  .wdvid <- which(tolower(names(data)) == "dvid")
  .wcmt <- which(tolower(names(data)) == "cmt")
  .env <- new.env(parent=emptyenv())
  .env$ui <- object$ui
  .env$data <- object$origData
  suppressMessages(.preProcessHooksRun(.env, "rxSolve"))
  .env$ui <- rxode2::rxUiDecompress(.env$ui)
  .info <- get("predDf", .env$ui)
  if (is.null(.info)) {
    return(invisible(data))
  }
  if (length(.info$cond) == 1L) {
    return(invisible(data))
  }
  .state0 <- rxode2::rxModelVars(.env$ui)$state
  .maxCmt <- max(.info$cmt)
  .cmtF <- character(max(length(.state0), .maxCmt))
  .dvidF <- character(length(.info$cmt))
  for (.i in seq_along(.state0)) {
    .cmtF[.i] <- .state0[.i]
  }
  for (.i in seq_along(.info$cmt)) {
    .cmtF[.info$cmt[.i]] <- paste(.info$cond[.i])
    .dvidF[.i] <- paste(.info$cond[.i])
  }
  if (length(.wdvid) == 1) {
    if (inherits(data[[.wdvid]], "numeric")) {
      data[[.wdvid]] <- as.integer(data[[.wdvid]])
    }
    if (!inherits(data[[.wdvid]], "factor") &&
          inherits(data[[.wdvid]], "integer")) {
      .tmp <- data[[.wdvid]]
      attr(.tmp, "levels") <- .dvidF
      attr(.tmp, "class") <- "factor"
      data[[.wdvid]] <- .tmp
    } else if (inherits(data[[.wdvid]], "character")) {
      data[[.wdvid]] <- factor(data[[.wdvid]], .dvidF)
    }
  }
  if (length(.wcmt) == 1) {
    if (inherits(data[[.wcmt]], "numeric")) {
      data[[.wcmt]] <- as.integer(data[[.wcmt]])
    }
    if (!inherits(data[[.wcmt]], "factor") &&
          inherits(data[[.wcmt]], "integer")) {
      .tmp <- data[[.wcmt]]
      attr(.tmp, "levels") <- .cmtF
      attr(.tmp, "class") <- "factor"
      data[[.wcmt]] <- .tmp
    } else if (inherits(data[[.wcmt]], "character")) {
      data[[.wcmt]] <- factor(data[[.wcmt]], .cmtF)
    }
  }
  data
}

#' Expand a VPC simulation
#'
#'
#' @param object nlmixr fit object
#' @param sim vpc simulation object
#' @param extra extra data from original fit to add
#' @param fullData is the full data (possibly modified); This is used
#'   for the vpc tad calculation
#' @return Expanded data frame with extra pieces added
#' @author Matthew L. Fidler
#' @export
#' @keywords internal
vpcSimExpand <- function(object, sim, extra, fullData=NULL) {
  if (is.null(extra)) return(sim)
  if (is.null(fullData)) {
    .fullData <- object$origData
  } else {
    .fullData <- fullData
  }
  .fullData$nlmixrRowNums <- seq_along(.fullData[, 1])
  .extra <- extra[extra %fin% names(.fullData)]
  .extra <- extra[!(extra %fin% names(sim))]
  if (length(.extra) == 0) return(sim)
  .wid <- which(tolower(names(.fullData)) == "id")
  names(.fullData)[.wid] <- "ID"
  .sim <- sim
  .wid <- which(tolower(names(.sim)) == "id")
  names(.sim)[.wid] <- "ID"
  .ret <- merge(.fullData, .sim, by=c("ID", "nlmixrRowNums"))
  .w <- which(names(.ret) == "nlmixrRowNums")
  vpcNameDataCmts(object, .ret[, -.w])
}

#' Get the least prediction simulation information for VPC
#'
#' @return The last prediction simulation from the `vpcSim` function (data.frame)
#'
#' @export
#' @keywords internal
.nlmixr2estLastPredSimulationInfo <- function() {
  nlmixr2global$lastPredSimulationInfo
}
