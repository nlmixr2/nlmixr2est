.lastPredSimulationInfo <- NULL # to get observation dataset with pred attached for pred_corr
#' VPC simulation
#'
#' @param object This is the nlmixr2 fit object
#' @param ... Other arguments sent to `rxSolve()`
#' @param keep Keep character vector
#' @param n Number of simulations
#' @param pred Should predictions be added to the simulation
#' @param seed Seed to set for the VPC simulation
#' @param nretry Number of times to retry the simulation if there is
#'   NA values in the simulation
#' @return data frame of the VPC simulation
#' @author Matthew L. Fidler
#' @export
#' @examples
#'
#' \donttest{
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
vpcSim <- function(object, ..., keep=NULL, n=300, pred=FALSE, seed=1009, nretry=50) {
  set.seed(seed)
  .si <- object$simInfo
  .si$object <- eval(.getSimModel(object, hideIpred=FALSE))
  .w <- which(names(.si) == "rx")
  .si <- .si[-.w]
  .si$nsim <- n
  .si <- c(.si, list(...))
  .pt <- proc.time()
  .si$keep <- unique(c(keep, "nlmixrRowNums"))
  .data <- .si$events
  .data$nlmixrRowNums <- seq_along(.data[, 1])
  .si$events <- .data
  .si$thetaMat <- NULL
  .si$dfSub <- NULL
  .si$dfObs <- NULL
  .si$returnType <- "data.frame.TBS"
  .sim <- do.call(rxode2::rxSolve, .si)
  # now look for how many have missing values
  .w <- which(is.na(.sim$ipred))
  .nretry <- 0
  while (length(.w) > 0 && .nretry < nretry) {
    .w <- which(is.na(.sim$ipred))
    .simIds <- unique(.sim$sim.id[.w])
    .sim <- .sim[!(.sim$sim.id %in% .simIds), ]
    .sim$sim.id <- as.integer(factor(.sim$sim.id))
    .mx <- max(.sim$sim.id)
    .si$n <- n - .mx
    .sim2 <- do.call(rxode2::rxSolve, .si)
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
    assignInMyNamespace(".lastPredSimulationInfo", .si2)
    .sim2 <- do.call(rxode2::rxSolve, .si2)
    .sim$pred <- .sim2$sim
  }
  .sim <- vpcNameDataCmts(object, .sim)
  class(.sim) <- c("nlmixr2vpcSim", class(.sim))
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
  .wdvid <- which(tolower(names(data)) == "dvid")
  .wcmt <- which(tolower(names(data)) == "cmt")
  .info <- get("predDf", object$ui)
  if (is.null(.info)) {
    return(invisible(data))
  }
  .state0 <- rxode2::rxModelVars(object$ui)$state
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
    if (inherits(data[[.wdvid]], "numeric")) data[[.wdvid]] <- as.integer(data[[.wdvid]])
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
    if (inherits(data[[.wcmt]], "numeric")) data[[.wcmt]] <- as.integer(data[[.wcmt]])
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
#' @return Expanded data frame with extra pieces added
#' @author Matthew L. Fidler
#' @export
#' @keywords internal
vpcSimExpand <- function(object, sim, extra) {
  if (is.null(extra)) return(sim)
  .fullData <- object$origData
  .fullData$nlmixrRowNums <- seq_along(.fullData[, 1])
  .extra <- extra[extra %in% names(.fullData)]
  .extra <- extra[!(extra %in% names(sim))]
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
  .lastPredSimulationInfo
}

#' Setup Observation data for VPC
#'
#' @param fit nlmixr2 fit
#' @param data replacement data
#' @return List with `namesObs`, `namesObsLower`, `obs` and `obsCols`
#' @author Matthew L. Fidler
#' @noRd
.vpcUiSetupObservationData <- function(fit, data=NULL) {
  if (!is.null(data)) {
    .obs <- data
  } else {
    .obs <- fit$origData
  }
  .no <- names(.obs)
  .nol <- tolower(.no)
  .wo <- which(.nol == "id")
  if (length(.wo) != 1) {
    stop("cannot find 'id' in original dataset",
         call.=FALSE)
  }
  .obsCols <- list(id=.no[.wo])
  .wo <- which(.nol == "dv")
  if (length(.wo) != 1) {
    stop("cannot find 'dv' in original dataset",
         call.=FALSE)
  }
  .obsCols <- c(.obsCols,
                list(dv=.no[.wo]))
  .wo <- which(.nol == "time")
  if (length(.wo) != 1) {
    stop("cannot find 'time' in original dataset",
         call.=FALSE)
  }
  .obsCols <- c(.obsCols,
                list(idv=.no[.wo]))
  list(namesObs=.no,
       namesObsLower=.nol,
       obs=.obs,
       obsCols=.obsCols)
}
