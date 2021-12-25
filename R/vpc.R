#' VPC simulation
#'
#' @param object This is the nlmixr2 fit object
#' @param ... Other arguments sent to `nlmixr2Sim()`
#' @param keep Keep character vector
#' @param n Number of simulations
#' @param pred Should predicitions be added to the simulation
#' @param seed Seed to set for the VPC simulation
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
vpcSim <- function(object, ..., keep=NULL, n=300, pred=FALSE, seed=1009) {
  set.seed(seed)
  .si <- object$simInfo
  .si$object <- object
  .si$nsim <- n
  .si <- c(.si, list(...))
  .si$modelName <- "VPC"
  .pt <- proc.time()
  .si$dfObs <- 0
  .si$dfSub <- 0
  .si$thetaMat <- NA
  .si$keep <- keep
  .si$returnType <- "data.frame"

  rxode2::.setWarnIdSort(FALSE)
  on.exit(rxode2::.setWarnIdSort(TRUE))
  .sim <- do.call("nlmixr2Sim", .si)
  if (pred) {
    .si2 <- .si
    .si2$modelName <- "Pred (for pcVpc)"
    .si2$params <- c(
      .si$params, setNames(rep(0, dim(.si$omega)[1]), dimnames(.si$omega)[[2]]),
      setNames(rep(0, dim(.si$sigma)[1]), dimnames(.si$sigma)[[2]])
    )
    .si2$omega <- NA
    .si2$sigma <- NA
    .si2$returnType <- "data.frame"
    .si2$nStud <- 1
    .si2$nsim <- NULL
    .sim2 <- do.call("nlmixr2Sim", .si2)
    .sim$pred <- .sim2$sim
  }
  return(.sim)
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

#' VPC based on ui model
#'
#' @param fit nlmixr2 fit object
#' @param data this is the data to use to augment the VPC fit.  By
#'     default is the fitted data, (can be retrieved by
#'     \code{\link[nlme]{getData}}), but it can be changed by specifying
#'     this argument.
#' @param n Number of VPC simulations.  By default 100
#' @inheritParams vpc::vpc
#' @inheritParams rxode2::rxSolve
#' @param ... Args sent to \code{\link[rxode2]{rxSolve}}
#' @return Simulated dataset (invisibly)
#' @author Matthew L. Fidler
#' @examples
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
#' vpcPlot(fit)
#'
#' }
#'
#' @export
vpcPlot <- function(fit, data = NULL, n = 300, bins = "jenks",
                  n_bins = "auto", bin_mid = "mean",
                  show = NULL, stratify = NULL, pred_corr = FALSE,
                  pred_corr_lower_bnd = 0, pi = c(0.05, 0.95), ci = c(0.05, 0.95),
                  uloq = NULL, lloq = NULL, log_y = FALSE, log_y_min = 0.001,
                  xlab = NULL, ylab = NULL, title = NULL, smooth = TRUE, vpc_theme = NULL,
                  facet = "wrap", labeller = NULL, vpcdb = FALSE, verbose = FALSE, ...,
                  seed=1009) {
  rxode2::rxReq("vpc")
  .ui <- fit$ui
  .obsLst <- .vpcUiSetupObservationData(fit, data)
  .no <- .obsLst$namesObs
  .nol <- .obsLst$namesObsLower
  .obs <- .obsLst$obs
  .obsCols <- .obsLst$obsCols
  # Setup stratify
  .wo <- which(.nol == "cmt")
  .multi <- length(fit$ui$predDf$cmt) > 1
  .w <- which(tolower(stratify) == "cmt")
  if (length(.w) == 0 && .multi && length(.wo) == 1) {
    stratify <- unique(c(stratify, .no[.wo]))
  } else {
    .wo <- which(.nol == "dvid")
    if (length(.wo) == 1 && .multi) {
      stratify <- unique(c(stratify, .no[.wo]))
    }
  }
  # Simulate with VPC
  .sim <- vpcSim(fit, ..., keep=stratify, n=n, pred=pred_corr, seed=seed)
  .simCols <- list(
    id="id",
    dv="sim",
    idv="time")
  if (pred_corr) {
    .simCols <- c(.simCols, list(pred="pred"))
  }

  vpc::vpc_vpc(sim=.sim, sim_cols=.simCols,
               obs=.obs, obs_cols=.obsCols,
               bins=bins, n_bins=n_bins, bin_mid=bin_mid,
               show = show, stratify = stratify, pred_corr = pred_corr,
               pred_corr_lower_bnd = pred_corr_lower_bnd, pi = pi, ci = ci,
               uloq = uloq, lloq = lloq, log_y = log_y, log_y_min = log_y_min,
               xlab = xlab, ylab = ylab, title = title, smooth = smooth, vpc_theme = vpc_theme,
               facet = facet, labeller = labeller, vpcdb = vpcdb, verbose = verbose)
}
