#' @export
rxUiGet.transUE <- function(x, ...) {
  .ui <- x[[1]]
  .iniDf <- .ui$iniDf
  .w <- which(.iniDf$neta1 == .iniDf$neta2)
  if (length(.w) == 0L) return(NULL)
  .n <- .iniDf$name[.w]
  .muRef <- .ui$muRefDataFrame
  vapply(.n, function(cur) {
    .w <- which(.muRef$eta == cur)
    if (length(.w) == 0L) return(cur)
    .muRef$theta[.w]
  }, character(1), USE.NAMES = TRUE)
}
attr(rxUiGet.transUE, "rstudio")  <- c(eta.ka="tka")

#' Get the parameter values for uninformative eta calculation
#'
#' @param name name of the eta parameter
#' @param ui rxode2 user interface object
#' @param pm plus or minus values
#' @param plus What type of deviate is this:
#'
#' - `plus=TRUE`: high estimate
#'
#' - `plus=FALSE`: low estimate
#'
#' - `plus=NA`: middle estimate
#'
#' @param saem is this a saem-style mu-referenced model?
#'
#' @return the named parameter value for plus/minus/mid calculation
#'
#' @noRd
#'
#' @author Matthew L. Fidler
.getMuValForUE <- function(name, trans, ui, pm, plus=TRUE, saem=TRUE, retName=FALSE) {
  if (!saem || name %fin% ui$nonMuEtas) {
    if (retName) return(name)
    if (is.na(plus)) {
      return(setNames(0, name))
    } else if (plus) {
      return(pm[name])
    } else {
      return(-pm[name])
    }
  }
  .n2 <- trans[name]
  .v0 <- ui$theta[.n2]
  if (retName) return(.n2)
  if (is.na(plus)) {
    setNames(.v0, .n2)
  } else if (plus) {
    setNames(.v0 + pm[name], .n2)
  } else {
    setNames(.v0 - pm[name], .n2)
  }
}
#' Get the low, middle, and hi subject values based on etas
#'
#' @param ui user interface function
#' @param data data.frame that will be used for fitting
#' @param alpha The alpha value to scale the interval from -1 to 1
#'   when choosing the quadrature points
#' @param saem boolean indicating if the model is a mu-referenced.
#' @param q plus or minus quadrature
#' @return a list with:
#'
#' - `trans` the translation of the parameters
#'   from etas to model pars (useful in mu-referenced modeling
#'   algorithms like `saem`)
#'
#' - `dat` rxode2 translated dataset (using etTrans)
#'
#' - `param` full list of parameters to solve for
#'
#' @noRd
#' @author Matthew L. Fidler
.uninformativeEtasExpand <- function(ui, data, trans, alpha=0.05,
                               saem=TRUE, q=sqrt(3/5)) {
  .trans <- rxode2::etTrans(data, ui)
  .lst <- attr(class(.trans), ".rxode2.lst")
  .n <- .lst$nid
  .iniDf <- ui$iniDf
  .eta <- .iniDf[which(.iniDf$neta1 == .iniDf$neta2), ]

  .pm <- setNames(qnorm(1 - alpha / 2) * sqrt(.eta$est) * q, .eta$name)

  .neta <- length(.eta$name)

  .nn <- trans

  .p <- do.call("rbind", lapply(seq_len(.n), function(i) {
    as.data.frame(t(vapply(names(.pm), .getMuValForUE, ui=ui, pm=.pm, plus=TRUE, saem=saem, trans=trans,
                           double(1), USE.NAMES=FALSE)))
  }))
  names(.p) <- .nn

  .m <- do.call("rbind", lapply(seq_len(.n), function(i) {
    as.data.frame(t(vapply(names(.pm), .getMuValForUE, ui=ui, pm=.pm, plus=FALSE, saem=saem, trans=trans,
                           double(1), USE.NAMES = FALSE)))
  }))
  names(.m) <- .nn

  .z <- do.call("rbind", lapply(seq_len(.n), function(i) {
    as.data.frame(t(vapply(names(.pm), .getMuValForUE, ui=ui, pm=.pm, plus=NA, saem=saem, trans=trans,
                           double(1), USE.NAMES = FALSE)))
  }))
  names(.z) <- .nn
  .env <- new.env(parent=emptyenv())
  .env$sim.id <- 1L
  .etas <- do.call("rbind", lapply(.nn, function(nm) {
    .df <- rbind(.m, .z, .p)
    for (cur in .nn) {
      if (cur == nm) next
      .df[[cur]] <- c(.z[[cur]], .z[[cur]], .z[[cur]])
    }
    .df$rxW <- which(nm == .nn)
    .df$rxPmz <- c(rep(-1L, .n), rep(0L, .n), rep(1L, .n))
    .df$id <- c(seq_len(.n), seq_len(.n), seq_len(.n))
    .df$sim.id <- .env$sim.id
    .tmp <- rep(.env$sim.id, .n)
    .env$sim.id <- .env$sim.id + 1L
    .tmp <- c(.tmp, rep(.env$sim.id, .n))
    .env$sim.id <-.env$sim.id  + 1L
    .tmp <- c(.tmp, rep(.env$sim.id, .n))
    .env$sim.id <- .env$sim.id + 1L
    .df$sim.id <- .tmp
    .df
  }))

  .fullN <- unique(c(names(.etas), names(ui$theta)))
  .full <- .etas
  for (n in .fullN) {
    if (n %fin% names(.full)) next
    .full[[n]] <- ui$theta[n]
  }

  list(trans=setNames(names(.pm), .nn), dat=.trans, param=.full, n=.n,
       neta=.neta)
}

.uninformativeEtas <- function(ui, handleUninformativeEtas=TRUE, data, model, alpha=0.05,
                               saem=TRUE, q=sqrt(3/5),
                               rxControl=NULL, tol=1e-7) {
  .rxControl <- rxControl
  if (is.null(rxControl)) .rxControl <- rxode2::rxControl()
  if (saem) {
    ui <- rxode2::assertRxUi(ui)
    .trans <- rxUiGet.transUE(list(ui))
    .pars <- .uninformativeEtasExpand(ui, data, trans=.trans, alpha=alpha, saem=TRUE, q=q)
    if (!handleUninformativeEtas) {
      .lst <- attr(class(.pars$dat), ".rxode2.lst")
      .n <- .lst$nid
      .mat <- matrix(rep(1L, .n*length(.pars$trans)),
                     nrow=.n,
                     ncol=length(.pars$trans),
                     dimnames=list(NULL, .pars$trans))
      return(.mat)
    }
    .minfo("calculate uninformed etas")
    .rxControl$returnType <- setNames(2L, "data.frame")
    # Get the predictions at +- etas
    .lst <- attr(class(.trans), ".rxode2.lst")

    .val <- do.call(rxode2::rxSolve, c(list(model, .pars$param, data), .rxControl))
    .val$id <- as.integer(.val$id)
    .ind <- .pars$param[,c("id", "sim.id", "rxW", "rxPmz")]
    .val <- merge(.val[, c("id", "sim.id", "rx_pred_")], .ind)
    .env <- new.env(parent=emptyenv())
    .env$nid <- .pars$n
    .env$neta <- .pars$neta
    .env$simId <- .val[["sim.id"]]
    .env$id <- .val[["id"]]
    .env$val <- .val[["rx_pred_"]]
    .env$w <- .val[["rxW"]]
    .env$pm <- .val[["rxPmz"]]
    .env$tol <- tol
    .mat <- .Call(`_nlmixr2est_uninformativeEta`, .env)
    dimnames(.mat) <- list(NULL, names(.pars$trans))
    .minfo("done")
    .mat
  }
}
