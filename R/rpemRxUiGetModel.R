# RPEM likelihood model (design/rpem/13-likelihood-model.md)
#
# The RPEM E-step needs p(Y_i | theta_i) evaluated for an externally supplied
# theta_i = mu-referenced population params + eta, with NO inner eta
# optimization.  This is a hybrid of two existing generators:
#
#   * nlmModel0  (R/nlm.R)   -- direct, no-optimization log-likelihood model
#                               (rxPredLlik=TRUE -> rx_pred_/rx_pred_f_/rx_r_),
#                               but maps every parameter to THETA[] (no eta).
#   * foceiModel0ll (R/focei.R) -- THETA[]+ETA[] parameterization
#                               (.uiGetThetaEta), but built to OPTIMIZE eta.
#
# rpemModel0 == nlmModel0 with the parameter prefix swapped to .uiGetThetaEta so
# etas appear as inputs (supplied per Monte Carlo sample, never optimized).  The
# rest of the pipeline (prune -> symengine -> assemble rxode2 model) mirrors the
# nlm population-log-likelihood pipeline.

#' @export
rxUiGet.rpemModel0 <- function(x, ...) {
  .ui <- rxode2::rxUiDecompress(x[[1]])
  nlmixr2global$rxPredLlik <- TRUE
  on.exit(nlmixr2global$rxPredLlik <- FALSE)
  .predDf <- .ui$predDf
  .save <- .predDf
  .predDf[.predDf$distribution == "norm", "distribution"] <- "dnorm"
  assign(".predDfFocei", .predDf, envir=.ui)
  on.exit(assign("predDf", .save, envir=.ui))
  .ret <- rxode2::rxCombineErrorLines(.ui, errLines=rxGetDistributionFoceiLines(.ui),
                                      prefixLines=.uiGetThetaEta(.ui),
                                      paramsLine=NA,
                                      modelVars=TRUE,
                                      cmtLines=FALSE,
                                      dvidLine=FALSE)
  .ret <- .ret[[2]]
  # Sign note (spec 13 OI-3): keep identical to the proven nlm llik model
  # (rx_pred_ negated) so this is a faithful THETA+ETA analog; the E-step decides
  # how to consume the summed value.  Revisit sign when wiring the E-step.
  .ret <- as.call(c(quote(`{`),
                    lapply(seq_along(.ret)[-1], function(i) {
                      .ret[[i]]
                    }),
                    list(str2lang("rx_pred_ <- -rx_pred_"))))
  as.call(c(list(quote(`rxModelVars`)), .ret))
}
attr(rxUiGet.rpemModel0, "rstudio") <- quote(rxModelVar({}))

#' Prune if/else branches of the RPEM log-likelihood model
#' @param x rxode2 UI object
#' @return pruned rxNorm string
#' @noRd
.rpemPrune <- function(x) {
  .x <- x[[1]]
  .x <- .x$rpemModel0[[-1]]
  .env <- new.env(parent = emptyenv())
  .env$.if <- NULL
  .env$.def1 <- NULL
  .malert("pruning branches ({.code if}/{.code else}) of RPEM log-likelihood model...")
  .ret <- rxode2::.rxPrune(.x, envir = .env,
                           strAssign = rxModelVars(x[[1]])$strAssign)
  .mv <- rxode2::rxModelVars(.ret)
  if (rxode2::.rxIsLinCmt() == 1L) {
    .vars <- c(.mv$params, .mv$lhs, .mv$slhs)
    .mv <- rxode2::.rxLinCmtGen(length(.mv$state), .vars)
  }
  .msuccess("done")
  rxode2::rxNorm(.mv)
}

#' @export
rxUiGet.loadPruneRpem <- function(x, ...) {
  .loadSymengine(.rpemPrune(x), promoteLinSens = FALSE)
}
attr(rxUiGet.loadPruneRpem, "rstudio") <- emptyenv()

#' @export
rxUiGet.rpemParams <- function(x, ...) {
  # THETA[] + ETA[] + covariates + DV.  Unlike nlmParams (THETA + DV), RPEM
  # carries the eta dimension so per-subject/per-sample etas are supplied inputs.
  .ui <- x[[1]]
  .str <- .uiGetThetaEtaParams(.ui, str=FALSE)
  .pars <- as.character(.str)[-1]
  .pars <- .pars[nzchar(.pars)]
  paste0("params(", paste(c(.pars, "DV"), collapse=","), ")")
}
attr(rxUiGet.rpemParams, "rstudio") <- "params()"

#' @export
rxUiGet.rpemRxModel <- function(x, ...) {
  .s <- rxUiGet.loadPruneRpem(x, ...)
  .isMatExp <- isTRUE(.rxInjectMatExpDdt(.s))
  .prd <- get("rx_pred_", envir = .s)
  .prd <- paste0("rx_pred_=", rxode2::rxFromSE(.prd))
  .ddt <- .s$..ddt
  if (is.null(.ddt)) .ddt <- ""
  .lhs <- character(0)
  if (.isMatExp) {
    .lhs <- .s$..lhs
    if (is.null(.lhs)) .lhs <- character(0)
  }
  .fr <- .nlmGetFRLines(.s)
  .ret <- paste(c(
    .lhs,
    .ddt,
    .prd,
    .fr$f_line,
    .fr$r_line,
    ""
  ), collapse = "\n")
  .sumProd <- rxode2::rxGetControl(x[[1]], "sumProd", FALSE)
  .optExpression <- rxode2::rxGetControl(x[[1]], "optExpression", TRUE)
  if (.sumProd) {
    .malert("stabilizing round off errors in RPEM log-likelihood model...")
    .ret <- rxode2::rxSumProdModel(.ret)
    .msuccess("done")
  }
  if (.optExpression) {
    .ret <- rxode2::rxOptExpr(.ret, "RPEM log-likelihood model")
    .msuccess("done")
  }
  .cmt <- rxUiGet.foceiCmtPreModel(x, ...)
  .interp <- rxUiGet.interpLinesStr(x, ...)
  if (.interp != "") {
    .cmt <- paste0(.cmt, "\n", .interp)
  }
  list(predOnly=rxode2::rxode2(paste(c(rxUiGet.rpemParams(x, ...), .cmt,
                                       .ret, .foceiToCmtLinesAndDvid(x[[1]])), collapse="\n")))
}
attr(rxUiGet.rpemRxModel, "rstudio") <- "params()"
