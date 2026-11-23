#' How the SAEM M-step updates each population theta
#'
#' A mu-referenced (phi1) theta is updated in closed form from the sufficient
#' statistics; a theta without an eta (phi0) is refined numerically under
#' `nonMuTheta="regress"` (and always for a general likelihood) or carried
#' stochastically otherwise.  A covariate coefficient follows its parent theta.
#' @param ui rxode2 ui
#' @param cfg `.configsaem()` configuration (0-based `i0`)
#' @param phi1Pooled `TRUE` when the general-likelihood phi1 refinement is set up
#' @param pseudoEtas saem's temporary etas (`.saemPseudoEtaRecord()`)
#' @return data frame with columns `theta` and `update` (`"closed-form"`,
#'   `"penalized"`, `"regress"`, `"stochastic"` or `"fixed"`)
#' @noRd
.saemMstepTable <- function(ui, cfg, phi1Pooled = FALSE, pseudoEtas = character(0)) {
  # phi columns in kernel order (covariate coefficients are not columns)
  .pars <- as.character(rxUiGet.saemParamsToEstimateCov(list(ui)))
  .general <- identical(as.integer(cfg$distribution), 4L)
  .nMix <- if (is.null(cfg$nMix)) 1L else as.integer(cfg$nMix)
  .phi0Upd <- if (.general || (identical(cfg$nonMuTheta, "regress") && .nMix == 1L)) {
    "regress"
  } else {
    "stochastic"
  }
  .upd <- rep(if (.general && phi1Pooled) "regress" else "closed-form", length(.pars))
  .i0 <- as.integer(cfg$i0) + 1L
  .upd[.i0] <- .phi0Upd
  # a population penalty refines the phi1 closed-form update (penMstep); the
  # kernel's penalty map covers every estimated theta
  .pen <- .likPenaltyInfo()
  if (.pen[["n"]] > 0L && .pen[["population"]] > 0L) {
    .upd[.upd == "closed-form"] <- "penalized"
  }
  .cov <- rxUiGet.saemMuRefCovariateDataFrame(list(ui))
  # a shared-eta mixture component theta (intercept only, its eta a phi1
  # column) follows its eta's recentering, as mixShiftResolve() decides
  .comp <- as.integer(cfg$mixShiftComp)
  .eta <- as.integer(cfg$mixShiftEta)
  if (
    .nMix > 1L &&
      identical(as.integer(cfg$mixSampleMethod), 0L) &&
      length(.comp) == length(.pars) &&
      length(.eta) == length(.pars)
  ) {
    .own <- .i0[
      .comp[.i0] >= 1L &
        .comp[.i0] <= .nMix &
        (.eta[.i0] - 1L) %in% as.integer(cfg$i1) &
        !(.pars[.i0] %in% .cov$theta)
    ]
    .upd[.own] <- "closed-form"
  }
  # a theta given a temporary eta is sampled like an eta
  .upd[.pars %in% .saemPseudoEtaThetaNames(ui, pseudoEtas)] <- "stochastic"
  .theta <- c(.pars, as.character(.cov$covariateParameter))
  .upd <- c(.upd, .upd[match(as.character(.cov$theta), .pars)])
  .fix <- ui$iniDf$fix[match(.theta, ui$iniDf$name)]
  .upd[!is.na(.fix) & .fix] <- "fixed"
  # free phi0 coordinates of the regress refinement (covariates are not)
  .nFree0 <- sum(.upd[seq_along(.pars)][.i0] == "regress")
  # an eta without a theta owns a phi column but is no population parameter
  .keep <- !(.theta %in% as.character(ui$nonMuEtas))
  .theta <- sub("^rxBoundedTr[.]", "", .theta[.keep])
  .upd <- .upd[.keep]
  # collapsed IOV: per-occasion rx.<theta>.<level> columns are one theta
  .iov <- .uiIovEnv$iovCollapsed
  if (!is.null(.iov)) {
    for (.t in as.character(.iov$pars$theta)) {
      .theta[.theta %in% paste0("rx.", .t, ".", .iov$levels)] <- .t
    }
    .d <- !duplicated(.theta)
    .theta <- .theta[.d]
    .upd <- .upd[.d]
  }
  .ret <- data.frame(theta = .theta, update = .upd, stringsAsFactors = FALSE)
  attr(.ret, "nFree0") <- .nFree0
  .ret
}

#' One-line note of the SAEM M-step split, or `NULL` when every theta is
#' closed form or fixed
#' @param tbl `.saemMstepTable()` result
#' @param cfg `.configsaem()` configuration
#' @return character note or `NULL`
#' @noRd
.saemMstepNote <- function(tbl, cfg) {
  .reg <- tbl$theta[tbl$update == "regress"]
  .sto <- tbl$theta[tbl$update == "stochastic"]
  .pen <- tbl$theta[tbl$update == "penalized"]
  if (length(.reg) == 0L && length(.sto) == 0L && length(.pen) == 0L) {
    return(NULL)
  }
  .cf <- tbl$theta[tbl$update == "closed-form"]
  .parts <- list()
  if (length(.cf) > 0L) {
    .parts <- c(.parts, list(c("closed form: ", .cf)))
  }
  if (length(.pen) > 0L) {
    .parts <- c(.parts, list(c("closed form+penalty: ", .pen)))
  }
  if (length(.reg) > 0L) {
    .opt <- if (identical(as.integer(cfg$distribution), 4L)) {
      "bobyqa"
    } else {
      # a single free coordinate always takes the 1-d search
      .nFree <- attr(tbl, "nFree0")
      if (is.null(.nFree)) {
        .nFree <- length(.reg)
      }
      .type <- if (.nFree <= 1L) 0L else as.integer(cfg$nonMuThetaOptType)
      c("optimize", "nelderMead", "newuoa")[.type + 1L]
    }
    .parts <- c(.parts, list(c(paste0("numerical (", .opt, "): "), .reg)))
  }
  if (length(.sto) > 0L) {
    .parts <- c(.parts, list(c("stochastic: ", .sto)))
  }
  # share the line budget left after the labels evenly between the classes
  .pre <- vapply(.parts, function(p) p[1], character(1))
  .share <- floor((74L - sum(nchar(.pre)) - 2L * (length(.parts) - 1L)) / length(.parts))
  paste(
    vapply(
      .parts,
      function(p) {
        paste0(p[1], .vaeTruncList(p[-1], budget = .share))
      },
      character(1)
    ),
    collapse = "; "
  )
}
