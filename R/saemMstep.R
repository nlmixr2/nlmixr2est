#' How the SAEM M-step updates each population theta
#'
#' A mu-referenced (phi1) theta is updated in closed form from the sufficient
#' statistics; a theta without an eta (phi0) is refined numerically under
#' `nonMuTheta="regress"` (and always for a general likelihood) or carried
#' stochastically otherwise.  A covariate coefficient follows its parent theta.
#' @param ui rxode2 ui
#' @param cfg `.configsaem()` configuration (0-based `i0`)
#' @param phi1Pooled `TRUE` when the general-likelihood phi1 refinement is set up
#' @return data frame with columns `theta` and `update` (`"closed-form"`,
#'   `"penalized"`, `"regress"`, `"stochastic"` or `"fixed"`)
#' @noRd
.saemMstepTable <- function(ui, cfg, phi1Pooled = FALSE) {
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
  .upd[as.integer(cfg$i0) + 1L] <- .phi0Upd
  .cov <- rxUiGet.saemMuRefCovariateDataFrame(list(ui))
  .cp <- as.character(.cov$covariateParameter)
  .theta <- c(.pars, .cp)
  .upd <- c(.upd, .upd[match(as.character(.cov$theta), .pars)])
  .fix <- ui$iniDf$fix[match(.theta, ui$iniDf$name)]
  # a population penalty refines the closed-form update numerically (penMstep)
  if (length(cfg$penMap) > 0L && .likPenaltyInfo()[["population"]] > 0L) {
    .upd[.upd == "closed-form"] <- "penalized"
  }
  .upd[!is.na(.fix) & .fix] <- "fixed"
  # an eta without a theta owns a phi column but is no population parameter
  .keep <- !(.theta %in% as.character(ui$nonMuEtas))
  data.frame(
    theta = sub("^rxBoundedTr[.]", "", .theta[.keep]),
    update = .upd[.keep],
    stringsAsFactors = FALSE
  )
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
      c("optimize", "nelderMead", "newuoa")[as.integer(cfg$nonMuThetaOptType) + 1L]
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
