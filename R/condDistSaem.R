#' Eta draws from a SAEM chain history
#'
#' @param saem the run's `saemFit` (`fit$saem`)
#' @param phiM the run's `N x nmc x niter x nphi` chain history
#' @param keep iterations of `phiM` to keep
#' @param etaTrans `ui$saemOmegaTrans`, the kernel eta column of each ui eta
#' @param etaNames the ui eta names
#' @return `N x neta x ndraw` array of eta draws (`NA` where SAEM does not
#'   sample the eta), or `NULL` when the shapes do not match
#' @noRd
.condDistSaemDraws <- function(saem, phiM, keep, etaTrans, etaNames) {
  .i1 <- attr(saem, "saem.cfg")$i1 + 1L
  .mprior <- saem$mprior_phi
  .ue <- saem$etaUe
  .d <- dim(phiM)
  if (
    length(.d) != 4L ||
      length(keep) == 0L ||
      min(keep) < 1L ||
      max(keep) > .d[3] ||
      length(.i1) == 0L ||
      !is.matrix(.mprior) ||
      !identical(dim(.mprior), .d[c(1L, 4L)]) ||
      !is.matrix(.ue) ||
      !identical(dim(.ue), c(.d[1], length(.i1))) ||
      length(etaTrans) != length(etaNames) ||
      anyNA(etaTrans) ||
      max(etaTrans) > length(.i1)
  ) {
    return(NULL)
  }
  # N x nphi1 x nmc x nkeep, so the per-subject prior mean and mask recycle
  .draw <- aperm(phiM[,, keep, .i1, drop = FALSE], c(1L, 4L, 2L, 3L))
  .draw <- .draw - as.vector(.mprior[, .i1, drop = FALSE])
  .draw[array(as.vector(.ue) == 0, dim(.draw))] <- NA_real_
  dim(.draw) <- c(.d[1], length(.i1), .d[2] * length(keep))
  .draw[, etaTrans, , drop = FALSE]
}

#' A held SAEM control: MCMC only, at the supplied population parameters
#'
#' Every population parameter is held (`.saemHoldCfg()`); with a chain state
#' the chains continue from it (`.saemWarmCfg()`).  No covariance, tables,
#' likelihood or conditional-mode search.
#' @param control a `saemControl()`
#' @param niter MCMC iterations
#' @param state `NULL`, or `.saemChainState()` of a SAEM fit
#' @return `saemControl()` object
#' @noRd
.condDistSaemControl <- function(control, niter, state = NULL) {
  .ctl <- control
  .ctl$mcmc$niter <- c(as.integer(niter), 0L)
  .ctl$covMethod <- ""
  .ctl$covMethodDeferred <- NA_character_
  .ctl$covFallback <- list()
  .ctl$calcTables <- FALSE
  .ctl$logLik <- FALSE
  .ctl$ebe <- "saem"
  .ctl$iterPrintControl <- saemControl(print = 0L)$iterPrintControl
  .ctl$saemHoldPar <- TRUE
  .ctl$saemWarmState <- state
  class(.ctl) <- "saemControl"
  .ctl
}

#' SAEM conditional-distribution draws of a fit
#'
#' A SAEM fit's own last iterations give the draws when the last half of its
#' second phase holds enough, so nothing is simulated.  Otherwise the SAEM MCMC kernels are run with every
#' population parameter held at the fit's estimates: for a SAEM fit continuing
#' its chains, for any other fit from a cold start, discarding `burn`
#' iterations.
#' @param fit nlmixr2 fit
#' @param nsamp draws per subject
#' @param burn discarded iterations of a held run
#' @param seed seed of a held run
#' @return list: `draws` (`N x neta x nsamp`) and `method`, `"saem"` for the
#'   fit's own chains or `"mcmc"` for a held run
#' @noRd
.condDistSaem <- function(fit, nsamp, burn, seed) {
  .etaNames <- names(fit$eta)[-1]
  .isSaem <- identical(tryCatch(fit$est, error = function(e) NULL), "saem")
  .lik <- .condDistLikControl(fit)
  .ctl <- if (.isSaem) {
    fit$saemControl
  } else {
    saemControl(seed = seed, addProp = .lik$addProp, rxControl = .lik$rxControl)
  }
  .nmc <- as.integer(.ctl$mcmc$nmc)
  .nit <- as.integer(ceiling(nsamp / .nmc))
  .take <- function(draws) draws[,, seq_len(nsamp), drop = FALSE]
  if (.isSaem) {
    .phiM <- tryCatch(fit$phiM, error = function(e) NULL)
    # the last half of the second phase, where the step size has fallen and
    # the parameters have nearly settled
    .nEm <- as.integer(.ctl$mcmc$niter[2])
    if (is.array(.phiM) && .nit <= .nEm %/% 2L) {
      .last <- dim(.phiM)[3]
      .d <- .condDistSaemDraws(fit$saem, .phiM, .last - rev(seq_len(.nit)) + 1L, fit$ui$saemOmegaTrans, .etaNames)
      if (!is.null(.d)) {
        return(list(draws = .take(.d), method = "saem"))
      }
    }
  }
  .ctl$seed <- seed
  .state <- if (.isSaem) .saemChainState(fit)
  .ctl <- .condDistSaemControl(.ctl, burn + .nit, .state)
  .fit2 <- .nativePinnedRefit(fit, "saem", .ctl, useEtaMat = FALSE)
  if (inherits(.fit2, "try-error")) {
    stop("held SAEM run at the fit's estimates failed; try method = \"imp\"", call. = FALSE)
  }
  if (length(.fit2$env$saemPseudoEtas) > 0L) {
    stop("SAEM needs temporary etas for this model; use method = \"imp\"", call. = FALSE)
  }
  .d <- .condDistSaemDraws(.fit2$saem, .fit2$phiM, burn + seq_len(.nit), .fit2$ui$saemOmegaTrans, .etaNames)
  if (is.null(.d)) {
    stop("held SAEM run has an unexpected chain shape", call. = FALSE)
  }
  list(draws = .take(.d), method = "mcmc")
}
