# vaeGrad.R -- analytic outer-gradient M-step for the VAE non-mu thetas
# (vaeControl(nonMuTheta="grad")).
#
# Replaces the bounded bobyqa regression (gVaeThetaObjR) with the exact FOCEi
# outer gradient: ONE complete augmented sensitivity solve per M-step, the same
# machinery foceiControl(fast=TRUE) drives.  bobyqa needs a full N-subject inner
# likelihood sweep per function evaluation; this needs one solve.
#
# The two differ in target as well as cost: bobyqa minimizes the JOINT likelihood
# at frozen encoder etas, whose optimum is displaced from the marginal one, while
# this differentiates the marginal (Laplace) objective.
#
# The augmented solve calls rxSolveFree() and replaces rxode2's global solve, so
# the C++ caller MUST restoreFitSolve_() before touching the inner problem again.

.vaeGradEnv <- new.env(parent = emptyenv())

#' Stash the per-fit context the M-step gradient needs.
#'
#' Called once from `.vaeTrain` before the C++ loop starts; `.vaeGradEval` then
#' takes only the values that move between M-steps.
#' @param ui rxode2 ui (post pre-processing hooks)
#' @param data estimation data (`dataSav`)
#' @param regNames names of the thetas the M-step regresses, in `regIdx` order
#' @return invisible NULL
#' @noRd
.vaeGradInit <- function(ui, data, regNames) {
  ## .vaeInnerSetup replaced the ui's control with the DERIVED focei control, so
  ## .analyticGradCaller (which rxUiGet.foceiOuter consults) would resolve to NA.
  ## Re-mark it so the augmented model builds for this caller.
  .ctl <- tryCatch(ui$control, error = function(e) NULL)
  if (!is.null(.ctl)) {
    .ctl$nonMuTheta <- "grad"
    assign("control", .ctl, envir = ui)
  }
  .vaeGradEnv$ui <- ui
  .vaeGradEnv$data <- data
  .vaeGradEnv$ids <- unique(data$ID)
  .vaeGradEnv$regNames <- regNames
  .vaeGradEnv$am <- NULL        # augmented model, built lazily on the first M-step
  .vaeGradEnv$failed <- FALSE   # sticky: a declined build never re-attempts
  invisible(NULL)
}

#' One outer-gradient evaluation for the VAE M-step.
#'
#' @param thVals full natural-scale theta vector, ntheta order
#' @param ebes current encoder etas, `N x neta` (already centered on the baseline)
#' @param omegaDiag current M-step omega diagonal
#' @return numeric gradient over `regNames` (same order), or `NULL` to make the
#'   caller fall back to the bobyqa regression for this M-step
#' @noRd
.vaeGradEval <- function(thVals, ebes, omegaDiag) {
  if (isTRUE(.vaeGradEnv$failed)) return(NULL)
  .ui <- .vaeGradEnv$ui
  .reg <- .vaeGradEnv$regNames
  tryCatch({
    .Om <- diag(as.numeric(omegaDiag), nrow = length(omegaDiag))
    .st <- .foceiAnalyticGradSetup(.ui, thVals, .Om, caller = "vae")
    if (is.null(.st)) { .vaeGradEnv$failed <- TRUE; return(NULL) }
    if (ncol(ebes) != .st$neta) { .vaeGradEnv$failed <- TRUE; return(NULL) }
    ## The augmented model depends only on the model + direction set, never on
    ## theta/eta/omega, so it is built once and reused for every M-step (the
    ## symbolic .rxSens pass dominates otherwise).
    if (is.null(.vaeGradEnv$am)) {
      .am <- tryCatch(.ui$foceiOuter, error = function(e) NULL)
      if (is.null(.am) || !inherits(.am$augMod, "rxode2")) {
        .vaeGradEnv$failed <- TRUE
        return(NULL)
      }
      .vaeGradEnv$am <- .am
    }
    .th <- setNames(as.numeric(thVals), paste0("THETA_", seq_along(thVals), "_"))
    .r <- .foceiAnalyticGradCore(.ui, .th, ebes, .vaeGradEnv$ids, .vaeGradEnv$data,
                                 .Om, .st$ef, .st$dir, .st$dOiEst, .st$tr28,
                                 .st$omNames, .foceiAnalyticSolveTol(.ui),
                                 interaction = .st$interaction,
                                 foceType = .st$foceType, am = .vaeGradEnv$am,
                                 nAGQ = .st$nAGQ)
    if (is.null(.r) || is.null(.r$g)) return(NULL)
    .g <- .r$g[.reg]
    ## a regressed theta the gradient does not carry (not in thStruct) means the
    ## direction set and the M-step disagree -- decline rather than step on NA
    if (anyNA(.g) || !all(is.finite(.g))) return(NULL)
    as.numeric(.g)
  }, error = function(e) NULL)
}
