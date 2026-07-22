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

#' 0-based lhs offsets of every column `vaeOuterSolve_` reads, resolved from the
#' augmented model's OWN lhs names.
#'
#' Resolving here (rather than re-deriving the naming scheme in C++) means a
#' renamed generated column fails loudly in R instead of silently reading the
#' wrong offset out of the lhs buffer.
#' @param am augmented model (`ui$foceiOuter`)
#' @return list of index vectors, or `NULL` when any expected column is missing
#' @noRd
.vaeOuterCols <- function(am) {
  .lhs <- as.character(rxode2::rxModelVars(am$augMod)$lhs)
  .dirs <- am$dirs
  .fDirs <- if (is.null(am$fDirs)) .dirs else am$fDirs
  .cm <- if (is.null(am$cols)) {
    .foceiAnalyticCols(.dirs, .fDirs, am$P2, if (is.null(am$P2r)) am$P2 else am$P2r, am$sigTh)
  } else am$cols
  .ix <- function(nm) { .i <- match(nm, .lhs); if (anyNA(.i)) NULL else as.integer(.i - 1L) }
  .hasR <- isTRUE(am$hasRvar)
  .hasT <- isTRUE(am$hasTrans)
  .predf <- .ix("rx_predf_")
  .f1 <- .ix(.cm$f1); .f2 <- .ix(.cm$f2)
  if (is.null(.predf) || is.null(.f1) || is.null(.f2)) return(NULL)
  .o <- list(predf = .predf, f1 = .f1, f2 = .f2,
             iiF = as.integer(.cm$iiF - 1L), jjF = as.integer(.cm$jjF - 1L),
             fDirIdx = as.integer(.cm$fDirIdx - 1L),
             nd = length(.dirs), hasR = .hasR, hasT = .hasT)
  if (.hasR) {
    .rvarf <- .ix("rx_rvarf_"); .rvar1 <- .ix(.cm$rvar1); .rvar2 <- .ix(.cm$rvar2)
    if (is.null(.rvarf) || is.null(.rvar1) || is.null(.rvar2)) return(NULL)
    .o$rvarf <- .rvarf; .o$rvar1 <- .rvar1; .o$rvar2 <- .rvar2
    .o$ii <- as.integer(.cm$ii - 1L); .o$jj <- as.integer(.cm$jj - 1L)
    .o$rsig <- .ix(.cm$rsig); .o$rsig2 <- .ix(.cm$rsig2)
    if (is.null(.o$rsig) || is.null(.o$rsig2)) return(NULL)
    .o$rsig1 <- lapply(.cm$rsig1, .ix)
    if (any(vapply(.o$rsig1, is.null, logical(1)))) return(NULL)
    .o$sigA <- if (is.null(.cm$sigP2)) integer(0) else as.integer(.cm$sigP2$a - 1L)
    .o$sigB <- if (is.null(.cm$sigP2)) integer(0) else as.integer(.cm$sigP2$b - 1L)
  }
  if (.hasT) {
    .t <- .ix(c("rx_tyj_", "rx_tlambda_", "rx_tlow_", "rx_thi_"))
    if (is.null(.t)) return(NULL)
    .o$trans <- .t
  }
  .o
}

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
  ## Pooled-solve wiring: resolve the lhs offsets once.  NULL leaves
  ## .foceiAnalyticSolveAll on the rxSolve path (correct, just slower).
  .vaeGradEnv$outerCols <- NULL
  .vaeGradEnv$cores <- 1L
  .am <- tryCatch(ui$foceiOuter, error = function(e) NULL)
  if (!is.null(.am) && inherits(.am$augMod, "rxode2")) {
    .vaeGradEnv$am <- .am
    ## outerCols stays NULL until the augmented model actually SIZES the shared
    ## pool (blocked on the parameter-name mismatch -- see .vaeInnerSetup).
    ## Enabling vaeOuterSolve_ against an inner-sized pool writes 26 states / 29
    ## lhs into buffers sized for 6 / 6: "double free or corruption", not an error.
    .vaeGradEnv$cores <- tryCatch({
      .c <- .am$cores
      if (is.null(.c) || is.na(.c) || .c < 1L) 1L else as.integer(.c)
    }, error = function(e) 1L)
  }
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
