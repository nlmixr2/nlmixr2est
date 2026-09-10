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
# Solve plumbing: the augmented model is solved IN THE SHARED FOCEi pool by
# vaeOuterSolve_ (function-pointer swap + ind->neqOverride), which frees nothing
# and leaves rxode2's global solve in place -- so there is NO rxSolveFree() and
# the C++ caller must NOT restoreFitSolve_(); the inner problem is still live for
# the next vaeInnerLikCore.  (An earlier revision did call rxSolveFree() and
# required a restore; that is no longer how this works.)

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
## Clear EVERY per-fit field.  .vaeGradEnv has session lifetime, so anything left
## here (notably `data` and `ids`) is retained until the next grad fit -- the
## dataset can be large.  .vaeGradInit does overwrite all of these, so this is
## memory hygiene rather than stale-state correctness, but a completed fit should
## not hold its data hostage for the rest of the session.
.vaeGradReset <- function() {
  .vaeGradEnv$outerCols <- NULL
  .vaeGradEnv$am <- NULL
  .vaeGradEnv$ui <- NULL
  .vaeGradEnv$data <- NULL
  .vaeGradEnv$ids <- NULL
  .vaeGradEnv$regNames <- NULL
  .vaeGradEnv$cores <- NULL
  .vaeGradEnv$failed <- NULL
  ## The C++ pooled setup is per-fit (it holds this model's lhs column maps); clearing
  ## the flag makes the next fit install its own rather than inherit this one's shape.
  .vaeGradEnv$pooledOk <- NULL
  .vaeGradEnv$order <- NULL
  .vaeGradEnv$dv <- NULL
  invisible(NULL)
}

#' Sensitivity ORDER the M-step's objective actually needs
#'
#' `mStepObjective="outer"` is the full FOCEi outer objective -- frozen-eta joint
#' PLUS the Laplace determinant, `0.5*log|Omega^-1|` and the DV-transform
#' Jacobian.  Differentiating the determinant needs `d(eta*)/d(p)` and the second
#' derivative of the prediction, hence `order = 2`.
#'
#' `mStepObjective="elbo"` is the plain variational bound: the frozen-eta joint
#' likelihood and nothing else.  With the etas held fixed there is no `eta*` to
#' differentiate and no determinant, so FIRST-order sensitivities are sufficient
#' and the whole second-order expansion is waste.  Measured on a two-latent
#' declared gamma model over the same 8 directions:
#'
#'   order = 2   146 ODE states, 36 second-order pairs, 76s to build
#'   order = 1    18 ODE states,  0 second-order pairs,  5s to build
#'
#' and one `order = 2` gradient evaluation costs ~285s, which is why a declared
#' `nonMuTheta="grad"` fit could not finish.  The full gradient is only built when
#' the full objective is what the M-step was asked to optimize.
#' @noRd
.vaeGradAugOrder <- function(ui) {
  if (identical(tryCatch(rxode2::rxGetControl(ui, "mStepObjective", "outer"),
                         error = function(e) "outer"), "elbo")) 1L else 2L
}

#' UNVALIDATED -- DO NOT ENABLE.  Gradient of the ELBO M-step objective.
#'
#' MEASURED WRONG, and the reason is recorded here so the next attempt does not
#' repeat it.  Rebuilding `0.5*sum[log(2pi) + log R + res^2/R]` in R from
#' `vaeOuterSolve_`'s `f`/`R` gives **2.9e18** where `sum(vaeInnerLik$obj)` --
#' the objective `mStepObjective="elbo"` actually optimizes -- is **1.02e8**, a
#' factor of 2.8e10 out, and the gradient is correspondingly wrong (analytic
#' 2.1e20 vs central difference 1.5e9, and `rxCor` even sign-flipped).
#'
#' It is NOT a convention error: the augmented solve's `f` is a genuine decaying
#' profile (4.31, 0.83, 0.040, 0.0016) against DV (10.9, 3.91, 0.32, 0.020), and
#' the solve's first-order sensitivities match finite differences of its own `f`
#' to 4-5 digits.  The failure is the PROPORTIONAL VARIANCE in the tail: with
#' `R = (prop.sd*f)^2` the last observation has `R ~ 2.6e-8` and `res^2/R ~ 1.4e4`,
#' so a naive `log R + res^2/R` is dominated by points `likInner0` evidently
#' safeguards (a variance floor, and/or a different variance construction).
#'
#' CONCLUSION for the next attempt: do NOT re-derive the residual/variance terms
#' in R from solve columns.  Assemble the ELBO gradient in C++ beside `likInner0`,
#' where `R` and the residual are already formed with the same safeguards the
#' objective uses, and reuse them.  Any R-side reimplementation has to reproduce
#' every guard exactly, and silently disagrees when it does not.
#'
#' Currently unreachable from a fit: `nlmixr2Est.vae` still downgrades
#' `nonMuTheta="grad"` to `"regress"` under `mStepObjective="elbo"`, so nothing
#' calls this. Leave that downgrade in place until this is fixed and validated.
#'
#' What IS correct and kept: the ORDER selection. The ELBO objective needs only
#' first-order sensitivities -- 18 ODE states against 146, and 59.5s against 285s
#' per evaluation on a declared model.
#'
#' Gradient of the ELBO M-step objective: the frozen-eta joint -2LL.
#'
#'   -2LL_joint = sum_obs [ log R + (y - f)^2 / R ]
#'   d/dtheta   = sum_obs [ (1/R - (y-f)^2/R^2) * dR/dtheta - 2(y-f)/R * df/dtheta ]
#'
#' The `eta' Omega^-1 eta` prior carries no theta dependence at frozen etas, and
#' the ELBO objective excludes the Laplace determinant and the DV-transform
#' Jacobian by definition, so neither appears.  Everything needed is FIRST order.
#'
#' Works in PARAMETER space, not direction space, and takes the map explicitly.
#' `dir$dirP = c(dirTh, dirSg)` is the codebase's own "every non-Omega param -> a
#' direction" map; a residual sigma enters the variance only (`df/dsigma == 0` by
#' construction, which is why sigma directions carry no prediction chain), so its
#' variance derivative comes from `RsigDir` rather than from `aR`.  Getting this
#' mapping wrong is silent -- the gradient stays finite and points somewhere else.
#' @param E per-subject list from `vaeOuterSolve_`
#' @param dv per-subject observed values, in solve order
#' @param dirTh direction index of each structural theta (1-based)
#' @param nsg number of residual-sigma parameters, appended after the thetas
#' @return numeric gradient over `c(thStruct, sgName)`, or NULL
#' @noRd
.vaeGradElboAssemble <- function(E, dv, dirTh, nsg) {
  .nth <- length(dirTh)
  .np <- .nth + nsg
  if (.np == 0L) return(NULL)
  .g <- numeric(.np)
  .ok <- attr(E, "ok")
  .used <- 0L
  for (.i in seq_along(E)) {
    .Ei <- E[[.i]]
    if (is.null(.Ei) || (!is.null(.ok) && !isTRUE(.ok[.i] == 1L))) next
    .f <- as.numeric(.Ei$f)
    .y <- if (.i <= length(dv)) as.numeric(dv[[.i]]) else NULL
    if (is.null(.y) || length(.y) != length(.f)) return(NULL)
    .r <- if (is.null(.Ei$R)) rep(1.0, length(.f)) else as.numeric(.Ei$R)
    if (any(!is.finite(.r)) || any(.r <= 0)) return(NULL)
    .res <- .y - .f
    .w1 <- 1 / .r - (.res * .res) / (.r * .r)      # multiplies dR/dtheta
    .w2 <- -2 * .res / .r                          # multiplies df/dtheta
    ## `Rsig` is dR/d(sigma_k), nobs x nsig -- that is the first-order variance
    ## derivative a residual sigma needs.  NOT `RsigDir`, which is nobs x ndir x
    ## nsig, the MIXED second derivative d2R/(d dir)(d sigma) and a 3-D array; it
    ## belongs to the determinant block of the full outer objective, which the
    ## ELBO does not have.  Indexing it as a matrix simply throws.
    .a <- .Ei$a; .aR <- .Ei$aR; .rs <- .Ei$Rsig
    for (.q in seq_len(.nth)) {
      .d <- dirTh[.q]
      if (is.na(.d) || .d < 1L) return(NULL)
      .dfk <- if (!is.null(.a) && .d <= ncol(.a)) .a[, .d] else rep(0, length(.f))
      .dRk <- if (!is.null(.aR) && .d <= ncol(.aR)) .aR[, .d] else rep(0, length(.f))
      .g[.q] <- .g[.q] + sum(.w1 * .dRk + .w2 * .dfk)
    }
    for (.j in seq_len(nsg)) {
      .dRj <- if (!is.null(.rs) && length(dim(.rs)) == 2L && .j <= ncol(.rs)) {
        as.numeric(.rs[, .j])
      } else rep(0, length(.f))
      .g[.nth + .j] <- .g[.nth + .j] + sum(.w1 * .dRj)
    }
    .used <- .used + 1L
  }
  if (.used == 0L || !all(is.finite(.g))) return(NULL)
  .g
}

.vaeGradInit <- function(ui, data, regNames, order = NULL) {
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
  ## Passed in, NOT read off the ui.  `.vaeInnerSetup` replaced `ui$control` with
  ## the derived FOCEi control, which does not carry `mStepObjective`, so reading
  ## it here silently resolves to the "outer" default -- the order then disagrees
  ## with the model `.vaeInnerSetup` actually registered and sized the pool for,
  ## and every gradient call declines.  (Same trap as `.analyticGradCaller`.)
  .vaeGradEnv$order <- if (is.null(order)) .vaeGradAugOrder(ui) else as.integer(order)
  ## Observed values per subject, in solve order -- the ELBO assembly differences
  ## (y - f) per observation and there is no other place to get `y`.
  .vaeGradEnv$dv <- tryCatch({
    .d0 <- data
    if (!is.null(.d0$EVID)) .d0 <- .d0[.d0$EVID == 0, , drop = FALSE]
    split(as.numeric(.d0$DV), factor(.d0$ID, levels = .vaeGradEnv$ids))
  }, error = function(e) NULL)
  .am <- if (identical(.vaeGradEnv$order, 1L)) {
    .d1 <- tryCatch(.foceiOuterDirs(ui, "vae"), error = function(e) NULL)
    if (is.null(.d1)) NULL else {
      tryCatch(.foceiAnalyticAugModelDirs(ui, .d1$dirs, order = 1L),
               error = function(e) NULL)
    }
  } else {
    tryCatch(ui$foceiOuter, error = function(e) NULL)
  }
  if (!is.null(.am) && inherits(.am$augMod, "rxode2")) {
    .vaeGradEnv$am <- .am
    ## Enables the pooled vaeOuterSolve_ path.  Valid ONLY because the augmented
    ## model also SIZED the pool (.vaeInnerSetup sets poolModel) -- against an
    ## inner-sized pool this writes 26 states / 29 lhs into 6 / 6 buffers and dies
    ## with "double free or corruption".  The two switches move together.
    .vaeGradEnv$outerCols <- tryCatch(.vaeOuterCols(.am), error = function(e) NULL)
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
#' @param omega current M-step omega: full matrix, or a vector taken as diagonal
#' @return numeric gradient over `regNames` (same order), or `NULL` to make the
#'   caller fall back to the bobyqa regression for this M-step
#' @noRd
.vaeGradEval <- function(thVals, ebes, omega) {
  if (isTRUE(.vaeGradEnv$failed)) return(NULL)
  .ui <- .vaeGradEnv$ui
  .reg <- .vaeGradEnv$regNames
  tryCatch({
    .Om <- if (is.matrix(omega)) omega else diag(as.numeric(omega), nrow = length(omega))
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
    ## The pooled setup describes the SHAPE (lhs column maps, direction indices, which
    ## kernel) and depends only on the model, so install it once and reuse it for every
    ## M-step; the point itself -- theta, the encoder etas, omega -- is passed per call.
    ## This is the same C++ core a focei fit's own gradient runs, which is the point:
    ## the R implementation this replaced was a second, drifting copy of it.
    if (!isTRUE(.vaeGradEnv$pooledOk)) {
      .ps <- .foceiGradPooledSetup(.ui)
      if (is.null(.ps) || !isTRUE(foceiGradPooledSetupLoad_(.ps))) {
        .vaeGradEnv$failed <- TRUE
        return(NULL)
      }
      .vaeGradEnv$pooledOk <- TRUE
    }
    ## ELBO objective: the frozen-eta joint, assembled from the FIRST-order solve.
    ## foceiGradPooledDirect_ differentiates the full outer objective (Laplace
    ## determinant included), which is a DIFFERENT functional -- stepping with it
    ## while scoring the ELBO would optimize one thing and report another, and it
    ## costs ~285s a call on a declared model against ~18 states here.
    if (identical(.vaeGradEnv$order, 1L)) {
      if (is.null(.vaeGradEnv$outerCols) || is.null(.vaeGradEnv$dv)) return(NULL)
      .E <- vaeOuterSolve_(as.numeric(thVals), as.matrix(ebes),
                           .vaeGradEnv$outerCols, .vaeGradEnv$cores)
      if (is.null(.E)) return(NULL)
      .ge <- .vaeGradElboAssemble(.E, .vaeGradEnv$dv, .st$dir$dirTh,
                                  length(.st$ef$sgName))
      if (is.null(.ge)) return(NULL)
      .nm <- c(.st$dir$thStruct, .st$ef$sgName)
      if (length(.ge) != length(.nm)) return(NULL)
      names(.ge) <- .nm
      .gv <- .ge[.reg]
      if (anyNA(.gv) || !all(is.finite(.gv))) return(NULL)
      return(as.numeric(.gv))
    }
    .g <- foceiGradPooledDirect_(as.numeric(thVals), as.matrix(ebes),
                                 solve(.Om), .st$dOiEst, as.numeric(.st$tr28),
                                 .vaeGradEnv$cores)
    if (is.null(.g)) return(NULL)
    names(.g) <- c(.st$dir$thStruct, .st$ef$sgName, .st$omNames)
    .g <- .g[.reg]
    ## a regressed theta the gradient does not carry (not in thStruct) means the
    ## direction set and the M-step disagree -- decline rather than step on NA
    if (anyNA(.g) || !all(is.finite(.g))) return(NULL)
    as.numeric(.g)
  }, error = function(e) NULL)
}
