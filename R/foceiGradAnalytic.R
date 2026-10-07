# Analytic FOCE/FOCEI outer (population) gradient from Almquist (2015) sensitivity
# equations -- the first-derivative precursor of the analytic observed-information
# R-matrix in foceiCovAnalytic.R.  Enabled by foceiControl(fast=TRUE).  The gradient
# is computed in C++ (analyticOuterGrad -> gradPooledCore); this file holds the R
# side: scope checks, the direction set, the augmented sensitivity model and the
# per-fit setup the C++ loads once.
#
# The outer gradient needs at most SECOND-order state sensitivities (Almquist Eqs
# 38-40) -- one order less than the covariance R-matrix -- so it reuses the
# direction set and error machinery from foceiCovAnalytic.R without the 3rd-order
# (Ath) tier.

#' Per-FIT constants for the all-C++ analytic outer gradient.
#'
#' Everything the gradient needs that does NOT change between outer iterations: the
#' augmented model's lhs column map, the direction bookkeeping and the problem
#' dimensions.  Computed ONCE per fit and cached in C++, so `analyticOuterGrad` can run
#' with no R interaction at all -- the omega derivatives it also needs come from the
#' fit's own `rxInv` handle, which C++ already holds for the inner problem.
#'
#' Deliberately returns plain atomic vectors: C++ copies them into a POD struct and keeps
#' no R objects alive across the fit.
#' @param ui model UI
#' @param e fit environment (for the live omega/rxInv reuse)
#' @return a plain list, or NULL when the pooled gradient is out of scope
#' @noRd
.foceiGradPooledSetup <- function(ui, e = NULL) {
  tryCatch(
    {
      if (.foceiAnalyticIsMixture(ui)) {
        return(NULL)
      }
      .thv <- tryCatch(
        {
          .ini <- ui$iniDf
          .r <- .ini[!is.na(.ini$ntheta), , drop = FALSE]
          .r <- .r[order(.r$ntheta), , drop = FALSE]
          setNames(as.numeric(.r$est), .r$name)
        },
        error = function(.) NULL
      )
      if (is.null(.thv)) {
        return(NULL)
      }
      ## Omega is only used here for the SHAPE of the estimation-scale derivative block
      ## (how many free omega parameters there are); the values themselves are recomputed
      ## in C++ from the fit's rxInv on every call.  So the initial Omega is fine, and it
      ## is what is available before the fit starts.
      .Om <- tryCatch(get("omega", e), error = function(.) NULL)
      if (is.null(.Om)) {
        .Om <- tryCatch(ui$omega, error = function(.) NULL)
      }
      if (is.null(.Om)) {
        return(NULL)
      }
      st <- .foceiAnalyticGradSetup(ui, .thv, .Om, e)
      if (is.null(st)) {
        return(NULL)
      }
      ## Shape.  `.foceiAnalyticGradSetup` reports interaction = 0 for an ll() model as well
      ## as for FOCE, so isLL has to be carried separately and tested first.
      .isLL <- isTRUE(st$ef$isLL)
      .interaction <- as.integer(st$interaction)
      .nAGQ <- as.integer(st$nAGQ)
      am <- .foceiAnalyticAugModelDirs(ui, st$dir$dirs)
      if (is.null(am)) {
        return(NULL)
      }
      .cols <- tryCatch(.vaeOuterCols(am), error = function(.) NULL)
      if (is.null(.cols)) {
        return(NULL)
      }
      ## The (f,R) kernels contract a variance model.  An ll() endpoint has none -- its
      ## rx_pred_ IS the per-observation log density -- so hasR is required for every shape
      ## except that one.
      if (!.isLL && !isTRUE(.cols$hasR)) {
        return(NULL)
      }
      ## AGQ solves its quadrature nodes through a 1st-order sibling model (26 ODE states
      ## down to 8 on a one-compartment model).  Prefer the one built and disk-cached at
      ## model setup; fall back to building it, and then to the order-2 model, which is
      ## what the nodes used before that optimization existed.
      .colsNode <- NULL
      if (.nAGQ > 1L) {
        .amN <- tryCatch(
          {
            .fmN <- ui$foceiModel
            if (inherits(.fmN$outerNode, "rxode2") && !is.null(.fmN$outerNodeMeta)) {
              c(list(augMod = .fmN$outerNode), .fmN$outerNodeMeta)
            } else {
              .foceiAnalyticAugModelDirs(ui, st$dir$dirs, order = 1L)
            }
          },
          error = function(.) NULL
        )
        if (!is.null(.amN) && identical(.amN$ndir, am$ndir)) {
          .colsNode <- tryCatch(.vaeOuterCols(.amN), error = function(.) NULL)
        }
        if (is.null(.colsNode)) .colsNode <- .cols
      }
      ## Map each OUTER-OPTIMIZER parameter to its slot in the kernel's output vector.
      ##
      ## The kernel emits nth theta directions, then nsg sigma, then nom omega.  That is NOT
      ## the optimizer's parameter vector, for two independent reasons:
      ##   * an estimated boxCox/yeoJohnson lambda appears in BOTH dir$thStruct (it is a
      ##     direction) and ef$sgName (it is a sigma theta of the augmented model), so the
      ##     kernel emits it twice;
      ##   * the mu-referenced (lin/irls) families profile some structural thetas out of the
      ##     outer problem entirely, so the optimizer has FEWER parameters than the kernel.
      ## The deleted R route hid both: it named the vector c(thStruct, sgName, omNames) and
      ## the caller subset it by name, which takes the first match and drops the rest.  C++
      ## is positional and cannot, so carry the map explicitly.  Getting this wrong is
      ## silent -- the arity guard in analyticOuterGradDirect() just declines to finite
      ## differences -- which is how both cases went unnoticed.
      .kernelNames <- c(st$dir$thStruct, st$ef$sgName, st$omNames)
      .thAll <- ui$iniDf$name[!is.na(ui$iniDf$ntheta)]
      .thRows <- ui$iniDf[!is.na(ui$iniDf$ntheta), , drop = FALSE]
      .thRows <- .thRows[order(.thRows$ntheta), , drop = FALSE]
      .thEst <- .thRows$name[!.thRows$fix]
      .parNames <- c(setdiff(.thEst, .foceiMuSkipThetaNames(ui, .thAll)), st$omNames)
      .gMap <- match(.parNames, .kernelNames)
      if (anyNA(.gMap)) {
        return(NULL)
      }
      .gMap <- as.integer(.gMap - 1L) # 0-based for C++
      ## ntheta position of each structural theta -- the ll() perturbation of a non-mu
      ## theta moves th[thPos[p]], which is not the direction index.
      .thPos <- tryCatch(as.integer(ui$iniDf$ntheta[match(st$dir$thStruct, ui$iniDf$name)]), error = function(.) {
        integer(0)
      })
      list(
        cols = .cols,
        colsNode = .colsNode,
        neta = as.integer(st$neta),
        nth = as.integer(st$dir$nth),
        nsg = as.integer(length(st$ef$sgName)),
        nom = as.integer(length(st$dOiEst)),
        dirTh = as.integer(st$dir$dirTh),
        sigCol = seq_along(st$ef$sgName),
        lamDir = as.integer(st$dir$lamDir),
        nLam = as.integer(length(st$dir$lamNames)),
        censOpt = as.integer(rxode2::rxGetControl(ui, "censOption", 0L)),
        isLL = .isLL,
        interaction = .interaction,
        foceType = as.integer(st$foceType),
        nAGQ = .nAGQ,
        ebeTol = {
          ## FOCE frozen-R0 Newton score tolerance -- a fixed 1e-9, NOT derived from
          ## sigdig.  foceiControl(foceEbeTol=) overrides it.
          .t <- suppressWarnings(as.numeric(rxode2::rxGetControl(ui, "foceEbeTol", NA_real_)))
          if (!is.finite(.t) || .t <= 0) {
            .t <- 1e-9
          }
          .t
        },
        ebeSkipTol = 1e-3, ## looser first-iteration "already stationary?" test
        dependsF0 = isTRUE(st$ef$dependsF0),
        canVanish = isTRUE(st$ef$canVanish),
        thPos = .thPos,
        gMap = .gMap
      )
    },
    error = function(e) NULL
  )
}

#' Subjects whose augmented solve failed, for the Phase 8D2 per-individual FD.
#'
#' Recorded here rather than acted on in the solve loop: that loop runs inside
#' OdeSwapEsBatch(odeSlotOuter), and the finite difference needs the INNER problem's
#' event sensitivities, which can only be installed at a batch boundary.
#'
#' `n` counts pooled solves abandoned because of a flag, and only ever grows.  A test
#' that wants to prove the POOLED result was used, rather than merely attempted, has to
#' check this too: `pooledSolveN` counts the attempt and rises either way.
#' @noRd
.foceiOuterFlagged <- new.env(parent = emptyenv())
.foceiOuterFlagged$ids <- integer(0)
.foceiOuterFlagged$n <- 0L

#' Threads for the pooled augmented solve, from the fit's `rxControl(cores=)`.
#'
#' `0` (rxControl's default) and `NA` mean "use rxode2's thread setting", the same
#' reading `rxSolve` gives them; anything >= 1 is taken literally.  Resolving 0 to a
#' literal 1 is what left the pooled route serial.  C++ still caps the result with
#' min2(cores, getOpCores(op)).
#' @noRd
.foceiPoolCores <- function(cores) {
  tryCatch(
    {
      .c <- suppressWarnings(as.integer(cores)[1L]) # a NULL/character cores -> NA
      if (is.na(.c)) {
        return(as.integer(rxode2::getRxThreads()))
      }
      if (.c < 1L) as.integer(rxode2::getRxThreads()) else .c
    },
    error = function(e) 1L
  )
}

.foceiAnalyticSolveAll <- function(am, thv, ebes, ids, data, obsTimes, tol) {
  ## Solve the augmented model IN THE SHARED FOCEi pool (which it sized) and take
  ## the per-subject E structures straight from C++, instead of routing through
  ## rxode2::rxSolve, which frees and rebuilds the global solve on every call.
  ##
  ## `tol` is the tolerance to solve at and applies to BOTH routes; pass NA to solve at
  ## the fit's instead.  A covariance caller wants its own (covSolveTol, else tightened
  ## from sigdig) because it differences these solves; a gradient caller wants the fit's,
  ## so that it differentiates the objective being minimized.  It has NO DEFAULT on
  ## purpose -- the two answers are different numbers and there is no value that is right
  ## for both, so the choice is the caller's to state.
  ##
  ## No session flag guards this any more.  vaeOuterSolve_ refuses unless the
  ## augmented model is registered AND the pool is at least its size
  ## (odeSwapCanPool -> odeDenyPoolNotSized otherwise), which is the structural
  ## form of what `.vaeGradEnv$active` was patching: a focei fast fit after a vae
  ## grad fit, running against a pool sized for its own inner model.  A NULL
  ## falls through to the rxSolve path below.
  ##
  ## DDE is IN scope for the pooled route.  The old exclusion assumed a delay model
  ## pins method="dop853"/dense per solve, which a shared pool cannot do -- but focei
  ## already forces that configuration at the FIT level (R/focei.R, the hasDelay block:
  ## method 0, stiff2 13, dense TRUE), so a DDE fit's pool is built that way to begin
  ## with and there is nothing to change per solve.
  ## (This used to be gated by .odeSwapNoPool, a verification-only opt-out that let a
  ## test evaluate the same fit through rxSolve instead of the pool.  Its only setter was
  ## the R gradient route, which is gone, so the gate could never fire.)
  .cols <- tryCatch(.vaeOuterCols(am), error = function(e) NULL)
  if (!is.null(.cols)) {
    .nc <- .foceiPoolCores(am$cores) # 0 means rxode2's threads, not one
    ## The pooled solve takes one tolerance for atol and rtol both, while the rxSolve
    ## fallback below reads a 2-vector as (atol, rtol).  Every caller passes a scalar;
    ## take the tighter of a pair rather than half the request.
    .tolP <- suppressWarnings(min(as.numeric(tol)))
    .Ec <- tryCatch(
      vaeOuterSolve_(
        as.numeric(thv),
        as.matrix(ebes),
        .cols,
        .nc,
        if (length(.tolP) != 1L || !is.finite(.tolP)) NA_real_ else .tolP
      ),
      error = function(e) NULL
    )
    ## vaeOuterSolve_ flags failed subjects per individual (attr "ok") rather than
    ## discarding the whole population.  Nothing here consumes the flags yet, so a
    ## flagged subject falls THROUGH to the rxSolve route below -- all or nothing.
    ##
    ## This branch used to zero-fill a flagged subject's E and return it, on the
    ## grounds that its column is replaced wholesale by the per-individual finite
    ## difference.  That was the R gradient route, which is gone; every
    ## caller now reads the E structures as they stand, so the zeros went straight
    ## into the covariance as a subject with no prediction and no sensitivity -- and
    ## a zero E is FINITE, so it did not even trip the callers' is.finite guards.
    ##
    ## The per-individual finite difference lives on the all-C++ route, which owns the
    ## per-subject gradient columns; this assembly has none to substitute into.
    if (!is.null(.Ec) && length(.Ec) > 0L) {
      .ok <- attr(.Ec, "ok")
      .foceiOuterFlagged$ids <- if (is.null(.ok)) integer(0) else which(.ok == 0L)
      if (length(.foceiOuterFlagged$ids) == 0L) {
        return(.Ec)
      }
      .foceiOuterFlagged$n <- .foceiOuterFlagged$n + 1L
    }
  }
  dirs <- am$dirs
  nd <- length(dirs)
  neta <- ncol(ebes)
  etav <- paste0("ETA_", seq_len(neta), "_")
  pars <- data.frame(ID = ids)
  for (k in seq_len(neta)) {
    pars[[etav[k]]] <- ebes[, k]
  }
  for (.nm in names(thv)) {
    pars[[.nm]] <- thv[[.nm]]
  }
  .ev <- .foceiAnalyticEvents(am, data) # reuse the pre-translated event table
  .nc <- if (is.null(am$cores)) 0L else am$cores # fit's rxControl thread count (parallel)
  # DDE: force pure dop853 (dense, no Jacobian) -- its 8th-order dense history reproduces the
  # delayed sensitivity solve exactly and needs no Jacobian, sidestepping the composite/ros4
  # on-the-fly Jacobian generation for this THETA/ETA-named augmented model.
  .ddeArgs <- if (isTRUE(rxode2::rxModelVars(am$augMod)$flags[["hasDelay"]] == 1L)) {
    list(method = "dop853", stiff2 = 0L, dense = TRUE)
  } else {
    list()
  }
  .sol <- tryCatch(
    withCallingHandlers(
      as.data.frame(do.call(
        rxode2::rxSolve,
        c(
          list(
            am$augMod,
            params = pars,
            events = .ev,
            cores = .nc,
            returnType = "data.frame",
            atol = tol[1],
            rtol = tol[length(tol)]
          ),
          .ddeArgs
        )
      )),
      warning = function(w) invokeRestart("muffleWarning")
    ),
    error = function(e) NULL
  )
  if (
    is.null(.sol) || !all(c("rx_predf_", paste0("rx_f1_", if (is.null(am$fDirs)) dirs else am$fDirs)) %in% names(.sol))
  ) {
    return(NULL)
  }
  # Extract every sensitivity column from the WHOLE solve as a matrix ONCE (the per-subject
  # data.frame [[ + paste0 dominated the gradient -- the ODE solve itself is ~4%); slice by row
  # per subject.  Column names + index maps are precomputed on `am$cols` at build time.
  .cm <- if (is.null(am$cols)) {
    .foceiAnalyticCols(
      dirs,
      if (is.null(am$fDirs)) dirs else am$fDirs,
      am$P2,
      if (is.null(am$P2r)) am$P2 else am$P2r,
      am$sigTh
    )
  } else {
    am$cols
  }
  np2 <- nrow(am$P2)
  np2r <- if (is.null(am$P2r)) np2 else nrow(am$P2r)
  .hasR <- isTRUE(am$hasRvar)
  .hasT <- isTRUE(am$hasTrans)
  # Every column subset below must be present, or `.sol[, cols]` throws "undefined columns
  # selected" (e.g. an rxode2 version/feature that did not emit an rvar/rsig/transform column).
  # Verify up front and cleanly return NULL so the caller falls back to finite differences.
  .need <- c("rx_predf_", "time", .cm$f1, .cm$f2)
  if (.hasR) {
    .need <- c(.need, "rx_rvarf_", .cm$rvar1, .cm$rvar2)
    if (length(am$sigTh) > 0L) .need <- c(.need, .cm$rsig, unlist(.cm$rsig1), .cm$rsig2)
  }
  if (.hasT) {
    .need <- c(.need, "rx_tyj_", "rx_tlambda_", "rx_tlow_", "rx_thi_")
  }
  if (!all(.need %in% names(.sol))) {
    return(NULL)
  }
  .M1 <- as.matrix(.sol[, .cm$f1, drop = FALSE])
  .M2 <- as.matrix(.sol[, .cm$f2, drop = FALSE])
  .fp <- .sol$rx_predf_
  .tm <- .sol$time
  if (.hasR) {
    .MR1 <- as.matrix(.sol[, .cm$rvar1, drop = FALSE])
    .MR2 <- as.matrix(.sol[, .cm$rvar2, drop = FALSE])
    .Rf <- .sol$rx_rvarf_
    .nsig <- length(am$sigTh)
    if (.nsig > 0L) {
      .MRs <- as.matrix(.sol[, .cm$rsig, drop = FALSE])
      .MRs1 <- lapply(.cm$rsig1, function(cc) as.matrix(.sol[, cc, drop = FALSE]))
      .MRs2 <- as.matrix(.sol[, .cm$rsig2, drop = FALSE])
    }
  }
  if (.hasT) {
    .tr <- list(yj = .sol$rx_tyj_, lambda = .sol$rx_tlambda_, low = .sol$rx_tlow_, hi = .sol$rx_thi_)
  }
  .idcol <- if ("id" %in% names(.sol)) .sol$id else .sol$ID
  .byIdSol <- split(seq_len(nrow(.sol)), as.character(.idcol))
  Es <- vector("list", length(ids))
  for (i in seq_along(ids)) {
    .ri <- .byIdSol[[as.character(ids[i])]]
    if (is.null(.ri)) {
      .ri <- .byIdSol[[as.character(i)]]
    }
    if (is.null(.ri)) {
      return(NULL)
    }
    .keep <- .ri[.tm[.ri] %in% obsTimes[[i]]]
    no <- length(.keep)
    if (no != length(obsTimes[[i]])) {
      return(NULL)
    }
    a <- matrix(0, no, nd)
    a[, .cm$fDirIdx] <- .M1[.keep, , drop = FALSE]
    A <- array(0, c(no, nd, nd))
    for (r in seq_len(np2)) {
      .v <- .M2[.keep, r]
      A[, .cm$iiF[r], .cm$jjF[r]] <- .v
      A[, .cm$jjF[r], .cm$iiF[r]] <- .v
    }
    .E <- list(f = .fp[.keep], a = a, A = A)
    if (.hasR) {
      aR <- .MR1[.keep, , drop = FALSE]
      AR <- array(0, c(no, nd, nd))
      for (r in seq_len(np2r)) {
        .v <- .MR2[.keep, r]
        AR[, .cm$ii[r], .cm$jj[r]] <- .v
        AR[, .cm$jj[r], .cm$ii[r]] <- .v
      }
      .E$R <- .Rf[.keep]
      .E$aR <- aR
      .E$AR <- AR
      if (.nsig > 0L) {
        .E$Rsig <- .MRs[.keep, , drop = FALSE]
        .E$RsigDir <- array(vapply(.MRs1, function(M) M[.keep, , drop = FALSE], matrix(0, no, nd)), c(no, nd, .nsig))
        .Rs2 <- array(0, c(no, .nsig, .nsig))
        for (r in seq_len(nrow(.cm$sigP2))) {
          .v <- .MRs2[.keep, r]
          .Rs2[, .cm$sigP2$a[r], .cm$sigP2$b[r]] <- .v
          .Rs2[, .cm$sigP2$b[r], .cm$sigP2$a[r]] <- .v
        }
        .E$Rsig2 <- .Rs2
      }
    }
    if (.hasT) {
      .E$trans <- lapply(.tr, `[`, .keep)
    } # both-sides transform: DV -> tbs(DV) scale
    Es[[i]] <- .foceiFloorRvar(.E)
  }
  Es
}

#' Floor a tiny residual variance the way the FOCEi objective does (#1132): a floored
#' (or zero -> 1) R is constant, so its sensitivities are zero too.
#' @noRd
.foceiFloorRvar <- function(E) {
  .fl <- sqrt(.Machine$double.eps)
  .w <- which(E$R < .fl)
  if (length(.w) == 0L) {
    return(E)
  }
  E$R[.w] <- ifelse(E$R[.w] <= 0, 1, .fl)
  E$aR[.w, ] <- 0
  E$AR[.w, , ] <- 0
  if (!is.null(E$Rsig)) {
    E$Rsig[.w, ] <- 0
    E$RsigDir[.w, , ] <- 0
    E$Rsig2[.w, , ] <- 0
  }
  E
}

.foceiAnalyticIsMixture <- function(ui) {
  isTRUE(tryCatch(length(ui$thetaMixIndex) > 0L, error = function(e) FALSE))
}

.foceiAnalyticGradSetup <- function(ui, thVals, Om, e = NULL, caller = .analyticGradCaller(ui)) {
  if (is.na(caller)) {
    return(NULL)
  }
  if (!.hasRxSens()) {
    return(NULL)
  }
  if (.foceiAnalyticIsMixture(ui)) {
    return(NULL)
  } # mixtures: weighted sum, no treatment yet
  if (.foceiUsesLinCmt(ui)) {
    return(NULL)
  } # linCmt(): no symbolic state sensitivities
  if (length(.foceiLaggedCalcVars(ui)) > 0L) {
    return(NULL)
  } # lag() of a calculated variable: no symbolic sensitivity through it
  if (!.analyticGradAllowsBoundedTr(ui, caller)) {
    return(NULL)
  }
  # tad/podo/tafd/tlast/tfirst/dosenum are functions of time and the dose record
  # only (no eta/theta dependence), so rxode2 treats them as zero-derivative
  # constants in the sensitivity expansion (.rxToSEDualVarFunction) -- they no
  # longer need to force the finite-difference fallback.
  if (isTRUE(as.logical(rxode2::rxGetControl(ui, "fo", FALSE)))) {
    return(NULL)
  }
  # ll()/generalized likelihood (needOptimHess -> interaction=0, EXACT inner
  # Hessian): rx_pred_ is the log-density, so skip the Gaussian ErrFull and set up
  # the direct-log-density core (gradPooledCoreLL in C++).
  if (.foceiLLGradInScope(ui, caller)) {
    .map <- .foceiEtaThetaMap(ui)
    neta <- length(.map$etaNames)
    .dir <- .foceiOuterDirsLL(ui)
    if (is.null(.dir)) {
      return(NULL)
    }
    .oe <- .foceiEstOmegaDeriv(ui, Om, e)
    if (is.null(.oe)) {
      return(NULL)
    }
    return(list(
      ef = list(isLL = TRUE),
      dir = .dir,
      dOiEst = .oe$dOi,
      tr28 = .oe$tr28,
      omNames = .oe$names,
      neta = neta,
      etaNames = .map$etaNames,
      interaction = 0L,
      foceType = 0L,
      nAGQ = 1L
    ))
  }
  interaction <- as.integer(rxode2::rxGetControl(ui, "interaction", 1L)) # 1 FOCEI / 0 FOCE
  foceType <- if (interaction == 0L) as.integer(rxode2::rxGetControl(ui, "foceType", 0L)) else 0L
  ## FOCE (interaction = 0) was declined here up front: its frozen-R0 EBE Newton could
  ## not reach the 1e-9 score target at the default solve, |S| flooring near 5e-3 at
  ## rtol = 1e-3 (nlmixr2/nlmixr2est#836).  That was measured BEFORE the shared ODE solve
  ## pool was fixed (#839), where a peer solve run under another slot's event-sensitivity
  ## shape corrupted the scratch the score is assembled from.  The gate is lifted so FOCE
  ## goes through gradPooledCore's isFoce/foceEbeNewton path like any other shape; a
  ## Newton that still cannot converge declines per fit at its own site rather than
  ## being refused for the whole method.
  nAGQ <- as.integer(rxode2::rxGetControl(ui, "nAGQ", 1L))
  # agqControl() forces interaction=TRUE, so only the FOCEI (f,R) kernel has a quadrature
  # form -- a FOCE-AGQ combination cannot arise.
  if (nAGQ > 1L && interaction != 1L) {
    return(NULL)
  }
  # the aqLow/aqHi clamp (inner.cpp) kinks the objective; both default to +/-Inf
  if (
    nAGQ > 1L &&
      (is.finite(as.numeric(rxode2::rxGetControl(ui, "agqLow", -Inf))) ||
        is.finite(as.numeric(rxode2::rxGetControl(ui, "agqHi", Inf))))
  ) {
    return(NULL)
  }
  # The grid is placed by Ht's Cholesky FACTOR (etaCur = etahat + chol(Ht)^-1 x), so we
  # must differentiate the exact factorization the objective used.  cholSEOpt forces the
  # generalized Cholesky, whose factor differs from chol() even for a PD matrix.  (The
  # runtime doChol flips are temporary and fire only when the plain Cholesky fails --
  # where chol(Ht) fails in R too, so those are already covered.)
  if (nAGQ > 1L && isTRUE(as.logical(rxode2::rxGetControl(ui, "cholSEOpt", FALSE)))) {
    return(NULL)
  }
  ef <- .foceiAnalyticErrFull(ui)
  if (is.null(ef)) {
    return(NULL)
  }
  ini <- ui$iniDf
  .map <- .foceiEtaThetaMap(ui)
  etaNames <- .map$etaNames
  neta <- length(etaNames)
  if (neta == 0L) {
    return(NULL)
  }
  thetaForEta <- .map$thetaForEta
  if (length(.uiIovEnv$iovVars) > 0L) {
    return(NULL)
  } # IOV out of Phase-1 scope
  .valc <- setNames(as.numeric(thVals[ef$sgName]), ef$sgVar)
  ef$ev <- local({
    v <- .valc
    function(expr, f, y, f0 = f) eval(expr, c(list(f = f, y = y, f0 = f0), as.list(v)))
  })
  if (any(.iniIsFixed(ini, thetaForEta))) {
    return(NULL)
  }
  keep <- !.iniIsFixed(ini, ef$sgName)
  ef$sgVar <- ef$sgVar[keep]
  ef$sgName <- ef$sgName[keep]
  .dir <- .foceiAnalyticDirections(ini, thetaForEta, ef$sgName, neta, sharedEta = unname(.foceiEtaOccurrence(ui) > 1L))
  if (is.null(.dir)) {
    return(NULL)
  }
  # multiple estimated lambdas (per-endpoint) need an endpoint->lambda DV mapping not yet
  # wired; keep those on FD.  A single estimated lambda is the ported case.
  if (length(.dir$lamNames) > 1L) {
    return(NULL)
  }
  .oe <- .foceiEstOmegaDeriv(ui, Om, e)
  if (is.null(.oe)) {
    return(NULL)
  }
  list(
    ef = ef,
    dir = .dir,
    dOiEst = .oe$dOi,
    tr28 = .oe$tr28,
    omNames = .oe$names,
    neta = neta,
    etaNames = etaNames,
    interaction = interaction,
    foceType = foceType,
    nAGQ = nAGQ
  )
}

#' Post-fit analytic natural-scale gradient, computed by the fit's OWN C++ path.
#'
#' Re-enters the estimation machinery at the fit's converged estimates and EBEs --
#' `est="none"` with zero inner/outer iterations and the fit's `etaMat`, the same
#' post-fit re-entry `setCov()` uses -- and reads back the gradient
#' `analyticOuterGrad()` stashed.  So this is the SHIPPING gradient, not a parallel
#' implementation of it: that distinction is the whole point, because a test that
#' validates a second implementation proves nothing about the code the fit runs.
#'
#' `NULL` when the fit is out of analytic scope -- the same signal
#' `foceiControl(fast=TRUE)` acts on when it falls back to finite differences -- and
#' also when the fit was not run with `fast=TRUE` at all.
#' @param fit nlmixr2 fit object
#' @return named natural-scale gradient (structural thetas, sigmas, om.chol), or `NULL`
#' @noRd
.foceiGradDirect <- function(fit) {
  tryCatch(
    {
      .env <- if (rxode2::rxIs(fit, "nlmixr2FitData")) fit$env else fit
      .est <- .env$est
      if (is.null(.est) || !nzchar(.est)) {
        return(NULL)
      }
      .control <- .env$foceiControl
      .control$maxInnerIterations <- 0L # evaluate at the fit's EBEs, do not re-optimize
      .control$maxOuterIterations <- 0L # no outer step: the gradient is at THIS theta
      .control$calcTables <- FALSE
      .control$covMethod <- 0L # no covariance step
      .control$skipCov <- fit$skipCov
      # `fast` is deliberately NOT forced on: this reports what the analytic gradient does
      # for THIS fit as configured, so a fast=FALSE fit correctly yields NULL rather than a
      # gradient it never used.
      # Re-run under the fit's OWN est, not est="none": the gradient SHAPE is the
      # estimation method (FOCE freezes the residual variance, AGQ adds quadrature), and
      # est="none" would silently evaluate every fit as plain FOCEI.
      .ui <- .uiPinTheta(fit$ui, fit)
      .etaMat <- .fitEtaMat(fit)
      if (!is.null(.etaMat)) {
        .control$etaMat <- .etaMat
      }
      # the nested re-fit resets mu-referencing global state (.muRefTrans$cur); restore it
      .savedMuRef <- .muRefTrans$cur
      on.exit(.muRefTrans$cur <- .savedMuRef, add = TRUE)
      .f2 <- suppressMessages(suppressWarnings(
        nlmixr2(.ui, data = getData(fit), est = .est, control = .control)
      ))
      .src <- tryCatch(.f2$env, error = function(.) NULL)
      if (is.null(.src) || !exists(".gradDirectFirst", .src, inherits = FALSE)) {
        return(NULL)
      }
      .g <- as.numeric(get(".gradDirectFirst", .src))
      # Name it the way the gradient assembly orders it: structural thetas, then sigmas,
      # then the estimation-scale omega (Cholesky) elements.
      .ini <- .ui$iniDf
      .thRows <- .ini[!is.na(.ini$ntheta), , drop = FALSE]
      .thRows <- .thRows[order(.thRows$ntheta), , drop = FALSE]
      .thv <- fit$theta[.thRows$name]
      if (anyNA(.thv)) {
        .thv <- .thRows$est
      }
      .st <- .foceiAnalyticGradSetup(.ui, stats::setNames(as.numeric(.thv), .thRows$name), fit$omega)
      # Everything from here on only NAMES `.g`, so a naming failure must not discard it.
      # Returning NULL here reported "no analytic gradient" for a gradient that had in fact
      # been computed -- exactly what a fix()ed theta did, since .foceiAnalyticGradSetup
      # declines for one.  That made the analytic path untestable on precisely the models
      # where full-theta and free-parameter indexing differ, which is the indexing the outer
      # FD fallback's step store gets wrong when it gets it wrong.  Hand back the unnamed
      # values; the same reasoning already applies to the gMap branch below.
      if (is.null(.st)) {
        return(.g)
      }
      # These names are in KERNEL space (nth + nsg + nom).  That is not the outer
      # optimizer's vector whenever a parameter occupies two kernel slots -- an estimated
      # boxCox/yeoJohnson lambda is both a theta direction and a sigma slot, so the kernel
      # names come out one longer than the gradient (9 vs 8 on a 1-cmt boxCox model) and
      # this used to bail to NULL, reporting "no analytic gradient" for a gradient that
      # had in fact been computed.  gMap is the same kernel -> outer gather the C++ uses
      # (analyticOuterGradDirect), so reuse it rather than re-deriving the correspondence.
      .nmKer <- c(.st$dir$thStruct, .st$ef$sgName, .st$omNames)
      .nm <- .nmKer
      if (length(.nmKer) != length(.g)) {
        .gp <- tryCatch(.foceiGradPooledSetup(.ui), error = function(e) NULL)
        .map <- if (is.null(.gp)) NULL else .gp$gMap
        if (is.null(.map) || length(.map) != length(.g) || any(.map < 0L) || any(.map >= length(.nmKer))) {
          return(.g)
        }
        .nm <- .nmKer[.map + 1L]
      }
      stats::setNames(.g, .nm)
    },
    error = function(e) NULL
  )
}

.analyticGradCaller <- function(ui) {
  if (isTRUE(as.logical(rxode2::rxGetControl(ui, "fast", FALSE)))) {
    return("focei")
  }
  if (identical(as.character(rxode2::rxGetControl(ui, "nonMuTheta", "")), "grad")) {
    return("vae")
  }
  NA_character_
}

#' Bounded-transform scope gate.
#'
#' `preProcessBoundedTransform` records the transforms on the ALREADY-REWRITTEN
#' ui, so by the time the gradient sees them the model is on the unconstrained
#' `rxBoundedTr.*` scale.  focei must still bail: it REPORTS a natural-scale
#' gradient to the outer optimizer, which would need a Jacobian correction that is
#' not applied.  The VAE consumes the gradient internally, on the same
#' unconstrained scale it takes its M-step on, so no correction arises.
#' @noRd
.analyticGradAllowsBoundedTr <- function(ui, caller) {
  if (identical(caller, "vae")) {
    return(TRUE)
  }
  is.null(ui$boundedTransforms) || length(ui$boundedTransforms) == 0L
}

#' Is the analytic outer gradient in scope for a VAE fit?
#'
#' Cheap direction-set probe -- no symengine/gcc pass -- covering every static
#' gate: `linCmt()`, `fo`, the distribution/error-model scope, IOV, and a model
#' with no eta.  A later build or solve failure still falls back at runtime.
#'
#' Two admissible shapes, the same pair `.foceiAnalyticGradSetup` dispatches on:
#' a conditionally Gaussian endpoint (the `(f,R)` direction set) or a single
#' non-Gaussian `ll()`/generalized endpoint (the direct-log-density set).  The
#' VAE consumes either through the SAME C++ gradient core, which
#' routes on `ef$isLL`.
#' @noRd
.vaeGradInScope <- function(ui) {
  if (!is.null(tryCatch(.foceiOuterDirs(ui, "vae"), error = function(e) NULL))) {
    return(TRUE)
  }
  isTRUE(.foceiLLGradInScope(ui, "vae")) &&
    !is.null(tryCatch(.foceiOuterDirsLL(ui), error = function(e) NULL))
}

#' Direction set for the augmented outer-gradient model, computed from the UI
#' alone (does not depend on theta/eta values): one direction per eta plus one per
#' non-mu-referenced structural theta.  `NULL` if out of analytic scope.
#' @noRd
.foceiOuterDirs <- function(ui, caller = .analyticGradCaller(ui)) {
  if (!.hasRxSens()) {
    return(NULL)
  }
  if (.foceiAnalyticIsMixture(ui)) {
    return(NULL)
  } # mixtures: weighted sum, no treatment yet
  if (.foceiUsesLinCmt(ui)) {
    return(NULL)
  } # linCmt(): no symbolic state sensitivities
  if (length(.foceiLaggedCalcVars(ui)) > 0L) {
    return(NULL)
  } # lag() of a calculated variable: no symbolic sensitivity through it
  if (!.analyticGradAllowsBoundedTr(ui, caller)) {
    return(NULL)
  }
  if (isTRUE(as.logical(rxode2::rxGetControl(ui, "fo", FALSE)))) {
    return(NULL)
  }
  ef <- .foceiAnalyticErrFull(ui)
  if (is.null(ef)) {
    return(NULL)
  }
  .map <- .foceiEtaThetaMap(ui)
  neta <- length(.map$etaNames)
  if (neta == 0L) {
    return(NULL)
  }
  if (length(.uiIovEnv$iovVars) > 0L) {
    return(NULL)
  }
  .foceiAnalyticDirections(
    ui$iniDf,
    .map$thetaForEta,
    ef$sgName,
    neta,
    sharedEta = unname(.foceiEtaOccurrence(ui) > 1L)
  )
}

#' Is a fit in scope for the ll()/generalized-likelihood analytic outer gradient?
#' The ll() objective uses the EXACT inner Hessian (needOptimHess), so `rx_pred_`
#' is the per-observation log-density and the gradient is assembled by
#' the log-density core (differentiating it directly) rather than
#' the Gaussian (f,R) path.  Scope: at least one non-Gaussian endpoint, no
#' linCmt/bounded transform/IOV/FO, at least one eta.  (Censoring and nAGQ are
#' handled by falling back to the finite-difference gradient.)
#'
#' `caller` only reaches the bounded-transform gate, which focei must fail and
#' the VAE need not -- see `.analyticGradAllowsBoundedTr`.  Defaulted, so the
#' focei callers keep their exact behavior.
#' @noRd
.foceiLLGradInScope <- function(ui, caller = .analyticGradCaller(ui)) {
  tryCatch(
    {
      if (!.hasRxSens()) {
        return(FALSE)
      }
      .pd <- ui$predDfFocei
      if (is.null(.pd) || nrow(.pd) < 1L) {
        return(FALSE)
      }
      ## Multiple endpoints ARE in scope.  They were gated off on the reading that the
      ## multi-endpoint gradient did not verify -- against central differences, add.pd was
      ## ~373x off and tka/tv/add.pk 4.2x/1.9x/2.5x off, while tcl (the only theta carrying
      ## an eta) was right.  The gradient was right and the OBJECTIVE it was differenced
      ## against was wrong: the endpoint's distribution was read one row early, so one
      ## observation per subject was scored as normal (nlmixr2/nlmixr2est#838, fixed in
      ## likInner0).  It looked like direction bookkeeping because the corrupted row is the
      ## subject's FIRST, which biases whichever endpoint that row belongs to.  With the
      ## objective fixed the analytic gradient matches central differences to 8e-3 relative
      ## on the 2-endpoint warfarin ll() model, the residual being the reference's own
      ## step noise.
      if (all(as.character(.pd$distribution) %in% c("norm", "dnorm"))) {
        return(FALSE)
      } # Gaussian -> (f,R) path
      # linCmt() anywhere (not only as the endpoint): no 2nd-order sensitivities, and
      # rxode2 >= 5.1.8 drops those terms silently instead of failing the build (#1103).
      if (.foceiUsesLinCmt(ui)) {
        return(FALSE)
      }
      if (!.analyticGradAllowsBoundedTr(ui, caller)) {
        return(FALSE)
      }
      if (isTRUE(as.logical(rxode2::rxGetControl(ui, "fo", FALSE)))) {
        return(FALSE)
      }
      if (as.integer(rxode2::rxGetControl(ui, "nAGQ", 1L)) > 1L) {
        return(FALSE)
      }
      if (length(.uiIovEnv$iovVars) > 0L) {
        return(FALSE)
      }
      length(.foceiEtaThetaMap(ui)$etaNames) > 0L
    },
    error = function(e) FALSE
  )
}

#' Direction set for the ll() analytic outer gradient: one direction per eta plus
#' one per non-mu-referenced structural theta.  For an ll() endpoint there is no
#' Gaussian residual-sigma set (add.sd et al. appear directly in the log-density),
#' so every non-mu structural theta gets its own direction (`sgName = character(0)`).
#' @noRd
.foceiOuterDirsLL <- function(ui) {
  .map <- .foceiEtaThetaMap(ui)
  neta <- length(.map$etaNames)
  if (neta == 0L) {
    return(NULL)
  }
  .d <- .foceiAnalyticDirections(
    ui$iniDf,
    .map$thetaForEta,
    character(0),
    neta,
    sharedEta = unname(.foceiEtaOccurrence(ui) > 1L)
  )
  if (is.null(.d) || is.null(.d$dirs)) {
    return(NULL)
  }
  .d
}

# Build the augmented outer-gradient sensitivity model (compiled model + `dirs` +
# `P2`) for a UI.  This is the persistent `..outer` sibling of the inner model:
# it depends only on the model + direction set (NOT theta/eta/omega), so it is
# built once during model setup (via `rxUiGet.foceiModel`/`foceModel`, which
# disk-cache the whole model list) and reused across every outer-gradient call.
# Callable independently as `ui$foceiOuter`.  `NULL` when out of analytic scope
# (the gradient then falls back to finite differences).
#' @export
rxUiGet.foceiOuter <- function(x, ...) {
  .ui <- x[[1]]
  .caller <- .analyticGradCaller(.ui)
  if (is.na(.caller)) {
    return(NULL)
  }
  interaction <- as.integer(rxode2::rxGetControl(.ui, "interaction", 1L))
  foceType <- if (interaction == 0L) as.integer(rxode2::rxGetControl(.ui, "foceType", 0L)) else 0L
  # nAGQ > 1 (adaptive Gaussian quadrature) uses the SAME augmented model at eta-hat: the
  # quadrature nodes are extra eta points on the same sensitivity solve, so the
  # direction set and the symbolic expansion are unchanged.  (The nodes themselves solve
  # a cheaper 1st-order model -- see rxUiGet.foceiOuterNode.)
  .dir <- .foceiOuterDirs(.ui, .caller)
  # ll()/generalized endpoint (needOptimHess, interaction=0): the Gaussian (f,R)
  # direction builder declines (ErrFull is norm-only), but rx_pred_ is the
  # log-density and the same augmented model supplies its 1st/2nd-order eta/theta
  # derivatives -- build over the ll() direction set instead.
  if (is.null(.dir) && .foceiLLGradInScope(.ui, .caller)) {
    .dir <- .foceiOuterDirsLL(.ui)
  }
  if (is.null(.dir)) {
    return(NULL)
  }
  .foceiAnalyticAugModelDirs(.ui, .dir$dirs)
}
attr(rxUiGet.foceiOuter, "rstudio") <- emptyenv()

#' Augmented model for the AGQ quadrature NODES: the same direction set as `foceiOuter`
#' but 1st order only.
#'
#' The nodes read only `f`/`R`/`a`/`aR`/`Rsig` -- they never touch `A`/`AR`/`RsigDir`,
#' because the 2nd-order block is used solely at eta-hat (the exact inner Hessian -> etaP,
#' and dHtD).  Dropping that tier takes the node solve from 26 ODE states to 8 on a
#' one-compartment model, and the nodes are `nAGQ^neta` solves per gradient -- 45-54% of
#' the gradient at neta=3/nAGQ=3 and 77-86% at neta=5.  Measured: 1.07x (neta=3, nAGQ=2)
#' to 1.91x (neta=5, nAGQ=3) on the whole gradient.
#'
#' Only built for nAGQ > 1; every other fast fit gets NULL and pays no extra build.  Like
#' `foceiOuter` this rides in the disk-cached `foceiModel` list, so the extra symengine+gcc
#' pass is paid once per model, not once per session.
#' @noRd
#' @export
rxUiGet.foceiOuterNode <- function(x, ...) {
  .ui <- x[[1]]
  if (!isTRUE(rxode2::rxGetControl(.ui, "fast", FALSE))) {
    return(NULL)
  }
  if (as.integer(rxode2::rxGetControl(.ui, "nAGQ", 1L)) <= 1L) {
    return(NULL)
  }
  .dir <- .foceiOuterDirs(.ui)
  if (is.null(.dir)) {
    return(NULL)
  }
  .foceiAnalyticAugModelDirs(.ui, .dir$dirs, order = 1L)
}
attr(rxUiGet.foceiOuterNode, "rstudio") <- emptyenv()

#' Estimation-scale (Cholesky) Omega-inverse derivatives for the outer gradient's
#' Omega block, from rxode2's `rxSymInvCholEnvCalculate` (rxSymInv.R d.omegaInv /
#' tr.28).  Returns `list(dOi = list(dOmega^{-1}/d theta_omega_k), tr28 = 0.5 *
#' tr(dOmega^{-1}_k Omega), names = <omega parameter names>)`, or `NULL`.
#' @noRd
.foceiEstOmegaDeriv <- function(ui, Om, e = NULL) {
  tryCatch(
    {
      # Build the rxSymInvChol env from the current Omega with the SAME diagonal
      # transform the optimizer uses, so the Cholesky parameter order/scale matches
      # op_focei's Omega slots (rxUiGet.focei builds env$rxInv the same way).
      .diagXform <- rxode2::rxGetControl(ui, "diagXform", "sqrt")
      # A fresh env pays ~60ms of one-time symbolic setup on its first `$d.omegaInv` (a
      # repeat read is free), which made this ~40% of each analytic gradient.  Reuse the
      # fit's persistent `env$rxInv` (C++ keeps it current via setOmegaTheta) -- but only
      # when it is already at this Omega, so a stale env falls back instead of silently
      # returning derivatives at the wrong one.
      .rxInv <- NULL
      if (!is.null(e)) {
        .cand <- tryCatch(get("rxInv", e), error = function(.) NULL)
        if (!is.null(.cand) && rxode2::rxIs(.cand, "rxSymInvCholEnv")) {
          .omChk <- tryCatch(as.matrix(.cand$omega), error = function(.) NULL)
          if (
            !is.null(.omChk) &&
              identical(dim(.omChk), dim(as.matrix(Om))) &&
              isTRUE(all.equal(unname(.omChk), unname(as.matrix(Om)), tolerance = 1e-10))
          ) {
            .rxInv <- .cand
          }
        }
      }
      if (is.null(.rxInv)) {
        .rxInv <- rxode2::rxSymInvCholCreate(mat = Om, diag.xform = .diagXform)
      }
      # `$.rxSymInvCholEnv` dispatches to the C rxSymInvCholEnvCalculate; d.omegaInv
      # is the list of dOmega^-1/d(chol theta_k), tr.28 = 0.5*tr(dOmega^-1_k Omega).
      .dOi <- .rxInv$d.omegaInv
      .tr28 <- .rxInv$tr.28
      if (is.null(.dOi) || is.null(.tr28)) {
        return(NULL)
      }
      list(dOi = .dOi, tr28 = as.numeric(.tr28), names = paste0("om.chol.", seq_along(.dOi)))
    },
    error = function(e) NULL
  )
}

#' Theta names excluded from the outer optimizer's free-parameter set by the
#' mu-referenced (lin/irls) regression -- mirrors inner.cpp isMuGroupSkip: the
#' mu-group thetas plus every mu-group covariate coefficient (bounded ones are
#' regression-updated with clamping too).  Index arrays are 0-based (see
#' `.muRefCppGroupSetup`).
#' @noRd
.foceiMuSkipThetaNames <- function(ui, thNames) {
  if (identical(rxode2::rxGetControl(ui, "muModel", "none"), "none")) {
    return(character(0))
  }
  .g <- as.integer(rxode2::rxGetControl(ui, "foceiMuGroupTheta", integer(0)))
  .ct <- as.integer(rxode2::rxGetControl(ui, "foceiMuGroupCovTheta", integer(0)))
  thNames[c(.g, .ct) + 1L]
}
