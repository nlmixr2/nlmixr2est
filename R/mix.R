#' Get the mixture probabilities from the estimated log-scale parameters
#'
#' @param val numeric vector of the full parameter set in focei
#'
#' @param idx integer vector of the indices of the mixture log-scale
#'   parameters
#'
#' @return A numeric vector of the mixture probabilities
#'
#' @noRd
#'
#' @author Matthew L. Fidler
.getMixFromLog <- function(val, idx) {
  v <- rxode2::mexpit(val[idx])
  c(v, 1-sum(v))
}
#' Get the mixture gradients of the estimated log-scale parameters
#'
#'
#' @param val numeric vector of the full parameter set in focei
#'
#' @param idx vector of the indices of the mixture log-scale
#'  parameters
#'
#' @return A numeric vector of the mixture probabilities
#'
#' @noRd
#'
#' @author Matthew L. Fidler
.getMixJacFromLog <- function(val, idx) {
  rxode2::dmexpit(val[idx])
}

#' Find all mix() calls in a parsed model expression list
#'
#' @param expr A single parsed expression (or sub-expression) from `ui$lstExpr`
#' @return A list of `mix()` call expressions found anywhere in `expr`
#' @noRd
#' @author Matthew L. Fidler
.findMixCalls <- function(expr) {
  if (is.call(expr)) {
    if (identical(expr[[1]], quote(mix))) {
      return(list(expr))
    }
    return(do.call(c, lapply(expr, .findMixCalls)))
  }
  return(NULL)
}

#' Extract ETA names referenced inside a mix() call's component expressions
#'
#' @param expr A parsed expression (or sub-expression) from inside a `mix()` call
#' @param etas Character vector of all known ETA names to match against
#' @return Character vector of ETA names found in `expr` (deduplicated)
#' @noRd
#' @author Matthew L. Fidler
.extractEtas <- function(expr, etas) {
  if (is.name(expr)) {
    .n <- as.character(expr)
    if (.n %in% etas) return(.n)
  } else if (is.call(expr)) {
    return(unique(unlist(lapply(expr[-1], .extractEtas, etas = etas))))
  }
  return(NULL)
}

#' Process mixture model information after a focei fit
#'
#' After the C++ focei fit, strips the MIXEST column from ranef, computes
#' posterior mixture probabilities from etaObfFull and theta priors, and
#' builds mixList (per-mixture ID/ETA/probability) and mixNum (best MIXNUM
#' per subject).
#'
#' @param env Fit environment (the C++ output environment)
#' @param ui rxode2 UI object
#' @return Nothing; modifies env in place for side effects
#' @noRd
#' @author Matthew L. Fidler
.mixFix <- function(env, ui) {
  .mixIdx <- try(get("mixIdx", envir=env), silent=TRUE)
  if (inherits(.mixIdx, "try-error")) return(invisible(NULL))
  if (length(.mixIdx) == 0L) return(invisible(NULL))
  if (!exists("etaObfFull", envir=env)) return(invisible(NULL))

  .etaFull <- get("etaObfFull", envir=env)
  .etaBest <- get("etaObf", envir=env)

  # Fix ranef: remove MIXEST column so nlmixr2Parameters() gets ID + ETAs only
  .ranef <- as.data.frame(get("ranef", envir=env))
  .wMix <- which(names(.ranef) == "MIXEST")
  if (length(.wMix) > 0L) {
    .ranef <- .ranef[, -.wMix, drop=FALSE]
    assign("ranef", .ranef, envir=env)
  }

  # Prefer env$mixProbabilities (back-transformed by preFinalParTableHook,
  # includes implicit last component); else compute from fixef (mlogit scale).
  if (exists("mixProbabilities", envir=env)) {
    .priorProbs <- get("mixProbabilities", envir=env)
  } else {
    .finalTheta <- get("fixef", envir=env)
    .priorProbs <- .getMixFromLog(.finalTheta, .mixIdx)
  }
  .nMix <- length(.priorProbs)
  .nSub <- nrow(.ranef)

  # etaObfFull: columns are "ID", "MIXEST"(1-indexed), "ETA[1]",...,"ETA[neta]", "OBJI"
  # Rows ordered by mixture then subject; sort to ensure consistency
  .etaFull <- .etaFull[order(.etaFull$MIXEST, as.integer(.etaFull$ID)), ]

  .etaCols <- grep("^ETA\\[", names(.etaFull), value=TRUE)
  .etaNames <- names(.ranef)[-1]  # eta names from fixed ranef (no ID, no MIXEST)

  # Build unnormalised posterior: exp(-OBJI/2) * prior_prob for each subject x mixture
  .llikMat <- matrix(NA_real_, nrow=.nSub, ncol=.nMix)
  for (k in seq_len(.nMix)) {
    .wk <- which(.etaFull$MIXEST == k)
    # Sort by subject ID to align rows correctly
    .wk <- .wk[order(as.integer(.etaFull$ID[.wk]))]
    .llikMat[, k] <- exp(-0.5 * .etaFull$OBJI[.wk]) * .priorProbs[k]
  }
  .rowTotals <- rowSums(.llikMat)
  .zeroRows <- which(.rowTotals <= 0 | !is.finite(.rowTotals))
  if (length(.zeroRows) > 0L) {
    warning(sprintf(
      "%d subject(s) had zero/underflowed mixture likelihood in all components; falling back to prior probabilities for those subjects",
      length(.zeroRows)), call. = FALSE)
    .rowTotals[.zeroRows] <- 1
    .llikMat[.zeroRows, ] <- matrix(.priorProbs, nrow = length(.zeroRows), ncol = .nMix, byrow = TRUE)
  }

  # Create mixList: one data frame per mixture component
  .mixList <- lapply(seq_len(.nMix), function(k) {
    .wk <- which(.etaFull$MIXEST == k)
    .wk <- .wk[order(as.integer(.etaFull$ID[.wk]))]
    .df <- .etaFull[.wk, .etaCols, drop=FALSE]
    .prob <- .llikMat[, k] / .rowTotals
    .ret <- cbind(data.frame(ID=.etaFull$ID[.wk]), .df, data.frame(prob=.prob))
    names(.ret) <- c("ID", .etaNames, "prob")
    row.names(.ret) <- NULL
    .ret
  })
  names(.mixList) <- paste0("mix", seq_len(.nMix))

  # Create mixNum: best mixture assignment per subject (1-indexed)
  # etaObf has columns: "ID", "MIXEST"(1-indexed best mix), eta names..., "OBJI"
  .wMix2 <- which(names(.etaBest) == "MIXEST")
  .mixNum <- data.frame(ID=.etaBest$ID,
                        mixnum=if (length(.wMix2) > 0L) as.integer(.etaBest[[.wMix2]]) else NA_integer_)
  row.names(.mixNum) <- NULL

  assign("mixList", .mixList, envir=env)
  assign("mixNum", .mixNum, envir=env)

  # Calculate Expected ETAs for shrinkage
  .etaExpected <- .etaFull[.etaFull$MIXEST == 1, .etaCols, drop=FALSE]
  names(.etaExpected) <- .etaNames
  for (.n in names(.etaExpected)) .etaExpected[[.n]] <- 0
  for (.m in .mixList) {
    for (.n in names(.etaExpected)) {
      if (.n %in% names(.m)) {
        .etaExpected[[.n]] <- .etaExpected[[.n]] + .m[[.n]] * .m$prob
      }
    }
  }
  .etaExpected <- cbind(data.frame(ID=.etaFull$ID[.etaFull$MIXEST == 1]), .etaExpected)
  assign("etaExpected", .etaExpected, envir=env)

  # iCov drives rxode2's per-subject mixture fixing during solve/table calc;
  # ID must be integer to match the data's ID column type.
  .iCov <- data.frame(ID=as.integer(.mixNum$ID), mixest=.mixNum$mixnum)
  assign("mixIcov", .iCov, envir=env)

  invisible(NULL)
}

#' Process mixture model information after a SAEM fit
#'
#' SAEM analogue of `.mixFix()`: builds `mixList` (per-mixture ID/ETA/
#' probability), `mixNum` (best mixture assignment), `mixIcov` (for rxode2's
#' mixture fixing during solve/table calc), and `mixProbabilities` (full
#' nMix-length vector for `.mixFixTable()`), all from the `mixWeights` matrix
#' already computed by the SAEM C++ engine (`env$saem$mixWeights`) -- unlike
#' `.mixFix()`, no `etaObfFull` is needed.
#'
#' @param env Fit environment (the SAEM output environment, before
#'   `nlmixr2CreateOutputFromUi`)
#' @param ui rxode2 UI object
#' @return Nothing; modifies `env` in place for side effects
#' @noRd
#' @author Matthew L. Fidler
.saemMixFix <- function(env, ui) {
  if (length(ui$mixProbs) == 0L) return(invisible(NULL))
  .saem <- env$saem
  if (is.null(.saem)) return(invisible(NULL))
  .mixWeights <- .saem$mixWeights  # N x nMix matrix of posterior weights
  if (is.null(.mixWeights) || nrow(.mixWeights) == 0L) return(invisible(NULL))
  .nMix <- ncol(.mixWeights)
  if (.nMix < 2L) return(invisible(NULL))

  # etaObf was populated by .getSaemOmega; columns: ID, eta names, OBJI
  .etaObf <- env$etaObf
  if (is.null(.etaObf) || nrow(.etaObf) == 0L) return(invisible(NULL))
  .nSub <- nrow(.etaObf)

  # eta column names (exclude ID and OBJI)
  .etaNames <- names(.etaObf)[!(names(.etaObf) %in% c("ID", "OBJI"))]

  # mixWeights rows correspond to subject order in etaObf
  # Ensure the matrix has a row for every subject
  if (nrow(.mixWeights) != .nSub) {
    warning("mixWeights row count doesn't match number of subjects; skipping SAEM mixFix",
            call.=FALSE)
    return(invisible(NULL))
  }

  .bestMix <- apply(.mixWeights, 1L, which.max)

  # Final mixture probabilities (full simplex, nMix elements)
  .mixProb <- .saem$mixProb
  if (length(.mixProb) == .nMix - 1L) {
    .mixProbabilities <- c(.mixProb, 1.0 - sum(.mixProb))
  } else if (length(.mixProb) == .nMix) {
    .mixProbabilities <- .mixProb
  } else {
    .mixProbabilities <- rep(1.0 / .nMix, .nMix)
  }
  env$mixProbabilities <- .mixProbabilities

  .allEtas <- ui$iniDf[!is.na(ui$iniDf$neta1), ]
  .allEtas <- .allEtas[.allEtas$neta1 == .allEtas$neta2, "name"]
  .mixCalls <- do.call(c, lapply(ui$lstExpr, .findMixCalls))
  
  .etaGroups <- list()
  for (.mc in .mixCalls) {
    .args <- as.list(.mc)[-1]
    .comps <- .args[seq(1, length(.args), by = 2)]
    .grpEtas <- unique(unlist(lapply(.comps, .extractEtas, etas = .allEtas)))
    if (length(.grpEtas) > 1L) {
      .etaGroups <- c(.etaGroups, list(.grpEtas))
    }
  }

  .omega <- env$omega
  .fixef <- env$fixef
  .muRef <- ui$muRefDataFrame
  
  if (length(.etaGroups) > 0L) {
    for (.grp in .etaGroups) {
      .rootName <- gsub("[0-9]+$", "", .grp[1])
      
      .sig02 <- .omega[.grp[1], .grp[1]]
      .thetas <- vapply(.grp, function(e) {
        .t <- .muRef$theta[.muRef$eta == e]
        if (length(.t) == 1L) .t else NA_character_
      }, character(1))
      
      .mus <- .fixef[.thetas]
      .mus[is.na(.mus)] <- 0.0
      
      .wGroup <- .mixProbabilities
      .meanMu <- sum(.wGroup * .mus)
      .overallVar <- .sig02
      
      .wIdx <- which(colnames(.omega) == .grp[1])
      if (length(.wIdx) == 1L) {
        colnames(.omega)[.wIdx] <- rownames(.omega)[.wIdx] <- .rootName
        .omega[.rootName, .rootName] <- .overallVar
      }
      
      .toRemove <- .grp[-1]
      .omega <- .omega[!(rownames(.omega) %in% .toRemove), !(colnames(.omega) %in% .toRemove), drop=FALSE]
      
      .etaObf[[.rootName]] <- vapply(seq_len(nrow(.etaObf)), function(i) {
        .etaObf[i, .grp[.bestMix[i]]]
      }, numeric(1))
      .etaObf <- .etaObf[, !(names(.etaObf) %in% .grp), drop=FALSE]

      .updateMat <- function(mat) {
        .dfMat <- as.data.frame(mat)
        .N <- nrow(.dfMat)
        .newCol <- vapply(seq_len(.N), function(i) {
          .subjIdx <- ((i - 1) %% .nSub) + 1
          .dfMat[i, .grp[.bestMix[.subjIdx]]]
        }, numeric(1))
        .dfMat[[.rootName]] <- .newCol
        .dfMat <- .dfMat[, !(names(.dfMat) %in% .grp), drop=FALSE]
        as.matrix(.dfMat)
      }
      if (exists(".etaMatBase", envir=env, inherits=FALSE) && !is.null(env$.etaMatBase)) {
        env$.etaMatBase <- .updateMat(env$.etaMatBase)
      }
      if (exists(".etaMat", envir=env, inherits=FALSE) && !is.null(env$.etaMat)) {
        env$.etaMat <- .updateMat(env$.etaMat)
      }
    }
    .funLines <- deparse(as.function(ui))
    for (.grp in .etaGroups) {
      .rootName <- gsub("[0-9]+$", "", .grp[1])
      for (.comp in .grp) {
        .funLines <- gsub(paste0("\\b", .comp, "\\b"), .rootName, .funLines)
      }
      .etaClLines <- grep(paste0("\\b", .rootName, "\\s*~"), .funLines)
      if (length(.etaClLines) > 1) {
        .funLines <- .funLines[-.etaClLines[-1]]
      }
    }
    .funText <- paste(.funLines, collapse="\n")
    .funNew <- eval(parse(text=.funText))
    .uiNew <- rxode2::rxode2(.funNew)
    if (exists("boundedTransforms", envir=ui$meta)) {
      assign("boundedTransforms", get("boundedTransforms", envir=ui$meta), envir=.uiNew$meta)
    }
    env$ui <- .uiNew
    env$omega <- .omega
    env$etaObf <- .etaObf
    .etaNames <- names(.etaObf)[!(names(.etaObf) %in% c("ID", "OBJI"))]
  }

  # Create mixList: one data frame per mixture component
  .mixList <- lapply(seq_len(.nMix), function(k) {
    .df <- as.data.frame(.etaObf[, .etaNames, drop=FALSE])
    .prob <- .mixWeights[, k]
    .ret <- cbind(data.frame(ID=.etaObf$ID), .df, data.frame(prob=.prob))
    names(.ret) <- c("ID", .etaNames, "prob")
    row.names(.ret) <- NULL
    .ret
  })
  names(.mixList) <- paste0("mix", seq_len(.nMix))

  # Create mixNum: best mixture assignment per subject (1-indexed)
  .mixNum <- data.frame(ID=.etaObf$ID,
                        mixnum=as.integer(.bestMix))
  row.names(.mixNum) <- NULL

  # Assign ranef:
  .ranef <- as.data.frame(.etaObf[, .etaNames, drop=FALSE])
  .ranef <- cbind(data.frame(ID=.etaObf$ID), .ranef)
  .ranef$mixnum <- as.integer(.bestMix)
  assign("ranef", .ranef, envir=env)

  assign("mixList", .mixList, envir=env)
  assign("mixNum", .mixNum, envir=env)

  # Calculate Expected ETAs for shrinkage
  .etaExpected <- as.data.frame(.etaObf[, .etaNames, drop=FALSE])
  for (.n in names(.etaExpected)) .etaExpected[[.n]] <- 0
  for (.m in .mixList) {
    for (.n in names(.etaExpected)) {
      if (.n %in% names(.m)) {
        .etaExpected[[.n]] <- .etaExpected[[.n]] + .m[[.n]] * .m$prob
      }
    }
  }
  .etaExpected <- cbind(data.frame(ID=.etaObf$ID), .etaExpected)
  assign("etaExpected", .etaExpected, envir=env)

  # Store iCov for the table/solve step: rxode2 reads 'mixest' from iCov to
  # fix each individual's mixture component during ODE solving.
  .iCov <- data.frame(ID=as.integer(.mixNum$ID), mixest=.mixNum$mixnum)
  assign("mixIcov", .iCov, envir=env)

  invisible(NULL)
}

#' Back-transform mixture probability columns in the "Back-Transformed" rows of parHistData
#'
#' C++ applies exp/expit for known transforms in "Back-Transformed" rows but
#' leaves mixture (mlogit-scale) parameters raw; this corrects just those
#' columns via \code{rxode2::mexpit()}. "Scaled"/"Unscaled" rows are left
#' untouched at their mlogit values.
#'
#' @param parHist data frame returned by C++ \code{parHistData()}
#' @param mixColNames character vector of column names to back-transform
#' @return \code{parHist} with mixture columns corrected in "Back-Transformed" rows
#' @noRd
#' @author Matthew L. Fidler
.backTransformParHistMix <- function(parHist, mixColNames) {
  .mixCols <- match(mixColNames, names(parHist))
  .mixCols <- .mixCols[!is.na(.mixCols)]
  if (length(.mixCols) == 0L) return(parHist)
  .btRows <- as.character(parHist$type) == "Back-Transformed"
  if (!any(.btRows)) return(parHist)
  .mlogitMat <- as.matrix(parHist[.btRows, .mixCols, drop=FALSE])
  .nBt <- sum(.btRows)
  .nMix <- length(.mixCols)
  # matrix() with explicit dims handles single-row/single-column edge cases
  # that apply()+t() alone mishandle.
  parHist[.btRows, .mixCols] <- matrix(
    t(apply(.mlogitMat, 1L, function(.row) rxode2::mexpit(.row)[seq_len(.nMix)])),
    nrow = .nBt, ncol = .nMix
  )
  parHist
}

#' Rotate a covariance's mixture-proportion block onto the probability scale
#'
#' The proportions are estimated as a multinomial logit, so the covariance a
#' covariance method produces is on that scale while the reported estimate is a
#' probability.  \code{p = mexpit(t)} is a softmax over the \code{K-1} free
#' coordinates, whose Jacobian is \code{dp_j/dt_l = p_j*(delta_jl - p_l)}, i.e.
#' \code{J = diag(p) - p \%*\% t(p)}.  The FULL Jacobian is used, not its
#' diagonal: the off-diagonal terms are what carry the proportions'
#' cross-covariances -- with each other and with the structural thetas -- onto
#' the reported scale.
#'
#' Same principle as \code{covFull}, which reports Omega on the natural
#' variance scale rather than the \code{chol(solve(omega))} estimation scale.
#'
#' @param cov covariance matrix with dimnames
#' @param mixNames mixture-proportion parameter names, in THETA-slot order
#' @param p free mixture probabilities, in \code{mixNames} order
#' @return \code{cov} with the mixture rows/columns on the probability scale;
#'   unchanged when there is no mixture block to rotate
#' @noRd
#' @author Matthew L. Fidler
.mixCovToProbScale <- function(cov, mixNames, p) {
  if (!is.matrix(cov) || length(mixNames) == 0L) return(cov)
  if (is.null(rownames(cov))) return(cov)
  .i <- match(mixNames, rownames(cov))
  if (anyNA(.i) || length(p) != length(.i) || !all(is.finite(p))) return(cov)
  .J <- diag(p, nrow = length(p)) - outer(p, p)
  .A <- diag(1, nrow(cov))
  .A[.i, .i] <- .J
  .out <- .A %*% cov %*% t(.A)
  dimnames(.out) <- dimnames(cov)
  .out
}

#' Rotate a fit's installed covariance onto the probability scale
#'
#' The pre-final table hook rotates \code{env$cov}, but the \code{covFull} and
#' analytic installers (\code{.foceiInstallFdFullCov} /
#' \code{.foceiInstallAnalyticCov}) then REPLACE it wholesale with a matrix
#' still on the mlogit estimation scale, so the hook's rotation is gone by the
#' time \code{.updateParFixed()} reads it.  This runs between the two and covers
#' whichever matrix ended up installed, plus the \code{covR}/\code{covS}/
#' \code{covRS} diagnostics so they stay on one scale.
#'
#' This is the ONLY fit-time rotation, so no matrix is rotated twice: the native
#' C++ covariance and both installers' output all arrive here on the mlogit
#' scale.  Post-fit \code{setCov()} installs rotate separately, in
#' \code{.covInstallResult()}, on a matrix the recompute engine produced.
#'
#' @param env fit environment
#' @return invisible \code{NULL}; called for its side effects on \code{env}
#' @noRd
#' @author Matthew L. Fidler
.mixInstallProbScaleCov <- function(env) {
  # a covariance handed in by the caller is already on the reported scale
  if (isTRUE(tryCatch(get(".mixCovPreRotated", envir = env, inherits = FALSE),
                      error = function(e) FALSE))) {
    return(invisible(NULL))
  }
  .mix <- .mixEnvPieces(env)
  if (is.null(.mix)) return(invisible(NULL))
  .mp <- .mix$names
  .p <- .mix$p
  for (.n in c("cov", "covR", "covS", "covRS")) {
    if (!exists(.n, envir = env, inherits = FALSE)) next
    .cur <- get(.n, envir = env)
    if (!is.matrix(.cur)) next
    assign(.n, .mixCovToProbScale(.cur, .mp, .p), envir = env)
  }
  .mixRefreshSeFromCov(env, .mp, .mix$idx)
  invisible(NULL)
}

#' Refresh the mixture rows of se/popDf from the rotated covariance
#'
#' \code{foceiFinalizeTables} fills \code{se}/\code{popDf} from the covariance
#' as it stood BEFORE the probability-scale rotation, and
#' \code{.updateParFixed()} reads \code{popDf} -- so without this the reported
#' SE stays on the mlogit estimation scale while the estimate beside it is a
#' probability.
#'
#' @param env fit environment
#' @param mixNames mixture-proportion parameter names, in THETA-slot order
#' @param mixIdx their positions in the theta vector, same order as
#'   \code{mixNames}
#' @return invisible \code{NULL}; called for its side effects on \code{env}
#' @noRd
#' @author Matthew L. Fidler
.mixRefreshSeFromCov <- function(env, mixNames, mixIdx) {
  .cov <- tryCatch(get("cov", envir = env, inherits = FALSE), error = function(e) NULL)
  .mixIdx <- mixIdx
  if (!is.matrix(.cov) || is.null(rownames(.cov)) ||
        is.null(.mixIdx) || length(.mixIdx) != length(mixNames)) {
    return(invisible(NULL))
  }
  .w <- match(mixNames, rownames(.cov))
  if (anyNA(.w)) return(invisible(NULL))
  .newSe <- sqrt(diag(.cov))[.w]
  if (exists("se", envir = env, inherits = FALSE)) {
    .se <- get("se", envir = env)
    if (length(.se) >= max(.mixIdx)) {
      .se[.mixIdx] <- .newSe
      assign("se", .se, envir = env)
    }
  }
  if (!exists("popDf", envir = env, inherits = FALSE)) return(invisible(NULL))
  .pd <- get("popDf", envir = env)
  if (!is.data.frame(.pd) || nrow(.pd) < max(.mixIdx) || !("SE" %in% names(.pd))) {
    return(invisible(NULL))
  }
  .pd[["SE"]][.mixIdx] <- .newSe
  if ("%RSE" %in% names(.pd)) {
    .e <- .pd[["Estimate"]][.mixIdx]
    .pd[["%RSE"]][.mixIdx] <-
      ifelse(is.finite(.e) & .e != 0, abs(.newSe / .e) * 100, NA_real_)
  }
  assign("popDf", .pd, envir = env)
  invisible(NULL)
}

#' Read a fit environment's mixture pieces, or NULL if it has none usable
#'
#' Both covariance consumers need the same three things off a fit env -- the
#' proportion parameter names, the free probabilities, and the per-subject
#' responsibility matrix -- with the same consistency checks between them.
#'
#' @param env fit environment
#' @param needResp when \code{TRUE} also require \code{$mixList} and return the
#'   responsibility matrix
#' @return list with \code{names}, \code{p} (free probabilities) and, when
#'   requested, \code{r} (subjects x components); \code{NULL} if unavailable
#' @noRd
#' @author Matthew L. Fidler
.mixEnvPieces <- function(env, needResp = FALSE) {
  .ui <- tryCatch(env$ui, error = function(e) NULL)
  if (is.null(.ui)) return(NULL)
  .idx <- tryCatch(.ui$thetaMixIndex, error = function(e) NULL)
  if (is.null(.idx) || length(.idx) == 0L) return(NULL)
  # Name the proportions by their THETA slot, not by ui$mixProbs.  The two are
  # the same set but NOT the same order: mixProbs follows the mix() call while
  # thetaMixIndex follows ini(), and everything downstream -- the covariance's
  # rows, op_focei.mixProb, $mixProbabilities, and the se/popDf rows -- is keyed
  # on the theta slot.  Using mixProbs to index a theta-ordered covariance puts
  # the Jacobian on the wrong rows whenever ini() lists them in a different
  # order than mix() uses them.
  .mp <- tryCatch(names(.ui$theta)[.idx], error = function(e) NULL)
  if (is.null(.mp) || length(.mp) != length(.idx) || anyNA(.mp)) return(NULL)
  .pi <- tryCatch(env$mixProbabilities, error = function(e) NULL)
  if (is.null(.pi) || length(.pi) != length(.mp) + 1L || !all(is.finite(.pi))) return(NULL)
  .ret <- list(names = .mp, idx = .idx, p = .pi[seq_along(.mp)], pi = .pi)
  if (!needResp) return(.ret)
  .ml <- tryCatch(env$mixList, error = function(e) NULL)
  if (is.null(.ml) || length(.ml) != length(.pi)) return(NULL)
  .r <- try(do.call(cbind, lapply(.ml, function(.z) .z$prob)), silent = TRUE)
  if (inherits(.r, "try-error") || !is.matrix(.r) || ncol(.r) != length(.pi)) return(NULL)
  .ret$r <- .r
  .ret
}

#' Probability-scale covariance of the mixture proportions from responsibilities
#'
#' NONMEM 7 Technical Guide eq. (7.51): the mixture parameters' information is
#' the outer product of the per-subject scores, \code{sum_i (r_i - p)(r_i - p)'}
#' on the mlogit scale.  Its inverse is rotated onto the probability scale with
#' the same full Jacobian every other method uses.
#'
#' @param r matrix of per-subject responsibilities, one column per FREE component
#' @param p free mixture probabilities
#' @return the probability-scale covariance block, or \code{NULL} if it is not
#'   invertible / not a usable covariance
#' @noRd
#' @author Matthew L. Fidler
.mixProbCovBlock <- function(r, p) {
  .d <- sweep(r, 2, p, "-")
  .blk <- try(solve(t(.d) %*% .d), silent = TRUE)
  if (inherits(.blk, "try-error") || !all(is.finite(.blk))) return(NULL)
  .j <- diag(p, nrow = length(p)) - outer(p, p)
  .blk <- .j %*% .blk %*% t(.j)
  if (!all(is.finite(.blk)) || any(diag(.blk) <= 0)) return(NULL)
  .blk
}

#' Append a mixture-proportion block to a covariance that has none
#'
#' \code{saem} excludes the mixture proportions from its kernel parameter vector
#' (they are updated by a separate EM step), so its Louis/linFim covariance has
#' no mixture rows at all and \code{p1} reports \code{SE = NA}.
#'
#' The block is NONMEM 7 Technical Guide eq. (7.51): the mixture parameters'
#' information is the outer product of the per-subject scores,
#' \code{sum_i (r_i - p)(r_i - p)'} on the mlogit scale, with \code{r_i} the
#' subject's posterior responsibilities -- which the fit already carries in
#' \code{$mixList}.  Its inverse is rotated onto the probability scale with the
#' same full Jacobian every other method uses.
#'
#' The cross terms (7.52)-(7.54) are NOT formed: they need per-subject scores for
#' the other parameters on the same footing, which the SAEM covariance does not
#' expose.  The appended block is therefore uncorrelated with the structural
#' parameters, so these SEs ignore that correlation and are mildly optimistic.
#' Measured against a focei fit of the same data, where the cross terms ARE
#' available, the difference is a couple of percent.
#'
#' @param env fit environment
#' @return invisible \code{NULL}; called for its side effects on \code{env}
#' @noRd
#' @author Matthew L. Fidler
.mixCovAppendBlock <- function(env) {
  .mix <- .mixEnvPieces(env, needResp = TRUE)
  if (is.null(.mix)) return(invisible(NULL))
  .mp <- .mix$names
  .pi <- .mix$pi
  .r <- .mix$r
  .cov <- tryCatch(get("cov", envir = env, inherits = FALSE), error = function(e) NULL)
  if (!is.matrix(.cov) || is.null(rownames(.cov))) return(invisible(NULL))
  if (any(.mp %in% rownames(.cov))) return(invisible(NULL))   # already covered
  .free <- seq_along(.mp)
  # An information matrix reports the precision of a MAXIMUM-likelihood estimate.
  # The mixture score is s_l = sum_i (r_il - p_l), so the fixed point is s == 0;
  # away from it the block is a confident-looking number attached to an estimate
  # that is not an MLE.  Refuse rather than report it, and say why -- saem can
  # land far off this (its proportions are updated by a separate EM step,
  # outside the kernel that converged everything else).
  #
  # Judge it by the SCORE STATISTIC s' solve(I) s, not by |mean(r) - p|: the
  # score is N*(mean(r) - p), so any absolute tolerance on the mean silently
  # loosens with the number of subjects (0.01 is a score of 1 at N=100 and 100
  # at N=10000).  s' I^-1 s is on a chi-square scale and does not drift with N.
  .d <- sweep(.r[, .free, drop = FALSE], 2, .pi[.free], "-")
  .s <- colSums(.d)
  .stat <- tryCatch(as.numeric(crossprod(.s, solve(crossprod(.d), .s))),
                    error = function(e) NA_real_)
  if (!is.finite(.stat) || .stat > 1e-3) {
    # warning(), not an assignment to runInfo: that is the channel the fit
    # collects run-time notes through, and a direct assignment here is
    # overwritten by the later table assembly.  Kept under 75 characters so it
    # renders on one line, and unprefixed -- the fit already reports its method.
    warning("mixture proportion SE skipped; not at the score-zero point",
            call. = FALSE)
    return(invisible(NULL))
  }
  .blk <- .mixProbCovBlock(.r[, .free, drop = FALSE], .pi[.free])
  if (is.null(.blk)) return(invisible(NULL))
  .n <- nrow(.cov)
  .out <- matrix(0, .n + length(.mp), .n + length(.mp))
  .out[seq_len(.n), seq_len(.n)] <- .cov
  .out[.n + .free, .n + .free] <- .blk
  .nm <- c(rownames(.cov), .mp)
  dimnames(.out) <- list(.nm, .nm)
  assign("cov", .out, envir = env)
  .updateParFixedRefreshSeFromCov(env, .out, onlyMissing = TRUE)
  invisible(NULL)
}

#' Pre-final parameter table hook: back-transform mixture probability parameters
#'
#' Registered via \code{preFinalParTableHooksAdd()}; converts mixture
#' probability parameters from mlogit scale to natural probability scale in
#' \code{env$theta$theta}. The full probability vector (including the implicit
#' last component) is stored in \code{env$mixProbabilities}, which
#' \code{.mixFix()} and \code{.mixInstallProbScaleCov()} both read.  The
#' covariance is rotated onto the probability scale later, by
#' \code{.mixInstallProbScaleCov()} -- not here, because the covFull/analytic
#' installers replace \code{env$cov} after this hook runs.
#'
#' @param env Fit environment containing \code{mixIdx} and \code{theta}
#' @return invisible \code{NULL}; called for its side effects on \code{env}
#' @noRd
#' @author Matthew L. Fidler
.aaaPostEstimationMixBacktransform <- function(env) {
  .mixIdx <- try(get("mixIdx", envir=env), silent=TRUE)
  if (inherits(.mixIdx, "try-error")) return(invisible(NULL))
  if (length(.mixIdx) == 0L) return(invisible(NULL))

  .thetaDf <- env$theta
  if (is.null(.thetaDf) || !is.data.frame(.thetaDf)) return(invisible(NULL))

  .mlogitVals <- .thetaDf$theta[.mixIdx]
  .probs <- rxode2::mexpit(.mlogitVals)

  .thetaDf$theta[.mixIdx] <- .probs[seq_along(.mixIdx)]
  env$theta <- .thetaDf
  # append implicit last component: mexpit only returns the n explicit
  # probabilities, and .mixFix uses length() as nMix.
  env$mixProbabilities <- c(.probs, 1 - sum(.probs))

  if (exists("parHistData", envir=env) && exists("thetaNames", envir=env)) {
    .phd <- env$parHistData
    if (is.data.frame(.phd) && nrow(.phd) > 0L) {
      env$parHistData <- .backTransformParHistMix(.phd, env$thetaNames[.mixIdx])
    }
  }

  invisible(NULL)
}
preFinalParTableHooksAdd(".aaaPostEstimationMixBacktransform", .aaaPostEstimationMixBacktransform)

#' Fix mixture LHS variables in the assembled fit table
#'
#' Safety fallback: replaces me/mn/mu with values from mixNum (me/mn) and
#' 1/nMix (mu), since older rxode2 versions silently reject iCov's mixest
#' and leave these columns as 0.
#'
#' @param fit nlmixr2FitData object (after addTable)
#' @param env fit environment
#' @param ui rxode2 UI object
#' @return modified fit (or fit unchanged for non-mixture models)
#' @noRd
#' @author Matthew L. Fidler
.mixFixTable <- function(fit, env, ui) {
  if (!inherits(fit, "nlmixr2FitData")) return(fit)
  if (!exists("mixNum", envir=env)) return(fit)
  .mn <- get("mixNum", envir=env)
  if (is.null(.mn) || nrow(.mn) == 0L) return(fit)
  .nMix <- length(ui$mixProbs) + 1L  # nMix = n_explicit_probs + 1
  # me and mn: best-fit mixture per individual (1-indexed)
  if ("me" %in% names(fit)) {
    .meMap <- setNames(as.integer(.mn$mixnum), as.integer(.mn$ID))
    fit[["me"]] <- .meMap[as.integer(fit[["ID"]])]
  }
  if ("mn" %in% names(fit)) {
    .meMap <- setNames(as.integer(.mn$mixnum), as.integer(.mn$ID))
    fit[["mn"]] <- .meMap[as.integer(fit[["ID"]])]
  }
  # mu: uniform mixture probability = 1/nMix (constant)
  if ("mu" %in% names(fit) && .nMix > 0L) {
    fit[["mu"]] <- 1.0 / .nMix
  }
  fit
}

#' @export
rxUiGet.thetaIniMix <- function(x, ...) {
  .ui <- x[[1]]
  .theta <- .ui$theta
  if (length(.ui$mixProbs) > 0) {
    .p <- .theta[.ui$mixProbs]
    if (all(.p >= 0 & .p <= 1) && sum(.p) <= 1.0) {
      .theta[.ui$mixProbs] <- rxode2::mlogit(.p)
    } else {
      # stop() (not warning()): an invalid value would silently reach
      # focei's initial parameter vector on the wrong scale with a
      # confusing downstream error, and warnings are dropped by
      # .collectWarn() when the fit ultimately errors out.
      stop("initial mixture probabilities are invalid (must each be in [0, 1] ",
           "and sum to no more than 1): ", paste(signif(.p, 3), collapse = ", "),
           call. = FALSE)
    }
  }
  .theta
}
attr(rxUiGet.thetaIniMix, "rstudio") <- stats::setNames(1, "a")

#' @export
rxUiGet.thetaMixIndex <- function(x, ...) {
  .ui <- x[[1]]
  .theta <- .ui$theta
  if (length(.ui$mixProbs) > 0) {
    which(names(.ui$theta) %in% .ui$mixProbs)
  } else {
    integer(0)
  }
}
attr(rxUiGet.thetaMixIndex, "rstudio") <- 1L
