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

#' Process mixture model information after a focei fit
#'
#' After the C++ focei fit completes, this function:
#'  1. Removes the MIXEST column from ranef so nlmixr2Parameters() works correctly
#'  2. Computes posterior mixture probabilities from etaObfFull and theta priors
#'  3. Creates mixList (one data frame per mixture with ID, ETAs, and probability)
#'  4. Creates mixNum (data frame with ID and best MIXNUM per subject)
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

  # Get final mixture prior probabilities.
  # If the preFinalParTableHook ran, env$mixProbabilities already holds the
  # back-transformed probability vector (including the implicit last component).
  # Otherwise fall back to computing from fixef (mlogit scale).
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

  # Store iCov for the table/solve step: rxode2 reads 'mixest' from iCov to
  # fix each individual's mixture component during ODE solving.
  # ID must be integer to match the data's ID column type.
  .iCov <- data.frame(ID=as.integer(.mixNum$ID), mixest=.mixNum$mixnum)
  assign("mixIcov", .iCov, envir=env)

  invisible(NULL)
}

#' Process mixture model information after a SAEM fit
#'
#' SAEM-specific analogue of `.mixFix()`.  After the SAEM fit completes and
#' `.getSaemOmega()` has populated `env$etaObf`, this function:
#'  1. Builds `mixList` (one data frame per mixture component with ID, ETAs, and
#'     posterior probability) from the `mixWeights` matrix returned by the C++
#'     SAEM engine.
#'  2. Builds `mixNum` (data frame with ID and most-likely mixture assignment).
#'  3. Builds `mixIcov` (data frame with ID and `mixest` column used by rxode2
#'     to fix the mixture component during ODE solving / table calculation).
#'  4. Stores `mixProbabilities` (full nMix-length probability vector including
#'     the implicit last component) so that `.mixFixTable()` can correct me/mn/mu.
#'
#' Unlike `.mixFix()`, this function does not require `etaObfFull` because the
#' per-subject posterior mixture weights are already computed by the SAEM
#' engine and returned in `env$saem$mixWeights`.
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

  .findMixCalls <- function(expr) {
    if (is.call(expr)) {
      if (identical(expr[[1]], quote(mix))) {
        return(list(expr))
      }
      return(do.call(c, lapply(expr, .findMixCalls)))
    }
    return(NULL)
  }
  
  .extractEtas <- function(expr, etas) {
    if (is.name(expr)) {
      .n <- as.character(expr)
      if (.n %in% etas) return(.n)
    } else if (is.call(expr)) {
      return(unique(unlist(lapply(expr[-1], .extractEtas, etas = etas))))
    }
    return(NULL)
  }

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
      
      .w_group <- .mixProbabilities
      .mean_mu <- sum(.w_group * .mus)
      .overall_var <- .sig02
      
      .w_idx <- which(colnames(.omega) == .grp[1])
      if (length(.w_idx) == 1L) {
        colnames(.omega)[.w_idx] <- rownames(.omega)[.w_idx] <- .rootName
        .omega[.rootName, .rootName] <- .overall_var
      }
      
      .to_remove <- .grp[-1]
      .omega <- .omega[!(rownames(.omega) %in% .to_remove), !(colnames(.omega) %in% .to_remove), drop=FALSE]
      
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
#' Each FOCEI iteration records three row types in \code{parHistData}: "Scaled",
#' "Unscaled" (shown in \code{fit$parHist}), and "Back-Transformed".  C++ applies
#' \code{exp}/\code{expit} for known transforms in the last row type but falls
#' through to the raw unscaled value for mixture parameters (mlogit scale).  This
#' function corrects only the "Back-Transformed" rows by applying
#' \code{rxode2::mexpit()} jointly to all mixture parameter columns.  "Scaled" and
#' "Unscaled" rows are intentionally left at their mlogit values.
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
  # apply() + t() transposes correctly only when nrow > 1 and ncol > 1.
  # Wrap in matrix() with explicit dimensions to handle the single-column
  # (one mixture parameter) and single-row edge cases uniformly.
  parHist[.btRows, .mixCols] <- matrix(
    t(apply(.mlogitMat, 1L, function(.row) rxode2::mexpit(.row)[seq_len(.nMix)])),
    nrow = .nBt, ncol = .nMix
  )
  parHist
}

#' Pre-final parameter table hook: back-transform mixture probability parameters
#'
#' Registered via \code{preFinalParTableHooksAdd()}. Runs after estimation
#' but before the final parameter table is built (called from C++
#' \code{foceiFinalizeTables} via \code{.preFinalParTableHooksRun}).
#'
#' Converts mixture probability parameters from the mlogit scale (used
#' internally during estimation) to the natural probability scale by modifying
#' \code{env$theta$theta} directly.  Because \code{skipCov} is already set to
#' \code{TRUE} for these indices (see \code{rxUiGet.foceiSkipCov}), the C++
#' table builder will leave SE and \%RSE as \code{NA} without any additional
#' patching.  The full probability vector (including the implicit last
#' component) is stored in \code{env$mixProbabilities} for use by
#' \code{.mixFix()}.
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
  # mexpit returns only the n explicit probabilities (not the implicit last component).
  # .mixFix uses length(mixProbabilities) as nMix, so we must append the last component.
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
#' With rxode2 >= 5.1.2 and the predOnlyModel built from the pruned mix() form,
#' iCov with mixest sets ind->mixest correctly so me/mn/mu come from the solve.
#' For older rxode2 versions (before the lName/liName typo fix in etTran.cpp),
#' iCov is silently rejected and me/mn/mu are 0.  This function replaces those
#' columns with correct values from mixNum (me/mn) and 1/nMix (mu) as a
#' safety fallback that works regardless of the rxode2 version.
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
      # A plain warning() here is not enough: this value feeds directly into
      # focei's initial parameter vector (env$thetaIni, R/focei.R), and
      # leaving it untransformed and on the wrong scale reliably crashes the
      # fit downstream with a confusing, unrelated error (e.g. "infinite
      # while evaluating initial objective function") rather than a clear
      # diagnostic -- confirmed empirically. Also, warnings raised during
      # estimation are collected by .collectWarn() (R/nlmixr2Est.R) and are
      # silently dropped whenever the fit ultimately errors out, so a
      # warning() would never even reach the user in exactly this failure
      # case. stop() early, before any of that, with a clear message.
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
