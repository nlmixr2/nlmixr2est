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

#' Pre-final parameter table hook: back-transform mixture probability parameters
#'
#' Registered via \code{preFinalParTableHooksAdd()}; converts mixture
#' probability parameters from mlogit scale to natural probability scale in
#' \code{env$theta$theta}. SE/\%RSE stay \code{NA} (skipCov already set for
#' these indices); the full probability vector (including the implicit last
#' component) is stored in \code{env$mixProbabilities} for \code{.mixFix()}.
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
