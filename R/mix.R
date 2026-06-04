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
  v <- do.call("rxode2::mexpit", list(val[idx]))
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
  do.call("rxode2::dmexpit", list(val[idx]))
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

  # Get final mixture prior probabilities from estimated theta
  .finalTheta <- get("fixef", envir=env)
  .priorProbs <- .getMixFromLog(.finalTheta, .mixIdx)
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
                        MIXNUM=if (length(.wMix2) > 0L) .etaBest[[.wMix2]] else NA_integer_)
  row.names(.mixNum) <- NULL

  assign("mixList", .mixList, envir=env)
  assign("mixNum", .mixNum, envir=env)
  invisible(NULL)
}

#' @export
rxUiGet.thetaIniMix <- function(x, ...) {
  .ui <- x[[1]]
  .theta <- .ui$theta
  if (length(.ui$mixProbs) > 0) {
    .theta[.ui$mixProbs] <- rxode2::mlogit(.theta[.ui$mixProbs])
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
