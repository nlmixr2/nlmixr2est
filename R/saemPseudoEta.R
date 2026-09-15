#' Transform line and internal starting value for one temporary eta
#'
#' @param t theta name
#' @param lo,hi range of the theta
#' @param est its natural-scale starting value
#' @return list with `type`, `init` and `line` (one or two model lines)
#' @noRd
.saemPseudoEtaLine <- function(t, lo, hi, est) {
  .core <- paste0("rxBoundedTr.", t, " + rx.eta.", t)
  if (is.finite(lo) && is.finite(hi)) {
    .eps <- (hi - lo) * 1e-6
    .e <- max(lo + .eps, min(hi - .eps, est))
    return(list(type = "logit", init = log((.e - lo) / (hi - .e)),
                line = paste0(t, " <- expit(", .core, ", ", lo, ", ", hi, ")")))
  }
  if (is.finite(lo)) {
    # rxode2 does not mu-reference `a + exp(theta + eta)`, so an offset takes a helper line
    .line <- if (lo == 0) {
      paste0(t, " <- exp(", .core, ")")
    } else {
      c(paste0("rx.l.", t, " <- exp(", .core, ")"), paste0(t, " <- ", lo, " + rx.l.", t))
    }
    return(list(type = "lower_exp", init = log(max(est - lo, 1e-6)), line = .line))
  }
  if (is.finite(hi)) {
    return(list(type = "upper_exp", init = log(max(hi - est, 1e-6)),
                line = paste0(t, " <- ", hi, " - exp(", .core, ")")))
  }
  list(type = "identity", init = est, line = paste0(t, " <- ", .core))
}

#' Give eta-less thetas a temporary mu-referenced eta on their range's scale
#'
#' Each theta is estimated through an internal `rxBoundedTr.<theta>` and a
#' prepended line on the scale of its range: `exp()` for a positive
#' parameter, `expit()` for a bounded one, additive when unbounded.  The
#' transforms are handed to the bounded-transform back-transform
#' (`.postEstimationBoundedTransform()`).
#'
#' @param ui rxode2 ui
#' @param spec data frame from `.saemPseudoEtaThetas()`
#' @param omega initial variance of each temporary eta
#' @return rewritten ui carrying `boundedTransforms`
#' @noRd
.saemAddPseudoEtas <- function(ui, spec, omega = 0.1) {
  .iniDf <- ui$iniDf
  .etaRows <- which(!is.na(.iniDf$neta1))
  .template <- .iniDf[if (length(.etaRows) > 0L) .etaRows[1] else 1L, , drop = FALSE]
  .maxEta <- max(c(0, .iniDf$neta1[.etaRows]))
  .newLines <- character(0)
  .transforms <- vector("list", nrow(spec))
  for (.k in seq_len(nrow(spec))) {
    .t <- spec$theta[.k]
    .w <- which(.iniDf$name == .t)
    .tr <- .saemPseudoEtaLine(.t, spec$lower[.k], spec$upper[.k], .iniDf$est[.w])
    .newLines <- c(.newLines, .tr$line)
    .transforms[[.k]] <- list(name = .t, internalName = paste0("rxBoundedTr.", .t),
                              type = .tr$type, lower = spec$lower[.k],
                              upper = spec$upper[.k], initTrans = .tr$init,
                              initOrig = .iniDf$est[.w], pseudoEta = TRUE)
    .iniDf[.w, c("name", "lower", "upper", "est", "err", "condition")] <-
      list(paste0("rxBoundedTr.", .t), -Inf, Inf, .tr$init, NA_character_, NA_character_)
    .maxEta <- .maxEta + 1
    .row <- .template
    .row[, c("ntheta", "neta1", "neta2", "name", "lower", "upper", "est", "fix",
             "label", "backTransform", "condition", "err")] <-
      list(NA_integer_, .maxEta, .maxEta, paste0("rx.eta.", .t), -Inf, Inf, omega, FALSE,
           NA_character_, NA_character_, "id", NA_character_)
    if (any(names(.row) == "prior")) .row$prior <- NA_character_
    .iniDf <- rbind(.iniDf, .row)
  }
  .model <- str2lang(paste0("model({",
                            paste(c(.newLines, vapply(ui$lstExpr, deparse1, character(1))),
                                  collapse = "\n"),
                            "})"))
  .ini <- as.expression(lotri::as.lotri(.iniDf))
  .ini[[1]] <- quote(`ini`)
  .fun <- .getUiFunFromIniAndModel(ui, .ini, .model)
  .newUi <- rxode2::rxUiDecompress(.fun())
  assign("modelName", ui$modelName, envir = .newUi)
  .newUi$boundedTransforms <- .transforms
  .newUi
}

#' Put saem's temporary-eta transforms back on its ui
#'
#' A later hook can drop them: an IOV rebuild loses the whole list, and the
#' bounded-transform hook replaces it with the user's own bounded thetas.  Add
#' back whichever specs are missing.
#'
#' @param ui rxode2 ui
#' @param stash the transform specs `.preProcessSaemModeledResid()` added
#' @return the ui carrying every spec
#' @noRd
.saemRestorePseudoTransforms <- function(ui, stash) {
  if (length(stash) == 0L) return(ui)
  .ui <- rxode2::rxUiDecompress(ui)
  .have <- vapply(.ui$boundedTransforms, function(tr) tr$internalName, character(1))
  .missing <- Filter(function(tr) !(tr$internalName %in% .have), stash)
  if (length(.missing) == 0L) return(ui)
  # assign(), not $<-: rxode2 refuses to replace an existing component on a compressed ui
  assign("boundedTransforms", c(.ui$boundedTransforms, .missing), envir = .ui)
  .ui
}

#' Fold saem's temporary etas back into their thetas
#'
#' An eta-less theta is fit as `transform(rxBoundedTr.theta + rx.eta.theta)`
#' (`.saemAddPseudoEtas()`).  Report it as one parameter: move the mean of the
#' temporary eta into the internal theta and centre the eta, which leaves every
#' individual prediction unchanged.
#'
#' @param env saem fit environment after `.getSaemTheta()`/`.getSaemOmega()`
#' @return Nothing, called for side effects
#' @noRd
.saemFoldPseudoEtas <- function(env) {
  .eta <- env$.etaMatBase
  if (is.null(.eta) || is.null(colnames(.eta))) return(invisible())
  .pseudo <- grep("^rx[.]eta[.]", colnames(.eta), value = TRUE)
  for (.e in .pseudo) {
    .t <- sub("^rx[.]eta[.]", "", .e)
    if (paste0("rxBoundedTr.", .t) %in% names(env$fullTheta)) .t <- paste0("rxBoundedTr.", .t)
    if (!(.t %in% names(env$fullTheta))) next
    .m <- mean(.eta[, .e])
    env$fullTheta[[.t]] <- env$fullTheta[[.t]] + .m
    env$.etaMatBase[, .e] <- env$.etaMatBase[, .e] - .m
    if (!is.null(env$.etaMat)) env$.etaMat[, .e] <- env$.etaMat[, .e] - .m
    if (!is.null(env$etaObf)) env$etaObf[[.e]] <- env$etaObf[[.e]] - .m
  }
  invisible()
}

#' Internal theta names that carry a temporary eta
#'
#' @param ui rxode2 ui
#' @return `rxBoundedTr.<name>` for each `rx.eta.<name>` in the model
#' @noRd
.saemPseudoEtaThetaNames <- function(ui) {
  .n <- ui$iniDf$name
  .t <- paste0("rxBoundedTr.", sub("^rx[.]eta[.]", "", grep("^rx[.]eta[.]", .n, value = TRUE)))
  .t[.t %in% .n]
}

#' Kernel phi1 columns (0-based) whose theta carries a temporary eta
#'
#' Uses the phi1 order of `.saemPhi1TargetMap()`; `ui$saemInit` carries no theta names.
#'
#' @param ui rxode2 ui
#' @return integer vector
#' @noRd
.saemPseudoPhi1Ix <- function(ui) {
  .pars <- rxUiGet.saemParamsToEstimateCov(list(ui))
  .phi1 <- .pars[.pars %in% ui$muRefDataFrame$theta]
  which(.phi1 %in% .saemPseudoEtaThetaNames(ui)) - 1L
}
