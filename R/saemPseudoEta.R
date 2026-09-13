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

#' Drop temporary-eta columns from the eta-shaped slots of a fit
#'
#' @param env fit environment
#' @param drop column names to drop
#' @return Nothing, called for side effects
#' @noRd
.saemPseudoEtaDropSlots <- function(env, drop) {
  for (.slot in c("ranef", "etaObf", "shrink")) {
    .d <- get0(.slot, envir = env, inherits = FALSE)
    if (is.data.frame(.d)) {
      assign(.slot, .d[, !(names(.d) %in% drop), drop = FALSE], envir = env)
    }
  }
}

#' Blank the variability cells of thetas that only had a temporary eta
#'
#' @param env fit environment
#' @param thetas theta names
#' @return Nothing, called for side effects
#' @noRd
.saemPseudoEtaBlankVariability <- function(env, thetas) {
  for (.slot in c("parFixedDf", "parFixed")) {
    .d <- get0(.slot, envir = env, inherits = FALSE)
    if (!is.data.frame(.d)) next
    .rows <- rownames(.d) %in% thetas
    .cols <- names(.d)[startsWith(names(.d), "BSV(") | names(.d) == "Shrink(SD)%"]
    .d[.rows, .cols] <- lapply(.d[.cols], function(x) if (is.numeric(x)) NA else "")
    assign(.slot, .d, envir = env)
  }
}

#' Drop columns from a fit data frame, keeping its class
#'
#' @param ret fit
#' @param drop column names to drop
#' @return the fit
#' @noRd
.saemDropFitColumns <- function(ret, drop) {
  .w <- which(names(ret) %in% drop)
  if (!inherits(ret, "data.frame") || length(.w) == 0L) return(ret)
  .cls <- class(ret)
  class(ret) <- "data.frame"
  ret <- ret[, -.w]
  class(ret) <- .cls
  ret
}

#' Drop temporary-eta rows the fit's model no longer uses
#'
#' The back-transform removes the temporary eta's model line; an `ini()` row
#' left behind makes the ui unparseable for a later rebuild (`.uiFinalizeIov()`).
#'
#' @param env fit environment
#' @return Nothing, called for side effects
#' @noRd
.saemPseudoEtaCleanUi <- function(env) {
  .ui <- get0("ui", envir = env, inherits = FALSE)
  if (!inherits(.ui, "rxUi")) return(invisible())
  .ui <- rxode2::rxUiDecompress(.ui)
  .iniDf <- .ui$iniDf
  .used <- unique(unlist(lapply(.ui$lstExpr, all.vars)))
  .rm <- grepl("^rx[.]eta[.]", .iniDf$name) & !(.iniDf$name %in% .used)
  if (!any(.rm)) return(invisible())
  .iniDf <- .iniDf[!.rm, , drop = FALSE]
  .e <- !is.na(.iniDf$neta1)
  .lev <- sort(unique(c(.iniDf$neta1[.e], .iniDf$neta2[.e])))
  .iniDf$neta1[.e] <- match(.iniDf$neta1[.e], .lev)
  .iniDf$neta2[.e] <- match(.iniDf$neta2[.e], .lev)
  .ini <- as.expression(lotri::as.lotri(.iniDf))
  .ini[[1]] <- quote(`ini`)
  .fun <- .getUiFunFromIniAndModel(.ui, .ini, rxode2::as.model(.ui$lstExpr))
  .new <- rxode2::rxUiDecompress(.fun())
  assign("modelName", .ui$modelName, envir = .new)
  assign("ui", .new, envir = env)
  invisible()
}

#' Remove saem's temporary etas from the finished fit
#'
#' The user's model and omega are restored by `.nlmixrEstUpdatesOrigModel()`.
#'
#' @param ret finished fit
#' @return the fit
#' @noRd
.saemPseudoEtaFinalize <- function(ret) {
  if (!is.environment(ret$env)) return(ret)
  .saemPseudoEtaCleanUi(ret$env)
  # the ui is already back-transformed here, so read the etas off the fit
  .ranef <- get0("ranef", envir = ret$env, inherits = FALSE)
  .pseudo <- unique(grep("^rx[.]eta[.]", c(names(.ranef), names(ret)), value = TRUE))
  if (length(.pseudo) == 0L) return(ret)
  .thetas <- sub("^rx[.]eta[.]", "", .pseudo)
  .drop <- c(.pseudo, paste0("rx.l.", .thetas), .thetas)
  .saemPseudoEtaDropSlots(ret$env, .drop)
  .saemPseudoEtaBlankVariability(ret$env, .thetas)
  .saemDropFitColumns(ret, .drop)
}

postFinalObjectHooksAdd(".saemPseudoEtaFinalize", .saemPseudoEtaFinalize)
