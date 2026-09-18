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
    if (!is.data.frame(.d)) {
      next
    }
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
  if (!inherits(ret, "data.frame") || length(.w) == 0L) {
    return(ret)
  }
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
  if (!inherits(.ui, "rxUi")) {
    return(invisible())
  }
  .ui <- rxode2::rxUiDecompress(.ui)
  .iniDf <- .ui$iniDf
  .used <- unique(unlist(lapply(.ui$lstExpr, all.vars)))
  .rm <- grepl("^rx[.]eta[.]", .iniDf$name) & !(.iniDf$name %in% .used)
  if (!any(.rm)) {
    return(invisible())
  }
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
  if (!is.environment(ret$env)) {
    return(ret)
  }
  .saemPseudoEtaCleanUi(ret$env)
  # the ui is already back-transformed here, so read the etas off the fit
  .ranef <- get0("ranef", envir = ret$env, inherits = FALSE)
  .pseudo <- unique(grep("^rx[.]eta[.]", c(names(.ranef), names(ret)), value = TRUE))
  if (length(.pseudo) == 0L) {
    return(ret)
  }
  .thetas <- sub("^rx[.]eta[.]", "", .pseudo)
  .drop <- c(.pseudo, paste0("rx.l.", .thetas), .thetas)
  .saemPseudoEtaDropSlots(ret$env, .drop)
  .saemPseudoEtaBlankVariability(ret$env, .thetas)
  .saemDropFitColumns(ret, .drop)
}

postFinalObjectHooksAdd(".saemPseudoEtaFinalize", .saemPseudoEtaFinalize)
