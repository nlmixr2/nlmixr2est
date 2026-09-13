#' Fold saem's temporary etas back into their thetas
#'
#' An eta-less residual theta is fit as `theta + rx.eta.theta`
#' (`.saemAddPseudoEtas()`).  Report it as one parameter: move the mean of
#' the temporary eta into the theta and centre the eta, which leaves every
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
    # a range-transformed theta is estimated on its internal scale
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

#' Remove saem's temporary etas from the finished fit
#'
#' Drops the `rx.eta.*`/`rx.mu.*` columns and blanks the variability cells of
#' the thetas they belonged to; the user's model and omega are restored by
#' `.nlmixrEstUpdatesOrigModel()`.
#'
#' @param ret finished fit
#' @return the fit
#' @noRd
.saemPseudoEtaFinalize <- function(ret) {
  if (!is.environment(ret$env)) return(ret)
  # the ui is already back-transformed here, so read the etas off the fit
  .ranef <- get0("ranef", envir = ret$env, inherits = FALSE)
  .pseudo <- unique(grep("^rx[.]eta[.]", c(names(.ranef), names(ret)), value = TRUE))
  if (length(.pseudo) == 0L) return(ret)
  .thetas <- sub("^rx[.]eta[.]", "", .pseudo)
  .drop <- c(.pseudo, paste0("rx.l.", .thetas), .thetas)
  for (.slot in c("ranef", "etaObf", "shrink")) {
    .d <- get0(.slot, envir = ret$env, inherits = FALSE)
    if (is.data.frame(.d)) {
      assign(.slot, .d[, !(names(.d) %in% .drop), drop = FALSE], envir = ret$env)
    }
  }
  for (.slot in c("parFixedDf", "parFixed")) {
    .d <- get0(.slot, envir = ret$env, inherits = FALSE)
    if (!is.data.frame(.d)) next
    .rows <- rownames(.d) %in% .thetas
    for (.c in names(.d)[startsWith(names(.d), "BSV(") | names(.d) == "Shrink(SD)%"]) {
      .d[.rows, .c] <- if (is.numeric(.d[[.c]])) NA else ""
    }
    assign(.slot, .d, envir = ret$env)
  }
  if (inherits(ret, "data.frame")) {
    .w <- which(names(ret) %in% .drop)
    if (length(.w) > 0L) {
      .cls <- class(ret)
      class(ret) <- "data.frame"
      ret <- ret[, -.w]
      class(ret) <- .cls
    }
  }
  ret
}

postFinalObjectHooksAdd(".saemPseudoEtaFinalize", .saemPseudoEtaFinalize)
