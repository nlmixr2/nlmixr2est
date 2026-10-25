#' Find the covariate factor a theta multiplies
#'
#' The factor is the whole product/quotient term holding `name` as a numerator
#' leaf, with `name` set to 1, so `b*WT/70` and `b*(WT/70)` both give `WT/70`.
#'
#' @param expr model expression to search
#' @param name theta name
#' @param covs covariate names
#' @return the covariate-only factor, or `NULL`
#' @noRd
.zeroThetaCovFactor <- function(expr, name, covs) {
  if (!is.call(expr)) {
    return(NULL)
  }
  .op <- as.character(expr[[1]])
  if (.op %in% c("*", "/") && length(expr) == 3L) {
    # numerator leaves of the */ chain (parentheses are transparent)
    .num <- function(e, inNum) {
      while (is.call(e) && identical(e[[1]], quote(`(`))) {
        e <- e[[2]]
      }
      if (is.call(e) && as.character(e[[1]]) %in% c("*", "/") && length(e) == 3L) {
        return(.num(e[[2]], inNum) || .num(e[[3]], inNum && identical(e[[1]], quote(`*`))))
      }
      inNum && identical(e, as.name(name))
    }
    if (.num(expr, TRUE)) {
      .f <- do.call(substitute, list(expr, stats::setNames(list(1), name)))
      .v <- all.vars(.f)
      if (length(.v) > 0L && all(.v %in% covs)) {
        return(.f)
      }
      return(NULL)
    }
  }
  for (.e in as.list(expr)[-1]) {
    .r <- .zeroThetaCovFactor(.e, name, covs)
    if (!is.null(.r)) {
      return(.r)
    }
  }
  NULL
}

#' Typical magnitude of a theta that starts at 0
#'
#' A coefficient of a covariate term (`theta * COV`, `theta * (WT - 70)`) gets
#' `1/sd()` of that term in the data, so a unit step moves the linear predictor by
#' one standard deviation of the covariate; anything else gets 1.
#'
#' @param ui rxode2 ui
#' @param name theta name
#' @param data estimation data (may be `NULL`)
#' @return positive typical magnitude
#' @noRd
.zeroThetaTypical <- function(ui, name, data) {
  if (!is.data.frame(data)) {
    return(1)
  }
  .covs <- ui$allCovs
  if (length(.covs) == 0L) {
    return(1)
  }
  .f <- NULL
  for (.e in ui$lstExpr) {
    .f <- .zeroThetaCovFactor(.e, name, .covs)
    if (!is.null(.f)) break
  }
  if (is.null(.f)) {
    return(1)
  }
  # covariate names match data columns without regard to case, as in estimation
  .v <- all.vars(.f)
  .w <- match(tolower(.v), tolower(names(data)))
  if (anyNA(.w)) {
    return(1)
  }
  .env <- setNames(lapply(.w, function(i) as.double(data[[i]])), .v)
  .x <- tryCatch(eval(.f, .env, baseenv()), error = function(e) NULL)
  if (!is.numeric(.x)) {
    return(1)
  }
  .sd <- stats::sd(.x[is.finite(.x)])
  if (length(.sd) != 1L || !is.finite(.sd) || .sd <= 0) {
    return(1)
  }
  1 / .sd
}

#' Handle theta initial estimates that are exactly zero
#'
#' A linear `theta` has no magnitude to be scaled by when it starts at `0`, so
#' this hook records a typical magnitude for each estimated, non-fixed,
#' non-residual theta that starts at `0` (`1/sd()` of the covariate term it
#' multiplies, else 1).  The scaling setup uses it as that parameter's `scaleC`.
#' For the FOCEi family (a `foceiControl`) the parameter is also moved to
#' `+zeroTheta` when that is within its bounds, otherwise `-zeroTheta`; if neither
#' is within the bounds it errors.  Runs before `.preProcessBoundedTransform` so
#' the moved value is what gets transformed.
#'
#' @param ui rxode2 ui model
#' @inheritParams nlmixr2
#' @return list with the ui (possibly modified)
#' @export
#' @author Matthew L. Fidler
.preProcessZeroTheta <- function(ui, est, data, control) {
  nlmixr2global$nlmixr2EstEnv$zeroThetaScale <- NULL
  .ui <- rxode2::rxUiDecompress(ui)
  .iniDf <- .ui$iniDf
  # Only structural thetas: residual-error parameters (is.na(err) FALSE) get their
  # own scaleC (0.5*|init|, or 1 for boxCox/yeoJohnson), and a 0-valued error sd is
  # a deliberate "disable this component" that must reduce to the smaller model.
  .w <- which(!is.na(.iniDf$ntheta) & is.na(.iniDf$err) & .iniDf$est == 0 & !.iniDf$fix)
  if (length(.w) == 0L) {
    return(NULL)
  }
  .mag <- if (inherits(control, "foceiControl")) control$zeroTheta else NULL
  .nudge <- is.numeric(.mag) && length(.mag) == 1L && is.finite(.mag) && .mag > 0
  .scale <- vapply(.w, function(i) .zeroThetaTypical(.ui, .iniDf$name[i], data), double(1))
  if (.nudge) {
    for (.i in .w) {
      .lower <- .iniDf$lower[.i]
      .upper <- .iniDf$upper[.i]
      if (.mag > .lower && .mag < .upper) {
        .new <- .mag
      } else if (-.mag > .lower && -.mag < .upper) {
        .new <- -.mag
      } else {
        stop(
          "cannot move the zero initial estimate of '",
          .iniDf$name[.i],
          "' off 0: neither ",
          .mag,
          " nor ",
          -.mag,
          " is within its bounds (",
          .lower,
          ", ",
          .upper,
          ")",
          call. = FALSE
        )
      }
      .iniDf$est[.i] <- .new
    }
  }
  # keyed by name AND start value, so a stale record never scales another model
  nlmixr2global$nlmixr2EstEnv$zeroThetaScale <- data.frame(
    name = .iniDf$name[.w],
    est = .iniDf$est[.w],
    scale = .scale,
    stringsAsFactors = FALSE
  )
  if (!.nudge) {
    return(NULL)
  }
  .ui$iniDf <- .iniDf
  .minfo(paste0(
    "moved zero initial estimate(s) off 0 (foceiControl(zeroTheta)): ",
    paste(.iniDf$name[.w], collapse = ", ")
  ))
  list(ui = .ui)
}
preProcessHooksAdd(".preProcessZeroTheta", .preProcessZeroTheta)
