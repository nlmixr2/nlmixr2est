#' Covariate factors of thetas, found in one walk of the model
#'
#' For each of `names`, the factor is the whole product/quotient term holding it
#' as a numerator leaf, with the theta set to 1, kept only when every variable
#' left is a covariate; so `b*WT/70` gives `1*WT/70`, and `WT/b` or `tcl*b*WT`
#' give nothing.  The first (outermost) match wins.
#'
#' @param lst list of model expressions
#' @param names theta names to look for
#' @param covs covariate names
#' @return named list of factor expressions (only the names found)
#' @noRd
.zeroThetaCovFactors <- function(lst, names, covs) {
  .found <- list()
  .isChain <- function(e) {
    is.call(e) && as.character(e[[1]]) %in% c("*", "/") && length(e) == 3L
  }
  .unparen <- function(e) {
    while (is.call(e) && identical(e[[1]], quote(`(`))) {
      e <- e[[2]]
    }
    e
  }
  # leaves of a */ chain, flagged by whether they sit in the numerator
  .leaves <- function(e, inNum) {
    e <- .unparen(e)
    if (.isChain(e)) {
      return(c(.leaves(e[[2]], inNum), .leaves(e[[3]], inNum && identical(e[[1]], quote(`*`)))))
    }
    list(list(e = e, num = inNum))
  }
  .walk <- function(e) {
    if (!is.call(e) || length(.found) == length(names)) {
      return(invisible())
    }
    if (!.isChain(e)) {
      for (.a in as.list(e)[-1]) {
        .walk(.a)
      }
      return(invisible())
    }
    .lv <- .leaves(e, TRUE)
    for (.l in .lv) {
      .n <- if (is.name(.l$e)) as.character(.l$e) else ""
      if (.l$num && .n %in% names && is.null(.found[[.n]])) {
        .f <- do.call(substitute, list(e, stats::setNames(list(1), .n)))
        .v <- all.vars(.f)
        if (length(.v) > 0L && all(.v %in% covs)) {
          .found[[.n]] <<- .f
        }
      }
    }
    # a theta in this chain is consumed here; only look inside its call leaves
    for (.l in .lv) {
      .walk(.l$e)
    }
    invisible()
  }
  for (.e in lst) {
    .walk(.e)
  }
  .found
}

#' 1/sd() of a covariate expression in the data, or 1
#'
#' Evaluated like mu2/mu3/mu4 covariates (`.mu2CovValue()`), so model variables
#' and string assignments resolve as they do for mu referencing.
#' @noRd
.zeroThetaInvSd <- function(expr, datEnv) {
  # mu4 comparisons evaluate to logical
  .x <- .mu2CovValue(datEnv, deparse1(expr))$value
  if (!is.numeric(.x) && !is.logical(.x)) {
    return(1)
  }
  .x <- as.double(.x)
  .sd <- stats::sd(.x[is.finite(.x)])
  if (length(.sd) != 1L || !is.finite(.sd) || .sd <= 0) {
    return(1)
  }
  1 / .sd
}

#' Typical magnitude of thetas that start at 0
#'
#' A covariate coefficient gets `1/sd()` of its covariate in the data, so a unit
#' step moves the linear predictor by one standard deviation; anything else gets
#' 1.  The covariate comes from rxode2's mu-reference tables
#' (`muRefCovariateDataFrame`, `mu2RefCovariateReplaceDataFrame`, which also
#' hold mu3/mu4 terms) when listed there; the rest are found in a single walk of
#' the model.
#'
#' @param ui rxode2 ui
#' @param names theta names
#' @param data estimation data (may be `NULL`)
#' @return named positive typical magnitudes
#' @noRd
.zeroThetaTypical <- function(ui, names, data) {
  .ret <- setNames(rep(1, length(names)), names)
  if (!is.data.frame(data) || length(names) == 0L) {
    return(.ret)
  }
  .cov <- list()
  for (.df in list(ui$muRefCovariateDataFrame, ui$mu2RefCovariateReplaceDataFrame)) {
    for (.n in setdiff(intersect(names, .df$covariateParameter), names(.cov))) {
      .cov[[.n]] <- str2lang(as.character(.df$covariate[match(.n, .df$covariateParameter)]))
    }
  }
  .rest <- setdiff(names, names(.cov))
  if (length(.rest) > 0L && length(ui$allCovs) > 0L) {
    .cov <- c(.cov, .zeroThetaCovFactors(ui$lstExpr, .rest, ui$allCovs))
  }
  if (length(.cov) == 0L) {
    return(.ret)
  }
  # covariate names match data columns without regard to case, as in estimation
  .covs <- ui$allCovs
  .w <- match(tolower(names(data)), tolower(.covs))
  names(data)[!is.na(.w)] <- .covs[.w[!is.na(.w)]]
  .datEnv <- .mu2CovEnv(ui, data)
  for (.n in names(.cov)) {
    .ret[.n] <- .zeroThetaInvSd(.cov[[.n]], .datEnv)
  }
  .ret
}

#' Handle theta initial estimates that are exactly zero
#'
#' A linear `theta` has no magnitude to be scaled by when it starts at `0`, so
#' this hook records a typical magnitude for each estimated, non-fixed,
#' non-residual theta that starts at `0` (`1/sd()` of the covariate it
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
  .scale <- unname(.zeroThetaTypical(.ui, .iniDf$name[.w], data))
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
