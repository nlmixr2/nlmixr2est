# A theta initialized at exactly 0 has no native magnitude, and FOCEi scales a
# linear parameter by 1/|init|.  .preProcessZeroTheta() therefore moves such a
# theta to +/-zeroTheta before estimation -- which makes the nudge do double
# duty, because it becomes the magnitude the outer search assumes.  A nudge of
# 0.001 tells the search the parameter is of size 0.001, it explores in steps
# that size, and it stops on the nudge: the reported "estimate" is the nudge.
#
# No single nudge fixes this, which is why this retry exists rather than a
# bigger default.  Measured, true value in brackets:
#
#   model                          nudge 0.001    nudge 0.1     truth
#   theo_sd, coefficient on        0.0009         0.9192        0.90
#     log(WT/70)
#   test-focei-zero-init-scale,    0.0311         0.0043        0.03
#     coefficient on WT (50-110)
#
# The two want scales 30x apart and a zero init carries no information about
# which.  What DOES separate them is the failure signature: a stalled fit comes
# back AT its nudge (0.0009 vs 0.001), while a successful one moves far off it
# (0.0311 is 31x its nudge).  So detect that, re-fit once from a larger nudge,
# and keep whichever objective function value is better -- the comparison is
# what makes this safe, since a retry that does not help is discarded.

#' Is a zero-nudged theta stalled on its nudge?
#'
#' @param fit focei fit
#' @param mag the `zeroTheta` magnitude used
#' @param tol multiple of `mag` within which an estimate counts as stalled
#' @return names of the stalled thetas (possibly empty)
#' @noRd
.foceiZeroThetaStalled <- function(fit, mag, tol) {
  if (!is.numeric(mag) || length(mag) != 1L || !is.finite(mag) || mag <= 0) {
    return(character(0))
  }
  .pf <- try(fit$parFixedDf, silent = TRUE)
  if (inherits(.pf, "try-error") || is.null(.pf) || nrow(.pf) == 0L) {
    return(character(0))
  }
  .ini <- try(rxode2::rxUiDecompress(fit$ui)$iniDf, silent = TRUE)
  if (inherits(.ini, "try-error") || is.null(.ini)) return(character(0))
  # the same set .preProcessZeroTheta() nudges: estimated structural thetas.
  # Residual-error parameters are excluded there (they carry their own scaleC)
  # and so are fixed ones.
  .cand <- .ini$name[!is.na(.ini$ntheta) & is.na(.ini$err) & !.ini$fix]
  .est <- setNames(.pf$Estimate, rownames(.pf))
  .cand <- intersect(.cand, names(.est))
  if (length(.cand) == 0L) return(character(0))
  .cand[is.finite(.est[.cand]) & abs(.est[.cand]) <= tol * mag]
}

#' Re-fit once from a larger zero-theta nudge when the first fit stalled on it
#'
#' Keeps whichever fit has the better objective function value, so a retry that
#' does not help changes nothing but time.
#'
#' @param env estimation environment (carries `data`)
#' @param fit the fit just produced
#' @param control the focei control used
#' @return `fit`, or the retry when it is better
#' @noRd
.foceiZeroThetaRetry <- function(env, fit, control) {
  .mag <- control$zeroTheta
  .fac <- control$zeroThetaRetry
  .tol <- control$zeroThetaRetryTol
  if (is.null(.fac) || !is.numeric(.fac) || length(.fac) != 1L ||
        !is.finite(.fac) || .fac <= 1) {
    return(fit)                         # retry disabled
  }
  if (is.null(.tol) || !is.numeric(.tol) || length(.tol) != 1L ||
        !is.finite(.tol) || .tol <= 0) {
    return(fit)
  }
  .stalled <- .foceiZeroThetaStalled(fit, .mag, .tol)
  if (length(.stalled) == 0L) return(fit)
  .obj0 <- try(fit$objf, silent = TRUE)
  if (inherits(.obj0, "try-error") || !is.numeric(.obj0) || !is.finite(.obj0)) {
    return(fit)
  }
  .data <- env$data
  if (is.null(.data)) return(fit)
  .ui <- try(rxode2::rxUiDecompress(fit$ui), silent = TRUE)
  if (inherits(.ui, "try-error")) return(fit)
  .iniDf <- .ui$iniDf
  .new <- .mag * .fac
  .moved <- character(0)
  for (.nm in .stalled) {
    .i <- which(.iniDf$name == .nm)
    if (length(.i) != 1L) next
    # keep the sign the first fit chose, and stay inside the bounds
    .sgn <- if (is.finite(.iniDf$est[.i]) && .iniDf$est[.i] < 0) -1 else 1
    .try <- .sgn * .new
    if (!(.try > .iniDf$lower[.i] && .try < .iniDf$upper[.i])) {
      .try <- -.try
      if (!(.try > .iniDf$lower[.i] && .try < .iniDf$upper[.i])) next
    }
    .iniDf$est[.i] <- .try
    .moved <- c(.moved, .nm)
  }
  if (length(.moved) == 0L) return(fit)
  .ui$iniDf <- .iniDf
  # zeroTheta is irrelevant on the retry (nothing is at 0 any more) and the
  # retry must not recurse
  .ctl <- control
  .ctl$zeroThetaRetry <- 1
  .fit2 <- try(suppressWarnings(nlmixr2(.ui, .data, est = "focei",
                                        control = .ctl)),
               silent = TRUE)
  if (inherits(.fit2, "try-error")) return(fit)
  .obj1 <- try(.fit2$objf, silent = TRUE)
  if (inherits(.obj1, "try-error") || !is.numeric(.obj1) || !is.finite(.obj1)) {
    return(fit)
  }
  if (.obj1 >= .obj0) {
    # the larger nudge did not help; say so, because the first fit's estimates
    # for these parameters are still sitting on the nudge
    warning("theta(s) ", paste(.moved, collapse = ", "),
            " came back at the foceiControl(zeroTheta=) nudge (", .mag,
            "); a re-fit from ", .new, " did not improve the objective (",
            signif(.obj1, 6), " vs ", signif(.obj0, 6),
            "), so the original fit is kept.  Treat those estimates as ",
            "un-estimated and set a non-zero initial value if that is wrong.",
            call. = FALSE)
    return(fit)
  }
  .minfo(paste0("theta(s) ", paste(.moved, collapse = ", "),
                " stalled on the zeroTheta nudge (", .mag,
                "); re-fit from ", .new, " improved the objective (",
                signif(.obj1, 6), " vs ", signif(.obj0, 6), ")"))
  .fit2
}
