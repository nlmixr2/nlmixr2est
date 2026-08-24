#' Warn when CENS/LIMIT data meets a distribution with no M2/M3/M4 support
#'
#' `censEst.h`'s M2/M3/M4 correction only understands normal (`add`/`prop`/
#' `pow`/combined, `lnorm`/`logitNorm`/`probitNorm` transforms, `dnorm`) and
#' `t()`/`cauchy()` (`doCensT1()`, #979/#992).  Every other generalized-
#' likelihood distribution (`pois`, `binom`, `beta`, `chisq`, `dexp`, `f`,
#' `geom`, `unif`, `weibull`, `dgamma`, `ordinal`, a general `ll()`) silently
#' scores a censored row with its ordinary (uncensored) density instead of
#' erroring or refusing (unlike `est="nls"`, which hard-stops).  This warns
#' instead of leaving that silent, so it is at least visible in the fit's
#' `$runInfo`.
#'
#' @param ui rxode2 ui
#' @inheritParams nlmixr2
#' @return list with the ui (unmodified; this hook only warns)
#' @noRd
NULL

#' Does this dataset carry any censoring information?
#'
#' A non-zero `CENS` (M3/M4) or a finite `LIMIT` (M2) on any row.  Shared with
#' `.foceiFamilyControl()`, which downgrades `fast = TRUE` for a censored
#' log-likelihood endpoint (#992).
#'
#' @inheritParams nlmixr2
#' @return logical
#' @noRd
.nlmixrDataHasCens <- function(data) {
  if (!is.data.frame(data)) return(FALSE)
  .nms <- tolower(names(data))
  .w <- which(.nms == "cens")
  if (length(.w) == 1L && isTRUE(any(data[[.w]] != 0, na.rm = TRUE))) return(TRUE)
  .w <- which(.nms == "limit")
  length(.w) == 1L && isTRUE(any(is.finite(data[[.w]]), na.rm = TRUE))
}

.preProcessCensDistWarn <- function(ui, est, data, control) {
  if (identical(est, "nls")) return(list(ui = ui))
  # t()/cauchy() M2/M3/M4 is wired for two kernels: src/nlm.cpp's
  # population-only one (#979) and src/inner.cpp's likInner0 (#992, which also
  # supplies the eta gradient).  Methods built on some OTHER kernel -- SAEM,
  # the importance-sampling family (imp/impmap/qrpem), nlme, the nonparametric
  # engines -- keep t()/cauchy() in the unsafe set.
  .nlmFamily <- c("nlm", "bobyqa", "lbfgsb3c", "n1qn1", "newuoa", "nlminb",
                  "optim", "uobyqa")
  # `fo`/`foi` are deliberately absent: they hard-stop on censoring rather
  # than ignoring it.
  .foceiFamily <- c("focei", "foce", "focep", "laplace", "agq", "posthoc",
                    "ifocei", "ifoce", "ifocep", "ilaplace", "iagq",
                    "mfocei", "mfoce", "mfocep", "mlaplace", "magq")
  .safe <- if (isTRUE(est %in% c(.nlmFamily, .foceiFamily))) {
    c("norm", "dnorm", "t", "cauchy")
  } else {
    c("norm", "dnorm")
  }
  .predDf <- tryCatch(ui$predDf, error = function(e) NULL)
  if (is.null(.predDf) || nrow(.predDf) == 0) return(list(ui = ui))
  .unsafeDist <- unique(.predDf$distribution[!(.predDf$distribution %in% .safe)])
  if (length(.unsafeDist) == 0) return(list(ui = ui))
  if (!.nlmixrDataHasCens(data)) return(list(ui = ui))
  for (.d in .unsafeDist) {
    warning(paste0("censoring ignored for '", .d, "' endpoint(s)"), call. = FALSE)
  }
  list(ui = ui)
}
preProcessHooksAdd(".preProcessCensDistWarn", .preProcessCensDistWarn)
