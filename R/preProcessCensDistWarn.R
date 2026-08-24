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
  .nms <- tolower(names(data))
  .hasCens <- FALSE
  .wCens <- which(.nms == "cens")
  if (length(.wCens) == 1L) {
    .hasCens <- isTRUE(any(data[[.wCens]] != 0, na.rm = TRUE))
  }
  .wLimit <- which(.nms == "limit")
  if (!.hasCens && length(.wLimit) == 1L) {
    .hasCens <- isTRUE(any(is.finite(data[[.wLimit]]), na.rm = TRUE))
  }
  if (!.hasCens) return(list(ui = ui))
  for (.d in .unsafeDist) {
    warning(paste0("censoring ignored for '", .d, "' endpoint(s)"), call. = FALSE)
  }
  list(ui = ui)
}
preProcessHooksAdd(".preProcessCensDistWarn", .preProcessCensDistWarn)
