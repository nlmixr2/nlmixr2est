# Issue 1114: focei with scaleType = "norm" and outerOpt = "nlminb" returned an
# objective of ~2e244 with every parameter left at its initial estimate and no
# diagnostic when a between-study variance started small (2.5e-6).
#
# Mechanism: "norm" scaling maps every parameter through one affine transform
# whose constant is the range of the internal parameter vector, and the omega
# enters that vector as omega^(-1/4) (sqrt of chol(omega^-1)) = 25 here, so a
# unit step in the scaled space moved `slope` from 0.0035 to 12.6.  The ~1e283
# objective at that trial point left the warm-started etas where the inner
# optimizer of 7.0.3 could not recover, so every later evaluation -- including
# the final re-evaluation at the untouched initial theta -- stayed near 1e250,
# nlminb reported "false convergence (8)", and the fit reported that value.
#
# 144 study-level rows from 12 studies (values perturbed).  DV is a percent
# change, COVARIATE another percent change, SE the row's standard error.
.d1114 <- structure(list(ID = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L,
2L, 2L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L,
3L, 3L, 4L, 4L, 4L, 4L, 4L, 4L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L,
5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 6L, 6L, 6L, 6L, 6L, 6L,
6L, 6L, 6L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 8L, 8L, 8L, 8L, 8L,
8L, 8L, 8L, 8L, 9L, 9L, 9L, 9L, 9L, 9L, 9L, 9L, 10L, 10L, 10L,
10L, 10L, 10L, 10L, 10L, 11L, 11L, 12L, 12L, 12L, 12L, 12L, 12L,
12L, 12L, 12L, 12L, 12L, 12L, 12L, 12L, 12L, 12L, 12L, 12L, 12L,
12L, 12L, 12L, 12L, 12L, 12L, 12L, 12L, 12L, 12L, 12L, 12L, 12L
), TIME = c(4, 4, 4, 8, 8, 8, 12, 12, 12, 16, 16, 16, 20, 20,
20, 24, 24, 24, 12, 12, 12, 12, 24, 24, 24, 24, 24, 36, 36, 36,
36, 36, 4, 4, 8, 8, 12, 12, 24, 24, 36, 36, 48, 48, 12, 12, 24,
24, 36, 36, 4, 4, 4, 12, 12, 12, 28, 28, 28, 52, 52, 52, 76,
76, 76, 100, 100, 100, 4, 4, 4, 12, 12, 12, 24, 24, 24, 12, 12,
24, 24, 36, 36, 48, 48, 24, 24, 24, 48, 48, 48, 72, 72, 72, 4,
4, 13, 13, 26, 26, 39, 39, 12, 24, 36, 48, 60, 72, 84, 96, 1.14285714285714,
1.14285714285714, 12, 12, 12, 12, 24, 24, 24, 24, 36, 36, 36,
36, 48, 48, 48, 48, 60, 60, 60, 60, 72, 72, 72, 72, 84, 84, 84,
84, 96, 96, 96, 96), DV = c(2.69, -4.23, -20.66, -7.84, -6.31,
-21.58, -1.15, -11.54, -33.56, 3.35, -32.31, -35.39, 13.15, -33.84,
-50.11, 1.48, -22.2, -54.08, -29.61, -16.22, -8.93, -24.73, -35.64,
-40.11, -26.34, -6.52, -28.78, -30.93, -41.72, -34.01, 2.87,
-32.98, 2.3, -7.05, 4.42, -15.27, 7.74, -21.89, 12.48, -31.43,
-2.32, -46.64, 8.91, -57.76, -28.3, -8.17, -39.14, -9.85, -46.78,
-7.26, -6.39, -7.11, -6.06, -31.66, -17.41, -19.12, -47.1, -32.03,
-29.81, -70.8, -48.7, -56.1, -68.73, -55.46, -54.8, -69.17, -65.66,
-60.41, -1.21, -9.51, -28.17, -39.81, -29.57, -44.35, -9.7, -30.4,
-60.62, -26.73, -35.29, -46.71, -53.95, -54.79, -53.81, -59.48,
-64.44, -26.61, -24.07, 26.67, -37.81, 6.65, -38.73, -48.3, -2.75,
27.1, 2.17, -16.63, -6.63, -39.02, -6.77, -53.66, -8.04, -49.3,
-25.67, -34.92, -37.69, -40.68, -43.45, -51.89, -52.76, -57.45,
-11.42, -6.19, -4.62, -21.78, -29.31, -14.15, -21.61, -30.44,
-37.24, 8.54, -12.93, -17.11, -22.74, -4.94, 10.03, -30.51, -24.91,
16.42, -14.6, -12.12, -46.77, 5, -16.44, -25.17, -31.04, 10.64,
-13.66, -30.24, -41.54, -3, -3.72, -16.24, -44.2, 4.69), COVARIATE = c(-13.08,
-17.74, -32.63, -10.91, -32.58, -37.51, -4.09, -36.62, -48.59,
-5.52, -36.05, -52.84, -6.69, -41.57, -52.4, 9.75, -44.9, -51.92,
-57.48, -39.38, -4.41, -57.22, -56.14, -29.39, -50.02, -3.41,
-62.08, -61.15, -31.87, -61.8, -4.23, -64.83, 1.66, -41.31, 1.32,
-54.28, -0.16, -56.22, 1.67, -64.18, 2.71, -64.4, 4.31, -75.1,
-54.81, 0.78, -63.04, 1.58, -62.54, 0.39, -41.67, -37.8, -36.08,
-56.28, -50.12, -55.48, -61.8, -62.05, -63.38, -68.23, -65.66,
-66.87, -72.68, -65.47, -67.18, -82.95, -66.4, -69.63, -0.28,
-29.09, -23.67, 3.52, -39.78, -41.12, 1.22, -50.33, -50.63, -62.37,
-57.14, -68.37, -59.62, -73.3, -59.65, -74.51, -71.61, -13.14,
-49.39, 4.08, -21.75, -51.61, -3.38, -22.13, -53.6, 9.63, -0.04,
-32.52, -0.38, -50.62, -1.33, -56.7, -2.23, -56.84, -47.42, -51.75,
-56.24, -62.96, -63.48, -66.9, -70.43, -68.04, -9.82, -12.78,
-32.91, -31.1, -23.57, -2.47, -28.52, -27.97, -27.56, -5.87,
-14.99, -19.1, -23.69, -1.08, -15.93, -19.51, -23.34, -1.44,
-5.4, -17.82, -22.21, 0.79, -10.37, -15.04, -23.63, 4.51, -8.6,
-6.83, -18.79, -7.04, -8.37, -13.85, -18.74, -0.13), SE = c(12.227,
11.543, 9.811, 13.08, 13.296, 11.938, 14.698, 12.271, 9.751,
15.391, 10.629, 11.23, 17.876, 10.608, 9.228, 19.037, 12.043,
10.625, 7.191, 8.092, 8.786, 6.154, 8.2, 11.775, 8.409, 9.824,
5.915, 8.393, 12.837, 7.813, 12.421, 6.081, 4.28, 3.874, 4.722,
3.705, 5.26, 4.223, 6.141, 4.4, 5.167, 3.395, 6.268, 3.213, 3.881,
5.637, 4.141, 6.154, 4.099, 5.509, 15.028, 8.664, 6.672, 14.273,
9.091, 7.966, 13.458, 10.069, 7.293, 12.316, 8.713, 7.01, 13.253,
9.684, 6.913, 16.164, 8.31, 7.333, 15.299, 15.502, 10.298, 13.302,
15.818, 11.373, 20.01, 17.779, 11.055, 7.128, 8.302, 7.554, 8.318,
6.713, 7.72, 7.706, 8.515, 26.434, 25.735, 34.474, 36.373, 37.021,
21.37, 27.585, 44.368, 52.38, 3.803, 2.834, 4.39, 2.815, 5.126,
2.767, 6.005, 2.932, 4.49, 4.557, 4.495, 4.738, 4.693, 5.441,
5.202, 5.237, 5.065, 6.051, 15.997, 14.633, 11.616, 14.36, 14.887,
14.156, 12.736, 20.446, 17.466, 17.325, 15.372, 19.276, 23.241,
16.9, 14.727, 23.315, 17.994, 18.939, 13.651, 21.896, 17.586,
17.928, 14.773, 23.628, 19.091, 16.725, 14.385, 18.757, 20.29,
19.702, 13.629, 20.998)), row.names = c(NA, -144L), class = "data.frame")

.m1114 <- function() {
  ini({
    b0 <- fix(0)
    slope <- 0.0035
    emax_scale <- 5
    let50 <- log(47)
    add.sd <- 0.9
    eta.slope ~ 2.5e-6
  })
  model({
    slope_study <- slope + eta.slope
    et50 <- exp(let50)
    lratio <- b0 + slope_study * (1 + emax_scale * TIME / (et50 + TIME)) * COVARIATE
    ratio <- exp(lratio)
    pct <- 100 * (ratio - 1)
    w <- add.sd * SE
    pct ~ add(w)
  })
}

# Fit and return list(fit=, w=) where `w` holds every note the run raised: the
# warnings nlmixr2Est() collected onto the fit's $runInfo (a warning() raised
# during a fit is routed there rather than re-emitted) plus any warning that
# still reached R.  An environment accumulator rather than `<<-` so the target
# is named at the assignment site.
.fitCollectWarnings <- function(...) {
  .acc <- new.env(parent = emptyenv())
  .acc$w <- character()
  .fit <- withCallingHandlers(
    suppressMessages(nlmixr(...)),
    warning = function(w) {
      .acc$w <- c(.acc$w, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  .runInfo <- tryCatch(.fit$runInfo, error = function(e) NULL)
  list(fit = .fit, w = c(as.character(unlist(.runInfo)), .acc$w))
}

.worseRegex <- "final objective function .* is worse than the initial objective function"

nmTest({
  test_that("issue 1114: scaleType='norm' with a 2.5e-6 starting omega fits", {
    .r <- .fitCollectWarnings(
      .m1114, .d1114, est = "focei",
      control = foceiControl(print = 0, scaleType = "norm", outerOpt = "nlminb")
    )
    fitNorm <- .r$fit
    # The fit is reported honestly: no worse-than-initial warning ...
    expect_false(any(grepl(.worseRegex, .r$w)))
    # ... and the objective improved on the initial 824.37
    expect_true(is.finite(fitNorm$objf))
    expect_lt(fitNorm$objf, 824.37)
    expect_lt(abs(fitNorm$objf - 813.66), 1)
    # Every estimated parameter moved off its initial estimate
    .theta <- fitNorm$theta
    expect_false(isTRUE(all.equal(unname(.theta[["slope"]]), 0.0035)))
    expect_false(isTRUE(all.equal(unname(.theta[["emax_scale"]]), 5)))
    expect_false(isTRUE(all.equal(unname(.theta[["let50"]]), log(47))))
    expect_false(isTRUE(all.equal(unname(.theta[["add.sd"]]), 0.9)))
    expect_equal(unname(.theta[["let50"]]), 3.9245, tolerance = 0.02)
    expect_equal(unname(.theta[["add.sd"]]), 0.9366, tolerance = 0.02)
    expect_false(grepl("false convergence", fitNorm$message))
  })

  test_that("issue 1114: the default scaling path is unchanged", {
    .r <- .fitCollectWarnings(
      .m1114, .d1114, est = "focei",
      control = foceiControl(print = 0, outerOpt = "nlminb")
    )
    fitDefault <- .r$fit
    expect_false(any(grepl(.worseRegex, .r$w)))
    expect_equal(fitDefault$objf, 822.9683, tolerance = 1e-4)
    .theta <- fitDefault$theta
    expect_equal(unname(.theta[["slope"]]), 0.003713, tolerance = 0.01)
    expect_equal(unname(.theta[["emax_scale"]]), 4.995, tolerance = 0.01)
    expect_equal(unname(.theta[["let50"]]), 4.0188, tolerance = 0.01)
    expect_equal(unname(.theta[["add.sd"]]), 0.9527, tolerance = 0.01)
  })

  test_that("a final objective worse than the initial one warns", {
    # A custom outer optimizer that returns a point it knows is worse: it
    # triples `slope` (par[1] on the default nlmixr2 scaling is slope's scaled
    # value 1; +2 is two scaleC = |init| steps).  The objective values it saw are
    # kept so the test asserts the premise (f1 > f0) and not just the warning.
    .seen <- new.env(parent = emptyenv())
    .badOuter <- function(par, fn, gr, lower = -Inf, upper = Inf, control = list(), ...) {
      .x <- par
      .x[1] <- par[1] + 2
      .seen$f0 <- fn(par)
      .seen$f1 <- fn(.x)
      list(x = .x, par = .x, objective = .seen$f1, convergence = 1L,
           message = "deliberately worse point")
    }
    .r <- .fitCollectWarnings(
      .m1114, .d1114, est = "focei",
      control = foceiControl(print = 0, outerOpt = .badOuter, covMethod = "",
                             calcTables = FALSE)
    )
    expect_gt(.seen$f1, .seen$f0)
    expect_equal(sum(grepl(.worseRegex, .r$w)), 1L)
    expect_true(is.finite(.r$fit$objf))
    expect_gt(.r$fit$objf, 824.37)
    # The warning carries both objective values
    .w <- .r$w[grepl(.worseRegex, .r$w)]
    expect_match(.w, "\\(824\\.", fixed = FALSE)

    # Negative control: the same custom optimizer returning its starting point
    # is a legitimate (if useless) run and must not warn.
    .sameOuter <- function(par, fn, gr, lower = -Inf, upper = Inf, control = list(), ...) {
      list(x = par, par = par, objective = fn(par), convergence = 0L,
           message = "returned the starting point")
    }
    .r0 <- .fitCollectWarnings(
      .m1114, .d1114, est = "focei",
      control = foceiControl(print = 0, outerOpt = .sameOuter, covMethod = "",
                             calcTables = FALSE)
    )
    expect_false(any(grepl(.worseRegex, .r0$w)))
    expect_equal(.r0$fit$objf, 824.37, tolerance = 1e-3)
  })

  test_that(".foceiFinalOfvWorse: band of 1% with a 0.1 floor, non-finite final", {
    # inside the band: not worse
    expect_false(.foceiFinalOfvWorse(824, 824))
    expect_false(.foceiFinalOfvWorse(824, 800))
    expect_false(.foceiFinalOfvWorse(824, 824 + 8.2))
    # just outside the 1% band
    expect_true(.foceiFinalOfvWorse(824, 824 + 8.3))
    expect_true(.foceiFinalOfvWorse(824, 2.1e244))
    # negative objectives: band is 1% of the magnitude
    expect_false(.foceiFinalOfvWorse(-500, -495.1))
    expect_true(.foceiFinalOfvWorse(-500, -494.9))
    # near zero the absolute floor of 0.1 applies
    expect_false(.foceiFinalOfvWorse(0, 0.09))
    expect_true(.foceiFinalOfvWorse(0, 0.11))
    expect_false(.foceiFinalOfvWorse(1, 1.09))
    expect_true(.foceiFinalOfvWorse(1, 1.11))
    # a non-finite final objective is always worse
    expect_true(.foceiFinalOfvWorse(824, Inf))
    expect_true(.foceiFinalOfvWorse(824, NaN))
    expect_true(.foceiFinalOfvWorse(824, NA_real_))
    # a non-finite initial objective cannot occur (the initial evaluation stops)
    # and is not reported as worse
    expect_false(.foceiFinalOfvWorse(NA_real_, 824))
    expect_false(.foceiFinalOfvWorse(Inf, 824))
  })
})
