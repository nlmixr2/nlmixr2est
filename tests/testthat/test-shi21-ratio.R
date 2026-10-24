test_that("shi21 ratio treatments of censored (r < 1) components (#1188)", {
  .types <- c(legacy = 0L, detected = 1L, substitute = 2L, lmomco = 3L)
  .r <- function(x, type) shi21RatioTest(x, .types[[type]])
  .hm <- function(x) length(x) / sum(1 / x)

  # a scalar ratio is returned unchanged by every treatment
  for (.t in names(.types)) {
    expect_equal(.r(4.6e-12, .t), 4.6e-12, info = .t)
    expect_equal(.r(3, .t), 3, info = .t)
  }

  # all detected: every new treatment equals the legacy harmonic mean
  .x <- c(1e6, 1.7e3, 2.5, 1)
  for (.t in names(.types)) {
    expect_equal(.r(.x, .t), .hm(.x), info = .t)
  }

  # one roundoff-level component (the #1179 column): only the legacy treatment
  # is pinned near zero
  .x <- c(6.7e6, 1.7e3, 4.6e-12)
  expect_lt(.r(.x, "legacy"), 1e-10)
  expect_equal(.r(.x, "detected"), .hm(.x[1:2]))
  expect_equal(.r(.x, "substitute"), .hm(c(.x[1:2], 1)))
  expect_equal(.r(.x, "lmomco"), .hm(.x[1:2]) * 2 / 3)

  # "detected" is invariant to adding components the coordinate does not touch
  for (.k in c(0, 1, 3, 10)) {
    .y <- c(.x, rep(0, .k), rep(1e-13, .k))
    expect_equal(.r(.y, "detected"), .hm(.x[1:2]), info = .k)
  }

  # all censored: "detected" returns the largest ratio, so the step still grows
  .x <- c(0, 0.3, 4.6e-12)
  expect_equal(.r(.x, "detected"), 0.3)
  expect_equal(.r(.x, "substitute"), 1)
  expect_equal(.r(.x, "lmomco"), 0)
})

test_that(".shi21RatioCensor() validates and sets the treatment", {
  .old <- .shi21RatioCensor("detected")
  on.exit(.shi21RatioCensor(.old))
  expect_equal(.shi21RatioCensor("substitute"), "detected")
  expect_equal(.shi21RatioCensor("lmomco"), "substitute")
  expect_error(.shi21RatioCensor("bogus"))
  withr::with_options(list(nlmixr2est.shi21RatioCensor = "legacy"), {
    expect_equal(.shi21RatioCensor(), "lmomco")
  })
  expect_error(shi21RatioCensorSet(4L))
})

test_that("shi21CentralWrap steps a column with an untouched component (#1188)", {
  # f'(t) = exp(4 t) has a large third derivative; the other two components do not
  # depend on t (exact zero and roundoff-level differences).
  .f <- function(t) c(exp(4 * t), 1e-3 * cos(t), 1 + 1e-17 * t)
  .df <- function(t) c(4 * exp(4 * t), -1e-3 * sin(t), 0)
  .t <- 0.3
  .old <- .shi21RatioCensor("legacy")
  on.exit(.shi21RatioCensor(.old))
  .ef <- 1e-12
  .cur <- shi21CentralWrap(.f, .t, .f(.t), 1L, .ef)
  .shi21RatioCensor("detected")
  .det <- shi21CentralWrap(.f, .t, .f(.t), 1L, .ef)
  .err <- function(s) abs(s$gr[1] - .df(.t)[1]) / .df(.t)[1]
  expect_lt(.det$h, .cur$h)
  expect_lt(.err(.det), 1e-6)
  expect_lt(.err(.det), .err(.cur))
})
