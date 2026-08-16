# The nlm population-likelihood C API (#953): value + analytic theta
# gradient over the nlmObjectiveSetup()-loaded problem, exposed for
# downstream samplers (nlmixr2stan tier 0).  Tested through the R shims
# over the SAME C entry points the external-pointer table hands out.

nmTest({

.popMod <- function() {
  ini({
    tcl <- 1
    tv <- 3
    add.sd <- c(0, 0.5)
  })
  model({
    cl <- exp(tcl)
    v <- exp(tv)
    cp <- 100 / v * exp(-cl / v * time)
    cp ~ add(add.sd)
  })
}

.popData <- local({
  set.seed(42)
  .tt <- c(0.5, 1, 2, 4, 8)
  do.call(rbind, lapply(1:4, function(id) {
    data.frame(ID = id, TIME = .tt,
               DV = 5 * exp(-0.05 * .tt) + stats::rnorm(5, 0, 0.5),
               AMT = 0, EVID = 0)
  }))
})

test_that("the table exists and reports the ABI version", {
  .p <- .nlmixr2estNlmPtrs()
  expect_type(.p, "list")
  expect_equal(names(.p), c("apiVersion", "dims", "eval"))
})

test_that("nlm eval: value is -logLik with constants; gradient is analytic (#953)", {
  nlmObjectiveSetup(.popMod, .popData, gradient = TRUE, scale = "natural")
  on.exit(.nlmFreeEnv(), add = TRUE)
  .d <- nlmLikDims_()
  expect_equal(.d$status, 0L)
  expect_equal(.d$ntheta, 3L)
  expect_equal(.d$nobs, 20L)
  # 0x01 gradient model loaded, 0x02 natural scale
  expect_equal(bitwAnd(.d$flags, 0x01), 0x01)
  expect_equal(bitwAnd(.d$flags, 0x02), 0x02)
  .th <- c(1, 3, 0.5)
  .e <- nlmLikEvalC_(.th)
  expect_equal(.e$nBad, 0L)
  # convention: value = -logLik INCLUDING normalization constants, so a
  # hand-written dnorm sum pins it exactly
  .hand <- sum(vapply(1:4, function(i) {
    .di <- .popData[.popData$ID == i, ]
    .f <- 100 / exp(3) * exp(-exp(1) / exp(3) * .di$TIME)
    sum(stats::dnorm(.di$DV, .f, 0.5, log = TRUE))
  }, numeric(1)))
  expect_equal(.e$value, -.hand, tolerance = 1e-8)
  # analytic gradient vs central differences of the value itself
  .h <- 1e-6
  .fd <- vapply(1:3, function(k) {
    .up <- .th
    .up[k] <- .up[k] + .h
    .dn <- .th
    .dn[k] <- .dn[k] - .h
    (nlmLikEvalC_(.up)$value - nlmLikEvalC_(.dn)$value) / (2 * .h)
  }, numeric(1))
  expect_equal(.e$grad, .fd, tolerance = 1e-6)
  # determinism: repeated evaluation is bitwise identical
  expect_identical(nlmLikEvalC_(.th), nlmLikEvalC_(.th))
})

test_that("gradient refused when only the pred model is loaded", {
  nlmObjectiveSetup(.popMod, .popData) # gradient = FALSE (predOnly)
  on.exit(.nlmFreeEnv(), add = TRUE)
  .d <- nlmLikDims_()
  expect_equal(bitwAnd(.d$flags, 0x01), 0L)
  expect_error(nlmLikEvalC_(c(1, 3, 0.5)), "-2")
})

test_that("not-loaded and bad-length report by return code", {
  .nlmFreeEnv()
  expect_equal(nlmLikDims_()$status, -1L)
  expect_error(nlmLikEvalC_(c(1, 3, 0.5)), "-1")
  nlmObjectiveSetup(.popMod, .popData, gradient = TRUE, scale = "natural")
  on.exit(.nlmFreeEnv(), add = TRUE)
  expect_error(nlmLikEvalC_(c(1, 3)), "-3")
})

})
