## The lik-contrib registry also fires from the population (nlm-family) C++
## objective in nlm.cpp (nlmSolveFid), so a plugin gets the EXACT per-obs
## error-model cotangent d(LL)/d(f) from the nlm solve -- not just from the
## focei-family inner engines.  nlmObjectiveSetup() loads the nlm predOnly
## problem once so repeated nlmSolveR() calls re-evaluate the objective (firing
## the hook) without re-optimizing -- the entry point nlmixr2nn's population
## weight fit uses to consume the C++ cotangent.

test_that("nlm population objective fires the lik-contrib hook + one-shot setup reuse", {
  skip_on_cran()

  pop.model <- function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka); cl <- exp(tcl); v <- exp(tv)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }

  .old <- rxode2::getRxThreads(); on.exit(rxode2::setRxThreads(.old), add = TRUE)
  rxode2::setRxThreads(1L)   # test contributor uses global accumulators

  .nObs <- sum(theo_sd$EVID == 0)
  .nsub <- length(unique(theo_sd$ID))

  ## a normal nlm fit fires the hook (the C++ cotangent capture path)
  .Call("_nlmixr2est_registerTestContrib", PACKAGE = "nlmixr2est")
  .fit <- suppressWarnings(suppressMessages(
    nlmixr2(pop.model, theo_sd, est = "nlm", control = nlmControl(print = 0L))))
  .resFit <- .Call("_nlmixr2est_getTestContrib", PACKAGE = "nlmixr2est")
  .Call("_nlmixr2est_removeTestContrib", PACKAGE = "nlmixr2est")
  names(.resFit) <- c("nObs", "sumDLLdf", "sumErr", "sumF", "nBegin", "nEnd")
  expect_gt(.resFit[["nObs"]], 0)                    # hook fired during the nlm fit
  expect_equal(.resFit[["nBegin"]], .resFit[["nEnd"]])
  expect_equal(.resFit[["nObs"]] / .resFit[["nBegin"]], .nObs / .nsub)  # obs/subject

  ## one-shot setup: repeated nlmSolveR reuse the single loaded problem
  .ui <- rxode2::rxUiDecompress(.fit$ui)
  .Call("_nlmixr2est_registerTestContrib", PACKAGE = "nlmixr2est")
  .parIni <- nlmObjectiveSetup(.ui, theo_sd)
  on.exit({.nlmFreeEnv(); .Call("_nlmixr2est_removeTestContrib", PACKAGE = "nlmixr2est")},
          add = TRUE)

  .o1 <- nlmSolveR(.parIni)
  .o2 <- nlmSolveR(.parIni)
  expect_true(is.finite(.o1))
  expect_equal(.o1, .o2)                             # deterministic reuse
  .res <- .Call("_nlmixr2est_getTestContrib", PACKAGE = "nlmixr2est")
  names(.res) <- c("nObs", "sumDLLdf", "sumErr", "sumF", "nBegin", "nEnd")
  expect_equal(.res[["nObs"]], 2L * .nObs)           # two evals fired the hook
  expect_equal(.res[["nObs"]] / .res[["nBegin"]], .nObs / .nsub)

  ## an extra constant LL per obs lowers the -LL objective by exactly c*nObs
  ## (nlmSolveFid folds _llAdd into ret(k)); reset counters via re-register
  .Call("_nlmixr2est_removeTestContrib", PACKAGE = "nlmixr2est")
  .Call("_nlmixr2est_registerTestContrib", PACKAGE = "nlmixr2est")
  .base <- nlmSolveR(.parIni)
  .cc <- 0.01
  .Call("_nlmixr2est_setTestContribAddLL", .cc, PACKAGE = "nlmixr2est")
  .shift <- nlmSolveR(.parIni)
  .Call("_nlmixr2est_setTestContribAddLL", 0.0, PACKAGE = "nlmixr2est")
  expect_equal(.shift - .base, -.cc * .nObs, tolerance = 1e-6)
})
