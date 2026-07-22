## Slow (weekly batch): full grad-vs-regress fits on theo_sd.
##
## nonMuTheta="regress" minimizes the JOINT likelihood at frozen encoder etas,
## whose optimum is displaced from the marginal one (the classic joint-vs-marginal
## inconsistency).  nonMuTheta="grad" differentiates the marginal (Laplace)
## objective instead, so it should land materially closer to the FOCEi MLE for the
## same non-mu theta.  Measured displacement on this model: regress ~-0.012,
## grad ~-0.0005 against a FOCEi MLE of tv = 3.4299.

nmTest({
  .mod <- function() {
    ini({ tka <- 0.45; tcl <- 1; tv <- c(2, 3.45, 5); add.sd <- c(0, 0.7, 5)
      eta.ka ~ 0.6; eta.cl ~ 0.3 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv)
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd) })
  }
  .ctl <- function(m) {
    vaeControl(nonMuTheta = m, print = 0L, calcTables = FALSE, returnVae = TRUE)
  }

  test_that("nonMuTheta='grad' lands closer to the FOCEi MLE than 'regress'", {
    skip_on_cran()
    .mle <- suppressMessages(
      nlmixr2(.mod(), nlmixr2data::theo_sd, est = "focei",
              control = foceiControl(print = 0L, covMethod = "", calcTables = FALSE)))$theta[["tv"]]

    .reg <- suppressWarnings(suppressMessages(rxode2::rxWithSeed(42,
      nlmixr2(.mod(), nlmixr2data::theo_sd, est = "vae", control = .ctl("regress")))))
    .grd <- suppressWarnings(suppressMessages(rxode2::rxWithSeed(42,
      nlmixr2(.mod(), nlmixr2data::theo_sd, est = "vae", control = .ctl("grad")))))

    ## the mechanism actually ran (a silent bobyqa fallback would still produce a
    ## plausible number, so assert the path, not just the value)
    expect_equal(.reg$nRegGrad, 0L)
    expect_true(.grd$nRegGrad > 0L)
    expect_equal(.grd$nRegFallback, 0L)

    .dReg <- abs(.reg$regressTheta[["tv"]] - .mle)
    .dGrd <- abs(.grd$regressTheta[["tv"]] - .mle)
    expect_true(.dGrd < .dReg)
    ## generous absolute guard: this is a stochastic training run, so pin the
    ## qualitative result (grad is close to the MLE) rather than a digit count
    expect_true(.dGrd < 0.005)
  })
})
