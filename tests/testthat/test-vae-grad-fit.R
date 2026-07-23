## Slow (weekly batch): full grad-vs-regress fits on theo_sd.
##
## Both modes drive the SAME (full outer) objective -- the Laplace determinant +
## 0.5*log|Omega^-1| + the DV-transform Jacobian -- with every mu-referenced theta
## held at its current M-step value.  They differ in the M-step update they take:
## "regress" fully re-solves the regressed thetas with a derivative-free bounded
## bobyqa search each M-step, while "grad" takes ONE Adam step along the exact
## analytic gradient each M-step.  So this is a comparison of two training
## dynamics, not two optimizers converging to one point.
##
## Measured on this model against a FOCEi MLE of tv = 3.4293: regress ~3.4291
## (d ~ 2e-4), grad ~3.4268 (d ~ 2.5e-3).  Both land within 5e-3 of the MLE, so
## the durable invariant is "the grad mechanism ran and produced a fit close to
## the MLE" -- NOT a strict ordering between the two (regress fully re-solving each
## M-step can, and here does, land marginally closer than grad's single step).

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

  test_that("nonMuTheta='grad' runs the analytic-gradient M-step and fits near the FOCEi MLE", {
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
    ## both training schemes land near the MLE; assert the durable invariant
    ## (each is close) rather than a fragile strict ordering between two heuristic
    ## M-step updates -- the analytic-gradient step is competitive with, not
    ## strictly better than, fully re-solving each M-step.
    expect_true(.dGrd < 0.005)
    expect_true(.dReg < 0.005)
  })
})
