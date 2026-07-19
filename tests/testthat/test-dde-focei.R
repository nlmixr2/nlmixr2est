# Delay differential equation (DDE) support in the FOCEI fast family, including
# the analytic covariance.  Guards the augmented-sensitivity-model pipeline
# (dense dop853 solve of the delayed sensitivity ODEs) end to end: a broken
# delay history or dropped sensitivity ODE shows up as non-finite / exploded SEs.

nmTest({
  .dde_focei_data <- function() {
    # method-of-steps reference DDE: y' = -k*delay(y, 1), constant history y = a.
    set.seed(42)
    .tk <- log(0.3); .ta <- log(2)
    .etak <- stats::rnorm(8, 0, sqrt(0.03))
    do.call(rbind, lapply(seq_len(8), function(id) {
      .k <- exp(.tk + .etak[id]); .a <- exp(.ta)
      .m <- rxode2::rxode2(sprintf(
        "k=%.10g\na=%.10g\ny(0)<-a\nd/dt(y)<- -k*delay(y,1)\npast(y,1)<-a\n", .k, .a))
      .s <- rxode2::rxSolve(.m, rxode2::et(seq(0.5, 5, by = 0.5)),
                            atol = 1e-9, rtol = 1e-9)
      data.frame(id = id, time = .s$time, dv = .s$y + stats::rnorm(nrow(.s), 0, 0.05))
    }))
  }

  .dde_focei_mod <- function() {
    ini({
      tk <- log(0.3); ta <- log(2)
      eta.k ~ 0.03
      add.sd <- 0.05
    })
    model({
      k <- exp(tk + eta.k); a <- exp(ta)
      y(0) <- a
      d/dt(y) <- -k * delay(y, 1)
      past(y, 1) <- a
      y ~ add(add.sd)
    })
  }

  test_that("DDE FOCEI (fast) fit produces a finite analytic covariance", {
    skip_on_cran()
    skip_on_ci()
    .dat <- .dde_focei_data()
    fit <- suppressWarnings(suppressMessages(
      nlmixr(.dde_focei_mod, .dat, "focei",
             foceiControl(print = 0L, fast = TRUE, covMethod = "analytic"))))
    # the augmented delayed-sensitivity solve fed the covariance without diverging
    expect_equal(fit$covMethod, "analytic")
    expect_true(all(is.finite(fit$parFixedDf$SE)))
    expect_true(all(fit$parFixedDf$SE > 0))
    # structural theta recovered near truth (tk=log(0.3), ta=log(2))
    expect_equal(unname(fit$theta[c("tk", "ta")]),
                 c(log(0.3), log(2)), tolerance = 0.15)
  })
})
