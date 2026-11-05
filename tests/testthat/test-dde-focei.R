# Delay differential equation (DDE) support in the FOCEI fast family, including
# the analytic covariance.  Guards the augmented-sensitivity-model pipeline
# (dense dop853 solve of the delayed sensitivity ODEs) end to end: a broken
# delay history or dropped sensitivity ODE shows up as non-finite / exploded SEs.

nmTest({
  .dde_focei_data <- function() {
    # method-of-steps reference DDE: y' = -k*delay(y, 1), constant history y = a.
    .testSeed(42)
    .tk <- log(0.3)
    .ta <- log(2)
    .etak <- stats::rnorm(8, 0, sqrt(0.03))
    do.call(
      rbind,
      lapply(seq_len(8), function(id) {
        .k <- exp(.tk + .etak[id])
        .a <- exp(.ta)
        .m <- rxode2::rxode2(sprintf(
        "k=%.10g\na=%.10g\ny(0)<-a\nd/dt(y)<- -k*delay(y,1)\npast(y,1)<-a\n", .k, .a))
        .s <- rxode2::rxSolve(.m, rxode2::et(seq(0.5, 5, by = 0.5)), atol = 1e-9, rtol = 1e-9)
        data.frame(id = id, time = .s$time, dv = .s$y + stats::rnorm(nrow(.s), 0, 0.05))
      })
    )
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
    # the delayed-sensitivity augmented solve needs an accurate tolerance to keep
    # the observed information PD; at the sigdig=4 default it diverges enough to
    # fall back to the "r" covariance, so pin sigdig=6 to exercise the analytic path.
    fit <- suppressWarnings(suppressMessages(
      nlmixr(.dde_focei_mod, .dat, "focei", foceiControl(print = 0L, fast = TRUE, covMethod = "analytic", sigdig = 6))
    ))
    # the augmented delayed-sensitivity solve fed the covariance without diverging
    expect_equal(.covBaseName(fit$covMethod), "analytic")
    expect_true(all(is.finite(fit$parFixedDf$SE)))
    expect_true(all(fit$parFixedDf$SE > 0))
    # structural theta recovered near truth (tk=log(0.3), ta=log(2))
    expect_equal(unname(fit$theta[c("tk", "ta")]), c(log(0.3), log(2)), tolerance = 0.15)
  })
})

# A parameter in a delay() duration is finite-differenced (#1300): the
# symbolic sensitivities miss the lag's own term and the dose breakpoint jumps.

nmTest({
  .ddeLagMod <- function() {
    ini({
      lkel <- log(0.3)
      lv <- log(5)
      lkin <- log(2)
      lTR <- log(6)
      ls <- log(0.5)
      eta.kin ~ 0.05
      eta.TR ~ 0.05
      pkadd <- 0.1
      pdprop <- 0.05
    })
    model({
      kel <- exp(lkel)
      v <- exp(lv)
      kin <- exp(lkin + eta.kin)
      TR <- exp(lTR + eta.TR)
      s <- exp(ls)
      C <- central / v
      R(0) <- kin * TR
      d/dt(central) <- -kel * central
      d/dt(R) <- kin * s * C - kin * s * delay(central, TR) / v
      C ~ add(pkadd) | central
      R ~ prop(pdprop) | R
    })
  }

  test_that("delay() duration parameters are found", {
    .ui <- rxode2::rxode2(.ddeLagMod)
    .dl <- .foceiDelayLagParams(.ui$loadPruneSens)
    expect_equal(.dl$eta, 2L)
    expect_equal(.dl$theta, 4L)
  })

  test_that("a constant delay() has no lag parameters", {
    .ui <- rxode2::rxode2(function() {
      ini({
        tk <- log(0.3)
        eta.k ~ 0.03
        add.sd <- 0.05
      })
      model({
        k <- exp(tk + eta.k)
        y(0) <- 2
        d/dt(y) <- -k * delay(y, 1)
        y ~ add(add.sd)
      })
    })
    .dl <- .foceiDelayLagParams(.ui$loadPruneSens)
    expect_length(.dl$eta, 0L)
    expect_length(.dl$theta, 0L)
  })

  test_that("an estimated delay() duration tightens the sigdig tolerances", {
    .rc <- .rxControlDelayLagTol(foceiControl()$rxControl, 3)
    expect_equal(.rc$rtol, 1e-9)
    expect_equal(.rc$atol, 1e-11)
    expect_equal(.rc$rtolSens, 1e-9)
    # a tolerance the user set is kept
    .rc <- .rxControlDelayLagTol(foceiControl(rxControl = list(rtol = 1e-5))$rxControl, 3)
    expect_equal(.rc$rtol, 1e-5)
    expect_equal(.rc$atol, 1e-11)
    .rc <- .rxControlDelayLagTol(foceiControl(rxControl = list(atol = 1e-14))$rxControl, 3)
    expect_equal(.rc$atol, 1e-14)
    expect_equal(.rc$rtol, 1e-9)
    # an explicit tolerance equal to the sigdig value, and a full rxControl(), are kept
    .rc <- .rxControlDelayLagTol(foceiControl(rxControl = list(rtol = 1e-6))$rxControl, 3)
    expect_equal(.rc$rtol, 1e-6)
    .rc <- .rxControlDelayLagTol(foceiControl(rxControl = rxode2::rxControl(rtol = 1e-6))$rxControl, 3)
    expect_equal(.rc$rtol, 1e-6)
    # a generated rxControl the user changed and passed back keeps the change
    .rc0 <- foceiControl()$rxControl
    .rc0$rtol <- 1e-12
    .rc <- .rxControlDelayLagTol(foceiControl(rxControl = .rc0)$rxControl, 3)
    expect_equal(.rc$rtol, 1e-12)
    .rc0 <- foceiControl()$rxControl
    .rc0$rtol <- 1e-12
    expect_equal(.rxControlDelayLagTol(.rc0, 3)$rtol, 1e-12)
  })
})
