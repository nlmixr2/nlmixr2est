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
})
