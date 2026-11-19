nmTest({
  test_that("Issue nlmixr2est#515: a constant likelihood rate gives a clear error, not 'Aborted calculation'", {
    # 'lamba' is a fixed population parameter, so the poisson likelihood does
    # not depend on any random effect; focei used to mask the real cause with
    # the generic "Aborted calculation" message.
    mod <- function() {
      ini({
        te0 <- log(10)
        eta.e0 ~ 0.9
        lamba <- 0.1
      })
      model({
        e0 = exp(te0 + eta.e0)
        kout = 1.1
        effect(0) = e0
        kin = e0 * kout
        d/dt(effect) = kin - kout * effect
        effect ~ dpois(lamba)
      })
    }
    d <- data.frame(ID = c(1, 1, 2, 2), TIME = c(0, 1, 0, 1), DV = c(1, 2, 1, 3), EVID = 0)
    .e <- expect_error(
      .nlmixr(mod, d, est = "focei", control = foceiControl(print = 0)),
      "none of the model predictions depend on a random effect"
    )
    # the informative error is no longer masked by the generic abort message
    expect_false(any(grepl("Aborted calculation", conditionMessage(.e))))
  })

  test_that("the d(prediction)/d(ETA) check warns only when some, not all, are zero", {
    # eta.z enters a variable the prediction does not use, so the prediction's
    # derivative with respect to it is identically zero
    .some <- function() {
      ini({
        tcl <- log(2)
        eta.cl ~ 0.1
        eta.z ~ 0.1
        add.sd <- 0.5
      })
      model({
        cl <- exp(tcl + eta.cl)
        z <- exp(eta.z)
        cp <- 10 * exp(-cl * t)
        cp ~ add(add.sd)
      })
    }
    .none <- function() {
      ini({
        tcl <- log(2)
        tv <- log(20)
        eta.cl ~ 0.1
        eta.v ~ 0.1
        add.sd <- 0.5
      })
      model({
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        cp <- 10 / v * exp(-cl / v * t)
        cp ~ add(add.sd)
      })
    }
    expect_warning(
      .s <- rxode2::assertRxUi(.some)$foceiHdEta,
      "some of the predictions do not depend on 'ETA'",
      fixed = TRUE
    )
    expect_identical(.s$..HdEta[2], "rx__sens_rx_pred__BY_ETA_2___=0")
    expect_no_warning(rxode2::assertRxUi(.none)$foceiHdEta)
  })
})
