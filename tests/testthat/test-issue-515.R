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
    d <- data.frame(ID = c(1, 1, 2, 2), TIME = c(0, 1, 0, 1),
                    DV = c(1, 2, 1, 3), EVID = 0)
    .e <- expect_error(
      .nlmixr(mod, d, est = "focei", control = foceiControl(print = 0)),
      "none of the model predictions depend on a random effect"
    )
    # the informative error is no longer masked by the generic abort message
    expect_false(any(grepl("Aborted calculation", conditionMessage(.e))))
  })
})
