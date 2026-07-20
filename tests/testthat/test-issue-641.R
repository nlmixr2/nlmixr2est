## Regression test for issue 641: FOCEI wasn't updating large-magnitude
## additive mu-referenced thetas because scaleC fell through to the C++
## default of 1/|init|, shrinking scaled-space steps to nothing for |init| >> 1.

nmTest({
  test_that("issue 641: FOCEI moves large-magnitude additive theta in algebraic model", {
    skip_on_cran()

    emax_fun <- function(x, e0, emax, ex50, hill = 1) {
      e0 + emax * x^hill / (ex50^hill + x^hill)
    }

    .testSeed(641)
    d_mod <- data.frame(
      WEEK = rep(0:12, 5),
      ID = rep(1:5, each = 13),
      COVV = 1
    )
    d_mod$emax_val <- rnorm(nrow(d_mod), mean = 20, sd = 5) + d_mod$ID
    d_mod$DV <- emax_fun(d_mod$WEEK, e0 = 0, emax = d_mod$emax_val, ex50 = 4)
    d_mod$TIME <- d_mod$WEEK

    mod <- function() {
      ini({
        e0 <- fixed(0)
        tvemax <- -40
        covValueEmax <- fixed(0)
        let50 <- log(3)
        covValueET50 <- fixed(0)
        iivemax ~ 10
        iivet50 ~ 0
        addSd <- 5
      })
      model({
        emax <- tvemax + iivemax + covValueEmax * COVV
        et50 <- exp(let50 + iivet50 + covValueET50 * COVV)
        eff <- e0 + emax * WEEK / (et50 + WEEK)
        eff ~ add(addSd)
      })
    }

    # The parallel inner loop is non-deterministic on poorly-conditioned models
    # (subject-level reductions reorder under cores > 1), so this test pins the
    # thread count to 1 to get a reproducible trajectory.  The bug being checked
    # is the scaleC computation, which is independent of threading.
    prevThreads <- rxode2::getRxThreads()
    on.exit(rxode2::setRxThreads(prevThreads), add = TRUE)
    rxode2::setRxThreads(1L)

    fit <- nlmixr(
      mod, data = d_mod, est = "focei",
      control = foceiControl(print = 0, covMethod = "", calcTables = FALSE)
    )

    tvemaxEst <- fit$parFixedDf["tvemax", "Estimate"]
    let50Est  <- fit$parFixedDf["let50",  "Estimate"]

    # Before the fix tvemax stayed within ~0.05 of the initial -40; after the
    # fix it moves tens of units away. Thresholds are loose to absorb optimizer
    # noise while still catching a regression to "no movement".
    expect_gt(abs(tvemaxEst - (-40)), 30)
    expect_lt(fit$objf, 1000)

    # The scaleCtheta accessor exposes what the fix computed; tvemax (init -40)
    # should now resolve to a sensible large positive scaleC.
    ui <- fit$ui
    scaleCs <- ui$scaleCtheta
    expect_true(any(scaleCs > 10),
                info = "expected at least one large-magnitude theta to get scaleC > 10")
  })
})
