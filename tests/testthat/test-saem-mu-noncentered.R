nmTest({
  test_that("the non-centered mu step runs and can be switched off", {
    one.cmt <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    .fit <- function(on, phase = "all") {
      suppressMessages(suppressWarnings(nlmixr2(
        one.cmt,
        theo_sd,
        "saem",
        saemControl(
          nBurn = 20,
          nEm = 20,
          seed = 7,
          print = 0,
          muNonCentered = on,
          omegaNonCentered = on,
          nonCenteredPhase = phase,
          covMethod = "",
          calcTables = FALSE
        )
      )))
    }
    f1 <- .fit(TRUE)
    f0 <- .fit(FALSE)
    expect_gt(f1$saem$muNonCenteredN, 0)
    expect_equal(f0$saem$muNonCenteredN, 0)
    expect_gt(f1$saem$omegaNonCenteredN, 0)
    # burn-in only: the steps stop when the gain drops below 1
    fb <- .fit(TRUE, "burnin")
    expect_lte(fb$saem$muNonCenteredN, 20)
    expect_gt(fb$saem$muNonCenteredN, 0)
    # the step changes the trajectory, not the model: both fits stay finite
    expect_true(all(is.finite(fixef(f1))))
    expect_true(all(is.finite(fixef(f0))))
  })

  test_that("the non-centered controls validate and are off by default", {
    expect_false(saemControl()$muNonCentered)
    expect_false(saemControl()$omegaNonCentered)
    expect_equal(saemControl()$nonCenteredPhase, "all")
    expect_equal(saemControl(nonCenteredPhase = "burnin")$nonCenteredPhase, "burnin")
    expect_error(saemControl(nonCenteredPhase = "sometimes"))
    expect_error(saemControl(muNonCentered = NA))
  })
})
