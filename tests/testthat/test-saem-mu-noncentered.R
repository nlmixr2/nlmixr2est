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
    .fit <- function(on) {
      suppressMessages(suppressWarnings(nlmixr2(
        one.cmt, theo_sd, "saem",
        saemControl(nBurn = 20, nEm = 20, seed = 7, print = 0,
                    muNonCentered = on, covMethod = "", calcTables = FALSE)
      )))
    }
    f1 <- .fit(TRUE)
    f0 <- .fit(FALSE)
    expect_gt(f1$saem$muNonCenteredN, 0)
    expect_equal(f0$saem$muNonCenteredN, 0)
    # the step changes the trajectory, not the model: both fits stay finite
    expect_true(all(is.finite(fixef(f1))))
    expect_true(all(is.finite(fixef(f0))))
  })

  test_that("mode 1B controls validate and default by declaration", {
    expect_error(saemControl(kappa1B = 1.5))
    expect_error(saemControl(kappa1B = -0.1))
    expect_null(saemControl()$nu1B)
    expect_equal(saemControl(nu1B = 3)$nu1B, 3L)
    expect_equal(saemControl(kappa1B = 0)$kappa1B, 0)
  })
})
