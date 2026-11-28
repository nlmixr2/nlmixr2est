nmTest({
  # rxSolve_() does not convert an ODE model to matExp() form, so
  # method = "indLin" on an ODE model failed every solve and the fit returned
  # the failed-solve penalty (2400 for theo) instead of an error
  one.cmt <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      eta.cl ~ 0.3
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }

  test_that("focei with method = 'indLin' on an ODE model errors at the start", {
    expect_error(
      suppressMessages(nlmixr(
        one.cmt,
        nlmixr2data::theo_sd,
        "focei",
        foceiControl(
          maxOuterIterations = 0,
          print = 0,
          covMethod = "",
          calcTables = FALSE,
          rxControl = rxode2::rxControl(method = "indLin")
        )
      )),
      "matExp"
    )
  })

  test_that("nlm with method = 'indLin' on an ODE model errors", {
    pk <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl)
        v <- exp(tv)
        d / dt(depot) <- -ka * depot
        d / dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }
    d <- nlmixr2data::theo_sd
    d <- d[d$ID == 1, ]
    expect_error(
      suppressMessages(nlmixr(
        pk,
        d,
        "nlm",
        nlmControl(rxControl = rxode2::rxControl(method = "indLin"), print = 0)
      )),
      "matExp"
    )
  })
})
