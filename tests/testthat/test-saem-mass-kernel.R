test_that("saemControl() validates the mode 1B options", {
  expect_equal(saemControl()$nu1B, 0L)
  expect_equal(saemControl(nu1B = 2)$nu1B, 2L)
  expect_error(saemControl(nu1B = -1))
  expect_error(saemControl(nb1B = 1))
  expect_error(saemControl(kappa1B = 0.2))
})

nmTest({
  .mod <- function() {
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
  .fit <- function(nu1B) {
    .nlmixr(
      .mod,
      nlmixr2data::theo_sd,
      est = "saem",
      control = saemControl(
        nBurn = 60,
        nEm = 60,
        seed = 42,
        print = 0,
        covMethod = "",
        calcTables = FALSE,
        nu1B = nu1B
      )
    )
  }

  test_that("mode 1B samples, is seeded and agrees with the default kernels", {
    .f0 <- .fit(0L)
    .f2 <- .fit(2L)
    expect_true(is.na(.f0$saem$accept1B))
    expect_gt(.f2$saem$accept1B, 0.1)
    expect_lt(.f2$saem$accept1B, 0.95)
    # the same seed reproduces the fit exactly
    expect_equal(fixef(.fit(2L)), fixef(.f2))
    expect_equal(fixef(.f2), fixef(.f0), tolerance = 0.1)
    # a wrong-signed proposal correction targets the wrong posterior (+5 OFV)
    expect_lt(abs(.f2$objf - .f0$objf), 2)
  })

  test_that("mode 1B proposes around an uninformative (masked) eta", {
    .modUe <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        tz <- 0
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        eta.z ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        z <- tz + eta.z
        linCmt() ~ add(add.sd)
      })
    }
    .f <- .nlmixr(
      .modUe,
      nlmixr2data::theo_sd,
      est = "saem",
      control = saemControl(
        nBurn = 40,
        nEm = 20,
        seed = 42,
        print = 0,
        covMethod = "",
        calcTables = FALSE,
        nu1B = 2L
      )
    )
    expect_gt(.f$saem$accept1B, 0.05)
    expect_true(all(is.finite(fixef(.f))))
  })
})
