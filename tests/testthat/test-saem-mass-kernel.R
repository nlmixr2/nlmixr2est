test_that("saemControl() validates the mass-matrix kernel options", {
  expect_equal(saemControl()$nuMass, 0L)
  expect_equal(saemControl(nuMass = 2)$nuMass, 2L)
  expect_error(saemControl(nuMass = -1))
  expect_error(saemControl(massStart = 1))
  expect_error(saemControl(massKappa = 0.2))
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
  .fit <- function(nuMass) {
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
        nuMass = nuMass
      )
    )
  }

  test_that("the mass-matrix kernel samples, is seeded and agrees with the default kernels", {
    .f0 <- .fit(0L)
    .f2 <- .fit(2L)
    expect_true(is.na(.f0$saem$massAccept))
    expect_gt(.f2$saem$massAccept, 0.1)
    expect_lt(.f2$saem$massAccept, 0.95)
    # the same seed reproduces the fit exactly
    expect_equal(fixef(.fit(2L)), fixef(.f2))
    expect_equal(fixef(.f2), fixef(.f0), tolerance = 0.1)
  })
})
