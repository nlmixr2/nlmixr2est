test_that("IOV model with a dotted occasion eta name fits with focei and saem (#1233)", {
  skip_on_cran()

  one.cmt <- function() {
    ini({
      tka <- log(1.5)
      tcl <- log(2.7)
      tv <- log(30)
      eta.cl ~ 0.1
      eta.v ~ 0.1
      eta.ka.occ ~ 0.1 | OCC
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka.occ)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }

  d <- nlmixr2data::theo_md
  d$OCC <- ave(d$EVID != 0, d$ID, FUN = cumsum)

  fitF <- suppressMessages(suppressWarnings(
    nlmixr2(one.cmt, d, est = "focei",
            control = foceiControl(maxOuterIterations = 0, covMethod = "", print = 0))
  ))
  expect_s3_class(fitF, "nlmixr2FitData")
  expect_true("OCC" %in% names(fitF$omega))

  fitS <- suppressMessages(suppressWarnings(
    nlmixr2(one.cmt, d, est = "saem", control = saemControl(print = 0, nBurn = 5, nEm = 5))
  ))
  expect_s3_class(fitS, "nlmixr2FitData")
  expect_true("OCC" %in% names(fitS$omega))
})
