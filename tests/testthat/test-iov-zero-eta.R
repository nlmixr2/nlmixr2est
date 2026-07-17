test_that("IOV model with a zero iiv eta fits with focei and saem (#627)", {
  skip_on_cran()

  one.cmt <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      eta.ka ~ 0
      eta.cl ~ 0.3
      eta.v ~ 0.1
      iov.cl ~ 0.1 | occ
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl + iov.cl)
      v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }

  theo_iov <- nlmixr2data::theo_md
  theo_iov$occ <- 1
  theo_iov$occ[theo_iov$TIME >= 144] <- 2

  # focei: previously errored with "initial 'omega' matrix inverse is
  # non-positive definite"; the removed zero eta is restored in the final ui
  fitF <- suppressMessages(suppressWarnings(
    nlmixr2(one.cmt, theo_iov, est = "focei",
            control = foceiControl(print = 0, maxOuterIterations = 0))))
  expect_s3_class(fitF, "nlmixr2FitData")
  expect_true("eta.ka" %in% fitF$ui$iniDf$name)

  # saem + cwres: previously errored while restoring the original model
  # because SAEM re-expresses the IOV eta into per-occasion id etas
  fitS <- suppressMessages(suppressWarnings(
    nlmixr2(one.cmt, theo_iov, est = "saem",
            control = saemControl(print = 0, nBurn = 5, nEm = 5),
            table = list(cwres = TRUE))))
  expect_s3_class(fitS, "nlmixr2FitData")
  expect_true("CWRES" %in% names(fitS))
})
