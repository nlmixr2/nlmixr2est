test_that("zero etas are restored in the final ui even with nested nlmixr2 calls (#741)", {
  skip_on_cran()

  .testSeed(42)
  d <- do.call(rbind, lapply(1:12, function(i) {
    b <- 2 + rnorm(1, 0, 0.3)
    data.frame(id = i, time = 1:5, DV = 1 + b * (1:5) + rnorm(5, 0, 0.2))
  }))

  mod <- function() {
    ini({
      a <- 1
      b <- 2
      bsva ~ 0
      bsvb ~ 0.1
      addSd <- 1
    })
    model({
      ai <- a + bsva
      bi <- b + bsvb
      ci <- ai + bi * time
      ci ~ add(addSd)
    })
  }

  # est="posthoc" adds the focei objective through a nested nlmixr2() call
  # (.setOfvFo), which used to reset the zero-eta restore info before the
  # original model was put back, dropping bsva from the fit ui
  fit <- suppressMessages(suppressWarnings(
    nlmixr2(mod, data = d, est = "posthoc",
            control = foceiControl(print = 0))))

  expect_true("bsva" %in% fit$ui$iniDf$name)
  expect_true(grepl("bsva", paste(deparse(body(fit$ui$fun)), collapse = " ")))

  # and the eta can be re-enabled by piping
  newUi <- suppressMessages(rxode2::ini(fit, bsva ~ 0.1))
  expect_equal(newUi$iniDf$est[newUi$iniDf$name == "bsva"], 0.1)
})

test_that("zero iiv eta is detected alongside an IOV eta (#627)", {
  # With IOV present ui$omega is a per-condition list, so the zero-eta
  # detector used to read dimnames() off a list, return NULL and leave the
  # zero eta in the omega -- producing a singular initial 'omega' matrix.
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
  ui <- rxode2::rxUiDecompress(rxode2::assertRxUi(one.cmt))

  expect_equal(.getZeroEtasFromModel(ui), "eta.ka")

  # downgrading drops the zero eta and leaves a non-singular id omega plus
  # the untouched occ (IOV) omega
  ui2 <- .downgradeEtas(ui, "eta.ka")
  expect_false("eta.ka" %in% ui2$iniDf$name)
  expect_true("iov.cl" %in% ui2$iniDf$name)
  expect_true(is.list(ui2$omega))
  expect_false("eta.ka" %in% dimnames(ui2$omega$id)[[1]])
  expect_true(all(eigen(ui2$omega$id, only.values = TRUE)$values > 0))
})
