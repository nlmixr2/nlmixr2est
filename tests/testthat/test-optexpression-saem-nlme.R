nmTest({
  # #864 for the two other estimators that build their final model after the
  # control was dropped from the ui: saem (predOnly) and nlme (EBE).
  .oneCmtEta <- function() {
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
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }

  test_that("optExpression is honored for the saem predOnly model", {
    .on <- .nlmixrMsg(.oneCmtEta, nlmixr2data::theo_sd, est = "saem",
                      control = saemControl(optExpression = TRUE, nBurn = 1, nEm = 1,
                                            print = 0),
                      table = tableControl(cwres = FALSE, npde = FALSE))
    expect_true(any(grepl("duplicate expressions in saem predOnly model", .on$msg)))
    expect_true(grepl("rx_expr", rxode2::rxNorm(.on$fit$env$saemModel$predOnly)))

    .off <- .nlmixrMsg(.oneCmtEta, nlmixr2data::theo_sd, est = "saem",
                       control = saemControl(optExpression = FALSE, nBurn = 1, nEm = 1,
                                             print = 0),
                       table = tableControl(cwres = FALSE, npde = FALSE))
    expect_false(any(grepl("duplicate expressions", .off$msg)))
    expect_false(grepl("rx_expr", rxode2::rxNorm(.off$fit$env$saemModel$predOnly)))
  })

  test_that("optExpression is honored for the nlme EBE model", {
    oneCmtLin <- function() {
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

    .on <- .nlmixrMsg(oneCmtLin, nlmixr2data::theo_sd, est = "nlme",
                      control = nlmeControl(optExpression = TRUE, verbose = FALSE),
                      table = tableControl(cwres = FALSE, npde = FALSE))
    expect_true(any(grepl("duplicate expressions in (Llik )?EBE model", .on$msg)))
    expect_true(grepl("rx_expr", rxode2::rxNorm(.on$fit$env$foceiModel$predOnly)))

    .off <- .nlmixrMsg(oneCmtLin, nlmixr2data::theo_sd, est = "nlme",
                       control = nlmeControl(optExpression = FALSE, verbose = FALSE),
                       table = tableControl(cwres = FALSE, npde = FALSE))
    expect_false(any(grepl("duplicate expressions", .off$msg)))
    expect_false(grepl("rx_expr", rxode2::rxNorm(.off$fit$env$foceiModel$predOnly)))
  })
})
