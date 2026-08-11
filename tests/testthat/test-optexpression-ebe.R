nmTest({
  # #864: the fit finalization removed the control from the ui before the
  # EBE/predOnly model was built, so optExpression (and sumProd) silently fell
  # back to the rxGetControl() defaults for that build.
  test_that("optExpression is honored for the nlm-family EBE model", {
    oneCmt <- function() {
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
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }

    .on <- .nlmixrMsg(oneCmt, nlmixr2data::theo_sd, est = "bobyqa",
                      control = bobyqaControl(optExpression = TRUE, maxfun = 10,
                                              print = 0),
                      table = tableControl(cwres = FALSE, npde = FALSE))
    expect_true(any(grepl("duplicate expressions in (Llik )?EBE model", .on$msg)))
    expect_true(grepl("rx_expr", rxode2::rxNorm(.on$fit$env$foceiModel$predOnly)))

    .off <- .nlmixrMsg(oneCmt, nlmixr2data::theo_sd, est = "bobyqa",
                       control = bobyqaControl(optExpression = FALSE, maxfun = 10,
                                               print = 0),
                       table = tableControl(cwres = FALSE, npde = FALSE))
    # the EBE model is still compiled ("compiling EBE model..."), just not rewritten
    expect_false(any(grepl("duplicate expressions", .off$msg)))
    expect_false(grepl("rx_expr", rxode2::rxNorm(.off$fit$env$foceiModel$predOnly)))
  })

  test_that("sumProd is honored for the nlm-family EBE model", {
    oneCmt <- function() {
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
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }

    .on <- .nlmixrMsg(oneCmt, nlmixr2data::theo_sd, est = "bobyqa",
                      control = bobyqaControl(sumProd = TRUE, optExpression = FALSE,
                                              maxfun = 10, print = 0),
                      table = tableControl(cwres = FALSE, npde = FALSE))
    expect_true(any(grepl("stabilizing round off errors in predictions or EBE model",
                          .on$msg)))
    .txt <- rxode2::rxNorm(.on$fit$env$foceiModel$predOnly)
    expect_true(grepl("prod(", .txt, fixed = TRUE))
    expect_true(grepl("sum(", .txt, fixed = TRUE))
  })
})
