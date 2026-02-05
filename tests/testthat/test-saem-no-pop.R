nmTest({
  test_that("all theta parameters are fixed", {

    one.compartment <- function() {
      ini({
        tka <- fixed(log(1.57))
        tcl <- fixed(log(2.72))
        tv <- fixed(log(31.5))
        eta.ka ~ 0.6
        add.sd <- 0.7
      })
      # and a model block with the error specification and model specification
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl)
        v <- exp(tv)
        cp <- linCmt()
        cp ~ add(add.sd)
      })
    }

    f <- .nlmixr(
      one.compartment, data = theo_sd, est="saem", control = saemControl(print=0, nEm=10, nBurn=10, literalFix=FALSE))

    expect_true(inherits(f, "nlmixr2FitData"))

    m2 <- suppressMessages(suppressMessages(rxode2::rxFixPop(one.compartment)))

    expect_error(suppressMessages(m2$saemModelPred), NA)

    f <-.nlmixr(
      one.compartment, data = theo_sd, est="saem", control = saemControl(print=0, nEm=10, nBurn=10, literalFix=TRUE))

    expect_true(inherits(f, "nlmixr2FitData"))

  })
})
