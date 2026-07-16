nmTest({
  test_that("covariance conversion refreshes condition numbers; laplace/agq objDf labels", {
    one.compartment <- function() {
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
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }

    fit <- .nlmixr(one.compartment, nlmixr2data::theo_sd, "focei",
                   foceiControl(print = 0))
    .cmBefore <- fit$covMethod
    .cnBefore <- fit$objDf[["Condition#(Cov)"]][1]
    expect_true(is.finite(.cnBefore))

    # converting to a different covariance must refresh the condition numbers
    invisible(nlmixr2est:::.setCov(fit, covMethod = "s", covType = "fd"))
    .cnAfter <- fit$objDf[["Condition#(Cov)"]][1]
    expect_true(is.finite(.cnAfter))
    expect_true(abs(.cnAfter - .cnBefore) > 1e-4)
    expect_equal(.cnAfter, fit$env$conditionNumberCov)
    expect_equal(fit$objDf[["Condition#(Cor)"]][1], fit$env$conditionNumberCor)
    .ev <- eigen(fit$cov, symmetric = TRUE, only.values = TRUE)$values
    expect_equal(.cnAfter, max(abs(.ev)) / min(abs(.ev)))

    # switching back through the covList keeps them in sync too
    invisible(setCov(fit, .cmBefore))
    expect_equal(fit$objDf[["Condition#(Cov)"]][1], .cnBefore, tolerance = 1e-6)

    # quadrature fits label their objDf row by method, not "FOCEi", and
    # ofvType matches the row so broom/setOfv row lookup works
    fitL <- .nlmixr(one.compartment, nlmixr2data::theo_sd, "laplace",
                    laplaceControl(print = 0))
    expect_equal(rownames(fitL$objDf)[1], "Laplace")
    expect_equal(fitL$ofvType, "laplace")
    expect_true(is.finite(fitL$objDf[["Condition#(Cov)"]][1]))

    fitA <- .nlmixr(one.compartment, nlmixr2data::theo_sd, "agq",
                    agqControl(print = 0, nAGQ = 2))
    expect_equal(rownames(fitA$objDf)[1], "AGQ2")
    expect_equal(fitA$ofvType, "agq2")
    expect_true(is.finite(fitA$objDf[["Condition#(Cov)"]][1]))

    # setOfv(fit, "focei") on a quadrature fit adds the real focei objective
    # (nAGQ is zeroed for the evaluation) and setOfv can switch back by row
    .agqObjf <- fitA$objDf["AGQ2", "OBJF"]
    fitA <- suppressMessages(setOfv(fitA, "focei"))
    expect_true("FOCEi" %in% rownames(fitA$objDf))
    expect_equal(fitA$ofvType, "FOCEi") # setOfv stores the type as passed
    expect_true(abs(fitA$objDf["FOCEi", "OBJF"] - .agqObjf) > 1e-4)
    fitA <- setOfv(fitA, "agq2")
    expect_equal(fitA$env$OBJF, .agqObjf)
  })
})
