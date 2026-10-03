nmTest({
  test_that("a saem quadrature objective falls back to the control's adjObf", {
    one.cmt <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv)
        linCmt() ~ add(add.sd)
      })
    }
    .ctl <- saemControl(
      print = 0,
      nBurn = 1,
      nEm = 1,
      nmc = 1,
      nu = c(1, 1, 1),
      adjObf = FALSE,
      calcTables = FALSE
    )
    fit <- .nlmixr(one.cmt, theo_sd, est = "saem", control = .ctl)
    # no adjObf on the fit environment, as in a fit saved before it was stored
    rm("adjObf", envir = fit$env)
    suppressMessages(setOfv(fit, "gauss3_1.6"))
    .row <- fit$objDf["gauss3_1.6", ]
    # adjObf = FALSE: OBJF keeps the normal constant, so it is -2 * log-likelihood
    expect_equal(.row$OBJF, -2 * .row[["Log-likelihood"]])
    expect_equal(fit$objf, .row$OBJF)
    # a name that is neither "laplace<nsd>" nor "gauss<nnodes>_<nsd>"
    expect_error(setOfv(fit, "gauss3"), "cannot switch objective function to 'gauss3' type")
  })
})
