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
    suppressMessages(setOfv(fit, "laplace1.6"))
    .row <- fit$objDf["laplace1.6", ]
    expect_equal(.row$OBJF, -2 * .row[["Log-likelihood"]])
    # a name that is neither "laplace<nsd>" nor "gauss<nnodes>_<nsd>"
    expect_error(setOfv(fit, "gauss3"), "cannot switch objective function to 'gauss3' type")
  })

  test_that("setOfv() switches a saem fit to an importance-sampling objective", {
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
    fit <- .nlmixr(
      one.cmt,
      theo_sd,
      est = "saem",
      control = saemControl(print = 0, nBurn = 10, nEm = 10, calcTables = FALSE)
    )
    # the first read of the objective computes the saem quadrature one
    expect_true(is.finite(fit$objf))
    .type0 <- fit$ofvType
    expect_identical(.type0, "gauss3_1.6")
    .extra0 <- fit$env$extra
    suppressWarnings(suppressMessages(setOfv(fit, "imp")))
    expect_identical(fit$ofvType, "IMP")
    expect_identical(fit$env$extra, crayon::silver$italic("OBJF by importance sampling (IMP)"))
    expect_identical(fit$objf, fit$objDf["IMP", "OBJF"])
    expect_identical(as.numeric(fit$logLik), fit$objDf["IMP", "Log-likelihood"])
    expect_identical(fit$AIC, fit$objDf["IMP", "AIC"])
    # and back to the objective the fit had
    setOfv(fit, .type0)
    expect_identical(fit$ofvType, .type0)
    expect_identical(fit$env$extra, .extra0)
    expect_identical(fit$objf, fit$objDf[.type0, "OBJF"])
    # an objective the saem fit cannot describe stops before anything is switched
    .state <- function(f) list(f$ofvType, f$env$extra, f$objf, as.numeric(f$logLik), f$AIC, f$BIC, f$objDf)
    .before <- .state(fit)
    .row <- data.frame(OBJF = 1, AIC = 2, BIC = 3, "Log-likelihood" = -0.5, check.names = FALSE)
    expect_error(
      nlmixrAddObjectiveFunctionDataFrame(fit, .row, "custom"),
      "the saem fit has no description of objective function 'custom'",
      fixed = TRUE
    )
    expect_identical(.state(fit), .before)
    # nothing was added, so the same call fails the same way again
    expect_error(
      nlmixrAddObjectiveFunctionDataFrame(fit, .row, "custom"),
      "the saem fit has no description of objective function 'custom'",
      fixed = TRUE
    )
  })
})
