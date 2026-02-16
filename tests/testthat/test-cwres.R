nmTest({
  test_that("cwres (and focei objective fun) is added to saem with addCwres", {

    one.cmt <- function() {
      ini({
        ## You may label each parameter with a comment
        tka <- 0.45 # Ka
        tcl <- log(c(0, 2.7, 100)) # Log Cl
        ## This works with interactive models
        ## You may also label the preceding line with label("label text")
        tv <- 3.45; label("log V")
        ## the label("Label name") works with all models
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

    fit <- .nlmixr(one.cmt, theo_sd, est = "saem", control = saemControlFast)


    expect_false(all(c("NPDE","EPRED","NPD","NPDE") %in% names(fit)))
    expect_warning(fit$etaSE)
    expect_warning(fit$etaRSE)
    expect_warning(fit$etaR)
    expect_false(any(names(fit$dataMergeInner) == "nlmixrLlikObs"))
    suppressMessages(expect_error(addCwres(fit), NA))
    expect_true(all(c("WRES","CPRED","CRES","CWRES") %in% names(fit)))
    expect_equal(row.names(fit$objDf), "FOCEi")
    expect_false(is.null(fit$etaSE))
    expect_false(is.null(fit$etaRSE))
    expect_false(is.null(fit$etaR))
    expect_true(any(names(fit$dataMergeInner) == "nlmixrLlikObs"))

    fit <- .nlmixr(one.cmt, theo_sd, est="saem", control = saemControlFast)

    expect_false(all(c("WRES","CPRED","CRES","CWRES") %in% names(fit)))
    suppressMessages(expect_error(addCwres(fit, focei=FALSE), NA))
    expect_true(all(c("WRES","CPRED","CRES","CWRES") %in% names(fit)))
    expect_equal(row.names(fit$objDf), "FOCE")

    fit <- .nlmixr(one.cmt, theo_sd, est="saem", control = saemControlFast)

    expect_false(all(c("WRES","CPRED","CRES","CWRES") %in% names(fit)))
    fit2 <- suppressMessages(addCwres(fit, updateObject=FALSE))
    expect_false(all(c("WRES","CPRED","CRES","CWRES") %in% names(fit)))
    expect_true(all(c("WRES","CPRED","CRES","CWRES") %in% names(fit2)))
    expect_equal(row.names(fit2$objDf), "FOCEi")
    expect_false(is.null(fit2$etaSE))

    fit <-
      .nlmixr(
        one.cmt, theo_sd, est = "saem",
        control = saemControlFast,
        table = tableControl(cwres = TRUE)
      )

    expect_true(all(c("WRES","CPRED","CRES","CWRES") %in% names(fit)))
    expect_false(is.null(fit$etaSE))
    expect_false(is.null(fit$etaRSE))
    expect_false(is.null(fit$etaR))
    expect_true(any(names(fit$dataMergeInner) == "nlmixrLlikObs"))
    expect_true(any(names(fit$fitMergeInner) == "nlmixrLlikObs"))
  })
})
