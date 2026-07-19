nmTest({
  test_that("cwres (and focei objective fun) is added to saem with addCwres", {
    .cloneFit <- function(fit) rlang::duplicate(fit, shallow = FALSE)

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

    baseFit <- .nlmixr(one.cmt, theo_sd, est = "saem", control = saemControlFast)

    fit <- .cloneFit(baseFit)


    expect_false(all(c("NPDE","EPRED","NPD","NPDE") %in% names(fit)))
    expect_warning(fit$etaSE)
    expect_warning(fit$etaRSE)
    expect_warning(fit$etaCI)
    expect_warning(fit$etaR)
    expect_false(any(names(fit$dataMergeInner) == "nlmixrLlikObs"))
    suppressMessages(expect_error(addCwres(fit), NA))
    expect_true(all(c("WRES","CPRED","CRES","CWRES") %in% names(fit)))
    expect_equal(row.names(fit$objDf)[1], "FOCEi")
    expect_false(is.null(fit$etaSE))
    expect_false(is.null(fit$etaRSE))
    expect_false(is.null(fit$etaCI))
    expect_false(is.null(fit$etaR))
    expect_true(any(names(fit$dataMergeInner) == "nlmixrLlikObs"))

    # etaCI brackets the eta estimate and matches etaSE at the fit ci level
    .etaCI <- fit$etaCI
    .etaSE <- fit$etaSE
    .eta <- fit$eta
    expect_equal(names(.etaCI),
                 c("ID", "eta.ka (2.5%)", "eta.ka (97.5%)",
                   "eta.cl (2.5%)", "eta.cl (97.5%)",
                   "eta.v (2.5%)", "eta.v (97.5%)"))
    .qn <- stats::qnorm(0.975)
    expect_equal(.etaCI[["eta.ka (2.5%)"]],
                 .eta$eta.ka - .qn * .etaSE[["se(eta.ka)"]])
    expect_equal(.etaCI[["eta.ka (97.5%)"]],
                 .eta$eta.ka + .qn * .etaSE[["se(eta.ka)"]])

    fit <- .cloneFit(baseFit)

    expect_false(all(c("WRES","CPRED","CRES","CWRES") %in% names(fit)))
    suppressMessages(expect_error(addCwres(fit, focei=FALSE), NA))
    expect_true(all(c("WRES","CPRED","CRES","CWRES") %in% names(fit)))
    expect_equal(row.names(fit$objDf)[1], "FOCE")

    fit <- .cloneFit(baseFit)

    expect_false(all(c("WRES","CPRED","CRES","CWRES") %in% names(fit)))
    fit2 <- suppressMessages(addCwres(fit, updateObject=FALSE))
    expect_false(all(c("WRES","CPRED","CRES","CWRES") %in% names(fit)))
    expect_true(all(c("WRES","CPRED","CRES","CWRES") %in% names(fit2)))
    expect_equal(row.names(fit2$objDf)[1], "FOCEi")
    expect_false(is.null(fit2$etaSE))

    fit <- one.compartment.fit.saem.cwres

    expect_true(all(c("WRES","CPRED","CRES","CWRES") %in% names(fit)))
    expect_false(is.null(fit$etaSE))
    expect_false(is.null(fit$etaRSE))
    expect_false(is.null(fit$etaR))
    expect_true(any(names(fit$dataMergeInner) == "nlmixrLlikObs"))
    expect_true(any(names(fit$fitMergeInner) == "nlmixrLlikObs"))
  })
})
