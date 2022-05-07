nmTest({

  test_that("cwres (and focei objective fun) is added to saem with addCwres", {

    .nlmixr <- function(...) suppressWarnings(suppressMessages(nlmixr(...)))

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

    fit <- .nlmixr(one.cmt, theo_sd, est="saem")

    expect_false(all(c("NPDE","EPRED","NPD","NPDE") %in% names(fit)))
    suppressMessages(expect_error(addCwres(fit), NA))
    expect_true(all(c("WRES","CPRED","CRES","CWRES") %in% names(fit)))
    expect_equal(row.names(fit$objDf), "FOCEi")

    fit <- .nlmixr(one.cmt, theo_sd, est="saem")

    expect_false(all(c("WRES","CPRED","CRES","CWRES") %in% names(fit)))
    suppressMessages(expect_error(addCwres(fit, focei=FALSE), NA))
    expect_true(all(c("WRES","CPRED","CRES","CWRES") %in% names(fit)))
    expect_equal(row.names(fit$objDf), "FOCE")

    fit <- .nlmixr(one.cmt, theo_sd, est="saem")

    expect_false(all(c("WRES","CPRED","CRES","CWRES") %in% names(fit)))
    fit2 <- suppressMessages(addCwres(fit, updateObject=FALSE))
    expect_false(all(c("WRES","CPRED","CRES","CWRES") %in% names(fit)))
    expect_true(all(c("WRES","CPRED","CRES","CWRES") %in% names(fit2)))
    expect_equal(row.names(fit2$objDf), "FOCEi")

    fit <- .nlmixr(one.cmt, theo_sd, est="saem",
                   table=tableControl(cwres=TRUE))

    expect_true(all(c("WRES","CPRED","CRES","CWRES") %in% names(fit)))

  })

})
