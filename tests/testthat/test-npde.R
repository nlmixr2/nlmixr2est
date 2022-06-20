nmTest({

  test_that("npde", {

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

    expect_false(all(c("EPRED","ERES","NPDE","NPD", "PDE", "PD") %in% names(fit)))
    suppressMessages(expect_error(addNpde(fit), NA))
    expect_true(all(c("EPRED","ERES","NPDE","NPD", "PDE", "PD") %in% names(fit)))
    .range1 <- range(fit$PDE)
    .range2 <- range(fit$PD)
    expect_true(.range1[1] > 0)
    expect_true(.range2[1] > 0)
    expect_true(.range1[2] < 1)
    expect_true(.range2[2] < 1)

    .range1 <- range(fit$NPDE)
    .range2 <- range(fit$NPD)
    .range3 <- range(fit$ERES)

    expect_true(.range1[1] < -2)
    expect_true(.range2[1] < -2)
    expect_true(.range3[1] < -2)
    expect_true(.range1[2] > 2)
    expect_true(.range2[2] > 2)
    expect_true(.range3[2] > 2)

    .range4 <- range(fit$EPRED)
    expect_true(.range4[1] > -0.1)
    expect_true(.range4[2] > 7)



    fit <- .nlmixr(one.cmt, theo_sd, est="saem")

    expect_false(all(c("EPRED","ERES","NPDE","NPD","PDE","PD") %in% names(fit)))
    fit2 <- suppressMessages(addNpde(fit, updateObject=FALSE))
    expect_false(all(c("EPRED","ERES","NPDE","NPD","PDE","PD") %in% names(fit)))
    expect_true(all(c("EPRED","ERES","NPDE","NPD","PDE","PD") %in% names(fit2)))

    fit <- .nlmixr(one.cmt, theo_sd, est="saem",
                   table=tableControl(npde=TRUE))

    expect_true(all(c("EPRED","ERES","NPDE","NPD", "PDE","PD") %in% names(fit)))

  })

})
