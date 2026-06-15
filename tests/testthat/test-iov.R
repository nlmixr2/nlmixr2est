nmTest({

  one.cmt <- function() {
    ini({
      ## You may label each parameter with a comment
      tka <- 0.45 # Log Ka
      tcl <- log(c(0, 2.7, 100)) # Log Cl
      ## This works with interactive models
      ## You may also label the preceding line with label("label text")
      tv <- 3.45; label("log V")
      ## the label("Label name") works with all models
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      iov.cl ~ 0.1 | occ
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl + iov.cl)
      v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }

  theo_iov <- nlmixr2data::theo_md
  theo_iov$occ <- 1
  theo_iov$occ[theo_iov$TIME >= 144] <- 2

  fit <- .nlmixr(one.cmt, theo_iov, est="focei",
                 control=foceiControlFast)

  .df1 <- as.data.frame(fit[,c("ID", "occ", "iov.cl")])
  .df1 <- .df1[!duplicated(paste0(.df1$ID, ";", .df1$occ)), ]
  row.names(.df1) <- NULL
  expect_equal(.df1,
               fit$iov$occ)

  expect_true(inherits(fit, "nlmixr2FitCore"))

  expect_false(fit$iniDf[fit$iniDf$name == "iov.cl", "fix"])

  expect_null(fit$control$etaMat)

  fit2 <- .nlmixr(fit, est="focei", control=foceiControlFast)

  expect_true(inherits(fit2, "nlmixr2FitCore"))

  expect_false(fit2$iniDf[fit2$iniDf$name == "iov.cl", "fix"])

  expect_equal(dimnames(fit2$control$etaMat)[[2]],
               c("eta.ka", "eta.cl", "eta.v", "rx.iov.cl.1", "rx.iov.cl.2"))

  fit <- .nlmixr(one.cmt, theo_iov, est="foce",
                 control=foceControlFast)

  expect_true(inherits(fit, "nlmixr2FitCore"))
  expect_false(fit$iniDf[fit$iniDf$name == "iov.cl", "fix"])

  fit <- .nlmixr(one.cmt, theo_iov, est="saem",
                 control=saemControlFast)

  expect_true(inherits(fit, "nlmixr2FitCore"))
  expect_false(fit$iniDf[fit$iniDf$name == "iov.cl", "fix"])


  one.cmt <- function() {
    ini({
      ## You may label each parameter with a comment
      tka <- 0.45 # Log Ka
      tcl <- log(c(0, 2.7, 100)) # Log Cl
      ## This works with interactive models
      ## You may also label the preceding line with label("label text")
      tv <- 3.45; label("log V")
      ## the label("Label name") works with all models
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      iov.cl ~ fix(0.1) | occ
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl + iov.cl)
      v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }

  .foceiCtl <- foceiControlFast
  fit <- .nlmixr(one.cmt, theo_iov, est="focei",
                 control=foceiControlFast)
  expect_true(inherits(fit, "nlmixr2FitCore"))
  expect_true( fit$iniDf[fit$iniDf$name == "iov.cl", "fix"])
  expect_true(any(names(fit) == "iov.cl"))

  .foceiCtl$iovXform <- "var"

  fit <- .nlmixr(one.cmt, theo_iov, est="focei",
                 control=.foceiCtl)
  expect_true(inherits(fit, "nlmixr2FitCore"))
  expect_true( fit$iniDf[fit$iniDf$name == "iov.cl", "fix"])
  expect_true(any(names(fit) == "iov.cl"))


  .foceiCtl$iovXform <- "logvar"
  fit <- .nlmixr(one.cmt, theo_iov, est="focei",
                 control=.foceiCtl)
  expect_true(inherits(fit, "nlmixr2FitCore"))
  expect_true( fit$iniDf[fit$iniDf$name == "iov.cl", "fix"])
  expect_true(any(names(fit) == "iov.cl"))


  .foceiCtl$iovXform <- "logsd"
  fit <- .nlmixr(one.cmt, theo_iov, est="focei",
                 control=.foceiCtl)
  expect_true(inherits(fit, "nlmixr2FitCore"))
  expect_true( fit$iniDf[fit$iniDf$name == "iov.cl", "fix"])
  expect_true(any(names(fit) == "iov.cl"))

  fit <- .nlmixr(one.cmt, theo_iov, est="foce",
                 control=foceControlFast)
  expect_true(inherits(fit, "nlmixr2FitCore"))
  expect_true(fit$iniDf[fit$iniDf$name == "iov.cl", "fix"])
  expect_true(any(names(fit) == "iov.cl"))


  fit <- .nlmixr(one.cmt, theo_iov, est="saem",
                 control=saemControlFast)
  expect_true(inherits(fit, "nlmixr2FitCore"))
  expect_true(fit$iniDf[fit$iniDf$name == "iov.cl", "fix"])
  expect_true(any(names(fit) == "iov.cl"))

  one.cmt <- function() {
    ini({
      ## You may label each parameter with a comment
      tka <- 0.45 # Log Ka
      tcl <- log(c(0, 2.7, 100)) # Log Cl
      ## This works with interactive models
      ## You may also label the preceding line with label("label text")
      tv <- 3.45; label("log V")
      ## the label("Label name") works with all models
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      iov.cl ~ fix(0.1) | whoa
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl + iov.cl)
      v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }

  expect_error(suppressWarnings(.nlmixr(one.cmt, theo_iov, est="focei",
                                        control=foceiControlFast)), "whoa")

  # Test that stale .uiIovEnv state is cleared when a non-IOV fit follows
  # an IOV fit in the same R session.
  one.cmt.no.iov <- function() {
    ini({
      tka <- 0.45
      tcl <- log(c(0, 2.7, 100))
      tv <- 3.45; label("log V")
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

  # First run an IOV fit to prime the .uiIovEnv state
  one.cmt.iov <- function() {
    ini({
      tka <- 0.45
      tcl <- log(c(0, 2.7, 100))
      tv <- 3.45; label("log V")
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      iov.cl ~ 0.1 | occ
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl + iov.cl)
      v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }

  fit_iov_stale <- .nlmixr(one.cmt.iov, theo_iov, est="focei",
                             control=foceiControlFast)
  expect_true(inherits(fit_iov_stale, "nlmixr2FitCore"))

  # Now run a non-IOV fit; stale .uiIovEnv state must not contaminate it
  fit_no_iov <- .nlmixr(one.cmt.no.iov, theo_iov, est="focei",
                         control=foceiControlFast)
  expect_true(inherits(fit_no_iov, "nlmixr2FitCore"))
  # No IOV table should be attached
  expect_null(fit_no_iov$iov)
  # No internal rx. IOV columns should leak into the output data frame
  expect_false(any(grepl("^rx[.]", names(fit_no_iov))))

  # Test IOV on two different conditioning variables (occ and occ2).
  # This exercises the fix for .uiIovEnv$lines saving ALL injected
  # expressions (not only the first conditioning variable's lines).
  one.cmt.two.cond <- function() {
    ini({
      tka <- 0.45
      tcl <- log(c(0, 2.7, 100))
      tv <- 3.45; label("log V")
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      iov.cl ~ 0.1 | occ
      iov.ka ~ 0.1 | occ2
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka + iov.ka)
      cl <- exp(tcl + eta.cl + iov.cl)
      v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }

  theo_iov2 <- theo_iov
  theo_iov2$occ2 <- 1L
  theo_iov2$occ2[theo_iov2$TIME >= 168] <- 2L

  fit_two_cond <- .nlmixr(one.cmt.two.cond, theo_iov2, est="focei",
                           control=foceiControlFast)
  expect_true(inherits(fit_two_cond, "nlmixr2FitCore"))
  # Both IOV parameters should be present (not fixed unless specified)
  expect_false(fit_two_cond$iniDf[fit_two_cond$iniDf$name == "iov.cl", "fix"])
  expect_false(fit_two_cond$iniDf[fit_two_cond$iniDf$name == "iov.ka", "fix"])
  # IOV tables for both conditioning variables should be attached
  expect_true(!is.null(fit_two_cond$iov$occ))
  expect_true(!is.null(fit_two_cond$iov$occ2))
  # No internal rx. injection lines should remain in the output columns
  expect_false(any(grepl("^rx[.]", names(fit_two_cond))))

})
