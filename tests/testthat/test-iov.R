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

  theo_iov <- theo_md
  theo_iov$occ <- 1
  theo_iov$occ[theo_iov$TIME >= 144] <- 2

  fit <- .nlmixr(one.cmt, theo_iov, est="focei",
                 control=foceiControlFast)
  expect_true(inherits(fit, "nlmixr2FitCore"))
  expect_false( fit$iniDf[fit$iniDf$name == "iov.cl", "fix"])

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

  fit <- .nlmixr(one.cmt, theo_iov, est="focei",
                 control=foceiControlFast)
  expect_true(inherits(fit, "nlmixr2FitCore"))
  expect_true( fit$iniDf[fit$iniDf$name == "iov.cl", "fix"])

  fit <- .nlmixr(one.cmt, theo_iov, est="foce",
                 control=foceControlFast)
  expect_true(inherits(fit, "nlmixr2FitCore"))
  expect_true(fit$iniDf[fit$iniDf$name == "iov.cl", "fix"])

  fit <- .nlmixr(one.cmt, theo_iov, est="saem",
                 control=saemControlFast)
  expect_true(inherits(fit, "nlmixr2FitCore"))
  expect_true(fit$iniDf[fit$iniDf$name == "iov.cl", "fix"])

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

})
