test_that("nmObject get tests", {

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
  skip_if_not(rxode2parse::.linCmtSens())

  fit <- suppressMessages(suppressWarnings(nlmixr(one.cmt, nlmixr2data::theo_sd, est="focei")))

  expect_equal(fit$modelName, "one.cmt")


  one.cmt <- function() {
    ini({
      ## You may label each parameter with a comment
      tka <- 0.45 # Log Ka
      tcl <- log(c(0, 2.7, 100)) # Log Cl
      ## This works with interactive models
      ## You may also label the preceding line with label("label text")
      tv <- 3.45; label("log V")
      ## the label("Label name") works with all models
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl)
      v <- exp(tv)
      linCmt() ~ add(add.sd)
    })
  }

  fit2 <- suppressMessages(suppressWarnings(nlmixr(one.cmt, nlmixr2data::theo_sd, est="focei")))

  expect_equal(fit2$modelName, "one.cmt")

})
