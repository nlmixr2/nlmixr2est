test_that("test focei llik", {
  # dnorm() works

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
      linCmt() ~ add(add.sd) + dnorm()
    })
  }

  f <- nlmixr(one.cmt, theo_sd, "focei")
  expect_true("CWRES" %in% names(f))
  
  f <- nlmixr(one.cmt, theo_sd, "foce")
  expect_true("CWRES" %in% names(f))

  expect_error(nlmixr(one.cmt, theo_sd, "fo"))

  one.cmt.ll <- function() {
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
      cp <- linCmt()
      ll <- -log(add.sd*sqrt(2*pi))-0.5*((DV-cp)/add.sd)^2
      n2ll(err) ~ -2*ll
    })
  }

  f <- try(nlmixr(one.cmt.ll, theo_sd, "focei"))
  
  expect_false(inherits(f, "try-error"))
  expect_equal(-2*f$ll, f$IPRED)
  expect_false("CWRES" %in% names(f))

  f <- try(nlmixr(one.cmt.ll, theo_sd, "foce"))
  
  expect_false(inherits(f, "try-error"))
  expect_equal(-2*f$ll, f$IPRED)
  expect_false("CWRES" %in% names(f))
  
  # no etas test

  one.cmt.noeta <- function() {
    ini({
      ## You may label each parameter with a comment
      tka <- 0.45 # Ka
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
      linCmt() ~ add(add.sd) + dnorm()
    })
  }

  f <- nlmixr(one.cmt.noeta, theo_sd, "focei")


  one.cmt.ll.noeta <- function() {
    ini({
      ## You may label each parameter with a comment
      tka <- 0.45 # Ka
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
      cp <- linCmt()
      ll <- -log(add.sd*sqrt(2*pi))-0.5*((DV-cp)/add.sd)^2
      n2ll(err) ~ -2*ll
    })
  }

  f <- nlmixr(one.cmt.ll.noeta, theo_sd, "focei")
  
  
})
