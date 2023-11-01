nmTest({

  dat <- Wang2007
  dat$DV <- dat$Y

  gg <- function(x, y) {
    x * y
  }

  f <- function() {
    ini({
      tke <- 0.5
      eta.ke ~ 0.04
      prop.sd <- sqrt(0.1)
    })
    model({
      ke <- gg(tke, exp(eta.ke))
      ipre <- gg(10, exp(-ke * t))
      lipre <- log(ipre)
      ipre ~ prop(prop.sd)
    })
  }

  .env <- new.env(parent=emptyenv())
  expect_error({.env$f <- nlmixr2(f, dat, "focei")}, NA)

  expect_equal(rxModelVars(.env$f$foceiModel$inner)$udf, c("gg"=2L))

  expect_equal(rxModelVars(.env$f$foceiModel$predOnly)$udf, c("gg"=2L))

  expect_equal(rxModelVars(.env$f$foceiModel$predNoLhs)$udf, c("gg"=2L))

  expect_error({.env$s <- nlmixr2(f, dat, "saem")}, NA)

  expect_error({.env$n <- nlmixr2(f, dat, "nlme")}, NA)

  g <- function() {
    ini({
      tke <- 0.5
      add.sd <- sqrt(0.1)
    })
    model({
      ke <- tke
      ipre <- gg(10, exp(-ke * t))
      lipre <- log(ipre)
      ipre ~ add(add.sd)
    })
  }

  expect_error({.env$nlm <- nlmixr2(g, dat, "nlm")}, NA)

  expect_error({.env$optim <- nlmixr2(g, dat, "optim")}, NA)

  #expect_error({.env$nls <- nlmixr2(g, dat, "nls")}, NA)

  expect_error({.env$nlminb <- nlmixr2(g, dat, "nlminb")}, NA)

  expect_error({.env$bobyqa <- nlmixr2(g, dat, "bobyqa")}, NA)

  expect_error({.env$lbfgsb3c <- nlmixr2(g, dat, "lbfgsb3c")}, NA)

  expect_error({.env$n1qn1 <- nlmixr2(g, dat, "n1qn1")}, NA)

  rxode2::rxFun(gg)

  rm(gg)

  rxClean()

  f <- nlmixr2(f, dat, "focei")


})
