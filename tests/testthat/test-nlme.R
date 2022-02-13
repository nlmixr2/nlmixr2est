.nlmixr <- function(...) suppressWarnings(suppressMessages(nlmixr(...)))

test_that("nlme one compartment theo_sd", {

  one.compartment <- function() {
    ini({
      tka <- 0.45 # Log Ka
      tcl <- 1 # Log Cl
      tv <- 3.45    # Log V
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      d/dt(depot) = -ka * depot
      d/dt(center) = ka * depot - cl / v * center
      cp = center / v
      cp ~ add(add.sd)
    })
  }

  nlme  <- .nlmixr(one.compartment, theo_sd, "nlme", control=nlmeControl(verbose=FALSE))

  expect_true(inherits(nlme, "nlmixr2FitData"))

  one.compartment <- function() {
    ini({
      tka <- 0.45 # Log Ka
      tcl <- 1 # Log Cl
      tv <- 3.45    # Log V
      eta.ka ~ 0.6
      eta.cl + eta.v ~ c(0.3,
                         0.001, 0.1)
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      d/dt(depot) = -ka * depot
      d/dt(center) = ka * depot - cl / v * center
      cp = center / v
      cp ~ add(add.sd)
    })
  }

  nlme  <- .nlmixr(one.compartment, theo_sd, "nlme", control=nlmeControl(maxIter=5, verbose=FALSE))

  expect_true(inherits(nlme, "nlmixr2FitData"))

  one.compartment <- function() {
    ini({
      tka <- exp(0.45) # Log Ka
      tcl <- 1 # Log Cl
      tv <- 3.45    # Log V
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- tka * exp(eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      d/dt(depot) = -ka * depot
      d/dt(center) = ka * depot - cl / v * center
      cp = center / v
      cp ~ add(add.sd)
    })
  }

  nlme  <- .nlmixr(one.compartment, theo_sd, "nlme", control=nlmeControl(maxIter=2, verbose=FALSE))

  expect_true(inherits(nlme, "nlmixr2FitData"))

})


test_that("Other error structures", {

  dat <- Wang2007
  dat$DV <- dat$Y

  mod <- function() {
    ini({
      tke <- 0.5
      eta.ke ~ 0.04
      prop.sd <- sqrt(0.1)
    })
    model({
      ke <- tke * exp(eta.ke)
      ipre <- 10 * exp(-ke * t)
      f2 <- ipre / (ipre + 5)
      f3 <- f2 * 3
      lipre <- log(ipre)
      ipre ~ prop(prop.sd)
    })
  }

  f <- mod()

  nlme  <- .nlmixr(f, dat, "nlme", control=nlmeControl(verbose=FALSE))

  expect_true(inherits(nlme, "nlmixr2FitData"))

  mod <- function() {
    ini({
      tke <- 0.5
      eta.ke ~ 0.04
      prop.sd <- sqrt(0.1)
      pw <- 4
    })
    model({
      ke <- tke * exp(eta.ke)
      ipre <- 10 * exp(-ke * t)
      f2 <- ipre / (ipre + 5)
      f3 <- f2 * 3
      lipre <- log(ipre)
      ipre ~ pow(prop.sd, pw)
    })
  }

  f <- mod()

  nlme  <- .nlmixr(mod, dat, "nlme", control=nlmeControl(verbose=FALSE))

  expect_true(inherits(nlme, "nlmixr2FitData"))

  mod <- function() {
    ini({
      tke <- 0.5
      eta.ke ~ 0.04
      add.sd <- sqrt(0.1)
      prop.sd <- sqrt(0.1)
    })
    model({
      ke <- tke * exp(eta.ke)
      ipre <- 10 * exp(-ke * t)
      f2 <- ipre / (ipre + 5)
      f3 <- f2 * 3
      lipre <- log(ipre)
      ipre ~ add(add.sd) + prop(prop.sd) + combined2()
    })
  }

  f <- mod()

  nlme  <- .nlmixr(mod, dat, "nlme", control=nlmeControl(msMaxIter=10000, verbose=FALSE))

  expect_true(inherits(nlme, "nlmixr2FitData"))

  mod <- function() {
    ini({
      tke <- 0.5
      eta.ke ~ 0.04
      add.sd <- sqrt(0.1)
      prop.sd <- sqrt(0.1)
    })
    model({
      ke <- tke * exp(eta.ke)
      ipre <- 10 * exp(-ke * t)
      f2 <- ipre / (ipre + 5)
      f3 <- f2 * 3
      lipre <- log(ipre)
      ipre ~ add(add.sd) + prop(prop.sd) + combined1()
    })
  }

  f <- mod()

  nlme  <- .nlmixr(mod, dat, "nlme", control=nlmeControl(msMaxIter=10000, verbose=FALSE))

  expect_true(inherits(nlme, "nlmixr2FitData"))

  mod <- function() {
    ini({
      tke <- 0.5
      eta.ke ~ 0.04
      add.sd <- sqrt(0.1)
      prop.sd <- sqrt(0.1)
      pw <- 1
    })
    model({
      ke <- tke * exp(eta.ke)
      ipre <- 10 * exp(-ke * t)
      f2 <- ipre / (ipre + 5)
      f3 <- f2 * 3
      lipre <- log(ipre)
      ipre ~ add(add.sd) + pow(prop.sd, pw) + combined1()
    })
  }

  f <- mod()

  nlme  <- .nlmixr(mod, dat, "nlme", control=nlmeControl(msMaxIter=10000, verbose=FALSE))

})

