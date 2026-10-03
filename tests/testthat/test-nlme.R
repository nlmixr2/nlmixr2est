nmTest({
  test_that("nlme will pick up interpolation", {
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
        cl <- exp(tcl + eta.cl) + wt
        v <- exp(tv + eta.v)
        d/dt(depot) = -ka * depot
        d/dt(center) = ka * depot - cl / v * center
        cp = center / v
        cp ~ add(add.sd)
      })
    }

    f <- one.compartment()

    expect_false(grepl("linear\\(wt\\)", suppressMessages(rxode2::rxNorm(f$nlmRxModel$predOnly))))

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
        linear(wt)
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl) + wt
        v <- exp(tv + eta.v)
        d/dt(depot) = -ka * depot
        d/dt(center) = ka * depot - cl / v * center
        cp = center / v
        cp ~ add(add.sd)
      })
    }

    f <- one.compartment()

    expect_true(grepl("linear\\(wt\\)", suppressMessages(rxode2::rxNorm(f$nlmRxModel$predOnly))))
  })

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

    nlme <- .nlmixr(one.compartment, theo_sd, "nlme", control = nlmeControl(verbose = FALSE, returnObject = TRUE))

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

    nlme <- .nlmixr(
      one.compartment,
      theo_sd,
      "nlme",
      control = nlmeControl(maxIter = 5, verbose = FALSE, returnObject = TRUE)
    )

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

    nlme <- .nlmixr(
      one.compartment,
      theo_sd,
      "nlme",
      control = nlmeControl(maxIter = 2, verbose = FALSE, returnObject = TRUE)
    )

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

    nlme <- .nlmixr(f, dat, "nlme", control = nlmeControl(verbose = FALSE, returnObject = TRUE))

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

    nlme <- .nlmixr(mod, dat, "nlme", control = nlmeControl(verbose = FALSE, returnObject = TRUE))

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

    nlme <- .nlmixr(mod, dat, "nlme", control = nlmeControl(msMaxIter = 10000, verbose = FALSE, returnObject = TRUE))

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

    nlme <- .nlmixr(mod, dat, "nlme", control = nlmeControl(msMaxIter = 10000, verbose = FALSE, returnObject = TRUE))

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

    nlme <- .nlmixr(mod, dat, "nlme", control = nlmeControl(msMaxIter = 10000, verbose = FALSE, returnObject = TRUE))
  })

  test_that("nlme random effects are returned", {
    one.cmt.all.mu.ref <- function() {
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
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }

    fit_nlme <- .nlmixr(
      one.cmt.all.mu.ref,
      theo_sd,
      est = "nlme",
      control = nlmeControl(verbose = FALSE, returnObject = TRUE)
    )
    expect_true(!all(ranef(fit_nlme)[[1]] == 0))

    one.cmt.one.mu.ref <- function() {
      ini({
        ## You may label each parameter with a comment
        tka <- 0.45 # Log Ka
        tcl <- log(c(0, 2.7, 100)) # Log Cl
        ## This works with interactive models
        ## You may also label the preceding line with label("label text")
        tv <- 3.45; label("log V")
        ## the label("Label name") works with all models
        eta.cl ~ 0.3
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv)
        linCmt() ~ add(add.sd)
      })
    }

    fit_nlme <- .nlmixr(
      one.cmt.one.mu.ref,
      theo_sd,
      est = "nlme",
      control = nlmeControl(verbose = FALSE, returnObject = TRUE)
    )
    expect_true(!all(ranef(fit_nlme)[[1]] == 0))

    one.cmt.non.mu.ref <- function() {
      ini({
        ## You may label each parameter with a comment
        tka <- 0.45 # Log Ka
        tcl <- c(0, 2.7, 100) # Log Cl
        ## This works with interactive models
        ## You may also label the preceding line with label("label text")
        tv <- 3.45; label("log V")
        ## the label("Label name") works with all models
        eta.cl ~ 0.3
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka)
        cl <- tcl*exp(eta.cl)
        v <- exp(tv)
        linCmt() ~ add(add.sd)
      })
    }

    fit_nlme <- .nlmixr(
      one.cmt.non.mu.ref,
      theo_sd,
      est = "nlme",
      control = nlmeControl(verbose = FALSE, returnObject = TRUE)
    )
    expect_true(!all(ranef(fit_nlme)[[1]] == 0))
  })

  test_that("a block omega reports nlme's estimate of the declared covariances (issue 1140)", {
    one.compartment <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        eta.ka ~ 0.6
        eta.cl + eta.v ~ c(0.3, 0.001, 0.1)
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }
    fit <- .nlmixr(one.compartment, theo_sd, "nlme", control = nlmeControl(verbose = FALSE, returnObject = TRUE))
    # nlme fits the declared structure: eta.ka alone, eta.cl and eta.v as a block
    .re <- fit$nlme$modelStruct$reStruct[[1]]
    expect_s3_class(.re, "pdBlocked")
    .eta <- c("eta.ka", "eta.cl", "eta.v")
    .est <- nlme::pdMatrix(.re)[.eta, .eta] * fit$nlme$sigma^2
    expect_equal(fit$omega, .est)
    expect_equal(fit$omega["eta.ka", c("eta.cl", "eta.v")], c(eta.cl = 0, eta.v = 0))
    # the estimate, not the ini() block
    expect_true(abs(fit$omega["eta.cl", "eta.v"] - 0.001) > 1e-3)
    expect_equal(diag(fit$omega), diag(.est))
    # refits, setCov() and setOfv() start from the fit's ui
    expect_equal(fit$ui$omega, fit$omega)
    .vc <- VarCorr(fit)
    expect_equal(as.numeric(.vc[.eta, "Variance"]), unname(diag(fit$omega)), tolerance = 1e-6)
    # ranef columns follow the ui's eta order
    expect_equal(names(fit$eta), c("ID", .eta))
    .re <- nlme::ranef(fit$nlme)
    .re <- .re[order(as.numeric(rownames(.re))), .eta]
    expect_equal(unname(as.matrix(fit$eta[, .eta])), unname(as.matrix(.re)))
  })

  test_that("the default nlme covariance is nlme's full fixed-effect covariance (issue 1140)", {
    one.compartment <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }
    fit <- .nlmixr(one.compartment, theo_sd, "nlme", control = nlmeControl(verbose = FALSE, returnObject = TRUE))
    expect_identical(fit$covMethod, "nlme")
    .th <- c("tka", "tcl", "tv")
    # the SEs are the ones summary() prints (ML: sigma adjusted to REML-like)
    .se <- summary(fit$nlme)$tTable[.th, "Std.Error"]
    expect_equal(sqrt(diag(fit$cov)), .se)
    expect_equal(fit$parFixedDf[.th, "SE"], .se)
    # and the correlations are nlme's; they used to be dropped
    expect_equal(stats::cov2cor(fit$cov), stats::cov2cor(vcov(fit$nlme)))
    expect_true(all(fit$cov[upper.tri(fit$cov)] != 0))
    # ML: vcov() is the same matrix before nlme's sigma adjustment
    .dims <- fit$nlme$dims
    expect_equal(fit$cov, vcov(fit$nlme) * .dims$N / (.dims$N - length(.th)))
  })

  test_that(".nlmeGetOmega returns nlme's matrix in the ui's eta order (issue 1140)", {
    # 4 etas: VarCorr()'s printed correlations used to be copied into the
    # wrong cells from the 4th eta on
    .eta <- c("eta.a", "eta.b", "eta.c", "eta.d")
    .cor <- matrix(
      c(
        1,
        0.1,
        0.2,
        0.3,
        0.1,
        1,
        0.4,
        0.5,
        0.2,
        0.4,
        1,
        0.6,
        0.3,
        0.5,
        0.6,
        1
      ),
      4,
      4
    )
    .sd <- c(0.5, 0.4, 0.3, 0.2)
    .m <- diag(.sd) %*% .cor %*% diag(.sd)
    dimnames(.m) <- list(.eta, .eta)
    # nlme holds the matrix relative to sigma^2, in its own (block) order
    .ord <- c("eta.b", "eta.a", "eta.d", "eta.c")
    .pd <- nlme::pdSymm(value = .m[.ord, .ord] / 0.25, form = eta.b + eta.a + eta.d + eta.c ~ 1)
    .fake <- list(modelStruct = list(reStruct = nlme::reStruct(list(ID = .pd))), sigma = 0.5)
    .ui <- list(omega = .m, muRefDataFrame = data.frame(theta = character(0), eta = character(0)))
    expect_equal(.nlmeGetOmega(.fake, .ui), .m)
  })
})
