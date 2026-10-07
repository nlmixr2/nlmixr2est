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
    # the off-diagonal from nlme's correlation and standard deviations, not from the
    # covariance matrix the fit is built from
    .cor <- nlme::corMatrix(fit$nlme$modelStruct$reStruct[[1]])
    .sd <- attr(.cor, "stdDev") * fit$nlme$sigma
    expect_equal(
      fit$omega["eta.cl", "eta.v"],
      .cor["eta.cl", "eta.v"] * .sd[["eta.cl"]] * .sd[["eta.v"]],
      # rebuilding a covariance from a correlation and two sds rounds at ~1e-8
      tolerance = 1e-6
    )
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
    # and the correlations are nlme's, not zero
    expect_equal(stats::cov2cor(fit$cov), stats::cov2cor(vcov(fit$nlme)))
    expect_true(all(fit$cov[upper.tri(fit$cov)] != 0))
    # ML: vcov() is the same matrix before nlme's sigma adjustment
    .dims <- fit$nlme$dims
    expect_equal(fit$cov, vcov(fit$nlme) * .dims$N / (.dims$N - length(.th)))
  })

  test_that("the residual parameters reproduce nlme's residual sd for every error model (issue 1140)", {
    # both error components matter: additive 0.3, proportional 0.1
    dat <- rxode2::rxWithSeed(1140, {
      .t <- c(0.5, 1, 2, 4, 6, 8)
      .d <- expand.grid(TIME = .t, ID = 1:20)
      .ke <- 0.3 * exp(rnorm(20, 0, 0.2))
      .ipre <- 10 * exp(-.ke[.d$ID] * .d$TIME)
      .d$DV <- .ipre + rnorm(nrow(.d), 0, sqrt(0.3^2 + (0.1 * .ipre)^2))
      .d[, c("ID", "TIME", "DV")]
    })
    base <- function() {
      ini({
        tke <- 0.3
        eta.ke ~ 0.04
        add.sd <- 0.3
        prop.sd <- 0.1
        pw <- 1
      })
      model({
        ke <- tke * exp(eta.ke)
        ipre <- 10 * exp(-ke * t)
        ipre ~ add(add.sd) + pow(prop.sd, pw)
      })
    }
    .cases <- list(
      add = list(quote(ipre ~ add(add.sd)), "combined2"),
      prop = list(quote(ipre ~ prop(prop.sd)), "combined2"),
      pow = list(quote(ipre ~ pow(prop.sd, pw)), "combined2"),
      combined1 = list(quote(ipre ~ add(add.sd) + prop(prop.sd)), "combined1"),
      combined2 = list(quote(ipre ~ add(add.sd) + prop(prop.sd)), "combined2"),
      combined1pow = list(quote(ipre ~ add(add.sd) + pow(prop.sd, pw)), "combined1")
    )
    for (.n in names(.cases)) {
      .m <- suppressMessages(eval(bquote(rxode2::model(base, .(.cases[[.n]][[1]])))))
      fit <- .nlmixr(
        .m,
        dat,
        "nlme",
        control = nlmeControl(verbose = FALSE, returnObject = TRUE, addProp = .cases[[.n]][[2]])
      )
      .nl <- fit$nlme
      .vs <- .nl$modelStruct$varStruct
      .th <- fit$theta
      # nlme's own residual sd of every observation (in its own row order)
      .sdNlme <- if (is.null(.vs)) .nl$sigma else unname(.nl$sigma / nlme::varWeights(.vs))
      .f <- if (is.null(.vs)) NULL else abs(attr(.vs, "covariate"))
      .sd <- unname(switch(
        .n,
        add = .th[["add.sd"]],
        prop = .th[["prop.sd"]] * .f,
        pow = .th[["prop.sd"]] * .f^.th[["pw"]],
        combined1 = .th[["add.sd"]] + .th[["prop.sd"]] * .f,
        combined2 = sqrt(.th[["add.sd"]]^2 + (.th[["prop.sd"]] * .f)^2),
        combined1pow = .th[["add.sd"]] + .th[["prop.sd"]] * .f^.th[["pw"]]
      ))
      expect_equal(.sd, .sdNlme, tolerance = 1e-6, info = .n)
      if (.n == "combined2") {
        # sigma is fixed at 1, so const and prop are the sds themselves
        expect_equal(.nl$sigma, 1)
      }
    }
  })

  test_that(".nlmeGetOmega returns nlme's matrix in the ui's eta order (issue 1140)", {
    # 4 etas: the correlations of every eta pair, the 4th eta's included, land
    # in their own cells
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

  test_that("nlme fits each declared omega block as its own pdSymm", {
    # eta.d and eta.f have no covariance of their own (a 0 row), but both share
    # one with eta.e, so the three are one block; eta.a/eta.b are a second block
    blockMod <- function() {
      ini({
        t1 <- 1
        t2 <- 1
        t3 <- 1
        t4 <- 1
        t5 <- 1
        t6 <- 1
        eta.a + eta.b ~ c(1, 0.2, 1)
        eta.c ~ 0.5
        eta.d + eta.e + eta.f ~ c(1, 0.1, 1, 0, 0.1, 1)
        add.sd <- 1
      })
      model({
        y <- t1 * exp(eta.a) + t2 * exp(eta.b) + t3 * exp(eta.c) + t4 * exp(eta.d) + t5 * exp(eta.e) +
          t6 * exp(eta.f)
        y ~ add(add.sd)
      })
    }
    .ui <- suppressWarnings(rxode2::rxode2(blockMod))
    .idf <- .ui$iniDf
    .etaName <- .idf$name[!is.na(.idf$neta1) & .idf$neta1 == .idf$neta2]
    .etaName <- .etaName[order(.idf$neta1[!is.na(.idf$neta1) & .idf$neta1 == .idf$neta2])]
    .blocks <- lapply(.nlmeOmegaBlocks(.idf), function(i) sort(.etaName[i]))
    expect_setequal(
      vapply(.blocks, paste, character(1), collapse = ","),
      c("eta.a,eta.b", "eta.c", "eta.d,eta.e,eta.f")
    )
    .pd <- .ui$nlmePdOmega
    expect_s3_class(.pd, "pdBlocked")
    .cls <- vapply(.pd, function(b) class(b)[1], character(1))
    .size <- vapply(.pd, function(b) nrow(as.matrix(b)), integer(1))
    expect_identical(sort(paste(.cls, .size)), c("pdDiag 1", "pdSymm 2", "pdSymm 3"))
  })

  test_that("one full omega block stays pdSymm and a diagonal omega pdDiag", {
    .mk <- function(om) {
      .f <- function() {
        ini({
          t1 <- 1
          t2 <- 1
          t3 <- 1
          add.sd <- 1
        })
        model({
          y <- t1 * exp(eta.a) + t2 * exp(eta.b) + t3 * exp(eta.c)
          y ~ add(add.sd)
        })
      }
      suppressWarnings(rxode2::rxode2(eval(bquote(rxode2::ini(.f, .(om))))))
    }
    .full <- .mk(quote(eta.a + eta.b + eta.c ~ c(1, 0.1, 1, 0.1, 0.1, 1)))
    expect_s3_class(.full$nlmePdOmega, "pdSymm")
    expect_false(inherits(.full$nlmePdOmega, "pdBlocked"))
    .diag <- .mk(quote({
      eta.a ~ 1
      eta.b ~ 1
      eta.c ~ 1
    }))
    expect_s3_class(.diag$nlmePdOmega, "pdDiag")
  })

  test_that("nlme refuses a fixed omega element", {
    .f <- function() {
      ini({
        tke <- 0.5
        eta.ke ~ fix(0.04)
        add.sd <- 0.1
      })
      model({
        ke <- tke * exp(eta.ke)
        ipre <- 10 * exp(-ke * t)
        ipre ~ add(add.sd)
      })
    }
    .d <- Wang2007
    .d$DV <- .d$Y
    expect_error(
      .nlmixr(.f, .d, "nlme", control = nlmeControl(verbose = FALSE)),
      "fixed omega elements are not supported"
    )
  })

  test_that("the nlme covariance of a single fixed effect is named and nlme's", {
    oneTheta <- function() {
      ini({
        tv <- 3.45
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- 1.5
        cl <- 2.7
        v <- exp(tv + eta.v)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }
    fit <- .nlmixr(oneTheta, theo_sd, "nlme", control = nlmeControl(verbose = FALSE, returnObject = TRUE))
    expect_identical(dimnames(fit$cov), list("tv", "tv"))
    expect_equal(fit$cov[1, 1], summary(fit$nlme)$tTable["tv", "Std.Error"]^2)
    expect_equal(unname(fit$parFixedDf["tv", "SE"]), summary(fit$nlme)$tTable["tv", "Std.Error"])
  })
})
