
nmTest({

  test_that("nlm makes sense", {

    dsn <- data.frame(i=1:1000)
    dsn$time <- exp(rnorm(1000))
    dsn$DV <- rbinom(1000,1,exp(-1+dsn$time)/(1+exp(-1+dsn$time)))

    mod <- function() {
      ini({
        E0 <- 0.5
        Em <- 0.5
        E50 <- 2
        g <- fix(2)
      })
      model({
        v <- E0+Em*time^g/(E50^g+time^g)
        p <- expit(v)
        ll(bin) ~ DV * v - log(1 + exp(v))
      })
    }

    fit2 <- nlmixr(mod, dsn, est="nlm")

    expect_true(inherits(fit2, "nlmixr2.nlm"))

    fit3 <- fit2 %>% ini(g=unfix) %>% nlmixr2(dsn, "nlm", nlmControl(solveType="grad"))

    expect_true(inherits(fit3, "nlmixr2.nlm"))

    fit4 <- fit2 %>% ini(g=unfix) %>% nlmixr2(dsn, "nlm", nlmControl(solveType="fun"))

    expect_true(inherits(fit4, "nlmixr2.nlm"))

    one.cmt <- function() {
      ini({
        tka <- 0.45
        tcl <- log(c(0, 2.7, 100))
        tv <- 3.45
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl)
        v <- exp(tv)
        linCmt() ~ add(add.sd)
      })
    }

    fit1 <- nlmixr(one.cmt, nlmixr2data::theo_sd, est="nlm")

    expect_true(inherits(fit1, "nlmixr2.nlm"))

  })

  test_that("nlm c++ level functions", {

    mod <- function() {
      ini({
        E0 <- 0.5
        Em <- 0.5
        E50 <- 2
        g <- fix(2)
      })
      model({
        v <- E0+Em*time^g/(E50^g+time^g)
        ll(bin) ~ DV * v - log(1 + exp(v))
      })
    }

    mod <- mod()

    dsn <- data.frame(i=1:1000)
    dsn$time <- exp(rnorm(1000))
    dsn$DV <- rbinom(1000,1,exp(-1+dsn$time)/(1+exp(-1+dsn$time)))

    f <- mod$nlmSensModel

    env <- new.env(parent=emptyenv())

    env$rxControl <- rxode2::rxControl()
    env$predOnly <- f$predOnly
    env$thetaGrad <- f$thetaGrad
    p <- mod$nlmParIni
    env$param <- setNames(p, sprintf("THETA[%d]", seq_along(p)))
    env$data <- dsn
    env$needFD <- f$eventTheta

    #.eventTypeIdx <- c("central" = 2L, "forward" = 3L)
    control <- list(stickyRecalcN=4,
                    maxOdeRecalc=5,
                    odeRecalcFactor=10^(0.5),

                    eventType=2L,
                    shi21maxFD=20L,
                    shiErr=(.Machine$double.eps)^(1/3),

                    optimHessType=2L,
                    shi21maxHess=20L,
                    hessErr=(.Machine$double.eps)^(1/3),

                    solveType=3L)

    env$control <- control

    .Call(`_nlmixr2est_nlmSetup`, env)

    ## .Call(`_nlmixr2est_optimFunC`, env$param, FALSE)

    ## .Call(`_nlmixr2est_optimFunC`, env$param, TRUE)

    ## .Call(`_nlmixr2est_optimFunC`, env$param * 0.5, FALSE)

    ## .Call(`_nlmixr2est_optimFunC`, env$param * 0.5, TRUE)

    ## .Call(`_nlmixr2est_nlminbFunC`, env$param, 1L)

    ## .Call(`_nlmixr2est_nlminbFunC`, env$param, 2L)

    ## .Call(`_nlmixr2est_nlminbFunC`, env$param, 3L)

    val <- .Call(`_nlmixr2est_nlmSolveR`, env$param)


    val1 <- val

    val2 <- .Call(`_nlmixr2est_nlmSolveGradR`, env$param)
    attr(val1, "gradient") <- attr(val2, "gradient")

    expect_equal(val1, val2)




    R <- deriv(~-(DV * (E0+Em*time^g/(E50^g+time^g)) - log(1 + exp(E0+Em*time^g/(E50^g+time^g)))),c("E0", "Em", "E50"))
    R <- R[[1]]
    R <- lapply(seq_along(R)[-1], function(i) R[[i]])

    g <- function(p, time, DV) {}
    body(g) <- as.call(c(list(quote(`{`),
                              str2lang("E0<-p[1]"),
                              str2lang("Em <- p[2]"),
                              str2lang("E50 <- p[3]"),
                              str2lang("g <- 2")),
                         R,
                         list(str2lang(".grad <- colSums(.grad)"),
                              str2lang(".value <- sum(.value)"),
                              str2lang("names(.grad) <- NULL"),
                              str2lang("attr(.value, 'gradient') <- .grad"),
                              str2lang(".value"))))

    expect_equal(.Call(`_nlmixr2est_nlmSolveGradR`, env$param),
                 .Call(`_nlmixr2est_nlmSolveGradR`, env$param))

    expect_equal(.Call(`_nlmixr2est_nlmSolveGradR`, env$param),
                 g(mod$nlmParIni, time=dsn$time, DV=dsn$DV))


    expect_equal(.Call(`_nlmixr2est_nlmSolveGradOnly`, env$param),
                 attr(val1, "gradient"))

    expect_equal(.Call(`_nlmixr2est_nlmSolveGradHess`, env$param),
                 .Call(`_nlmixr2est_nlmSolveGradHess`, env$param))

    R <- deriv(~-(DV * (E0+Em*time^g/(E50^g+time^g)) - log(1 + exp(E0+Em*time^g/(E50^g+time^g)))),c("E0", "Em", "E50"), hessian=TRUE)

    R <- R[[1]]

    R <- lapply(seq_along(R)[-1], function(i) R[[i]])

    g2 <- function(p, time, DV) {}
    body(g2) <- as.call(c(list(quote(`{`),
                              str2lang("E0<-p[1]"),
                              str2lang("Em <- p[2]"),
                              str2lang("E50 <- p[3]"),
                              str2lang("g <- 2")),
                         R))

    i <- seq_along(dsn$time)

    envH <- new.env(parent=emptyenv())
    envH$hess <- NULL

    hess <- lapply(seq_along(dsn$time), function(i) {
      .ret <- attr(g2(p, dsn$time[i], dsn$DV[i]), "hessian")
      .ret0 <- matrix(rep(0.0, 3 * 3), 3, 3,
                      dimnames=list(c("E0", "Em", "E50"),
                                    c("E0", "Em", "E50")))
      for (i in c("E0", "Em", "E50")) {
        for (j in c("E0", "Em", "E50")) {
          .ret0[i, j] <- .ret[,, i][j]
        }
      }
      if (is.null(envH$hess)) {
        envH$hess <- .ret0
      } else {
        envH$hess <- envH$hess + .ret0
      }
      .ret0
    })

    val3 <- .Call(`_nlmixr2est_nlmSolveGradHess`, env$param)

    attr(val2, "hessian") <- attr(val3, "hessian")

    expect_equal(val2, val3)

    expect_equal(.Call(`_nlmixr2est_nlmSolveSwitch`, env$param), val2)

    env$control$solveType <- 2L

    .Call(`_nlmixr2est_nlmSetup`, env)

    expect_equal(.Call(`_nlmixr2est_nlmSolveSwitch`, env$param),
                 .Call(`_nlmixr2est_nlmSolveGradR`, env$param))

    env$control$solveType <- 1L
    .Call(`_nlmixr2est_nlmSetup`, env)

    expect_equal(.Call(`_nlmixr2est_nlmSolveSwitch`, env$param),
                 .Call(`_nlmixr2est_nlmSolveR`, env$param))

    # now test forward differences
    env$control$eventType <- 1L
    env$control$optimHessType <- 1L

    .Call(`_nlmixr2est_nlmSetup`, env)

    expect_equal(.Call(`_nlmixr2est_nlmSolveGradHess`, env$param),
                 .Call(`_nlmixr2est_nlmSolveGradHess`, env$param))

    # now get a solved system with a finite difference event

    one.compartment <- function() {
      ini({
        tka <- 0.45 # Log Ka
        tcl <- 1 # Log Cl
        tv <- 3.45    # Log V
        tlag <- log(0.5)
        add.sd <- 1.2
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl)
        v <- exp(tv)
        d/dt(depot) = -ka * depot
        lag(depot) <- exp(tlag)
        d/dt(center) = ka * depot - cl / v * center
        cp = center / v
        cp ~ add(add.sd)
      })
    }

    mod <- one.compartment()

    f <- mod$nlmSensModel

    env <- new.env(parent=emptyenv())

    env$rxControl <- rxode2::rxControl()
    env$predOnly <- f$predOnly
    env$thetaGrad <- f$thetaGrad

    rxode2::rxDynLoad(env$predOnly)
    rxode2::rxDynLoad(env$thetaGrad)

    p <- mod$nlmParIni
    env$param <- setNames(p, sprintf("THETA[%d]", seq_along(p)))
    d <- nlmixr2data::theo_sd
    env$data <- d
    env$needFD <- f$eventTheta

    #.eventTypeIdx <- c("central" = 2L, "forward" = 3L)
    control <- list(stickyRecalcN=4,
                    maxOdeRecalc=5,
                    odeRecalcFactor=10^(0.5),
                    shi21maxInner = 20L,
                    eventType=2L,
                    shi21maxFD=20L,
                    shiErr=(.Machine$double.eps)^(1/3),
                    optimHessType=2L,
                    shi21maxHess=20L,
                    hessErr=(.Machine$double.eps)^(1/3),
                    solveType=3L)

    env$control <- control
    .Call(`_nlmixr2est_nlmSetup`, env)

    .Call(`_nlmixr2est_nlmSolveGradHess`, env$param)

    .Call(`_nlmixr2est_nlmWarnings`)


  })

})
