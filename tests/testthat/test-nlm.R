nmTest({

  library(dplyr)


  test_that("nlm makes sense", {

    dsn <- data.frame(i=1:1000) %>%
      mutate(time=exp(rnorm(n())),DV=rbinom(n(),1,exp(-1+time)/(1+exp(-1+time)))) %>%
      mutate(id=1)

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

    fit2 <- nlmixr(mod, dsn, est="nlm")

    expect_true(inherits(fit2, "nlmixr2.nlm"))

    fit3 <- fit2 %>% ini(g=unfix) %>% nlmixr2(dsn, "nlm", nlmControl(covMethod="nlm"))

    expect_true(inherits(fit3, "nlmixr2.nlm"))

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


    library(dplyr)

    dsn <- data.frame(i=1:1000) %>%
      mutate(time=exp(rnorm(n())),DV=rbinom(n(),1,exp(-1+time)/(1+exp(-1+time)))) %>%
      mutate(id=1)

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

    val <- .Call(`_nlmixr2est_nlmSolveR`, env$param)

    expect_equal(.Call(`_nlmixr2est_nlmSolveGradR`, env$param),
                 .Call(`_nlmixr2est_nlmSolveGradR`, env$param))

    val1 <- val

    val2 <- .Call(`_nlmixr2est_nlmSolveGradR`, env$param)
    attr(val1, "gradient") <- attr(val2, "gradient")

    expect_equal(val1, val2)

    expect_equal(.Call(`_nlmixr2est_nlmSolveGradHess`, env$param),
                 .Call(`_nlmixr2est_nlmSolveGradHess`, env$param))

    val3 <- .Call(`_nlmixr2est_nlmSolveGradHess`, env$param)

    attr(val2, "hessian") <- attr(val3, "hessian")

    expect_equal(val2, val3)

    expect_equal(.Call(`_nlmixr2est_nlmSolveSwitch`, env$param), val2)

    env$control$solveType <- 2L

    .Call(`_nlmixr2est_nlmSetup`, env)

    exepct_equal(.Call(`_nlmixr2est_nlmSolveSwitch`, env$param),
                 .Call(`_nlmixr2est_nlmSolveGradR`, env$param))

    env$control$solveType <- 1L
    .Call(`_nlmixr2est_nlmSetup`, env)

    expect_equal(.Call(`_nlmixr2est_nlmSolveSwitch`, env$param),
                 .Call(`_nlmixr2est_nlmSolveGradR`, env$param))

  })

})
