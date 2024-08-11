nmTest({
  test_that("saem building works; Issue nlmixr#281", {

    Lesion7 <- function() {
      ini({
        temaxD <- -0.05 ; label("typical value of drug emax")
        tec50 <- -0.04 ; label("typical value of ec50")
        temaxT <- -0.9 ; label("typical value of emax for placebo effect as a function of time")
        tet50 <- 6.64 ; label("typical value of et50")
        tkin <- -7.8
        tslope <- -1 ; label("typical value of growth parameter")
        TRXslope <- -0.1
        tdelay <- 7.15
        eta.emaxD ~ 0.3
        eta.emaxT ~ 0.3
        eta.slope ~ 0.3
        eta.delay ~ 0.3
        add.err <- .1 ; label("add. residual variability")
      })
      model({
        Resp(0) <- 1 # default Resp(0) = 0
        emaxD <- exp(temaxD + eta.emaxD)
        ec50 <- exp(tec50)
        emaxT <- exp(temaxT + eta.emaxT)
        et50 <- exp(tet50)
        slope <- exp(tslope + TRX * TRXslope + eta.slope)
        delay <- exp(tdelay + eta.delay)
        kin <- exp(tkin)
        kout <- kin
        GAM <- exp(1) / 2 # Hill coef
        C2 <- centr / V
        CONC <- (C2)^2
        Stim1 <- emaxT * (time) / (time + et50)
        Stim2 <- emaxD * (CONC^GAM) / (CONC^GAM + ec50^GAM)
        Stim <- Stim1 * (1 + TRX * Stim2)
        Delta <- 1 / (1 + exp(-20 * (time - delay)))
        d / dt(depot) <- -KA * depot
        d / dt(centr) <- KA * depot - CL * C2
        d / dt(Resp) <- kin * (1 + Delta * slope) - kout * (1 + Stim) * Resp
        Resp ~ add(add.err)
      })
    }

    tmp <- nlmixr(Lesion7)

    expect_error(tmp$saemModel, NA)

    expect_false(grepl("linear\\(TRX\\)", tmp$saemModel))

    expect_error(tmp$saemModelPred, NA)
    expect_false(grepl("linear\\(TRX\\)", rxode2::rxNorm(tmp$saemModelPred$predOnly)))

    m <- tmp$foceiModel

    expect_false(grepl("linear\\(TRX\\)", rxode2::rxNorm(m$inner)))

    expect_false(grepl("linear\\(TRX\\)", rxode2::rxNorm(m$predOnly)))

    expect_false(grepl("linear\\(TRX\\)", rxode2::rxNorm(m$predNoLhs)))
  })

  test_that("locf and other issues", {

    Lesion7 <- function() {
      ini({
        temaxD <- -0.05 ; label("typical value of drug emax")
        tec50 <- -0.04 ; label("typical value of ec50")
        temaxT <- -0.9 ; label("typical value of emax for placebo effect as a function of time")
        tet50 <- 6.64 ; label("typical value of et50")
        tkin <- -7.8
        tslope <- -1 ; label("typical value of growth parameter")
        TRXslope <- -0.1
        tdelay <- 7.15
        eta.emaxD ~ 0.3
        eta.emaxT ~ 0.3
        eta.slope ~ 0.3
        eta.delay ~ 0.3
        add.err <- .1 ; label("add. residual variability")
      })
      model({
        linear(TRX)
        Resp(0) <- 1 # default Resp(0) = 0
        emaxD <- exp(temaxD + eta.emaxD)
        ec50 <- exp(tec50)
        emaxT <- exp(temaxT + eta.emaxT)
        et50 <- exp(tet50)
        slope <- exp(tslope + TRX * TRXslope + eta.slope)
        delay <- exp(tdelay + eta.delay)
        kin <- exp(tkin)
        kout <- kin
        GAM <- exp(1) / 2 # Hill coef
        C2 <- centr / V
        CONC <- (C2)^2
        Stim1 <- emaxT * (time) / (time + et50)
        Stim2 <- emaxD * (CONC^GAM) / (CONC^GAM + ec50^GAM)
        Stim <- Stim1 * (1 + TRX * Stim2)
        Delta <- 1 / (1 + exp(-20 * (time - delay)))
        d / dt(depot) <- -KA * depot
        d / dt(centr) <- KA * depot - CL * C2
        d / dt(Resp) <- kin * (1 + Delta * slope) - kout * (1 + Stim) * Resp
        Resp ~ add(add.err)
      })
    }

    tmp <- nlmixr(Lesion7)

    expect_error(tmp$saemModel, NA)
    expect_true(grepl("linear\\(TRX\\)", tmp$saemModel))

    expect_error(tmp$saemModelPred, NA)
    expect_true(grepl("linear\\(TRX\\)", rxode2::rxNorm(tmp$saemModelPred$predOnly)))

    m <- tmp$foceiModel

    expect_true(grepl("linear\\(TRX\\)", rxode2::rxNorm(m$inner)))

    expect_true(grepl("linear\\(TRX\\)", rxode2::rxNorm(m$predOnly)))

    expect_true(grepl("linear\\(TRX\\)", rxode2::rxNorm(m$predNoLhs)))
  })

})
