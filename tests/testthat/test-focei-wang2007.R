nmTest({

  .nlmixr <- function(...) {
    rxode2::rxUnloadAll()
    suppressMessages(suppressWarnings(nlmixr(...)))
  }
  #.nlmixr <- function(...) nlmixr(...)


  dat <- Wang2007
  dat$DV <- dat$Y

  ## Enhance data frame to include dosing records.
  dat2 <- dat[dat$Time == 0, ]
  dat2$EVID <- 101
  dat2$AMT <- 10
  dat2 <- rbind(dat2, data.frame(dat, EVID = 0, AMT = 0))
  dat2 <- dat2[(order(dat2$ID, -dat2$EVID, dat2$Time)), ]

  datl <- dat
  datl$DV <- log(datl$DV)

  datl2 <- dat2
  datl2$DV <- log(datl2$DV)

  f <- function() {
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

  f <- nlmixr(f)

  fit.prop <- .nlmixr(f, dat, "focei", foceiControl(maxOuterIterations = 0, covMethod = ""))

  fo <- function() {
    ini({
      tke <- 0.5
      eta.ke ~ 0.04
      prop.sd <- sqrt(0.1)
    })
    model({
      ke <- tke * exp(eta.ke)
      d/dt(ipre) <- -ke * ipre
      f2 <- ipre / (ipre + 5)
      f3 <- f2 * 3
      lipre <- log(ipre)
      ipre ~ prop(prop.sd)
    })
  }

  fo <- nlmixr(fo)

  testErr <- function(type, fun, val = rep(NA_real_, 6), addProp = 2, log=FALSE) {
    .f <- fun(f)
    .fo <- fun(fo)
    .dode <- dat2
    .d <- dat
    if (log) {
      .d <- datl
      .dode <- datl2
    }
    fit1 <- .nlmixr(.fo, .dode, "focei",
                    control = foceiControl(
                      maxOuterIterations = 0, covMethod = "",
                      addProp = paste0("combined", addProp)))

    fit2 <- .nlmixr(.f, .d, "focei",
                    control = foceiControl(
                      maxOuterIterations = 0, covMethod = "",
                      addProp = paste0("combined", addProp)))

    fit3 <- .nlmixr(.fo, .dode, "foce",
                    control = foceiControl(
                      maxOuterIterations = 0, covMethod = "",
                      addProp = paste0("combined", addProp)))

    fit4 <- .nlmixr(.f, .d, "foce",
                    control = foceiControl(
                      maxOuterIterations = 0, covMethod = "",
                      addProp = paste0("combined", addProp)))

    fit5 <- .nlmixr(.fo, .dode, "fo",
                    control = foceiControl(
                      maxOuterIterations = 0, covMethod = "",
                      addProp = paste0("combined", addProp)))
    setOfv(fit5, "fo")

    fit6 <- .nlmixr(.f, .dode, "fo", control = foceiControl(
      maxOuterIterations = 0, covMethod = "",
      addProp = paste0("combined", addProp)))
    setOfv(fit6, "fo")


    .n <- paste(type, c("focei ode", "focei", "foce ode", "foce", "fo ode", "fo"), paste0("combined", addProp))
    ret <- c(fit1$objective, fit2$objective, fit3$objective, fit4$objective, fit5$objective, fit6$objective)
    ret <- setNames(ret, .n)
    val <- setNames(val, .n)
    ## Now test
    if (!all(is.na(val))) {
      test_that(
        type,{
          expect_equal(setNames(ret, NULL), setNames(val, NULL), tolerance=1e-3)
        })
      test_that(paste0(type, " print"),{
        utils::capture.output(suppressMessages({
          withr::with_options(list(cli.unicode=TRUE),expect_error(print(fit1), NA))
          withr::with_options(list(cli.unicode=FALSE),expect_error(print(fit1), NA))
          withr::with_options(list(cli.unicode=TRUE),expect_error(print(fit2), NA))
          withr::with_options(list(cli.unicode=FALSE),expect_error(print(fit2), NA))
          withr::with_options(list(cli.unicode=TRUE),expect_error(print(fit3), NA))
          withr::with_options(list(cli.unicode=FALSE),expect_error(print(fit3), NA))
          withr::with_options(list(cli.unicode=TRUE),expect_error(print(fit4), NA))
          withr::with_options(list(cli.unicode=FALSE),expect_error(print(fit4), NA))
          withr::with_options(list(cli.unicode=TRUE),expect_error(print(fit5), NA))
          withr::with_options(list(cli.unicode=FALSE),expect_error(print(fit5), NA))
          withr::with_options(list(cli.unicode=TRUE),expect_error(print(fit6), NA))
          withr::with_options(list(cli.unicode=FALSE),expect_error(print(fit6), NA))
        }))
      })
    }
    return(invisible(ret))
  }


  fit.prop2 <- .nlmixr(fo, dat2, "focei", foceiControl(maxOuterIterations = 0, covMethod = ""))

  out.focei.prop <- qs::qread("out.focei.prop.qs")

  test_that("Matches NONMEM objective proportional function; (Based on Wang2007)", {
    expect_equal(fit.prop$objective, 39.458, tolerance=1e-3) # Matches Table 2 Prop FOCEI for NONMEM
    expect_equal(fit.prop$eta.ke, out.focei.prop$ETA1, tolerance=1e-3) # match NONMEM output
    ## Individual properties
    expect_equal(fit.prop$IPRED, out.focei.prop$IPRE, tolerance=1e-3)
    expect_equal(fit.prop$IRES, out.focei.prop$IRES, tolerance=1e-3)
    expect_equal(fit.prop$IWRES, out.focei.prop$IWRES, tolerance=1e-3)
    ## WRES variants
    expect_equal(fit.prop$PRED, out.focei.prop$NPRED, tolerance=1e-3) # matches output of PRED from NONMEM
    expect_equal(fit.prop$PRED, out.focei.prop$PRED, tolerance=1e-3) # matches output of PRED from NONMEM
    expect_equal(fit.prop$RES, out.focei.prop$RES, tolerance=1e-3) # match NONMEM output
    expect_equal(fit.prop$RES, out.focei.prop$NRES, tolerance=1e-3) # match NONMEM output
    ## FOI equivalents
    expect_equal(fit.prop$PRED, out.focei.prop$PREDI, tolerance=1e-3) # matches output of PRED from NONMEM
    ## CWRES variants
    expect_equal(fit.prop$CRES, out.focei.prop$CRES, tolerance=1e-3) # match NONMEM output
    expect_equal(fit.prop$CPRED, out.focei.prop$CPRED, tolerance=1e-3) # match NONMEM output
    expect_equal(fit.prop$CWRES, out.focei.prop$CWRES, tolerance=1e-3) # match NONMEM output
    ## Note that E[x] for CPRED and CPREDI are equal
    expect_equal(fit.prop$CRES, out.focei.prop$CRESI, tolerance=1e-3) # match NONMEM output
    expect_equal(fit.prop$CPRED, out.focei.prop$CPREDI, tolerance=1e-3) # match NONMEM output
  })


  test_that("Matches NONMEM objective proportional function; (Based on Wang2007; unoptimized)", {
    # Check unoptimized expression
    expect_equal(fit.prop2$objective, 39.458, tolerance=1e-3) # Matches Table 2 Prop FOCEI for NONMEM
    expect_equal(fit.prop2$eta.ke, out.focei.prop$ETA1, tolerance=1e-3) # match NONMEM output
    ## Individual properties
    expect_equal(fit.prop2$IPRED, out.focei.prop$IPRE, tolerance=1e-3)
    expect_equal(fit.prop2$IRES, out.focei.prop$IRES, tolerance=1e-3)
    expect_equal(fit.prop2$IWRES, out.focei.prop$IWRES, tolerance=1e-3)
    ## WRES variants
    expect_equal(fit.prop2$PRED, out.focei.prop$NPRED, tolerance=1e-3) # matches output of PRED from NONMEM
    expect_equal(fit.prop2$PRED, out.focei.prop$PRED, tolerance=1e-3) # matches output of PRED from NONMEM
    expect_equal(fit.prop2$RES, out.focei.prop$RES, tolerance=1e-3) # match NONMEM output
    expect_equal(fit.prop2$RES, out.focei.prop$NRES, tolerance=1e-3) # match NONMEM output
    ## FOI equivalents
    expect_equal(fit.prop2$PRED, out.focei.prop$PREDI, tolerance=1e-3) # matches output of PRED from NONMEM
    ## CWRES variants
    expect_equal(fit.prop2$CRES, out.focei.prop$CRES, tolerance=1e-3) # match NONMEM output
    expect_equal(fit.prop2$CPRED, out.focei.prop$CPRED, tolerance=1e-3) # match NONMEM output
    expect_equal(fit.prop2$CWRES, out.focei.prop$CWRES, tolerance=1e-3) # match NONMEM output
    ## Note that E[x] for CPRED and CPREDI are equal
    expect_equal(fit.prop2$CRES, out.focei.prop$CRESI, tolerance=1e-3) # match NONMEM output
    expect_equal(fit.prop2$CPRED, out.focei.prop$CPREDI, tolerance=1e-3) # match NONMEM output
  })


  .propVals <- c(39.458, 39.458, 39.275, 39.207, 39.213, 39.207)
  .propModVals <- c(63.353, 63.353, 63.001, 63.063, 63.063, 63.063)
  .propFVals <- c(6.496, 6.496, 6.488, 6.275, 9.262, 9.262)
  .propFModVals <- c(19.177, 19.177, 19.07, 18.202, 18.333, 18.333)

  .addVals <- c(-2.059, -2.059, -2.059, -2.059, 0.026, 0.026)

  .powVals <- c(9.966, 9.966, 9.948, 9.331, 9.651, 9.651)
  .powFModVals <- c(0.776, 0.776, 0.772, 0.58, 3.2, 3.2)
  .powF1Vals <- c(27.301, 27.301, 27.147, 26.854, 26.888, 26.888)
  .powF2Vals <- c(17.831, 17.831, 17.785, 16.617, 16.852, 16.852)
  .powF3Vals <- c(79.733, 79.733, 79.371, 79.448, 79.448, 79.448)
  .powFMod1Vals <- c(24.877, 24.877, 24.773, 24.518, 24.554, 24.554)
  .powFMod2Vals <- c(49.312, 49.312, 49.311, 49.296, 52.603, 52.603)
  .powFMod3Vals <- c(10.848, 10.848, 10.784, 9.446, 9.96, 9.96)
  .addProp2 <- c(39.735, 39.735, 39.562, 39.499, 39.505, 39.499)
  .addModPropMod2 <- c(106.308, 106.308, 106.013, 106.079, 106.079, 106.079)
  .addProp1 <- c(43.554, 43.554, 43.416, 43.394, 43.398, 43.394)
  .addPow1 <- c(16.231, 16.231, 16.219, 16.008, 16.093, 16.093)
  .addPow2 <- c(10.886, 10.886, 10.868, 10.417, 10.662, 10.662)
  .addModVals <- c(3.238, 3.238, 3.207, 2.438, 3.311, 3.311)
  .addModPropMod1 <- c(106.308, 106.308, 106.013, 106.079, 106.079, 106.079)
  .addModPropFModVals2 <- c(54.317, 54.317, 54.14, 54.165, 54.166, 54.166)
  .addPropFVals2 <- c(-2.321, -2.321, -2.322, -2.454, -0.65, -0.65)

  ################################################################################
  # Propotional tests
  ################################################################################
  testErr("prop", function(f) {
    f %>% model(ipre ~ prop(prop.sd)) %>% ini(prop.sd=sqrt(0.1))
  }, .propVals)

  testErr("propMod", function(f) {
    f %>% model(ipre ~ prop(f2))
  }, .propModVals)

  # In this case propT = prop
  testErr("propT", function(f) {
    f %>% model(ipre ~ propT(prop.sd)) %>% ini(prop.sd=sqrt(0.1))
  }, .propVals)

  testErr("propTMod", function(f) {
    f %>% model(ipre ~ propT(f2))
  }, .propModVals)

  testErr("propF", function(f) {
    f %>% model(ipre ~ propF(prop.sd, f2)) %>% ini(prop.sd=sqrt(0.1))
  }, .propFVals)

  testErr("propFMod", function(f) {
    f %>% model(ipre ~ propF(lipre, f2))
  }, .propFModVals)

  ################################################################################
  # Additive Model Tests
  ################################################################################
  testErr("add", function(f) {
    f %>% model(ipre ~ add(add.sd)) %>% ini(add.sd=sqrt(0.1))
  }, .addVals)

  testErr("addMod", function(f) {
    f %>% model(ipre ~ add(f2))
  }, .addModVals)

  skip_on_cran()
  rxode2::rxUnloadAll() # don't do too much on windows because of dll overloading
  skip_on_os("windows")

  ################################################################################
  # Power Model Tests
  ################################################################################
  testErr("pow1=prop", function(f) {
    f %>% model(ipre ~ pow(pow.sd, pw)) %>% ini(pow.sd=sqrt(0.1), pw=1)
  }, .propVals)

  testErr("powT1=propT=prop", function(f) {
    f %>% model(ipre ~ powT(pow.sd, pw)) %>% ini(pow.sd=sqrt(0.1), pw=1)
  }, .propVals)

  testErr("powF1=propF", function(f) {
    f %>% model(ipre ~ powF(pow.sd, pw, f2)) %>% ini(pow.sd=sqrt(0.1), pw=1)
  }, .propFVals)

  testErr("pow1Mod=propMod", function(f) {
    f %>% model(ipre ~ pow(f2, pw)) %>% ini(pw=1)
  }, .propModVals)

  testErr("pow", function(f) {
    f %>% model(ipre ~ pow(pow.sd, pw)) %>% ini(pow.sd=sqrt(0.1), pw=0.5)
  }, .powVals)

  testErr("powF1", function(f) {
    f %>% model(ipre ~ pow(f2, pw)) %>% ini(pw=0.5)
  }, .powF1Vals)

  testErr("powF2", function(f) {
    f %>% model(ipre ~ pow(pow.sd, f2)) %>% ini(pow.sd=sqrt(0.1))
  }, .powF2Vals)

  testErr("powF3", function(f) {
    f %>% model(ipre ~ pow(lipre, f2))
  }, .powF3Vals)

  testErr("powT", function(f) {
    f %>% model(ipre ~ powT(pow.sd, pw)) %>% ini(pow.sd=sqrt(0.1), pw=0.5)
  }, .powVals)

  testErr("powTF1", function(f) {
    f %>% model(ipre ~ powT(f2, pw)) %>% ini(pw=0.5)
  }, .powF1Vals)

  testErr("powFT2", function(f) {
    f %>% model(ipre ~ powT(pow.sd, f2)) %>% ini(pow.sd=sqrt(0.1))
  }, .powF2Vals)

  testErr("powTF3", function(f) {
    f %>% model(ipre ~ powT(lipre, f2))
  }, .powF3Vals)

  testErr("powFMod", function(f) {
    f %>% model(ipre ~ powF(pow.sd, pw, f2)) %>% ini(pow.sd=sqrt(0.1), pw=0.5)
  }, .powFModVals)

  testErr("powFMod1", function(f) {
    f %>% model(ipre ~ powF(lipre, pw, f2)) %>% ini(pw=0.5)
  }, .powFMod1Vals)

  testErr("powFMod2", function(f) {
    f %>% model(ipre ~ powF(pow.sd, lipre, f2)) %>% ini(pow.sd=sqrt(0.1))
  }, .powFMod2Vals)

  testErr("powFMod3", function(f) {
    f %>% model(ipre ~ powF(lipre, f3, f2))
  }, .powFMod3Vals)

  ################################################################################
  # Add+Proportional tests (combined 2)
  ################################################################################
  testErr("add+prop, combined 2->add", function(f) {
    f %>% model(ipre ~ add(add.sd) + prop(prop.sd)) %>% ini(add.sd=sqrt(0.1), prop.sd=0)
  }, .addVals, addProp = 2)

  testErr("addMod+prop, combined 2->addMod", function(f) {
    f %>% model(ipre ~ add(f2) + prop(prop.sd)) %>% ini(prop.sd=0)
  }, .addModVals, addProp = 2)

  testErr("add+prop, combined 2->prop", function(f) {
    f %>% model(ipre ~ add(add.sd) + prop(prop.sd)) %>% ini(add.sd=0, prop.sd=sqrt(0.1))
  }, .propVals, addProp = 2)

  testErr("add+propMod, combined 2->propMod", function(f) {
    f %>% model(ipre ~ add(add.sd) + prop(f2)) %>% ini(add.sd=0)
  }, .propModVals, addProp = 2)

  testErr("addMod+propMod, combined 2->propMod", function(f) {
    f %>% model(ipre ~ add(f2) + prop(f3))
  },  .addModPropMod2, addProp = 2)

  testErr("add+prop, combined 2", function(f) {
    f %>% model(ipre ~ add(add.sd) + prop(prop.sd)) %>% ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1))
  }, .addProp2, addProp = 2)

  testErr("add+prop, combined 2", function(f) {
    f %>% model(ipre ~ add(add.sd) + prop(prop.sd) + combined2()) %>% ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1))
  }, .addProp2, addProp = 1)

  ## propT
  testErr("add+propT, combined 2->add", function(f) {
    f %>% model(ipre ~ add(add.sd) + propT(prop.sd)) %>% ini(add.sd=sqrt(0.1), prop.sd=0)
  }, .addVals, addProp = 2)

  testErr("addMod+propT, combined 2->addMod", function(f) {
    f %>% model(ipre ~ add(f2) + propT(prop.sd)) %>% ini(prop.sd=0)
  }, .addModVals, addProp = 2)

  testErr("add+propT, combined 2->prop", function(f) {
    f %>% model(ipre ~ add(add.sd) + propT(prop.sd)) %>% ini(add.sd=0, prop.sd=sqrt(0.1))
  }, .propVals, addProp = 2)

  testErr("add+propTMod, combined 2->propMod", function(f) {
    f %>% model(ipre ~ add(add.sd) + propT(f2)) %>% ini(add.sd=0)
  }, .propModVals, addProp = 2)

  testErr("addMod+propTMod, combined 2->propMod", function(f) {
    f %>% model(ipre ~ add(f2) + propT(f3))
  },  .addModPropMod2, addProp = 2)

  testErr("add+propT, combined 2", function(f) {
    f %>% model(ipre ~ add(add.sd) + propT(prop.sd)) %>% ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1))
  }, .addProp2, addProp = 2)

  testErr("add+propT, combined 2 (specified)", function(f) {
    f %>% model(ipre ~ add(add.sd) + propT(prop.sd) + combined2()) %>% ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1))
  }, .addProp2, addProp = 1)

  # propF
  testErr("add+propF, combined 2->add", function(f) {
    f %>% model(ipre ~ add(add.sd) + propF(prop.sd, f2)) %>% ini(add.sd=sqrt(0.1), prop.sd=0)
  }, .addVals, addProp = 2)

  testErr("addMod+propF, combined 2->addMod", function(f) {
    f %>% model(ipre ~ add(f2) + propF(prop.sd, f2)) %>% ini(prop.sd=0)
  }, .addModVals, addProp = 2)

  testErr("add+propF, combined 2->prop", function(f) {
    f %>% model(ipre ~ add(add.sd) + propF(prop.sd, f2)) %>% ini(add.sd=0, prop.sd=sqrt(0.1))
  }, .propFVals, addProp = 2)

  testErr("add+propFMod, combined 2->propMod", function(f) {
    f %>% model(ipre ~ add(add.sd) + propF(lipre, f2)) %>% ini(add.sd=0)
  }, .propFModVals, addProp = 2)

  testErr("addMod+propFMod, combined 2->propMod", function(f) {
    f %>% model(ipre ~ add(f2) + propF(lipre, f3))
  }, .addModPropFModVals2, addProp = 2)

  testErr("add+propF, combined 2", function(f) {
    f %>% model(ipre ~ add(add.sd) + propF(prop.sd, f2)) %>% ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1))
  },  .addPropFVals2, addProp = 2)

  testErr("add+propF, combined 2 (specified)", function(f) {
    f %>% model(ipre ~ add(add.sd) + propF(prop.sd, f2) + combined2()) %>% ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1))
  }, .addPropFVals2, addProp = 1)

  ################################################################################
  # Add+Proportional tests (combined 1)
  ################################################################################
  testErr("add+prop, combined 1->add", function(f) {
    f %>% model(ipre ~ add(add.sd) + prop(prop.sd)) %>% ini(add.sd=sqrt(0.1), prop.sd=0)
  }, .addVals, addProp = 1)

  testErr("addMod+prop, combined 1->addMod", function(f) {
    f %>% model(ipre ~ add(f2) + prop(prop.sd)) %>% ini(prop.sd=0)
  }, .addModVals, addProp = 1)

  testErr("add+prop, combined 1->prop", function(f) {
    f %>% model(ipre ~ add(add.sd) + prop(prop.sd)) %>% ini(add.sd=0, prop.sd=sqrt(0.1))
  }, .propVals, addProp = 1)

  testErr("add+propMod, combined 1->propMod", function(f) {
    f %>% model(ipre ~ add(add.sd) + prop(f2)) %>% ini(add.sd=0)
  }, .propModVals, addProp = 1)

  testErr("addMod+propMod, combined 1->propMod", function(f) {
    f %>% model(ipre ~ add(f2) + prop(f3))
  },  .addModPropMod1, addProp = 2)

  testErr("add+prop, combined 2", function(f) {
    f %>% model(ipre ~ add(add.sd) + prop(prop.sd)) %>% ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1))
  }, .addProp2, addProp = 2)

  testErr("add+prop, combined 2", function(f) {
    f %>% model(ipre ~ add(add.sd) + prop(prop.sd) + combined2()) %>% ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1))
  }, .addProp2, addProp = 1)

  ## propT
  testErr("add+propT, combined 2->add", function(f) {
    f %>% model(ipre ~ add(add.sd) + propT(prop.sd)) %>% ini(add.sd=sqrt(0.1), prop.sd=0)
  }, .addVals, addProp = 2)

  testErr("addMod+propT, combined 2->addMod", function(f) {
    f %>% model(ipre ~ add(f2) + propT(prop.sd)) %>% ini(prop.sd=0)
  }, .addModVals, addProp = 2)

  testErr("add+propT, combined 2->prop", function(f) {
    f %>% model(ipre ~ add(add.sd) + propT(prop.sd)) %>% ini(add.sd=0, prop.sd=sqrt(0.1))
  }, .propVals, addProp = 2)

  testErr("add+propTMod, combined 2->propMod", function(f) {
    f %>% model(ipre ~ add(add.sd) + propT(f2)) %>% ini(add.sd=0)
  }, .propModVals, addProp = 2)

  testErr("addMod+propTMod, combined 2->propMod", function(f) {
    f %>% model(ipre ~ add(f2) + propT(f3))
  },  .addModPropMod2, addProp = 2)

  testErr("add+propT, combined 2", function(f) {
    f %>% model(ipre ~ add(add.sd) + propT(prop.sd)) %>% ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1))
  }, .addProp2, addProp = 2)

  testErr("add+propT, combined 2 (specified)", function(f) {
    f %>% model(ipre ~ add(add.sd) + propT(prop.sd) + combined2()) %>% ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1))
  }, .addProp2, addProp = 1)

  # propF
  testErr("add+propF, combined 2->add", function(f) {
    f %>% model(ipre ~ add(add.sd) + propF(prop.sd, f2)) %>% ini(add.sd=sqrt(0.1), prop.sd=0)
  }, .addVals, addProp = 2)

  testErr("addMod+propF, combined 2->addMod", function(f) {
    f %>% model(ipre ~ add(f2) + propF(prop.sd, f2)) %>% ini(prop.sd=0)
  }, .addModVals, addProp = 2)

  testErr("add+propF, combined 2->prop", function(f) {
    f %>% model(ipre ~ add(add.sd) + propF(prop.sd, f2)) %>% ini(add.sd=0, prop.sd=sqrt(0.1))
  }, .propFVals, addProp = 2)

  testErr("add+propFMod, combined 2->propMod", function(f) {
    f %>% model(ipre ~ add(add.sd) + propF(lipre, f2)) %>% ini(add.sd=0)
  }, .propFModVals, addProp = 2)

  testErr("addMod+propFMod, combined 2->propMod", function(f) {
    f %>% model(ipre ~ add(f2) + propF(lipre, f3))
  }, .addModPropFModVals2, addProp = 2)

  testErr("add+propF, combined 2", function(f) {
    f %>% model(ipre ~ add(add.sd) + propF(prop.sd, f2)) %>% ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1))
  },  .addPropFVals2, addProp = 2)

  testErr("add+propF, combined 2 (specified)", function(f) {
    f %>% model(ipre ~ add(add.sd) + propF(prop.sd, f2) + combined2()) %>% ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1))
  }, .addPropFVals2, addProp = 1)

  #################################################################################################
  testErr("add+prop, combined 1", function(f) {
    f %>% model(ipre ~ add(add.sd) + prop(prop.sd)) %>% ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1))
  }, .addProp1, addProp = 1)

  testErr("add+prop, combined 1->add", function(f) {
    f %>% model(ipre ~ add(add.sd) + prop(prop.sd)) %>% ini(add.sd=sqrt(0.1), prop.sd=0)
  }, .addVals, addProp = 1)

  testErr("add+prop, combined 1->prop", function(f) {
    f %>% model(ipre ~ add(add.sd) + prop(prop.sd)) %>% ini(add.sd=0, prop.sd=sqrt(0.1))
  }, .propVals, addProp = 1)

  testErr("add+prop, combined 1 (override)", function(f) {
    f %>% model(ipre ~ add(add.sd) + prop(prop.sd) + combined1()) %>% ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1))
  }, .addProp1, addProp = 2)

  testErr("add+prop, combined 2 (override)", function(f) {
    f %>% model(ipre ~ add(add.sd) + prop(prop.sd) + combined2()) %>% ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1))
  }, .addProp2, addProp = 1)

  ################################################################################
  # Add+Pow tests (combined 2)
  ################################################################################

  testErr("add+pow combined 2 -> add+prop combined2", function(f) {
    f %>% model(ipre ~ add(add.sd) + pow(prop.sd, pw)) %>% ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=1)
  }, .addProp2, addProp = 2)

  testErr("add+pow combined 2", function(f) {
    f %>% model(ipre ~ add(add.sd) + pow(prop.sd, pw)) %>% ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5)
  }, c(10.886, 10.886, 10.868, 10.417, 10.662, 10.662), addProp = 2)

  testErr("add+pow combined 1->add+prop combined1", function(f) {
    f %>% model(ipre ~ add(add.sd) + pow(prop.sd, pw)) %>% ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=1)
  }, .addProp1, addProp = 1)

  testErr("add+pow combined 1", function(f) {
    f %>% model(ipre ~ add(add.sd) + pow(prop.sd, pw)) %>% ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5)
  }, .addPow1, addProp = 1)

  testErr("add+pow combined 1 (override)", function(f) {
    f %>% model(ipre ~ add(add.sd) + pow(prop.sd, pw) + combined1()) %>% ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5)
  }, .addPow1, addProp = 2)


  ## start looking at transformations

  ################################################################################
  ## lognormal error structures
  ################################################################################

  .lnorm <- c(40.039, 40.039, 40.039, 40.039, 40.055, 40.039)

  .lnormProp <- c(118.419, 118.419, 118.279, 118.311, 118.311, 118.311)
  .lnormPropT <- c(65.886, 65.886, 65.827, 65.832, 65.832, 65.832)
  .lnormPropF <- c(27.293, 27.293, 27.274, 26.643, 27.035, 27.035)

  .lnormPow <- c(77.87, 77.87, 77.826, 77.836, 77.836, 77.836)
  .lnormPowT <- c(52.616, 52.616, 52.601, 52.584, 52.587, 52.587)
  .lnormPowF <- c(32.834, 32.834, 32.834, 32.697, 32.784, 32.784)

  .lnormProp1 <- c(123.318, 123.318, 123.219, 123.24, 123.24, 123.24)
  .lnormPropT1 <- c(81.152, 81.152, 81.13, 81.135, 81.135, 81.135)
  .lnormPropF1 <- c(56.785, 56.785, 56.78, 56.774, 56.775, 56.775)

  .lnormPowT1 <- c(72.356, 72.356, 72.351, 72.351, 72.351, 72.351)
  .lnormPowF1 <- c(60.494, 60.494, 60.492, 60.489, 60.49, 60.49)
  .lnormPow1 <- c(89.824, 89.824, 89.803, 89.809, 89.809, 89.809)

  .lnormProp2 <- c(118.777, 118.777, 118.646, 118.676, 118.676, 118.676)
  .lnormPropT2 <- c(69.981, 69.981, 69.947, 69.951, 69.951, 69.951)
  .lnormPropF2 <- c(45.498, 45.498, 45.495, 45.482, 45.49, 45.49)

  .lnormPow2 <- c(80.244, 80.244, 80.212, 80.219, 80.219, 80.219)
  .lnormPowF2 <- c(48.202, 48.202, 48.2, 48.193, 48.198, 48.198)
  .lnormPowT2 <- c(59.837, 59.837, 59.831, 59.826, 59.827, 59.827)

  testErr("lnorm", function(f) {
    f %>% model(ipre ~ lnorm(lnorm.sd)) %>% ini(lnorm.sd=sqrt(0.1))
  }, .lnorm, addProp = 1)

  test_that("lognormal likelihood can be determined by data too", {
    suppressWarnings({
      expect_equal(setNames(round(testErr("add lnorm", function(f) {
        f %>% model(lipre ~ add(add.sd)) %>% ini(add.sd=sqrt(0.1))
      },  log=TRUE) + 2 * sum(datl$DV), 3), NULL),
      .lnorm, tolerance=5e-3)
    })
  })


  ################################################################################
  ## lnorm(NA) tests
  ################################################################################

  testErr("lnorm(NA)+prop", function(f) {
    f %>% model(ipre ~ lnorm(NA) + prop(prop.sd)) %>% ini(prop.sd=sqrt(0.1))
  }, .lnormProp, addProp = 1)

  testErr("lnorm(NA)+propT", function(f) {
    f %>% model(ipre ~ lnorm(NA) + propT(prop.sd)) %>% ini(prop.sd=sqrt(0.1))
  }, .lnormPropT, addProp = 1)

  testErr("lnorm(NA)+propF", function(f) {
    f %>% model(ipre ~ lnorm(NA) + propF(prop.sd, f2)) %>% ini(prop.sd=sqrt(0.1))
  }, .lnormPropF, addProp = 1)

  testErr("lnorm(NA)+pow->lnorm(NA)+prop", function(f) {
    f %>% model(ipre ~ lnorm(NA) + pow(prop.sd, pw)) %>% ini(prop.sd=sqrt(0.1), pw=1)
  }, .lnormProp, addProp = 1)

  testErr("lnorm(NA)+powT->lnorm(NA)+propT", function(f) {
    f %>% model(ipre ~ lnorm(NA) + powT(prop.sd, pw)) %>% ini(prop.sd=sqrt(0.1), pw=1)
  }, .lnormPropT, addProp = 1)

  testErr("lnorm(NA)+powF->lnorm(NA)+propF", function(f) {
    f %>% model(ipre ~ lnorm(NA) + powF(prop.sd, pw, f2)) %>% ini(prop.sd=sqrt(0.1), pw=1)
  }, .lnormPropF, addProp = 1)

  testErr("lnorm(NA)+pow", function(f) {
    f %>% model(ipre ~ lnorm(NA) + pow(prop.sd, pw)) %>% ini(prop.sd=sqrt(0.1), pw=0.5)
  }, .lnormPow, addProp = 1)

  testErr("lnorm(NA)+powT", function(f) {
    f %>% model(ipre ~ lnorm(NA) + powT(prop.sd, pw)) %>% ini(prop.sd=sqrt(0.1), pw=0.5)
  }, .lnormPowT, addProp = 1)

  testErr("lnorm(NA)+powF", function(f) {
    f %>% model(ipre ~ lnorm(NA) + powF(prop.sd, pw, f2)) %>% ini(prop.sd=sqrt(0.1), pw=0.5)
  }, .lnormPowF, addProp = 1)


  ################################################################################
  ## lnorm combined1
  ################################################################################

  testErr("lnorm+prop combined1->lnorm(NA)+prop", function(f) {
    f %>% model(ipre ~ lnorm(lnorm.sd) + prop(prop.sd)) %>% ini(lnorm.sd=0, prop.sd=sqrt(0.1))
  }, .lnormProp, addProp = 1)

  testErr("lnorm+propT combined1->lnorm(NA)+propT", function(f) {
    f %>% model(ipre ~ lnorm(lnorm.sd) + propT(prop.sd)) %>% ini(lnorm.sd=0, prop.sd=sqrt(0.1))
  }, .lnormPropT, addProp = 1)

  testErr("lnorm+propF combined1->lnorm(NA)+propF", function(f) {
    f %>% model(ipre ~ lnorm(lnorm.sd) + propF(prop.sd, f2)) %>% ini(lnorm.sd=0, prop.sd=sqrt(0.1))
  }, .lnormPropF, addProp = 1)

  testErr("lnorm+prop combined1->lnorm(NA)", function(f) {
    f %>% model(ipre ~ lnorm(lnorm.sd) + prop(prop.sd)) %>% ini(lnorm.sd=sqrt(0.1), prop.sd=0)
  }, .lnorm, addProp = 1)

  testErr("lnorm+propT combined1->lnorm(NA)", function(f) {
    f %>% model(ipre ~ lnorm(lnorm.sd) + propT(prop.sd)) %>% ini(lnorm.sd=sqrt(0.1), prop.sd=0)
  }, .lnorm, addProp = 1)

  testErr("lnorm+propF combined1->lnorm(NA)", function(f) {
    f %>% model(ipre ~ lnorm(lnorm.sd) + propF(prop.sd, f2)) %>% ini(lnorm.sd=sqrt(0.1), prop.sd=0)
  }, .lnorm, addProp = 1)

  testErr("lnorm+prop combined1", function(f) {
    f %>% model(ipre ~ lnorm(lnorm.sd) + prop(prop.sd)) %>% ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1))
  }, .lnormProp1, addProp = 1)

  testErr("lnorm+propT combined1", function(f) {
    f %>% model(ipre ~ lnorm(lnorm.sd) + propT(prop.sd)) %>% ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1))
  }, .lnormPropT1, addProp = 1)

  testErr("lnorm+propF combined1", function(f) {
    f %>% model(ipre ~ lnorm(lnorm.sd) + propF(prop.sd, f2)) %>% ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1))
  }, .lnormPropF1, addProp = 1)

  testErr("lnorm+powF->lnorm+propF combined1", function(f) {
    f %>% model(ipre ~ lnorm(lnorm.sd) + powF(prop.sd, pw, f2)) %>% ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=1)
  }, .lnormPropF1, addProp = 1)

  testErr("lnorm+pow->lnorm+prop combined1", function(f) {
    f %>% model(ipre ~ lnorm(lnorm.sd) + pow(prop.sd, pw)) %>% ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=1)
  }, .lnormProp1, addProp = 1)

  testErr("lnorm+powT->lnorm+propT combined1", function(f) {
    f %>% model(ipre ~ lnorm(lnorm.sd) + powT(prop.sd, pw)) %>% ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=1)
  }, .lnormPropT1, addProp = 1)

  testErr("lnorm+powF combined1", function(f) {
    f %>% model(ipre ~ lnorm(lnorm.sd) + powF(prop.sd, pw, f2)) %>% ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5)
  }, .lnormPowF1, addProp = 1)

  testErr("lnorm+pow combined1", function(f) {
    f %>% model(ipre ~ lnorm(lnorm.sd) + pow(prop.sd, pw)) %>% ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5)
  }, .lnormPow1, addProp = 1)

  testErr("lnorm+powT combined1", function(f) {
    f %>% model(ipre ~ lnorm(lnorm.sd) + powT(prop.sd, pw)) %>% ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5)
  }, .lnormPowT1, addProp = 1)

  ################################################################################
  ## lnorm combined2
  ################################################################################

  testErr("lnorm+propT combined2->lnorm(NA)+propT", function(f) {
    f %>% model(ipre ~ lnorm(lnorm.sd) + propT(prop.sd)) %>% ini(lnorm.sd=0, prop.sd=sqrt(0.1))
  }, .lnormPropT, addProp = 2)

  testErr("lnorm+propF combined2->lnorm(NA)+propF", function(f) {
    f %>% model(ipre ~ lnorm(lnorm.sd) + propF(prop.sd, f2)) %>% ini(lnorm.sd=0, prop.sd=sqrt(0.1))
  }, .lnormPropF, addProp = 2)

  testErr("lnorm+prop combined2->lnorm(NA)", function(f) {
    f %>% model(ipre ~ lnorm(lnorm.sd) + prop(prop.sd)) %>% ini(lnorm.sd=sqrt(0.1), prop.sd=0)
  }, .lnorm, addProp = 2)

  testErr("lnorm+propT combined2->lnorm(NA)", function(f) {
    f %>% model(ipre ~ lnorm(lnorm.sd) + propT(prop.sd)) %>% ini(lnorm.sd=sqrt(0.1), prop.sd=0)
  }, .lnorm, addProp = 2)

  testErr("lnorm+propF combined2->lnorm(NA)", function(f) {
    f %>% model(ipre ~ lnorm(lnorm.sd) + propF(prop.sd, f2)) %>% ini(lnorm.sd=sqrt(0.1), prop.sd=0)
  }, .lnorm, addProp = 2)

  testErr("lnorm+pow->lnorm+prop combined2", function(f) {
    f %>% model(ipre ~ lnorm(lnorm.sd) + pow(prop.sd, pw)) %>% ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=1)
  }, .lnormProp2, addProp = 2)

  testErr("lnorm+prop combined2", function(f) {
    f %>% model(ipre ~ lnorm(lnorm.sd) + prop(prop.sd)) %>% ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1))
  }, .lnormProp2 ,addProp = 2)

  testErr("lnorm+propT combined2", function(f) {
    f %>% model(ipre ~ lnorm(lnorm.sd) + propT(prop.sd)) %>% ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1))
  }, .lnormPropT2, addProp = 2)

  testErr("lnorm+propF combined1", function(f) {
    f %>% model(ipre ~ lnorm(lnorm.sd) + propF(prop.sd, f2)) %>% ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1))
  }, .lnormPropF2, addProp = 2)

  testErr("lnorm+powF->lnorm+propF combined2", function(f) {
    f %>% model(ipre ~ lnorm(lnorm.sd) + powF(prop.sd, pw, f2)) %>% ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=1)
  }, .lnormPropF2, addProp = 2)

  testErr("lnorm+pow->lnorm+prop combined2", function(f) {
    f %>% model(ipre ~ lnorm(lnorm.sd) + pow(prop.sd, pw)) %>% ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=1)
  }, .lnormProp2, addProp = 2)

  testErr("lnorm+powT->lnorm+propT combined2", function(f) {
    f %>% model(ipre ~ lnorm(lnorm.sd) + powT(prop.sd, pw)) %>% ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=1)
  }, .lnormPropT2, addProp = 2)

  testErr("lnorm+powF combined2", function(f) {
    f %>% model(ipre ~ lnorm(lnorm.sd) + powF(prop.sd, pw, f2)) %>% ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5)
  }, .lnormPowF2, addProp = 2)

  testErr("lnorm+pow combined2", function(f) {
    f %>% model(ipre ~ lnorm(lnorm.sd) + pow(prop.sd, pw)) %>% ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5)
  }, .lnormPow2, addProp = 2)

  testErr("lnorm+powT combined2", function(f) {
    f %>% model(ipre ~ lnorm(lnorm.sd) + powT(prop.sd, pw)) %>% ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5)
  }, .lnormPowT2, addProp = 2)


  #################################################################################
  # Box-Cox Regression
  #################################################################################

  .boxCoxPropTVals <- c(35.132, 35.132, 34.879, 34.694, 34.709, 34.709)
  .boxCoxPropTModVals <- c(58.21, 58.21, 57.741, 57.802, 57.803, 57.803)
  .boxCoxPowTVals <- c(9.11, 9.11, 9.095, 8.215, 8.644, 8.644)
  .boxCoxPowTF1Vals <- c(25.57, 25.57, 25.388, 24.951, 25.002, 25.002)
  .boxCoxPowTF2Vals <- c(16.507, 16.507, 16.423, 14.981, 15.314, 15.314)
  .boxCoxPowTF3Vals <- c(76.496, 76.496, 76.097, 76.179, 76.179, 76.179)
  .boxCoxAddModPropTMod2 <- c(100.737, 100.737, 100.351, 100.437, 100.437, 100.437)
  .boxCoxAddPropT2Vals <- c(35.457, 35.457, 35.224, 35.056, 35.07, 35.07)
  .boxCoxPropT2ModVals <- c(58.21, 58.21, 57.741, 57.802, 57.803, 57.803)
  .boxCoxPropAddPropT2Vals <- c(100.737, 100.737, 100.351, 100.437, 100.437, 100.437)

  testErr("prop+boxCox->prop", function(f) {
    f %>% model(ipre ~ prop(prop.sd) + boxCox(lambda)) %>% ini(prop.sd=sqrt(0.1), lambda=1)
  }, .propVals)

  testErr("prop+boxCox->propMod", function(f) {
    f %>% model(ipre ~ prop(f2) + boxCox(lambda)) %>% ini(lambda=1)
  }, .propModVals)

  # In this case propT = prop
  testErr("propT+boxCox->propT", function(f) {
    f %>% model(ipre ~ propT(prop.sd) + boxCox(lambda)) %>% ini(prop.sd=sqrt(0.1), lambda=1)
  }, .boxCoxPropTVals)

  testErr("propTMod+boxCox->propTMod", function(f) {
    f %>% model(ipre ~ propT(f2) + boxCox(lambda)) %>% ini(lambda=1)
  }, .boxCoxPropTModVals)

  testErr("propF+boxCox->propF", function(f) {
    f %>% model(ipre ~ propF(prop.sd, f2) + boxCox(lambda)) %>% ini(prop.sd=sqrt(0.1), lambda=1)
  }, .propFVals)

  testErr("propFMod+boxCox->propFMod", function(f) {
    f %>% model(ipre ~ propF(lipre, f2) + boxCox(lambda)) %>% ini(lambda=1)
  }, .propFModVals)

  ################################################################################
  # boxCox -> Additive Model Tests
  ################################################################################
  testErr("add+boxCox->add", function(f) {
    f %>% model(ipre ~ add(add.sd) + boxCox(lambda)) %>% ini(add.sd=sqrt(0.1), lambda=1)
  }, .addVals)

  testErr("addMod+boxCox->addMod", function(f) {
    f %>% model(ipre ~ add(f2) + boxCox(lambda)) %>% ini(lambda=1)
  }, .addModVals)

  ################################################################################
  # boxCox -> Power Model Tests
  ################################################################################
  testErr("boxCox+pow1->prop", function(f) {
    f %>% model(ipre ~ pow(pow.sd, pw) + boxCox(lambda)) %>% ini(pow.sd=sqrt(0.1), pw=1, lambda=1)
  }, .propVals)

  testErr("boxCox+powT1=propT=prop", function(f) {
    f %>% model(ipre ~ powT(pow.sd, pw) + boxCox(lambda)) %>% ini(pow.sd=sqrt(0.1), pw=1, lambda=1)
  }, .boxCoxPropTVals)

  testErr("boxCox+powF1=propF", function(f) {
    f %>% model(ipre ~ powF(pow.sd, pw, f2) + boxCox(lambda)) %>% ini(pow.sd=sqrt(0.1), pw=1, lambda=1)
  }, .propFVals)

  testErr("boxCox+pow1Mod=propMod", function(f) {
    f %>% model(ipre ~ pow(f2, pw) + boxCox(lambda)) %>% ini(pw=1, lambda=1)
  }, .propModVals)

  testErr("boxCox+pow->pow", function(f) {
    f %>% model(ipre ~ pow(pow.sd, pw) + boxCox(lambda)) %>% ini(pow.sd=sqrt(0.1), pw=0.5, lambda=1)
  }, .powVals)

  testErr("boxCox+powF1->powF1", function(f) {
    f %>% model(ipre ~ pow(f2, pw) + boxCox(lambda)) %>% ini(pw=0.5, lambda=1)
  }, .powF1Vals)

  testErr("boxCox+powF2->powF2", function(f) {
    f %>% model(ipre ~ pow(pow.sd, f2) + boxCox(lambda)) %>% ini(pow.sd=sqrt(0.1), lambda=1)
  }, .powF2Vals)

  testErr("boxCox+powF3->powF3", function(f) {
    f %>% model(ipre ~ pow(lipre, f2) + boxCox(lambda)) %>% ini(lambda=1)
  }, .powF3Vals)

  testErr("boxCox+powT->powT", function(f) {
    f %>% model(ipre ~ powT(pow.sd, pw) + boxCox(lambda)) %>% ini(pow.sd=sqrt(0.1), pw=0.5, lambda=1)
  }, .boxCoxPowTVals)

  testErr("boxCox+powTF1->powTF1", function(f) {
    f %>% model(ipre ~ powT(f2, pw) + boxCox(lambda)) %>% ini(pw=0.5, lambda=1)
  }, .boxCoxPowTF1Vals)

  testErr("boxCox+powFT2->powFT2", function(f) {
    f %>% model(ipre ~ powT(pow.sd, f2) + boxCox(lambda)) %>% ini(pow.sd=sqrt(0.1), lambda=1)
  }, .boxCoxPowTF2Vals)

  testErr("boxCox+powTF3->powTF3", function(f) {
    f %>% model(ipre ~ powT(lipre, f2) + boxCox(lambda)) %>% ini(lambda=1)
  }, .boxCoxPowTF3Vals)

  testErr("boxCox+powFMod->powFMod", function(f) {
    f %>% model(ipre ~ powF(pow.sd, pw, f2) + boxCox(lambda)) %>% ini(pow.sd=sqrt(0.1), pw=0.5, lambda=1)
  }, .powFModVals)

  testErr("boxCox+powFMod1->powFMod1", function(f) {
    f %>% model(ipre ~ powF(lipre, pw, f2) + boxCox(lambda)) %>% ini(pw=0.5, lambda=1)
  }, .powFMod1Vals)

  testErr("boxCox+powFMod2->powFMod2", function(f) {
    f %>% model(ipre ~ powF(pow.sd, lipre, f2) + boxCox(lambda)) %>% ini(pow.sd=sqrt(0.1), lambda=1)
  }, .powFMod2Vals)

  testErr("boxCox+powFMod3->powFMod3", function(f) {
    f %>% model(ipre ~ powF(lipre, f3, f2) + boxCox(lambda)) %>% ini(lambda=1)
  }, .powFMod3Vals)

  ################################################################################
  # Box-Cox Add+Proportional tests (combined 2)
  ################################################################################
  testErr("boxCox+add+prop, combined 2->add", function(f) {
    f %>% model(ipre ~ add(add.sd) + prop(prop.sd) + boxCox(lambda)) %>% ini(add.sd=sqrt(0.1), prop.sd=0, lambda=1)
  }, .addVals, addProp = 2)

  testErr("boxCox+addMod+prop, combined 2->addMod", function(f) {
    f %>% model(ipre ~ add(f2) + prop(prop.sd) + boxCox(lambda)) %>% ini(prop.sd=0, lambda=1)
  }, .addModVals, addProp = 2)

  testErr("boxCox+add+prop, combined 2->prop", function(f) {
    f %>% model(ipre ~ add(add.sd) + prop(prop.sd) + boxCox(lambda)) %>% ini(add.sd=0, prop.sd=sqrt(0.1), lambda=1)
  }, .propVals, addProp = 2)

  testErr("boxCox+add+propMod, combined 2->propMod", function(f) {
    f %>% model(ipre ~ add(add.sd) + prop(f2) + boxCox(lambda)) %>% ini(add.sd=0, lambda=1)
  }, .propModVals, addProp = 2)

  testErr("boxCox+addMod+propMod, combined 2->propMod", function(f) {
    f %>% model(ipre ~ add(f2) + prop(f3) + boxCox(lambda)) %>% ini(lambda=1)
  }, .addModPropMod2, addProp = 2)

  testErr("boxCox+add+prop, combined 2->add+prop, combined 2", function(f) {
    f %>% model(ipre ~ add(add.sd) + prop(prop.sd) + boxCox(lambda)) %>% ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), lambda=1)
  }, .addProp2, addProp = 2)

  testErr("boxCox+add+prop, combined 2->add+prop, combined 2", function(f) {
    f %>% model(ipre ~ add(add.sd) + prop(prop.sd) + boxCox(lambda) + combined2()) %>%
      ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), lambda=1)
  }, .addProp2, addProp = 1)

  ## propT
  testErr("boxCox+add+propT, combined 2->add", function(f) {
    f %>% model(ipre ~ add(add.sd) + propT(prop.sd) + boxCox(lambda)) %>%
      ini(add.sd=sqrt(0.1), prop.sd=0, lambda=1)
  }, .addVals, addProp = 2)

  testErr("boxCox+addMod+propT, combined 2->addMod", function(f) {
    f %>% model(ipre ~ add(f2) + propT(prop.sd) + boxCox(lambda)) %>%
      ini(prop.sd=0, lambda=1)
  }, .addModVals, addProp = 2)

  testErr("boxCox+add+propT, combined 2->propT", function(f) {
    f %>% model(ipre ~ add(add.sd) + propT(prop.sd) + boxCox(lambda)) %>%
      ini(add.sd=0, prop.sd=sqrt(0.1), lambda=1)
  }, addProp = 2)

  testErr("boxCox+add+propTMod, combined 2->propMod", function(f) {
    f %>% model(ipre ~ add(add.sd) + propT(f2) + boxCox(lambda)) %>%
      ini(add.sd=0, lambda=1)
  }, .boxCoxPropTModVals, addProp = 2)

  testErr("boxCox+addMod+propTMod, combined 2->propTMod", function(f) {
    f %>% model(ipre ~ add(f2) + propT(f3) + boxCox(lambda)) %>%
      ini(lambda=1)
  }, .boxCoxAddModPropTMod2, addProp = 2)

  .boxCoxAddPropTc2Vals <- c(35.457, 35.457, 35.224, 35.056, 35.07, 35.07)

  testErr("boxCox+add+propT, combined 2", function(f) {
    f %>% model(ipre ~ add(add.sd) + propT(prop.sd) + boxCox(lambda)) %>%
      ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), lambda=1)
  }, .boxCoxAddPropTc2Vals, addProp = 2)

  testErr("boxCox+add+propT, combined 2 (specified)", function(f) {
    f %>% model(ipre ~ add(add.sd) + propT(prop.sd) + boxCox(lambda) + combined2()) %>%
      ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), lambda=1)
  }, .boxCoxAddPropTc2Vals, addProp = 1)

  # propF
  testErr("boxCox+add+propF, combined 2->add", function(f) {
    f %>% model(ipre ~ add(add.sd) + propF(prop.sd, f2) + boxCox(lambda)) %>% ini(add.sd=sqrt(0.1), prop.sd=0, lambda=1)
  }, .addVals, addProp = 2)

  testErr("boxCox+addMod+propF, combined 2->addMod", function(f) {
    f %>% model(ipre ~ add(f2) + propF(prop.sd, f2) + boxCox(lambda)) %>% ini(prop.sd=0, lambda=1)
  }, .addModVals, addProp = 2)

  testErr("boxCox+add+propFMod, combined 2->propMod", function(f) {
    f %>% model(ipre ~ add(add.sd) + propF(lipre, f2) + boxCox(lambda)) %>% ini(add.sd=0, lambda=1)
  }, .propFModVals, addProp = 2)

  testErr("boxCox+addMod+propFMod, combined 2->propMod", function(f) {
    f %>% model(ipre ~ add(f2) + propF(lipre, f3) + boxCox(lambda)) %>%
      ini(lambda=1)
  }, .addModPropFModVals2, addProp = 2)

  testErr("boxCox+add+propF, combined 2", function(f) {
    f %>% model(ipre ~ add(add.sd) + propF(prop.sd, f2) + boxCox(lambda)) %>%
      ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), lambda=1)
  },  .addPropFVals2, addProp = 2)

  testErr("boxCox+add+propF, combined 2 (specified)", function(f) {
    f %>% model(ipre ~ add(add.sd) + propF(prop.sd, f2) + boxCox(lambda) + combined2()) %>%
      ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), lambda=1)
  }, .addPropFVals2, addProp = 1)

  ################################################################################
  # Box-Box Add+Proportional tests (combined 1)
  ################################################################################
  testErr("boxCox+add+prop, combined 1->add", function(f) {
    f %>% model(ipre ~ add(add.sd) + prop(prop.sd) + boxCox(lambda)) %>%
      ini(add.sd=sqrt(0.1), prop.sd=0, lambda=1)
  }, .addVals, addProp = 1)

  testErr("boxCox+addMod+prop, combined 1->addMod", function(f) {
    f %>% model(ipre ~ add(f2) + prop(prop.sd) + boxCox(lambda)) %>%
      ini(prop.sd=0, lambda=1)
  }, .addModVals, addProp = 1)

  testErr("boxCox+add+prop, combined 1->prop", function(f) {
    f %>% model(ipre ~ add(add.sd) + prop(prop.sd) + boxCox(lambda)) %>%
      ini(add.sd=0, prop.sd=sqrt(0.1), lambda=1)
  }, .propVals, addProp = 1)

  testErr("boxCox+add+propMod, combined 1->propMod", function(f) {
    f %>% model(ipre ~ add(add.sd) + prop(f2) + boxCox(lambda)) %>%
      ini(add.sd=0, lambda=1)
  }, .propModVals, addProp = 1)

  testErr("boxCox+addMod+propMod, combined 1->propMod", function(f) {
    f %>% model(ipre ~ add(f2) + prop(f3) + boxCox(lambda)) %>%
      ini(lambda=1)
  },  .addModPropMod1, addProp = 2)

  testErr("boxCox+add+prop, combined 2->add+prop, combined 2", function(f) {
    f %>% model(ipre ~ add(add.sd) + prop(prop.sd) + boxCox(lambda)) %>%
      ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), lambda=1)
  }, .addProp2, addProp = 2)

  testErr("boxCox+add+prop, combined 2", function(f) {
    f %>% model(ipre ~ add(add.sd) + prop(prop.sd) + boxCox(lambda) + combined2()) %>%
      ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), lambda=1)
  }, .addProp2, addProp = 1)

  ## propT
  testErr("boxCox+add+propT, combined 2->add", function(f) {
    f %>% model(ipre ~ add(add.sd) + propT(prop.sd) + boxCox(lambda)) %>%
      ini(add.sd=sqrt(0.1), prop.sd=0, lambda=1)
  }, .addVals, addProp = 2)

  testErr("boxCox+addMod+propT, combined 2->addMod", function(f) {
    f %>% model(ipre ~ add(f2) + propT(prop.sd) + boxCox(lambda)) %>%
      ini(prop.sd=0, lambda=1)
  }, .addModVals, addProp = 2)

  .boxCoxPropT <- c(35.132, 35.132, 34.879, 34.694, 34.709, 34.709)

  testErr("boxCox+add+propT, combined 2->prop", function(f) {
    f %>% model(ipre ~ add(add.sd) + propT(prop.sd) + boxCox(lambda)) %>%
      ini(add.sd=0, prop.sd=sqrt(0.1), lambda=1)
  }, .boxCoxPropT, addProp = 2)

  testErr("boxCox+add+propTMod, combined 2->propMod", function(f) {
    f %>% model(ipre ~ add(add.sd) + propT(f2) + boxCox(lambda)) %>%
      ini(add.sd=0, lambda=1)
  }, .boxCoxPropT2ModVals, addProp = 2)

  .boxCoxAddPropT3Vals <- c(100.737, 100.737, 100.351, 100.437, 100.437, 100.437)

  testErr("boxCox+addMod+propTMod, combined 2->propMod", function(f) {
    f %>% model(ipre ~ add(f2) + propT(f3) + boxCox(lambda)) %>%
      ini(lambda=1)
  },  .boxCoxAddPropT3Vals, addProp = 2)

  .boxCoxAddPropT2Vals <- c(35.457, 35.457, 35.224, 35.056, 35.07, 35.07)

  testErr("boxCox+add+propT, combined 2 (specified)", function(f) {
    f %>% model(ipre ~ add(add.sd) + propT(prop.sd) + boxCox(lambda) + combined2()) %>%
      ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), lambda=1)
  }, .boxCoxAddPropT2Vals, addProp = 1)

  # propF
  testErr("boxCox+add+propF, combined 2->add", function(f) {
    f %>% model(ipre ~ add(add.sd) + propF(prop.sd, f2) + boxCox(lambda)) %>%
      ini(add.sd=sqrt(0.1), prop.sd=0, lambda=1)
  }, .addVals, addProp = 2)

  testErr("boxCox+addMod+propF, combined 2->addMod", function(f) {
    f %>% model(ipre ~ add(f2) + propF(prop.sd, f2) + boxCox(lambda)) %>%
      ini(prop.sd=0, lambda=1)
  }, .addModVals, addProp = 2)

  testErr("boxCox+add+propF, combined 2->prop", function(f) {
    f %>% model(ipre ~ add(add.sd) + propF(prop.sd, f2) + boxCox(lambda)) %>%
      ini(add.sd=0, prop.sd=sqrt(0.1), lambda=1)
  }, .propFVals, addProp = 2)

  testErr("boxCox+add+propFMod, combined 2->propMod", function(f) {
    f %>% model(ipre ~ add(add.sd) + propF(lipre, f2) + boxCox(lambda)) %>%
      ini(add.sd=0, lambda=1)
  }, .propFModVals, addProp = 2)

  testErr("boxCox+addMod+propFMod, combined 2->propMod", function(f) {
    f %>% model(ipre ~ add(f2) + propF(lipre, f3) + boxCox(lambda)) %>%
      ini(lambda=1)
  }, .addModPropFModVals2, addProp = 2)

  testErr("boxCox+add+propF, combined 2", function(f) {
    f %>% model(ipre ~ add(add.sd) + propF(prop.sd, f2) + boxCox(lambda)) %>%
      ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), lambda=1)
  },  .addPropFVals2, addProp = 2)

  testErr("boxCox+add+propF, combined 2 (specified)", function(f) {
    f %>% model(ipre ~ add(add.sd) + propF(prop.sd, f2) + boxCox(lambda) + combined2()) %>%
      ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), lambda=1)
  }, .addPropFVals2, addProp = 1)

  #################################################################################################
  testErr("boxCox+add+prop, combined 1", function(f) {
    f %>% model(ipre ~ add(add.sd) + prop(prop.sd) + boxCox(lambda)) %>%
      ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), lambda=1)
  }, .addProp1, addProp = 1)

  testErr("boxCox+add+prop, combined 1->add", function(f) {
    f %>% model(ipre ~ add(add.sd) + prop(prop.sd) + boxCox(lambda)) %>%
      ini(add.sd=sqrt(0.1), prop.sd=0, lambda=1)
  }, .addVals, addProp = 1)

  testErr("boxCox+add+prop, combined 1->prop", function(f) {
    f %>% model(ipre ~ add(add.sd) + prop(prop.sd) + boxCox(lambda)) %>%
      ini(add.sd=0, prop.sd=sqrt(0.1), lambda=1)
  }, .propVals, addProp = 1)

  testErr("boxCox+add+prop, combined 1 (override)", function(f) {
    f %>% model(ipre ~ add(add.sd) + prop(prop.sd) + boxCox(lambda) + combined1()) %>%
      ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), lambda=1)
  }, .addProp1, addProp = 2)

  testErr("boxCox+add+prop, combined 2 (override)", function(f) {
    f %>% model(ipre ~ add(add.sd) + prop(prop.sd) + boxCox(lambda) + combined2()) %>%
      ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), lambda=1)
  }, .addProp2, addProp = 1)

  ################################################################################
  # Add+Pow tests (combined 2)
  ################################################################################
  testErr("boxCox+add+pow combined 2 -> add+prop combined2", function(f) {
    f %>% model(ipre ~ add(add.sd) + pow(prop.sd, pw) + boxCox(lambda)) %>%
      ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=1, lambda=1)
  }, .addProp2, addProp = 2)

  testErr("boxCox+add+pow combined 2 -> add+pow combined 2", function(f) {
    f %>% model(ipre ~ add(add.sd) + pow(prop.sd, pw) + boxCox(lambda)) %>%
      ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5, lambda=1)
  }, c(10.886, 10.886, 10.868, 10.417, 10.662, 10.662), addProp = 2)

  testErr("boxCox+add+pow combined 1 -> add+prop combined1", function(f) {
    f %>% model(ipre ~ add(add.sd) + pow(prop.sd, pw) + boxCox(lambda)) %>%
      ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=1, lambda=1)
  }, .addProp1, addProp = 1)

  testErr("boxCox+add+pow combined 1", function(f) {
    f %>% model(ipre ~ add(add.sd) + pow(prop.sd, pw) + boxCox(lambda)) %>%
      ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5, lambda=1)
  }, .addPow1, addProp = 1)

  testErr("boxCox+add+pow combined 1 (override)->add+pow combined 1 (override)", function(f) {
    f %>% model(ipre ~ add(add.sd) + pow(prop.sd, pw) + boxCox(lambda) + combined1()) %>%
      ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5, lambda=1)
  }, .addPow1, addProp = 2)


  ################################################################################
  ## BoxCox(0) and lnorm equivalence
  ################################################################################

  testErr("boxCox(0)+add-> lnorm", function(f) {
    f %>% model(ipre ~ add(add.sd) + boxCox(lambda)) %>%
      ini(add.sd=sqrt(0.1), lambda=0)
  }, .lnorm, addProp = 1)

  ################################################################################
  ## boxCox(0)+prop -> lnorm(NA) tests
  ################################################################################

  testErr("boxCox(0)+prop->lnorm(NA)+prop", function(f) {
    f %>% model(ipre ~ boxCox(lambda) + prop(prop.sd)) %>% ini(prop.sd=sqrt(0.1), lambda=0)
  }, .lnormProp, addProp = 1)

  testErr("boxCox(0)+propT->lnorm(NA)+propT", function(f) {
    f %>% model(ipre ~ boxCox(lambda) + propT(prop.sd)) %>%
      ini(prop.sd=sqrt(0.1), lambda=0)
  }, .lnormPropT, addProp = 1)

  testErr("boxCox(0)+propF->lnorm(NA)+propF", function(f) {
    f %>% model(ipre ~ boxCox(lm) + propF(prop.sd, f2)) %>%
      ini(prop.sd=sqrt(0.1), lm=0)
  }, .lnormPropF, addProp = 1)

  testErr("boxCox(0)+pow->lnorm(NA)+prop", function(f) {
    f %>% model(ipre ~ boxCox(lm) + pow(prop.sd, pw)) %>%
      ini(prop.sd=sqrt(0.1), pw=1, lm=0)
  }, .lnormProp, addProp = 1)

  testErr("boxCox(0)+powT->lnorm(NA)+propT", function(f) {
    f %>% model(ipre ~ boxCox(lm) + powT(prop.sd, pw)) %>%
      ini(prop.sd=sqrt(0.1), pw=1, lm=0)
  }, .lnormPropT, addProp = 1)

  testErr("boxCox(0)+powF->lnorm(NA)+propF", function(f) {
    f %>% model(ipre ~ boxCox(lm) + powF(prop.sd, pw, f2)) %>%
      ini(prop.sd=sqrt(0.1), pw=1, lm=0)
  }, .lnormPropF, addProp = 1)

  testErr("boxCox(0)+pow->lnorm(NA)+pow", function(f) {
    f %>% model(ipre ~ boxCox(lm) + pow(prop.sd, pw)) %>%
      ini(prop.sd=sqrt(0.1), pw=0.5, lm=0)
  }, .lnormPow, addProp = 1)

  testErr("boxCox(0)+powT->lnorm(NA)+powT", function(f) {
    f %>% model(ipre ~ boxCox(lm) + powT(prop.sd, pw)) %>%
      ini(prop.sd=sqrt(0.1), pw=0.5, lm=0)
  }, .lnormPowT, addProp = 1)

  testErr("boxCox(0)+powF->lnorm(NA)+powF", function(f) {
    f %>% model(ipre ~ boxCox(lm) + powF(prop.sd, pw, f2)) %>%
      ini(prop.sd=sqrt(0.1), pw=0.5, lm=0)
  }, .lnormPowF, addProp = 1)


  ################################################################################
  ## lnorm combined1
  ################################################################################

  testErr("boxCox(0)+add+prop combined1->lnorm(NA)+prop", function(f) {
    f %>% model(ipre ~ add(lnorm.sd) + prop(prop.sd) + boxCox(lm)) %>%
      ini(lnorm.sd=0, prop.sd=sqrt(0.1), lm=0)
  }, .lnormProp, addProp = 1)

  testErr("boxCox(0)+add+propT combined1->lnorm(NA)+propT", function(f) {
    f %>% model(ipre ~ add(lnorm.sd) + propT(prop.sd) + boxCox(lm)) %>%
      ini(lnorm.sd=0, prop.sd=sqrt(0.1), lm=0)
  }, .lnormPropT, addProp = 1)

  testErr("boxCox(0)+add+propF combined1->lnorm(NA)+propF", function(f) {
    f %>% model(ipre ~ add(lnorm.sd) + propF(prop.sd, f2) + boxCox(lm)) %>%
      ini(lnorm.sd=0, prop.sd=sqrt(0.1), lm=0)
  }, .lnormPropF, addProp = 1)

  testErr("boxCox(0)+add+prop combined1->lnorm(NA)", function(f) {
    f %>% model(ipre ~ add(lnorm.sd) + prop(prop.sd) + boxCox(lm)) %>%
      ini(lnorm.sd=sqrt(0.1), prop.sd=0, lm=0)
  }, .lnorm, addProp = 1)

  testErr("boxCox(0)+add+propT combined1->lnorm(NA)", function(f) {
    f %>% model(ipre ~ add(lnorm.sd) + propT(prop.sd) + boxCox(lm)) %>%
      ini(lnorm.sd=sqrt(0.1), prop.sd=0, lm=0)
  }, .lnorm, addProp = 1)

  testErr("boxCox(0)+add+propF combined1->lnorm(NA)", function(f) {
    f %>% model(ipre ~ add(lnorm.sd) + propF(prop.sd, f2) + boxCox(lm)) %>%
      ini(lnorm.sd=sqrt(0.1), prop.sd=0, lm=0)
  }, .lnorm, addProp = 1)

  testErr("boxCox(0)+add+prop combined1", function(f) {
    f %>% model(ipre ~ add(lnorm.sd) + prop(prop.sd) + boxCox(lm)) %>%
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), lm=0)
  }, .lnormProp1, addProp = 1)

  testErr("boxCox(0)+add+propT combined1", function(f) {
    f %>% model(ipre ~ add(lnorm.sd) + propT(prop.sd) + boxCox(lm)) %>%
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), lm=0)
  }, .lnormPropT1, addProp = 1)

  testErr("boxCox(0)+add+propF combined1", function(f) {
    f %>% model(ipre ~ add(lnorm.sd) + propF(prop.sd, f2) + boxCox(lm)) %>%
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), lm=0)
  }, .lnormPropF1, addProp = 1)

  testErr("boxCox(0)+add+powF->lnorm+propF combined1", function(f) {
    f %>% model(ipre ~ add(lnorm.sd) + powF(prop.sd, pw, f2) + boxCox(lm)) %>%
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=1, lm=0)
  }, .lnormPropF1, addProp = 1)

  testErr("boxCox(0)+add+pow->lnorm+prop combined1", function(f) {
    f %>% model(ipre ~ add(lnorm.sd) + pow(prop.sd, pw) + boxCox(lm)) %>%
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=1, lm=0)
  }, .lnormProp1, addProp = 1)

  testErr("boxCox(0)+add+powT->lnorm+propT combined1", function(f) {
    f %>% model(ipre ~ add(lnorm.sd) + powT(prop.sd, pw) + boxCox(lm)) %>%
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=1, lm=0)
  }, .lnormPropT1, addProp = 1)

  testErr("boxCox(0)+add+powF combined1", function(f) {
    f %>% model(ipre ~ add(lnorm.sd) + powF(prop.sd, pw, f2) + boxCox(lm)) %>%
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5, lm=0)
  }, .lnormPowF1, addProp = 1)

  testErr("boxCox(0)+add+pow combined1", function(f) {
    f %>% model(ipre ~ add(lnorm.sd) + pow(prop.sd, pw) + boxCox(lm)) %>%
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5, lm=0)
  }, .lnormPow1, addProp = 1)

  testErr("boxCox(0)+add+powT combined1", function(f) {
    f %>% model(ipre ~ add(lnorm.sd) + powT(prop.sd, pw) + boxCox(lm)) %>%
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5, lm=0)
  }, .lnormPowT1, addProp = 1)

  ################################################################################
  ## lnorm combined2
  ################################################################################

  testErr("boxCox(0)+add+propT combined2->lnorm(NA)+propT", function(f) {
    f %>% model(ipre ~ add(add.sd) + propT(prop.sd) + boxCox(lambda)) %>%
      ini(add.sd=0, prop.sd=sqrt(0.1), lambda=0)
  }, .lnormPropT, addProp = 2)

  testErr("boxCox(0)+add+propF combined2->lnorm(NA)+propF", function(f) {
    f %>% model(ipre ~ add(lnorm.sd) + boxCox(lm)+ propF(prop.sd, f2)) %>%
      ini(lnorm.sd=0, prop.sd=sqrt(0.1), lm=0)
  }, .lnormPropF, addProp = 2)

  testErr("boxCox(0)+add+prop combined2->lnorm(NA)", function(f) {
    f %>% model(ipre ~ add(lnorm.sd) + prop(prop.sd) + boxCox(lm)) %>%
      ini(lnorm.sd=sqrt(0.1), prop.sd=0, lm=0)
  }, .lnorm, addProp = 2)

  testErr("boxCox(0)+add+propT combined2->lnorm(NA)", function(f) {
    f %>% model(ipre ~ boxCox(lm) + add(lnorm.sd) + propT(prop.sd)) %>%
      ini(lnorm.sd=sqrt(0.1), prop.sd=0, lm=0)
  }, .lnorm, addProp = 2)

  testErr("boxCox(0)+add+propF combined2->lnorm(NA)", function(f) {
    f %>% model(ipre ~ add(lnorm.sd) + boxCox(lm) + propF(prop.sd, f2)) %>%
      ini(lnorm.sd=sqrt(0.1), prop.sd=0, lm=0)
  }, .lnorm, addProp = 2)

  testErr("boxCox(0)+add+pow->lnorm+prop combined2", function(f) {
    f %>% model(ipre ~ boxCox(lm) + add(lnorm.sd) + pow(prop.sd, pw)) %>%
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=1, lm=0)
  }, .lnormProp2, addProp = 2)

  testErr("boxCox(0)+add+prop combined2", function(f) {
    f %>% model(ipre ~ add(lnorm.sd) + prop(prop.sd) + boxCox(lm)) %>%
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), lm=0)
  }, .lnormProp2 ,addProp = 2)

  testErr("boxCox(0)+add+propT combined2", function(f) {
    f %>% model(ipre ~ add(lnorm.sd) + propT(prop.sd) + boxCox(lm)) %>%
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), lm=0)
  }, .lnormPropT2, addProp = 2)

  testErr("boxCox(0)+add+propF combined1", function(f) {
    f %>% model(ipre ~ add(lnorm.sd) + boxCox(lm)+ propF(prop.sd, f2)) %>%
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), lm=0)
  }, .lnormPropF2, addProp = 2)

  testErr("boxCox(0)+add+powF->lnorm+propF combined2", function(f) {
    f %>% model(ipre ~ add(lnorm.sd) + powF(prop.sd, pw, f2) + boxCox(lm)) %>%
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=1, lm=0)
  }, .lnormPropF2, addProp = 2)


  testErr("boxCox(0)+add+pow->lnorm+prop combined2", function(f) {
    f %>% model(ipre ~ add(lnorm.sd) + boxCox(lm)+ pow(prop.sd, pw)) %>%
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=1, lm=0)
  }, .lnormProp2, addProp = 2)

  testErr("boxCox(0)+add+powT->lnorm+propT combined2", function(f) {
    f %>% model(ipre ~ add(lnorm.sd) + powT(prop.sd, pw) + boxCox(lm)) %>%
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=1, lm=0)
  }, .lnormPropT2, addProp = 2)

  testErr("boxCox(0)+add+powF combined2", function(f) {
    f %>% model(ipre ~ boxCox(lm) + add(lnorm.sd) + powF(prop.sd, pw, f2)) %>%
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5, lm=0)
  }, .lnormPowF2, addProp = 2)

  testErr("boxCox(0)+add+pow combined2", function(f) {
    f %>% model(ipre ~ boxCox(lm) + add(lnorm.sd) + pow(prop.sd, pw)) %>%
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5, lm=0)
  }, .lnormPow2, addProp = 2)

  testErr("boxCox(0)+add+powT combined2", function(f) {
    f %>% model(ipre ~ boxCox(lm) + add(lnorm.sd) + powT(prop.sd, pw)) %>%
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5, lm=0)
  }, .lnormPowT2, addProp = 2)


  ## Test lambda=0.5

  ################################################################################
  ## BoxCox(0.5)
  ################################################################################

  .boxCoxAdd <- c(10.402, 10.402, 10.402, 10.402, 11.078, 11.078)

  testErr("boxCox(0.5)+add", function(f) {
    f %>% model(ipre ~ add(add.sd) + boxCox(lambda)) %>%
      ini(add.sd=sqrt(0.1), lambda=0.5)
  }, .boxCoxAdd, addProp = 1)

  ################################################################################
  ## boxCox(0.5)+prop
  ################################################################################

  .boxCoxProp <- c(77.857, 77.857, 77.698, 77.732, 77.732, 77.732)
  testErr("boxCox(0.5)+prop", function(f) {
    f %>% model(ipre ~ boxCox(lambda) + prop(prop.sd)) %>% ini(prop.sd=sqrt(0.1), lambda=0.5)
  }, .boxCoxProp, addProp = 1)

  .boxCoxPropT <- c(48.455, 48.455, 48.321, 48.287, 48.291, 48.291)
  testErr("boxCox(0.5)+propT", function(f) {
    f %>% model(ipre ~ boxCox(lambda) + propT(prop.sd)) %>%
      ini(prop.sd=sqrt(0.1), lambda=0.5)
  }, .boxCoxPropT, addProp = 1)

  .boxCoxPropF <- c(4.243, 4.243, 4.216, 3.537, 6.584, 6.584)
  testErr("boxCox(0.5)+propF", function(f) {
    f %>% model(ipre ~ boxCox(lm) + propF(prop.sd, f2)) %>%
      ini(prop.sd=sqrt(0.1), lm=0.5)
  }, .boxCoxPropF, addProp = 1)

  testErr("boxCox(0.5)+pow->boxCox(0.5)+prop", function(f) {
    f %>% model(ipre ~ boxCox(lm) + pow(prop.sd, pw)) %>%
      ini(prop.sd=sqrt(0.1), pw=1, lm=0.5)
  }, .boxCoxProp, addProp = 1)

  testErr("boxCox(0.5)+powT->boxCox(0.5)+propT", function(f) {
    f %>% model(ipre ~ boxCox(lm) + powT(prop.sd, pw)) %>%
      ini(prop.sd=sqrt(0.1), pw=1, lm=0.5)
  }, .boxCoxPropT, addProp = 1)

  testErr("boxCox(0.5)+powF->boxCox(0.5)+propF", function(f) {
    f %>% model(ipre ~ boxCox(lm) + powF(prop.sd, pw, f2)) %>%
      ini(prop.sd=sqrt(0.1), pw=1, lm=0.5)
  }, .boxCoxPropF, addProp = 1)

  .boxCoxPow <- c(39.674, 39.674, 39.629, 39.563, 39.573, 39.573)
  testErr("boxCox(0.5)+pow", function(f) {
    f %>% model(ipre ~ boxCox(lm) + pow(prop.sd, pw)) %>%
      ini(prop.sd=sqrt(0.1), pw=0.5, lm=0.5)
  }, .boxCoxPow, addProp = 1)

  .boxCoxPowT <- c(27.311, 27.311, 27.296, 27.082, 27.148, 27.148)
  testErr("boxCox(0.5)+powT", function(f) {
    f %>% model(ipre ~ boxCox(lm) + powT(prop.sd, pw)) %>%
      ini(prop.sd=sqrt(0.1), pw=0.5, lm=0.5)
  }, .boxCoxPowT, addProp = 1)

  .boxCoxPowF <- c(7.046, 7.046, 7.032, 6.617, 8.287, 8.287)
  testErr("boxCox(0.5)+powF", function(f) {
    f %>% model(ipre ~ boxCox(lm) + powF(prop.sd, pw, f2)) %>%
      ini(prop.sd=sqrt(0.1), pw=0.5, lm=0.5)
  }, .boxCoxPowF, addProp = 1)


  ################################################################################
  ## lnorm combined1
  ################################################################################

  testErr("boxCox(0.5)+add+prop combined1->boxCox(0.5)+prop", function(f) {
    f %>% model(ipre ~ add(lnorm.sd) + prop(prop.sd) + boxCox(lm)) %>%
      ini(lnorm.sd=0, prop.sd=sqrt(0.1), lm=0.5)
  }, .boxCoxProp, addProp = 1)

  testErr("boxCox(0.5)+add+propT combined1->boxCox(0.5)+propT", function(f) {
    f %>% model(ipre ~ add(lnorm.sd) + propT(prop.sd) + boxCox(lm)) %>%
      ini(lnorm.sd=0, prop.sd=sqrt(0.1), lm=0.5)
  }, .boxCoxPropT, addProp = 1)

  testErr("boxCox(0.5)+add+propF combined1->boxCox(0.5)+propF", function(f) {
    f %>% model(ipre ~ add(lnorm.sd) + propF(prop.sd, f2) + boxCox(lm)) %>%
      ini(lnorm.sd=0, prop.sd=sqrt(0.1), lm=0.5)
  }, .boxCoxPropF, addProp = 1)

  testErr("boxCox(0.5)+add+prop combined1->boxCox(0.5)+add", function(f) {
    f %>% model(ipre ~ add(lnorm.sd) + prop(prop.sd) + boxCox(lm)) %>%
      ini(lnorm.sd=sqrt(0.1), prop.sd=0, lm=0.5)
  }, .boxCoxAdd, addProp = 1)

  testErr("boxCox(0.5)+add+propT combined1->boxCox(0.5)+add", function(f) {
    f %>% model(ipre ~ add(lnorm.sd) + propT(prop.sd) + boxCox(lm)) %>%
      ini(lnorm.sd=sqrt(0.1), prop.sd=0, lm=0.5)
  }, .boxCoxAdd, addProp = 1)

  testErr("boxCox(0.5)+add+propF combined1->boxCox(0.5)+add", function(f) {
    f %>% model(ipre ~ add(lnorm.sd) + propF(prop.sd, f2) + boxCox(lm)) %>%
      ini(lnorm.sd=sqrt(0.1), prop.sd=0, lm=0.5)
  }, .boxCoxAdd, addProp = 1)

  .boxCoxAddProp1 <- c(82.622, 82.622, 82.508, 82.534, 82.534, 82.534)
  testErr("boxCox(0.5)+add+prop combined1", function(f) {
    f %>% model(ipre ~ add(lnorm.sd) + prop(prop.sd) + boxCox(lm)) %>%
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), lm=0.5)
  }, .boxCoxAddProp1, addProp = 1)

  .boxCoxAddPropT1 <- c(57.327, 57.327, 57.252, 57.256, 57.257, 57.257)
  testErr("boxCox(0.5)+add+propT combined1", function(f) {
    f %>% model(ipre ~ add(lnorm.sd) + propT(prop.sd) + boxCox(lm)) %>%
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), lm=0.5)
  }, .boxCoxAddPropT1, addProp = 1)

  .boxCoxAddPropF1 <- c(22.065, 22.065, 22.064, 21.96, 22.067, 22.067)
  testErr("boxCox(0.5)+add+propF combined1", function(f) {
    f %>% model(ipre ~ add(lnorm.sd) + propF(prop.sd, f2) + boxCox(lm)) %>%
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), lm=0.5)
  }, .boxCoxAddPropF1, addProp = 1)

  testErr("boxCox(0.5)+add+powF->boxCox(0.5)+add+propF combined1", function(f) {
    f %>% model(ipre ~ add(lnorm.sd) + powF(prop.sd, pw, f2) + boxCox(lm)) %>%
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=1, lm=0.5)
  }, .boxCoxAddPropF1, addProp = 1)

  testErr("boxCox(0.5)+add+pow->boxCox(0.5)+add+prop combined1", function(f) {
    f %>% model(ipre ~ add(lnorm.sd) + pow(prop.sd, pw) + boxCox(lm)) %>%
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=1, lm=0.5)
  }, .boxCoxAddProp1, addProp = 1)

  testErr("boxCox(0.5)+add+powT->boxCox(0.5)+add+propT combined1", function(f) {
    f %>% model(ipre ~ add(lnorm.sd) + powT(prop.sd, pw) + boxCox(lm)) %>%
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=1, lm=0.5)
  }, .boxCoxAddPropT1, addProp = 1)

  .boxCoxAddPowF1 <- c(24.676, 24.676, 24.676, 24.631, 24.692, 24.692)
  testErr("boxCox(0.5)+add+powF combined1", function(f) {
    f %>% model(ipre ~ add(lnorm.sd) + powF(prop.sd, pw, f2) + boxCox(lm)) %>%
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5, lm=0.5)
  }, .boxCoxAddPowF1, addProp = 1)


  .boxCoxAddPow1 <- c(50.265, 50.265, 50.241, 50.23, 50.232, 50.232)
  testErr("boxCox(0.5)+add+pow combined1", function(f) {
    f %>% model(ipre ~ add(lnorm.sd) + pow(prop.sd, pw) + boxCox(lm)) %>%
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5, lm=0.5)
  }, .boxCoxAddPow1, addProp = 1)

  .boxCoxAddPowT1 <- c(40.45, 40.45, 40.436, 40.41, 40.416, 40.416)
  testErr("boxCox(0.5)+add+powT combined1", function(f) {
    f %>% model(ipre ~ add(lnorm.sd) + powT(prop.sd, pw) + boxCox(lm)) %>%
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5, lm=0.5)
  }, .boxCoxAddPowT1, addProp = 1)

  ################################################################################
  ## boxCox(0.5) combined2
  ################################################################################

  testErr("boxCox(0.5)+add+propT combined2->boxCox(0.5)+propT", function(f) {
    f %>% model(ipre ~ add(add.sd) + propT(prop.sd) + boxCox(lambda)) %>%
      ini(add.sd=0, prop.sd=sqrt(0.1), lambda=0.5)
  }, .boxCoxPropT, addProp = 2)

  testErr("boxCox(0.5)+add+propF combined2->boxCox(0.5)+propF", function(f) {
    f %>% model(ipre ~ add(lnorm.sd) + boxCox(lm)+ propF(prop.sd, f2)) %>%
      ini(lnorm.sd=0, prop.sd=sqrt(0.1), lm=0.5)
  }, .boxCoxPropF, addProp = 2)

  testErr("boxCox(0.5)+add+prop combined2->boxCox(0.5)+add", function(f) {
    f %>% model(ipre ~ add(lnorm.sd) + prop(prop.sd) + boxCox(lm)) %>%
      ini(lnorm.sd=sqrt(0.1), prop.sd=0, lm=0.5)
  }, .boxCoxAdd, addProp = 2)

  testErr("boxCox(0.5)+add+propT combined2->boxCox(0.5)+add", function(f) {
    f %>% model(ipre ~ boxCox(lm) + add(lnorm.sd) + propT(prop.sd)) %>%
      ini(lnorm.sd=sqrt(0.1), prop.sd=0, lm=0.5)
  }, .boxCoxAdd, addProp = 2)

  testErr("boxCox(0.5)+add+propF combined2->boxCox(0.5)+add", function(f) {
    f %>% model(ipre ~ add(lnorm.sd) + boxCox(lm) + propF(prop.sd, f2)) %>%
      ini(lnorm.sd=sqrt(0.1), prop.sd=0, lm=0.5)
  }, .boxCoxAdd, addProp = 2)

  .boxCoxAddProp2 <- c(78.202, 78.202, 78.052, 78.084, 78.084, 78.084)
  testErr("boxCox(0.5)+add+pow->boxCox(0.5)+add+prop combined2", function(f) {
    f %>% model(ipre ~ boxCox(lm) + add(lnorm.sd) + pow(prop.sd, pw)) %>%
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=1, lm=0.5)
  }, .boxCoxAddProp2, addProp = 2)

  testErr("boxCox(0.5)+add+prop combined2", function(f) {
    f %>% model(ipre ~ add(lnorm.sd) + prop(prop.sd) + boxCox(lm)) %>%
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), lm=0.5)
  }, .boxCoxAddProp2 ,addProp = 2)


  .boxCoxAddPropT2 <- c(49.798, 49.798, 49.689, 49.668, 49.671, 49.671)
  testErr("boxCox(0.5)+add+propT combined2", function(f) {
    f %>% model(ipre ~ add(lnorm.sd) + propT(prop.sd) + boxCox(lm)) %>%
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), lm=0.5)
  }, .boxCoxAddPropT2, addProp = 2)

  .boxCoxAddPropF2 <- c(14.213, 14.213, 14.21, 14.085, 14.498, 14.498)
  testErr("boxCox(0.5)+add+propF combined1", function(f) {
    f %>% model(ipre ~ add(lnorm.sd) + boxCox(lm)+ propF(prop.sd, f2)) %>%
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), lm=0.5)
  }, .boxCoxAddPropF2, addProp = 2)

  testErr("boxCox(0.5)+add+powF->lnorm+propF combined2", function(f) {
    f %>% model(ipre ~ add(lnorm.sd) + powF(prop.sd, pw, f2) + boxCox(lm)) %>%
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=1, lm=0.5)
  }, .boxCoxAddPropF2, addProp = 2)

  .boxCoxAddPropF2 <- c(78.202, 78.202, 78.052, 78.084, 78.084, 78.084)
  testErr("boxCox(0.5)+add+pow->lnorm+prop combined2", function(f) {
    f %>% model(ipre ~ add(lnorm.sd) + boxCox(lm)+ pow(prop.sd, pw)) %>%
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=1, lm=0.5)
  }, .boxCoxAddPropF2, addProp = 2)

  testErr("boxCox(0.5)+add+powT->boxCox(0.5)+add+propT combined2", function(f) {
    f %>% model(ipre ~ add(lnorm.sd) + powT(prop.sd, pw) + boxCox(lm)) %>%
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=1, lm=0.5)
  }, .boxCoxAddPropT2, addProp = 2)

  .boxCoxAddPowF2 <- c(15.857, 15.857, 15.856, 15.771, 16.06, 16.06)
  testErr("boxCox(0.5)+add+powF combined2", function(f) {
    f %>% model(ipre ~ boxCox(lm) + add(lnorm.sd) + powF(prop.sd, pw, f2)) %>%
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5, lm=0.5)
  }, .boxCoxAddPowF2, addProp = 2)

  .boxCoxAddPow2 <- c(41.665, 41.665, 41.631, 41.589, 41.596, 41.596)
  testErr("boxCox(0.5)+add+pow combined2", function(f) {
    f %>% model(ipre ~ boxCox(lm) + add(lnorm.sd) + pow(prop.sd, pw)) %>%
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5, lm=0.5)
  }, .boxCoxAddPow2, addProp = 2)

  .boxCoxAddPowT2 <- c(30.736, 30.736, 30.722, 30.625, 30.657, 30.657)
  testErr("boxCox(0.5)+add+powT combined2", function(f) {
    f %>% model(ipre ~ boxCox(lm) + add(lnorm.sd) + powT(prop.sd, pw)) %>%
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5, lm=0.5)
  }, .boxCoxAddPowT2, addProp = 2)


  # Now yeoJohnson
  .yeoJohnsonAdd <- c(11.279, 11.279, 11.279, 11.279, 11.738, 11.738)
  testErr("add+yeoJohnson", function(f) {
    f %>% model(ipre ~ add(add.sd) + yeoJohnson(lm)) %>% ini(add.sd=sqrt(0.1), lm=0.5)
  })

  .yeoJohnsonProp <- c(80.257, 80.257, 80.102, 80.136, 80.136, 80.136)
  testErr("prop+yeoJohnson", function(f) {
    f %>% model(ipre ~ prop(prop.sd) + yeoJohnson(lm)) %>% ini(prop.sd=sqrt(0.1), lm=0.5)
  })

  .yeoJohnsonPow <- c(41.644, 41.644, 41.598, 41.552, 41.559, 41.559)
  testErr("pow+yeoJohnson", function(f) {
    f %>% model(ipre ~ pow(prop.sd, pw) + yeoJohnson(lm)) %>% ini(prop.sd=sqrt(0.1), lm=0.5, pw=0.5)
  }, c(41.644, 41.644, 41.598, 41.552, 41.559, 41.559))

  .yeoJohnsonAddProp1 <- c(85.048, 85.048, 84.937, 84.961, 84.962, 84.962)
  testErr("add+prop+yeoJohnson", function(f) {
    f %>% model(ipre ~ add(add.sd) + prop(prop.sd) + yeoJohnson(lm)) %>%
      ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), lm=0.5)
  }, .yeoJohnsonAddProp1, addProp = 1)

  .yeoJohnsonAddProp2 <- c(80.605, 80.605, 80.458, 80.49, 80.491, 80.491)
  testErr("add+prop+yeoJohnson", function(f) {
    f %>% model(ipre ~ add(add.sd) + prop(prop.sd) + yeoJohnson(lm)) %>%
      ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), lm=0.5)
  }, .yeoJohnsonAddProp2, addProp = 2)

  .yeoJohnsonAddPow1 <- c(52.481, 52.481, 52.456, 52.451, 52.452, 52.452)
  testErr("add+pow+yeoJohnson", function(f) {
    f %>% model(ipre ~ add(add.sd) + pow(prop.sd, pw) + yeoJohnson(lm)) %>%
      ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5, lm=0.5)
  }, .yeoJohnsonAddPow1, addProp = 1)

  .yeoJohnsonAddPow2 <- c(43.704, 43.704, 43.669, 43.641, 43.645, 43.645)
  testErr("add+pow+yeoJohnson", function(f) {
    f %>% model(ipre ~ add(add.sd) + pow(prop.sd, pw) + yeoJohnson(lm)) %>%
      ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5, lm=0.5)
  }, .yeoJohnsonAddPow2, addProp = 2)

  ## logitNorm
  .logitNormAdd <- c(0.612, 0.612, 0.612, 0.612, 0.786, 0.786)
  testErr("logitNorm", function(f) {
    f %>% model(ipre ~ logitNorm(logit.sd, 0, 12)) %>%
      ini(logit.sd=sqrt(0.1))
  }, .logitNormAdd)

  .logitNormProp <- c(67.882, 67.882, 67.731, 67.765, 67.765, 67.765)
  testErr("logitNorm(NA)+prop", function(f) {
    f %>% model(ipre ~ logitNorm(NA, 0, 12) + prop(prop.sd)) %>%
      ini(prop.sd=sqrt(0.1))
  }, .logitNormProp)

  testErr("logitNorm(NA)+pow->logitNorm(NA)+prop", function(f) {
    f %>% model(ipre ~ logitNorm(NA, 0, 12) + pow(prop.sd, pw)) %>%
      ini(prop.sd=sqrt(0.1), pw=1)
  }, .logitNormProp)

  .logitNormPow <- c(29.055, 29.055, 29.007, 28.987, 28.989, 28.989)
  testErr("logitNorm(NA)+pow->logitNorm(NA)+prop", function(f) {
    f %>% model(ipre ~ logitNorm(NA, 0, 12) + pow(prop.sd, pw)) %>%
      ini(prop.sd=sqrt(0.1), pw=0.5)
  }, .logitNormPow)

  .logitNormAddProp1 <- c(72.699, 72.699, 72.591, 72.615, 72.615, 72.615)
  testErr("logitNorm+prop", function(f) {
    f %>% model(ipre ~ logitNorm(logit.sd, 0, 12) + prop(prop.sd)) %>%
      ini(logit.sd=sqrt(0.1), prop.sd=sqrt(0.1))
  }, .logitNormAddProp1, addProp = 1)

  testErr("logitNorm+pow->logitNorm+prop", function(f) {
    f %>% model(ipre ~ logitNorm(logit.sd, 0, 12) + pow(prop.sd, pw)) %>%
      ini(logit.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=1)
  }, .logitNormAddProp1, addProp = 1)

  .logitNormAddPow1 <- c(40.053, 40.053, 40.028, 40.028, 40.029, 40.029)
  testErr("logitNorm+prop", function(f) {
    f %>% model(ipre ~ logitNorm(logit.sd, 0, 12) + pow(prop.sd, pw)) %>%
      ini(logit.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5)
  }, .logitNormAddPow1, addProp = 1)

  ## logitNorm + yeoJohnson
  .logitNormAddYeoJohnson <- c(9.019, 9.019, 9.019, 9.019, 9.576, 9.576)
  testErr("logitNorm+yeoJohnson", function(f) {
    f %>% model(ipre ~ logitNorm(logit.sd, 0, 12) + yeoJohnson(lm)) %>%
      ini(logit.sd=sqrt(0.1), lm=0.5)
  }, .logitNormAddYeoJohnson)

  .logitNormPropAddYeoJohnson <- c(78.136, 78.136, 77.984, 78.017, 78.017, 78.017)
  testErr("logitNorm(NA)+prop+yeoJohnson", function(f) {
    f %>% model(ipre ~ logitNorm(NA, 0, 12) + prop(prop.sd) + yeoJohnson(lm)) %>%
      ini(prop.sd=sqrt(0.1), lm=0.5)
  }, .logitNormPropAddYeoJohnson)

  testErr("logitNorm(NA)+pow+yeoJohnson->logitNorm(NA)+prop+yeoJohnson", function(f) {
    f %>% model(ipre ~ logitNorm(NA, 0, 12) + pow(prop.sd, pw) + yeoJohnson(lm)) %>%
      ini(prop.sd=sqrt(0.1), lm=0.5, pw=1)
  }, .logitNormPropAddYeoJohnson)

  .logitNormPowYeoJohnson <- c(39.334, 39.334, 39.287, 39.251, 39.256, 39.256)
  testErr("logitNorm(NA)+pow+yeoJohnson", function(f) {
    f %>% model(ipre ~ logitNorm(NA, 0, 12) + pow(prop.sd, pw) + yeoJohnson(lm)) %>%
      ini(prop.sd=sqrt(0.1), lm=0.5, pw=0.5)
  }, .logitNormPowYeoJohnson)

  .logitNormAddPropAddYeoJohnson1 <- c(82.941, 82.941, 82.833, 82.857, 82.857, 82.857)
  testErr("logitNorm+add+prop+yeoJohnson combined 1", function(f) {
    f %>% model(ipre ~ logitNorm(logit.sd, 0, 12) + prop(prop.sd) + yeoJohnson(lm)) %>%
      ini(logit.sd=sqrt(0.1), prop.sd=sqrt(0.1), lm=0.5)
  }, .logitNormAddPropAddYeoJohnson1, addProp=1)

  .logitNormAddPropAddYeoJohnson2 <- c(78.485, 78.485, 78.341, 78.373, 78.373, 78.373)
  testErr("logitNorm+add+prop+yeoJohnson combined2", function(f) {
    f %>% model(ipre ~ logitNorm(logit.sd, 0, 12) + prop(prop.sd) + yeoJohnson(lm)) %>%
      ini(logit.sd=sqrt(0.1), prop.sd=sqrt(0.1), lm=0.5)
  }, .logitNormAddPropAddYeoJohnson2, addProp = 2)


  .logitNormAddPowAddYeoJohnson1 <- c(82.941, 82.941, 82.833, 82.857, 82.857, 82.857)
  testErr("logitNorm+pow+yeoJohnson combined2", function(f) {
    f %>% model(ipre ~ logitNorm(logit.sd, 0, 12) + pow(prop.sd, pw) + yeoJohnson(lm)) %>%
      ini(logit.sd=sqrt(0.1), prop.sd=sqrt(0.1), lm=0.5)
  }, .logitNormAddPowAddYeoJohnson1, addProp = 1)

  .logitNormAddPowAddYeoJohnson2 <- c(78.485, 78.485, 78.341, 78.373, 78.373, 78.373)
  testErr("logitNorm+pow+yeoJohnson combined2", function(f) {
    f %>% model(ipre ~ logitNorm(logit.sd, 0, 12) + pow(prop.sd, pw) + yeoJohnson(lm)) %>%
      ini(logit.sd=sqrt(0.1), prop.sd=sqrt(0.1), lm=0.5)
  }, .logitNormAddPowAddYeoJohnson2, addProp = 2)


  .probitNormAdd <- c(12.827, 12.827, 12.827, 12.827, 12.847, 12.847)
  testErr("probitNorm", function(f) {
    f %>% model(ipre ~ probitNorm(logit.sd, 0, 12)) %>%
      ini(logit.sd=sqrt(0.1))
  }, .probitNormAdd)

  .probitNormProp <- c(88.875, 88.875, 88.733, 88.766, 88.766, 88.766)
  testErr("probitNorm(NA)+prop", function(f) {
    f %>% model(ipre ~ probitNorm(NA, 0, 12) + prop(prop.sd)) %>%
      ini(prop.sd=sqrt(0.1))
  }, .probitNormProp)

  testErr("probitNorm(NA)+pow->probitNorm(NA)+prop", function(f) {
    f %>% model(ipre ~ probitNorm(NA, 0, 12) + pow(prop.sd, pw)) %>%
      ini(prop.sd=sqrt(0.1), pw=1)
  }, .probitNormProp)

  .probitNormPow <- c(48.625, 48.625, 48.579, 48.587, 48.587, 48.587)
  testErr("probitNorm(NA)+pow->probitNorm(NA)+prop", function(f) {
    f %>% model(ipre ~ probitNorm(NA, 0, 12) + pow(prop.sd, pw)) %>%
      ini(prop.sd=sqrt(0.1), pw=0.5)
  }, .probitNormPow)

  .probitNormAddProp1 <- c(93.761, 93.761, 93.661, 93.682, 93.682, 93.682)
  testErr("probitNorm+prop", function(f) {
    f %>% model(ipre ~ probitNorm(probit.sd, 0, 12) + prop(prop.sd)) %>%
      ini(probit.sd=sqrt(0.1), prop.sd=sqrt(0.1))
  }, .probitNormAddProp1, addProp = 1)

  testErr("probitNorm+pow->probitNorm+prop", function(f) {
    f %>% model(ipre ~ probitNorm(probit.sd, 0, 12) + pow(prop.sd, pw)) %>%
      ini(probit.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=1)
  }, .probitNormAddProp1, addProp = 1)

  .probitNormAddPow1 <- c(60.418, 60.418, 60.396, 60.401, 60.401, 60.401)
  testErr("probitNorm+pow, combined1", function(f) {
    f %>% model(ipre ~ probitNorm(probit.sd, 0, 12) + pow(prop.sd, pw)) %>%
      ini(probit.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5)
  }, .probitNormAddPow1, addProp = 1)

  ## probitNorm + yeoJohnson
  .probitNormAddYeoJohnson <- c(19.69, 19.69, 19.69, 19.69, 19.729, 19.729)
  testErr("probitNorm+yeoJohnson", function(f) {
    f %>% model(ipre ~ probitNorm(probit.sd, 0, 12) + yeoJohnson(lm)) %>%
      ini(probit.sd=sqrt(0.1), lm=0.5)
  }, .probitNormAddYeoJohnson)

  .probitNormPropAddYeoJohnson <- c(96.041, 96.041, 95.899, 95.931, 95.931, 95.931)
  testErr("probitNorm(NA)+prop+yeoJohnson", function(f) {
    f %>% model(ipre ~ probitNorm(NA, 0, 12) + prop(prop.sd) + yeoJohnson(lm)) %>%
      ini(prop.sd=sqrt(0.1), lm=0.5)
  }, .probitNormPropAddYeoJohnson)

  testErr("probitNorm(NA)+pow+yeoJohnson->probitNorm(NA)+prop+yeoJohnson", function(f) {
    f %>% model(ipre ~ probitNorm(NA, 0, 12) + pow(prop.sd, pw) + yeoJohnson(lm)) %>%
      ini(prop.sd=sqrt(0.1), lm=0.5, pw=1)
  }, .probitNormPropAddYeoJohnson)

  .probitNormPowYeoJohnson <- c(55.798, 55.798, 55.751, 55.758, 55.758, 55.758)
  testErr("probitNorm(NA)+pow+yeoJohnson", function(f) {
    f %>% model(ipre ~ probitNorm(NA, 0, 12) + pow(prop.sd, pw) + yeoJohnson(lm)) %>%
      ini(prop.sd=sqrt(0.1), lm=0.5, pw=0.5)
  }, .probitNormPowYeoJohnson)

  .probitNormAddPropAddYeoJohnson1 <- c(100.925, 100.925, 100.824, 100.846, 100.846, 100.846)
  testErr("probitNorm+add+prop+yeoJohnson combined 1", function(f) {
    f %>% model(ipre ~ probitNorm(probit.sd, 0, 12) + prop(prop.sd) + yeoJohnson(lm)) %>%
      ini(probit.sd=sqrt(0.1), prop.sd=sqrt(0.1), lm=0.5)
  }, .probitNormAddPropAddYeoJohnson1, addProp=1)

  .probitNormAddPropAddYeoJohnson2 <- c(96.398, 96.398, 96.264, 96.295, 96.295, 96.295)
  testErr("probitNorm+add+prop+yeoJohnson combined2", function(f) {
    f %>% model(ipre ~ probitNorm(probit.sd, 0, 12) + prop(prop.sd) + yeoJohnson(lm)) %>%
      ini(probit.sd=sqrt(0.1), prop.sd=sqrt(0.1), lm=0.5)
  }, .probitNormAddPropAddYeoJohnson2, addProp = 2)

  .probitNormAddPowAddYeoJohnson1 <- c(100.925, 100.925, 100.824, 100.846, 100.846, 100.846)
  testErr("probitNorm+pow+yeoJohnson combined1", function(f) {
    f %>% model(ipre ~ probitNorm(probit.sd, 0, 12) + pow(prop.sd, pw) + yeoJohnson(lm)) %>%
      ini(probit.sd=sqrt(0.1), prop.sd=sqrt(0.1), lm=0.5)
  }, .probitNormAddPowAddYeoJohnson1, addProp = 1)

  .probitNormAddPowAddYeoJohnson2 <- c(96.398, 96.398, 96.264, 96.295, 96.295, 96.295)
  testErr("probitNorm+pow+yeoJohnson combined2", function(f) {
    f %>% model(ipre ~ probitNorm(probit.sd, 0, 12) + pow(prop.sd, pw) + yeoJohnson(lm)) %>%
      ini(probit.sd=sqrt(0.1), prop.sd=sqrt(0.1), lm=0.5)
  }, .probitNormAddPowAddYeoJohnson2, addProp = 2)

  rxode2::rxUnloadAll()
})
