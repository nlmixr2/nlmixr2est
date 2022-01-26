.nlmixr <- function(...) suppressMessages(suppressWarnings(nlmixr(...)))
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
  expect_equal(setNames(round(expect_warning(expect_warning(testErr("add lnorm", function(f) {
    f %>% model(lipre ~ add(add.sd)) %>% ini(add.sd=sqrt(0.1))
  },  log=TRUE) + 2 * sum(datl$DV), "changed to"), "changed to"), 3), NULL),
  .lnorm, tolerance=5e-3)
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
testErr("add+propT, combined 2->add", function(f) {
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

testErr("boxCox+add+propT, combined 2", function(f) {
  f %>% model(ipre ~ add(add.sd) + propT(prop.sd) + boxCox(lambda)) %>%
    ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), lambda=1)
}, .boxCoxAddPropT2Vals, addProp = 2)

testErr("boxCox+add+propT, combined 2 (specified)", function(f) {
  f %>% model(ipre ~ add(add.sd) + propT(prop.sd) + boxCox(lambda) + combined2()) %>%
    ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), lambda=1)
}, .boxCoxAddPropT2Vals, addProp = 1)

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

testErr("boxCox+add+propT, combined 2", function(f) {
  f %>% model(ipre ~ add(add.sd) + propT(prop.sd) + boxCox(lambda)) %>%
    ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), lambda=1)
}, .addProp2, addProp = 2)

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
  f %>% model(ipre ~ add(add.sd) + pow(prop.sd, pw) + combined1()) %>%
    ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5, lambda=1)
}, .addPow1, addProp = 2)

skip()
## Box Cox
testErr("prop+boxCox", function() {
  return(prop(.1) + boxCox(.5))
}, c(61.473, 61.473, 61.298, 61.324, 61.325, 61.324))

testErr("pow+boxCox", function() {
  return(pow(.1, 0.5) + boxCox(.5))
}, c(39.567, 39.567, 39.477, 39.402, 39.411, 39.402))

testErr("add+prop+boxCox", function() {
  return(add(0.1) + prop(.1) + boxCox(.5))
}, c(66.075, 66.075, 65.949, 65.973, 65.974, 65.973), addProp = 1)

testErr("add+prop+boxCox", function() {
  return(add(0.1) + prop(.1) + boxCox(.5))
}, c(61.802, 61.802, 61.636, 61.661, 61.662, 61.661), addProp = 2)

testErr("add+pow+boxCox", function() {
  return(add(0.1) + pow(.1, 0.5) + boxCox(.5))
}, c(46.768, 46.768, 46.709, 46.693, 46.695, 46.693), addProp = 1)

testErr("add+pow+boxCox", function() {
  return(add(0.1) + pow(.1, 0.5) + boxCox(.5))
}, c(40.451, 40.451, 40.372, 40.313, 40.32, 40.313), addProp = 2)

# Now yeoJohnson
testErr("add+yeoJohnson", function() {
  return(add(.1) + yeoJohnson(.5))
}, c(2.339, 2.339, 2.339, 2.339, 3.575, 2.339))

testErr("prop+yeoJohnson", function() {
  return(prop(.1) + yeoJohnson(.5))
}, c(62.821, 62.821, 62.647, 62.676, 62.677, 62.676))

testErr("pow+yeoJohnson", function() {
  return(pow(.1, 0.5) + yeoJohnson(.5))
}, c(40.724, 40.724, 40.632, 40.575, 40.581, 40.575))

testErr("add+prop+yeoJohnson", function() {
  return(add(0.1) + prop(.1) + yeoJohnson(.5))
}, c(67.453, 67.453, 67.329, 67.353, 67.354, 67.353), addProp = 1)

testErr("add+prop+yeoJohnson", function() {
  return(add(0.1) + prop(.1) + yeoJohnson(.5))
}, c(63.152, 63.152, 62.988, 63.016, 63.017, 63.016), addProp = 2)

testErr("add+pow+yeoJohnson", function() {
  return(add(0.1) + pow(.1, 0.5) + yeoJohnson(.5))
}, c(48.036, 48.036, 47.978, 47.967, 47.969, 47.967), addProp = 1)

testErr("add+pow+yeoJohnson", function() {
  return(add(0.1) + pow(.1, 0.5) + yeoJohnson(.5))
}, c(41.628, 41.628, 41.548, 41.503, 41.509, 41.503), addProp = 2)

## logitNorm
testErr("logitNorm", function() {
  return(logitNorm(.1, 0, 12))
}, c(0.612, 0.612, 0.612, 0.612, 0.786, 0.612))

testErr("logitNorm(NA)+prop", function() {
  return(logitNorm(NA, 0, 12) + prop(0.1))
}, c(67.882, 67.882, 67.731, 67.765, 67.765, 67.765))

testErr("logitNorm(NA)+pow", function() {
  return(logitNorm(NA, 0, 12) + pow(0.1, 0.5))
}, c(44.632, 44.632, 44.542, 44.556, 44.556, 44.556))

testErr("logitNorm+prop", function() {
  return(logitNorm(.1, 0, 12) + prop(0.1))
}, c(72.699, 72.699, 72.591, 72.615, 72.615, 72.615), addProp = 1)

testErr("logitNorm+prop", function() {
  return(logitNorm(.1, 0, 12) + prop(0.1))
}, c(68.233, 68.233, 68.09, 68.122, 68.123, 68.122), addProp = 2)

testErr("logitNorm+pow", function() {
  return(logitNorm(.1, 0, 12) + pow(0.1, 0.5))
}, c(52.641, 52.641, 52.589, 52.6, 52.6, 52.6), addProp = 1)

testErr("logitNorm+pow", function() {
  return(logitNorm(.1, 0, 12) + pow(0.1, 0.5))
}, c(45.668, 45.668, 45.59, 45.603, 45.603, 45.603), addProp = 2)

## logitNorm + yeoJohnson
testErr("logitNorm+yeoJohnson", function() {
  return(logitNorm(.1, 0, 12) + yeoJohnson(0.5))
}, c(5.127, 5.127, 5.127, 5.127, 5.484, 5.127))

testErr("logitNorm(NA)+prop+yeoJohnson", function() {
  return(logitNorm(NA, 0, 12) + prop(0.1) + yeoJohnson(0.5))
}, c(73.881, 73.881, 73.73, 73.764, 73.764, 73.764))

testErr("logitNorm(NA)+pow+yeoJohnson", function() {
  return(logitNorm(NA, 0, 12) + pow(0.1, 0.5) + yeoJohnson(0.5))
}, c(50.631, 50.631, 50.54, 50.551, 50.552, 50.551))

testErr("logitNorm+prop+yeoJohnson", function() {
  return(logitNorm(.1, 0, 12) + prop(0.1) + yeoJohnson(0.5))
}, c(78.693, 78.693, 78.585, 78.609, 78.609, 78.609), addProp = 1)

testErr("logitNorm+prop+yeoJohnson", function() {
  return(logitNorm(.1, 0, 12) + prop(0.1) + yeoJohnson(0.5))
}, c(74.231, 74.231, 74.088, 74.12, 74.12, 74.12), addProp = 2)

testErr("logitNorm+pow+yeoJohnson", function() {
  return(logitNorm(.1, 0, 12) + pow(0.1, 0.5) + yeoJohnson(0.5))
}, c(58.628, 58.628, 58.575, 58.586, 58.586, 58.586), addProp = 1)

testErr("logitNorm+pow+yeoJohnson", function() {
  return(logitNorm(.1, 0, 12) + pow(0.1, 0.5) + yeoJohnson(0.5))
}, c(51.662, 51.662, 51.584, 51.595, 51.595, 51.595), addProp = 2)

## probitNorm
testErr("probitNorm", function() {
  return(probitNorm(.1, 0, 12))
}, c(12.827, 12.827, 12.827, 12.827, 12.847, 12.827))

testErr("probitNorm(NA)+prop", function() {
  return(probitNorm(NA, 0, 12) + prop(0.1))
}, c(88.875, 88.875, 88.733, 88.766, 88.766, 88.766))

testErr("probitNorm(NA)+pow", function() {
  return(probitNorm(NA, 0, 12) + pow(0.1, 0.5))
}, c(65.098, 65.098, 65.02, 65.037, 65.037, 65.037))

testErr("probitNorm+prop", function() {
  return(probitNorm(0.1, 0, 12) + prop(0.1))
}, c(93.761, 93.761, 93.661, 93.682, 93.682, 93.682), addProp = 1)

testErr("probitNorm+prop", function() {
  return(probitNorm(0.1, 0, 12) + prop(0.1))
}, c(89.232, 89.232, 89.098, 89.129, 89.129, 89.129), addProp = 2)

testErr("probitNorm+pow", function() {
  return(probitNorm(0.1, 0, 12) + pow(0.1, 0.5))
}, c(73.405, 73.405, 73.359, 73.37, 73.37, 73.37), addProp = 1)

testErr("probitNorm+pow", function() {
  return(probitNorm(0.1, 0, 12) + pow(0.1, 0.5))
}, c(66.187, 66.187, 66.12, 66.135, 66.135, 66.135), addProp = 2)

## probitNorm + yeoJohnson

testErr("probitNorm+yeoJohnson", function() {
  return(probitNorm(.1, 0, 12) + yeoJohnson(0.5))
}, c(16.768, 16.768, 16.768, 16.768, 16.799, 16.768))

testErr("probitNorm(NA)+prop+yeoJohnson", function() {
  return(probitNorm(NA, 0, 12) + prop(0.1) + yeoJohnson(0.5))
}, c(93.071, 93.071, 92.929, 92.962, 92.962, 92.962))

testErr("probitNorm(NA)+pow+yeoJohnson", function() {
  return(probitNorm(NA, 0, 12) + pow(0.1, 0.5) + yeoJohnson(0.5))
}, c(69.295, 69.295, 69.217, 69.234, 69.234, 69.234))

testErr("probitNorm(0.1)+prop+yeoJohnson", function() {
  return(probitNorm(0.1, 0, 12) + prop(0.1) + yeoJohnson(0.5))
}, c(97.957, 97.957, 97.856, 97.878, 97.878, 97.878), addProp = 1)

testErr("probitNorm(0.1)+prop+yeoJohnson", function() {
  return(probitNorm(0.1, 0, 12) + prop(0.1) + yeoJohnson(0.5))
}, c(93.429, 93.429, 93.295, 93.326, 93.326, 93.326), addProp = 2)

testErr("probitNorm(0.1)+pow+yeoJohnson", function() {
  return(probitNorm(0.1, 0, 12) + pow(0.1, 0.5) + yeoJohnson(0.5))
}, c(77.599, 77.599, 77.553, 77.564, 77.564, 77.564), addProp = 1)

testErr("probitNorm(0.1)+pow+yeoJohnson", function() {
  return(probitNorm(0.1, 0, 12) + pow(0.1, 0.5) + yeoJohnson(0.5))
}, c(70.383, 70.383, 70.316, 70.331, 70.331, 70.331), addProp = 2)

## lognormal -- equivalent to add on log-space and back-transformed.

## Next run on the log-transformed space
datl <- dat
datl$DV <- log(datl$DV)
datl2 <- dat2
datl2$DV <- log(datl2$DV)
predl <- function() log(ipre)

fit.lnorm <- .foceiFit(dat, inits, mypar1, mod, pred, function() {
  return(lnorm(.1))
},
control = foceiControl(maxOuterIterations = 0, covMethod = "")
)

test_that("Matches NONMEM objective lognormal function; (Based on Wang2007)", {
  expect_equal(fit.lnorm$objective, 40.039, tol=1e-3)
})

fit.lnorm0 <- .foceiFit(datl, inits, mypar1, mod, predl, function() {
  return(add(.1))
},
control = foceiControl(maxOuterIterations = 0, covMethod = "")
)

test_that("Matches NONMEM objective lognormal function; (Based on Wang2007)", {
  expect_equal(fit.lnorm$objective, 40.039, tol=1e-3)
  expect_equal(fit.lnorm0$objective + 2 * sum(datl$DV), fit.lnorm$objective)
  expect_equal(fit.lnorm0$objective, -42.106, tol=1e-3)
})

fit.lnorm2 <- .foceiFit(dat2, inits, mypar1, m1, pred, function() {
  return(lnorm(.1))
},
control = foceiControl(maxOuterIterations = 0, covMethod = "")
)

fit.lnorm20 <- .foceiFit(datl2, inits, mypar1, m1, predl, function() {
  return(add(.1))
},
control = foceiControl(maxOuterIterations = 0, covMethod = "")
)

test_that("Matches NONMEM objective lognormal function; ODE (Based on Wang2007)", {
  expect_equal(fit.lnorm2$objective, 40.039, tol=1e-3)
  expect_equal(fit.lnorm20$objective + 2 * sum(datl$DV), fit.lnorm2$objective)
  expect_equal(fit.lnorm20$objective, -42.106, tol=1e-3)
})

fit.lnorm <- .foceiFit(dat, inits, mypar1, mod, pred, function() {
  return(lnorm(.1))
},
control = foceiControl(maxOuterIterations = 0, covMethod = "", interaction = FALSE)
)

fit.lnorm0 <- .foceiFit(datl, inits, mypar1, mod, predl, function() {
  return(add(.1))
},
control = foceiControl(maxOuterIterations = 0, covMethod = "", interaction = FALSE)
)

test_that("Matches NONMEM objective lognormal error FOCE (Based on Wang2007)", {
  expect_equal(fit.lnorm$objective, 40.039, tol=1e-3)
  expect_equal(fit.lnorm0$objective + 2 * sum(datl$DV), fit.lnorm$objective)
  expect_equal(fit.lnorm0$objective, -42.106, tol=1e-3)
})

fit.lnorm2 <- .foceiFit(dat2, inits, mypar1, m1, pred, function() {
  return(lnorm(.1))
},
control = foceiControl(maxOuterIterations = 0, covMethod = "", interaction = FALSE)
)

fit.lnorm20 <- .foceiFit(datl2, inits, mypar1, m1, predl, function() {
  return(add(.1))
},
control = foceiControl(maxOuterIterations = 0, covMethod = "", interaction = FALSE)
)

test_that("Matches NONMEM objective lognormal error FOCE; ODE (Based on Wang2007)", {
  expect_equal(fit.lnorm2$objective, 40.039, tol=1e-3)
  expect_equal(fit.lnorm20$objective + 2 * sum(datl$DV), fit.lnorm2$objective)
  expect_equal(fit.lnorm20$objective, -42.106, tol=1e-3)
})

fit.lnorm <- .foceiFit(dat, inits, mypar1, mod, pred, function() {
  return(lnorm(.1))
},
control = foceiControl(maxOuterIterations = 0, covMethod = "", interaction = FALSE, fo = TRUE)
)

fit.lnorm0 <- .foceiFit(datl, inits, mypar1, mod, predl, function() {
  return(add(.1))
},
control = foceiControl(maxOuterIterations = 0, covMethod = "", interaction = FALSE, fo = TRUE)
)

test_that("Matches NONMEM objective lognormal error FO (Based on Wang2007)", {
  expect_equal(fit.lnorm$objective, 40.055, tol=1e-3)
  expect_equal(fit.lnorm0$objective + 2 * sum(datl$DV), fit.lnorm$objective)
  expect_equal(fit.lnorm0$objective, -42.09, tol=1e-3)
})

fit.lnorm2 <- .foceiFit(dat2, inits, mypar1, m1, pred, function() {
  return(lnorm(.1))
},
control = foceiControl(maxOuterIterations = 0, covMethod = "", interaction = FALSE, fo = TRUE)
)

fit.lnorm20 <- .foceiFit(datl2, inits, mypar1, m1, predl, function() {
  return(add(.1))
},
control = foceiControl(maxOuterIterations = 0, covMethod = "", interaction = FALSE, fo = TRUE)
)

test_that("Matches NONMEM objective lognormal error FO; ODE (Based on Wang2007)", {
  expect_equal(fit.lnorm2$objective, 40.055, tol=1e-3)
  expect_equal(fit.lnorm20$objective + 2 * sum(datl$DV), fit.lnorm2$objective)
  expect_equal(fit.lnorm20$objective, -42.09, tol=1e-3)
})
