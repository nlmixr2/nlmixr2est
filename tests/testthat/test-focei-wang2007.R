nlmixr2Test(
  {

    ## For some reason the ODE and solved FOCE proportional models
    ## give quite different results.  However, FOCE doesn't work as
    ## well with prop models.

    .nlmixr <- function(...) suppressWarnings(nlmixr(...))

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
      return(ret)
    }

    .propVals <- c(39.458, 39.458, 39.275, 39.207, 39.213, 39.207)

    .propModVals <- c(63.353, 63.353, 63.001, 63.063, 63.063, 63.063)

    .propFVals <- c(6.496, 6.496, 6.488, 6.275, 9.262, 9.262)

    .propFModVals <- c(19.177, 19.177, 19.07, 18.202, 18.333, 18.333)

    .addVals <- c(-2.059, -2.059, -2.059, -2.059, 0.026, 0.026)

    .addProp2 <- c(39.735, 39.735, 39.562, 39.499, 39.505, 39.499)

    .addProp1 <- c(43.554, 43.554, 43.416, 43.394, 43.398, 43.394)

    .addPow1 <- c(16.231, 16.231, 16.219, 16.008, 16.093, 16.093)

    .addPow2 <- c(10.886, 10.886, 10.868, 10.417, 10.662, 10.662)

    .lnorm <- c(40.039, 40.039, 40.039, 40.039, 40.055, 40.039)

    .addModVals <- c(6.081, 6.081, 6.073, 5.783, 6.162, 6.162)

    testErr("prop", function(f) {
      f %>% model(ipre ~ prop(prop.sd)) %>% ini(prop.sd=sqrt(0.1))
    }, .propVals)

    # In this case propT = prop
    testErr("propT", function(f) {
      f %>% model(ipre ~ propT(prop.sd)) %>% ini(prop.sd=sqrt(0.1))
    }, .propVals)

    testErr("propMod", function(f) {
      f %>% model(ipre ~ prop(f2))
    }, .propModVals)

    testErr("propF", function(f) {
      f %>% model(ipre ~ propF(prop.sd, f2)) %>% ini(prop.sd=sqrt(0.1))
    }, .propFVals)

    testErr("propFMod", function(f) {
      f %>% model(ipre ~ propF(lipre, f2))
    }, .propFModVals)

    testErr("add", function(f) {
      f %>% model(ipre ~ add(add.sd)) %>% ini(add.sd=sqrt(0.1))
    }, .addVals)

    testErr("addMod", function(f) {
      f %>% model(ipre ~ add(f2))
    }, .addModVals)

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
    }, c(9.966, 9.966, 9.948, 9.331, 9.651, 9.651))

    testErr("add+prop, combined 2->add", function(f) {
      f %>% model(ipre ~ add(add.sd) + prop(prop.sd)) %>% ini(add.sd=sqrt(0.1), prop.sd=0)
    }, .addVals, addProp = 2)

    testErr("add+prop, combined 2->prop", function(f) {
      f %>% model(ipre ~ add(add.sd) + prop(prop.sd)) %>% ini(add.sd=0, prop.sd=sqrt(0.1))
    }, .propVals, addProp = 2)

    testErr("add+prop, combined 2", function(f) {
      f %>% model(ipre ~ add(add.sd) + prop(prop.sd)) %>% ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1))
    }, .addProp2, addProp = 2)

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

    testErr("add+pow combined 2 (override)", function(f) {
      f %>% model(ipre ~ add(add.sd) + pow(prop.sd, pw) + combined2()) %>% ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5)
    }, .addPow2, addProp = 1)

    testErr("add+pow combined 1 (override)", function(f) {
      f %>% model(ipre ~ add(add.sd) + pow(prop.sd, pw) + combined1()) %>% ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5)
    }, .addPow1, addProp = 2)


    ## start looking at transformations

    testErr("lnorm", function(f) {
      f %>% model(ipre ~ lnorm(lnorm.sd)) %>% ini(lnorm.sd=sqrt(0.1))
    }, .lnorm, addProp = 1)

    testErr("add lnorm", function(f) {
      f %>% model(log(ipre) ~ add(add.sd)) %>% ini(add.sd=sqrt(0.1))
    },  log=TRUE) + 2 * sum(datl$DV)

    .lnorm


    testErr("lnorm", function(f) {
      f %>% model(ipre ~ lnorm(lnorm.sd)) %>% ini(lnorm.sd=sqrt(0.1))
    }, c(40.039, 40.039, 40.039, 40.039, 40.055, 40.039), addProp = 1)

    testErr("lnorm(NA)+prop", function(f) {
      f %>% model(ipre ~ lnorm(NA) + prop(prop.sd)) %>% ini(prop.sd=sqrt(0.1))
    }, c(118.419, 118.419, 118.279, 118.311, 118.311, 118.311), addProp = 1)

    testErr("lnorm(NA)+pow", function(f) {
      f %>% model(ipre ~ lnorm(NA) + pow(prop.sd, pw)) %>% ini(prop.sd=sqrt(0.1), pw=0.5)
    }, c(94.535, 94.535, 94.461, 94.478, 94.478, 94.478), addProp = 1)

    testErr("lnorm+prop combined1", function(f) {
      f %>% model(ipre ~ lnorm(lnorm.sd) + prop(prop.sd)) %>% ini(lnorm.sd=sqrt(0.1), prop.sd=0.5)
    }, c(123.318, 123.318, 123.219, 123.24, 123.24, 123.24), addProp = 1)

    testErr("lnorm+prop combined2", function(f) {
      f %>% model(ipre ~ lnorm(lnorm.sd) + prop(prop.sd)) %>% ini(lnorm.sd=sqrt(0.1), prop.sd=0.5)
    }, c(118.777, 118.777, 118.646, 118.676, 118.676, 118.676), addProp = 2)

    testErr("lnorm+pow", function() {
      lnorm(0.1) + pow(0.1, 0.5)
    }, c(102.899, 102.899, 102.855, 102.865, 102.865, 102.865), addProp = 1)

    testErr("lnorm+pow", function() {
      lnorm(0.1) + pow(0.1, 0.5)
    }, c(95.634, 95.634, 95.57, 95.585, 95.585, 95.585), addProp = 2)

    ## Box Cox
    testErr("add+boxCox", function() {
      return(add(.1) + boxCox(.5))
    }, c(2.06, 2.06, 2.06, 2.06, 3.529, 2.06))

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
  },
  test = "focei"
)
