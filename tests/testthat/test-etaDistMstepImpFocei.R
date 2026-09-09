# Every way a declared-distribution M-step can decline is SILENT, and a fit
# whose estimates look reasonable is no evidence at all that it engaged -- that
# is the lesson `test-etaDistMstepEngages.R` already records for saem.  imp and
# focei had no such coverage: no test ran a declared model through est="imp" at
# all, and the focei tests checked parameter-history bookkeeping rather than
# whether the M-step ran.
#
# Both estimators publish a counter for exactly this reason:
#   imp    fit$env$impEtaDistN   (src/imp.cpp, alongside a printed warning at 0)
#   focei  foceiEtaDistN_()      (src/inner.cpp)
#
# So these tests assert the mechanism, not the numbers.

nmTest({

  .edGamma <- function() {
    function() {
      ini({
        tka   <- 0.45
        lclm  <- log(1.1)
        lclrv <- log(0.09)
        tv    <- 3.45
        eta.ka ~ 0.6
        add.sd <- 0.7
      })
      model({
        dist(cl) ~ dgamma(shape = 1/exp(lclrv),
                          rate  = 1/(exp(lclrv)*exp(lclm)))
        ka <- exp(tka + eta.ka)
        v  <- exp(tv)
        d/dt(depot)  <- -ka*depot
        d/dt(center) <-  ka*depot - cl/v*center
        cp <- center/v
        cp ~ add(add.sd)
      })
    }
  }

  test_that("imp's declared-distribution M-step actually fires", {
    skip_on_cran()
    .d <- nlmixr2data::theo_sd
    .f <- suppressWarnings(
      nlmixr2(.edGamma(), .d, est = "imp",
              control = impControl(nIter = 5L, isample = 50L, print = 0L,
                                   covMethod = "", calcTables = FALSE,
                                   etaDistMstep = TRUE,
                                   etaDistWarmStart = FALSE,
                                   mceta = 0L)))
    expect_true(!is.null(.f$env$impEtaDistN))
    expect_gt(.f$env$impEtaDistN, 0)
  })

  test_that("imp reports it when the M-step is off", {
    skip_on_cran()
    .d <- nlmixr2data::theo_sd
    .f <- suppressWarnings(
      nlmixr2(.edGamma(), .d, est = "imp",
              control = impControl(nIter = 5L, isample = 50L, print = 0L,
                                   covMethod = "", calcTables = FALSE,
                                   etaDistMstep = FALSE,
                                   etaDistWarmStart = FALSE,
                                   mceta = 0L)))
    # off means zero, not absent -- the counter is always published
    expect_identical(as.numeric(.f$env$impEtaDistN), 0)
  })

  test_that("focei's declared-distribution M-step actually fires", {
    skip_on_cran()
    .d <- nlmixr2data::theo_sd
    .f <- suppressWarnings(
      nlmixr2(.edGamma(), .d, est = "focei",
              control = foceiControl(maxOuterIterations = 3L, print = 0L,
                                     covMethod = "", calcTables = FALSE,
                                     etaDistMstep = TRUE,
                                     etaDistWarmStart = FALSE,
                                     mceta = 0L)))
    expect_gt(nlmixr2est:::foceiEtaDistN_(), 0)
    # and it must not have silently declined into the inert path
    expect_false(any(grepl("etaDistMstep=TRUE was requested",
                           as.character(.f$runInfo), fixed = TRUE)))
  })

  test_that("mceta is a real impmapControl argument", {
    # imp has always USED mceta -- it goes through the shared FOCEi inner MAP and
    # reaches the C++ where focei's does -- but it was reachable only through
    # `...`, so it was neither validated nor documented.  It matters here: at
    # mceta=10 a single Newton M-step moved log relative variance from -2.46 to
    # +60.3 on Bauer's gamma model, which is what forced the Levenberg-Marquardt
    # damping in src/imp.cpp.
    expect_true("mceta" %in% names(formals(impmapControl)))
    expect_identical(impmapControl()$mceta, -2L)
    expect_identical(impmapControl(mceta = 10L)$mceta, 10L)
    expect_error(impmapControl(mceta = -5L), "mceta")
    # and it survives the control round-trip
    .ct <- impmapControl(mceta = 4L)
    expect_identical(do.call(impmapControl, .ct)$mceta, 4L)
  })

})
