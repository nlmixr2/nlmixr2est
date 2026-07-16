## est="npag" error-structure support.  The residual enters likInner0 as
## f + sqrt(r)*eps with eps ~ N(0,1), so any transformably-linear structure
## carried in r (proportional, combined) flows through the conditional
## likelihood; transform-both-sides links (lnorm/boxCox/yeoJohnson) transform the
## DV and add the per-observation Jacobian so lambda-type parameters are
## estimable.  Real fits -> weekly slow batch.

nmTest({

  test_that("est='npag' fits proportional error and certifies optimality (gamma)", {
    .m <- function() {
      ini({ tka<-log(1.5); tv<-log(31.5); tke<-log(0.08); prop.sd<-0.1
        eta.ka~0.3; eta.ke~0.1 })
      model({ ka<-exp(tka+eta.ka); v<-exp(tv); ke<-exp(tke+eta.ke)
        d/dt(depot)<- -ka*depot; d/dt(center)<-ka*depot-ke*center
        cp<-center/v; cp~prop(prop.sd) })
    }
    f <- nlmixr2(.m, nlmixr2data::theo_sd, est="npag",
                 control=npagControl(points=500L, cycles=30L, gammaOptimize=TRUE))
    expect_s3_class(f, "nlmixr2FitData")
    ## with the gamma-consistent certificate the NPML is reached (D(F) ~ 0)
    expect_true(is.finite(f$env$npagDF) && abs(f$env$npagDF) < 1e-2)
  })

  test_that("est='npag' estimates the residual: gamma is folded into add.sd", {
    ## start add.sd deliberately too large; with gamma optimization the fitted
    ## assay-error multiplier is folded back into add.sd, so the REPORTED add.sd
    ## recovers the true residual magnitude (near the FOCEI estimate ~0.78)
    .m <- function() {
      ini({ tka<-log(1.5); tv<-log(31.5); tke<-log(0.08); add.sd<-2.0
        eta.ka~0.3; eta.ke~0.1 })
      model({ ka<-exp(tka+eta.ka); v<-exp(tv); ke<-exp(tke+eta.ke)
        d/dt(depot)<- -ka*depot; d/dt(center)<-ka*depot-ke*center
        cp<-center/v; cp~add(add.sd) })
    }
    f <- nlmixr2(.m, nlmixr2data::theo_sd, est="npag",
                 control=npagControl(points=500L, cycles=40L, gammaOptimize=TRUE))
    expect_s3_class(f, "nlmixr2FitData")
    ## recovered residual is far from the (wrong) ini 2.0 and near the truth
    expect_true(f$parFixedDf["add.sd", "Estimate"] < 1.2)
    expect_equal(unname(f$parFixedDf["add.sd", "Estimate"]), 0.78, tolerance = 0.15)
  })

  test_that("est='npag' with residOptimize='none' holds the residual at ini", {
    ## residOptimize='none' (and no gamma warm-start) holds residual params fixed
    .m <- function() {
      ini({ tka<-log(1.5); tv<-log(31.5); tke<-log(0.08); add.sd<-2.0
        eta.ka~0.3; eta.ke~0.1 })
      model({ ka<-exp(tka+eta.ka); v<-exp(tv); ke<-exp(tke+eta.ke)
        d/dt(depot)<- -ka*depot; d/dt(center)<-ka*depot-ke*center
        cp<-center/v; cp~add(add.sd) })
    }
    f <- nlmixr2(.m, nlmixr2data::theo_sd, est="npag",
                 control=npagControl(points=300L, cycles=10L,
                                     gammaOptimize=FALSE, residOptimize="none"))
    expect_equal(unname(f$parFixedDf["add.sd", "Estimate"]), 2.0, tolerance = 1e-8)
  })

  test_that("est='npag' fits combined additive+proportional error", {
    .m <- function() {
      ini({ tka<-log(1.5); tv<-log(31.5); tke<-log(0.08); add.sd<-0.3; prop.sd<-0.1
        eta.ka~0.3; eta.ke~0.1 })
      model({ ka<-exp(tka+eta.ka); v<-exp(tv); ke<-exp(tke+eta.ke)
        d/dt(depot)<- -ka*depot; d/dt(center)<-ka*depot-ke*center
        cp<-center/v; cp~add(add.sd)+prop(prop.sd) })
    }
    f <- nlmixr2(.m, nlmixr2data::theo_sd, est="npag",
                 control=npagControl(points=300L, cycles=12L, gammaOptimize=FALSE))
    expect_s3_class(f, "nlmixr2FitData")
    expect_true(is.finite(as.numeric(f$objf)))
  })

  test_that("est='npag' estimates the boxCox transform lambda", {
    ## gamma cannot represent a transform shape parameter, so npag optimizes lambda
    ## directly.  Start lambda deliberately high (1.5); it should move well away
    ## toward the data-supported value (FOCEI ~ 0.44 on theo).
    .m <- function() {
      ini({ tka<-log(1.5); tv<-log(31.5); tke<-log(0.08); add.sd<-0.7; lambda<-1.5
        eta.ka~0.3; eta.ke~0.1 })
      model({ ka<-exp(tka+eta.ka); v<-exp(tv); ke<-exp(tke+eta.ke)
        d/dt(depot)<- -ka*depot; d/dt(center)<-ka*depot-ke*center
        cp<-center/v; cp~add(add.sd)+boxCox(lambda) })
    }
    f <- nlmixr2(.m, nlmixr2data::theo_sd, est="npag",
                 control=npagControl(points=400L, cycles=40L, gammaOptimize=TRUE))
    expect_true(f$parFixedDf["lambda", "Estimate"] < 0.9)   # moved far from ini 1.5
    expect_equal(unname(f$parFixedDf["lambda", "Estimate"]), 0.44, tolerance = 0.2)
  })

  test_that("est='npag' estimates the AR(1) correlation", {
    ## theo has essentially no residual autocorrelation; started at 0.7 the AR
    ## coordinate search should drive ar1.cor back toward 0.
    .m <- function() {
      ini({ tka<-log(1.5); tv<-log(31.5); tke<-log(0.08); add.sd<-0.7; ar1.cor<-0.7
        eta.ka~0.3; eta.ke~0.1 })
      model({ ka<-exp(tka+eta.ka); v<-exp(tv); ke<-exp(tke+eta.ke)
        d/dt(depot)<- -ka*depot; d/dt(center)<-ka*depot-ke*center
        cp<-center/v; cp~add(add.sd)+ar(ar1.cor) })
    }
    f <- nlmixr2(.m, nlmixr2data::theo_sd, est="npag",
                 control=npagControl(points=400L, cycles=40L, gammaOptimize=TRUE))
    expect_true(f$parFixedDf["ar1.cor", "Estimate"] < 0.2)  # moved from ini 0.7 -> ~0
  })

  test_that("est='npag' fits a boxCox transform-both-sides model (Jacobian)", {
    .m <- function() {
      ini({ tka<-log(1.5); tv<-log(31.5); tke<-log(0.08); add.sd<-0.7; lambda<-0.5
        eta.ka~0.3; eta.ke~0.1 })
      model({ ka<-exp(tka+eta.ka); v<-exp(tv); ke<-exp(tke+eta.ke)
        d/dt(depot)<- -ka*depot; d/dt(center)<-ka*depot-ke*center
        cp<-center/v; cp~add(add.sd)+boxCox(lambda) })
    }
    f <- nlmixr2(.m, nlmixr2data::theo_sd, est="npag",
                 control=npagControl(points=300L, cycles=12L, gammaOptimize=FALSE))
    expect_s3_class(f, "nlmixr2FitData")
    expect_true(is.finite(as.numeric(f$objf)))
  })

  test_that("est='npag' fits a lognormal residual and certifies optimality", {
    ## drop the t=0 observations: the structural prediction is 0 there, and the
    ## lnorm (log) link is undefined at 0 -- see the clear-error test below
    .d <- nlmixr2data::theo_sd
    .d <- .d[!(.d$EVID==0 & .d$TIME==0), ]
    .m <- function() {
      ini({ tka<-log(1.5); tv<-log(31.5); tke<-log(0.08); add.sd<-0.5
        eta.ka~0.3; eta.ke~0.1 })
      model({ ka<-exp(tka+eta.ka); v<-exp(tv); ke<-exp(tke+eta.ke)
        d/dt(depot)<- -ka*depot; d/dt(center)<-ka*depot-ke*center
        cp<-center/v; cp~lnorm(add.sd) })
    }
    f <- nlmixr2(.m, .d, est="npag",
                 control=npagControl(points=300L, cycles=15L, gammaOptimize=FALSE))
    expect_s3_class(f, "nlmixr2FitData")
    expect_true(is.finite(f$env$npagDF) && abs(f$env$npagDF) < 5e-2)
  })

  test_that("est='npag' reports a clear error when a transform sees a zero prediction", {
    .m <- function() {
      ini({ tka<-log(1.5); tv<-log(31.5); tke<-log(0.08); add.sd<-0.5
        eta.ka~0.3; eta.ke~0.1 })
      model({ ka<-exp(tka+eta.ka); v<-exp(tv); ke<-exp(tke+eta.ke)
        d/dt(depot)<- -ka*depot; d/dt(center)<-ka*depot-ke*center
        cp<-center/v; cp~lnorm(add.sd) })
    }
    ## theo_sd keeps the t=0 observation where cp is structurally 0 -> lnorm(0)
    expect_error(
      nlmixr2(.m, nlmixr2data::theo_sd, est="npag",
              control=npagControl(points=100L, cycles=3L, gammaOptimize=FALSE)),
      "zero conditional density at every support point")
  })

})
