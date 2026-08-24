# End-to-end numeric validation of M2/M3/M4 censoring for a llik-forced
# t()/cauchy()/dnorm() endpoint under FOCEi (#992).  Three independent
# references, because a self-consistent-but-wrong wiring (reading the wrong
# lhs column, say) would satisfy any one of them alone:
#
#  1. population-only (neta == 0): objf is exactly -2*(sum ll - nobs*log(sqrt(2pi))),
#     so it is compared against the censored t() log-likelihood recomputed in R
#     from the fit's own IPRED.
#  2. with etas: the same censored likelihood written out by hand as a general
#     ll() endpoint (pcauchy via atan) must give the SAME objective -- both go
#     through the identical FOCEi log-likelihood machinery, so log|H| and the
#     adjLik bookkeeping cancel exactly.
#  3. the analytic inner eta-gradient must match central differences of the
#     inner objective (focei_tCensDll vs focei_tCensLl).
#
# Multiple fits per test -- weekly-batched via .slowBatches in tests/testthat.R,
# so do NOT add skip_on_ci().

.censTheo <- function(n = 4L) {
  .d <- nlmixr2data::theo_sd
  .d[.d$ID <= n, ]
}
# rows censored below the LLOQ, with DV replaced by the LLOQ itself
.censRows <- function(d, lloq = 4) {
  which(d$EVID == 0 & d$DV < lloq & d$DV > 0)
}
# one dataset per censoring mode: "none", "m3" (CENS=1 at the LLOQ),
# "m4" (M3 plus a finite LIMIT), "m2" (LIMIT only, nothing censored)
.censDat <- function(d, mode, lloq = 4, limit = 0.5) {
  .lo <- .censRows(d, lloq)
  d$CENS <- 0L
  d$LIMIT <- NA_real_
  if (mode %in% c("m3", "m4")) {
    d$CENS[.lo] <- 1L
    d$DV[.lo] <- lloq
  }
  if (mode == "m4") d$LIMIT[.lo] <- limit
  if (mode == "m2") d$LIMIT[d$EVID == 0] <- limit
  d
}
# the same rows, as plain covariate columns for the hand-written ll() model
.manDat <- function(d, mode, lloq = 4, limit = 0.5) {
  .lo <- .censRows(d, lloq)
  d$CENSF <- 0
  d$LOQ <- lloq
  d$LIM <- limit
  if (mode %in% c("m3", "m4")) {
    d$CENSF[.lo] <- 1
    d$DV[.lo] <- lloq
  }
  d
}

# censored Student-t log-likelihood, straight out of Beal (2001) -- the R
# counterpart of censEst.h's doCensT1()
.refLl <- function(dv, cens, lim, f, sd, nu) {
  .z <- function(x) (x - f) / sd
  if (cens != 0) {
    if (is.finite(lim)) {
      return(log(stats::pt(cens * .z(dv), nu) - stats::pt(cens * .z(lim), nu)) -
               stats::pt(cens * .z(lim), nu, lower.tail = FALSE, log.p = TRUE))
    }
    return(stats::pt(cens * .z(dv), nu, log.p = TRUE))
  }
  .ll <- stats::dt(.z(dv), nu, log = TRUE) - log(sd)
  if (is.finite(lim)) {
    .zz <- ifelse(lim < f, 1, -1) * (lim - f) / sd
    .ll <- .ll - stats::pt(.zz, nu, lower.tail = FALSE, log.p = TRUE)
  }
  .ll
}

nmTest({

  ## ---------------------------------------------------------------- 1. value
  .popT <- function() {
    ini({
      tka <- fix(0.45)
      tcl <- fix(1)
      tv <- fix(3.45)
      add.sd <- 0.7
      nu <- fix(5)
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl)
      v <- exp(tv)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd) + t(nu)
    })
  }

  test_that("population-only FOCEi t() censoring matches the closed form (#992)", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    .d <- .censTheo()
    .lo <- .censRows(.d)
    # -2*(sum ll - nObs*log(sqrt(2*pi))) IS the neta == 0 objective, with the
    # fit's own IPRED standing in for f so no ODE tolerance enters
    .refObjf <- function(fit, dat) {
      .obs <- dat[dat$EVID == 0, ]
      .f <- as.data.frame(fit)$IPRED
      expect_equal(length(.f), nrow(.obs))
      .ll <- vapply(seq_along(.f), function(k) {
        .refLl(.obs$DV[k], .obs$CENS[k], .obs$LIMIT[k], .f[k], 0.7, 5)
      }, numeric(1))
      -2 * (sum(.ll) - length(.ll) * log(sqrt(2 * pi)))
    }
    for (.nm in c("none", "m3", "m4", "m2")) {
      .dat <- .censDat(.d, .nm, limit = 0.25)
      .fit <- .nlmixr(.popT, .dat, est = "focei",
                      foceiControl(print = 0L, covMethod = "",
                                   maxOuterIterations = 0L))
      expect_equal(as.numeric(.fit$objf), .refObjf(.fit, .dat), tolerance = 1e-8,
                   info = paste("case", .nm))
    }
  })

  ## ------------------------------------------------- 2. value, with etas
  .etaCauchy <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      add.sd <- c(0, 0.7)
      eta.ka ~ 0.2
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl)
      v <- exp(tv)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd) + cauchy()
    })
  }
  # pcauchy(q, loc, scale) = 0.5 + atan((q-loc)/scale)/pi -- so the censored
  # cauchy likelihood is writable as a plain general ll() endpoint
  .manM3 <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      add.sd <- c(0, 0.7)
      eta.ka ~ 0.2
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl)
      v <- exp(tv)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      sdv <- abs(add.sd)
      zz <- (DV - cp) / sdv
      lcens <- log(0.5 + atan((LOQ - cp) / sdv) / pi)
      lunc <- -log(pi * sdv) - log(1 + zz * zz)
      ll(err) ~ CENSF * lcens + (1 - CENSF) * lunc
    })
  }
  .manM4 <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      add.sd <- c(0, 0.7)
      eta.ka ~ 0.2
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl)
      v <- exp(tv)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      sdv <- abs(add.sd)
      zz <- (DV - cp) / sdv
      pHi <- 0.5 + atan((LOQ - cp) / sdv) / pi
      pLo <- 0.5 + atan((LIM - cp) / sdv) / pi
      lcens <- log(pHi - pLo) - log(1 - pLo)
      lunc <- -log(pi * sdv) - log(1 + zz * zz)
      ll(err) ~ CENSF * lcens + (1 - CENSF) * lunc
    })
  }
  .manM2 <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      add.sd <- c(0, 0.7)
      eta.ka ~ 0.2
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl)
      v <- exp(tv)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      sdv <- abs(add.sd)
      zz <- (DV - cp) / sdv
      zl <- (2 * (LIM < cp) - 1) * (LIM - cp) / sdv
      ll(err) ~ -log(pi * sdv) - log(1 + zz * zz) - log(0.5 - atan(zl) / pi)
    })
  }

  test_that("FOCEi t()/cauchy() censoring matches a hand-written ll() (#992)", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    .d <- .censTheo()
    .ctl <- foceiControl(print = 0L, covMethod = "", calcTables = FALSE,
                         maxOuterIterations = 0L)
    .run <- function(m, dat) as.numeric(.nlmixr(m, dat, est = "focei", .ctl)$objf)

    expect_equal(.run(.etaCauchy, .censDat(.d, "m3")),
                 .run(.manM3, .manDat(.d, "m3")), tolerance = 1e-6)
    # est="foce" reaches the same place by a different route -- it is the
    # explicit form of the interaction=0 build .foceiFitInternal picks anyway
    .runFoce <- function(m, dat) {
      as.numeric(.nlmixr(m, dat, est = "foce", foceiControl(
        print = 0L, covMethod = "", calcTables = FALSE,
        maxOuterIterations = 0L))$objf)
    }
    expect_equal(.runFoce(.etaCauchy, .censDat(.d, "m3")),
                 .runFoce(.manM3, .manDat(.d, "m3")), tolerance = 1e-6)
    expect_equal(.run(.etaCauchy, .censDat(.d, "m4")),
                 .run(.manM4, .manDat(.d, "m4")), tolerance = 1e-5)
    expect_equal(.run(.etaCauchy, .censDat(.d, "m2", limit = 0.25)),
                 .run(.manM2, .manDat(.d, "m2", limit = 0.25)), tolerance = 1e-6)
  })

  ## ------------------------------------------------ mechanism-used evidence
  test_that("FOCEi reports the censoring method it applied (#992)", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    .d <- .censTheo(2L)
    .ctl <- foceiControl(print = 0L, covMethod = "", calcTables = FALSE,
                         maxOuterIterations = 0L)
    .m3 <- .censDat(.d, "m3", lloq = 3)
    .naive <- .m3
    .naive$CENS <- 0L
    .fitCens <- .nlmixr(.etaCauchy, .m3, est = "focei", .ctl)
    .fitNaive <- .nlmixr(.etaCauchy, .naive, est = "focei", .ctl)
    # the FOCEi family appends the inner-Hessian censoring treatment
    expect_equal(as.character(.fitCens$censInformation), "M3 censoring (gauss)")
    expect_equal(as.character(.fitNaive$censInformation), "No censoring")
    # the correction has to MOVE the objective; that it ran is not enough
    expect_true(abs(.fitCens$objf - .fitNaive$objf) > 1e-6)
    # ...and it must no longer warn that censoring was ignored
    expect_no_warning(suppressMessages(
      nlmixr2(.etaCauchy, .m3, est = "focei", .ctl)))
  })

  test_that("a llik-forced dnorm() endpoint honors censoring too (#992)", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    # dnorm() takes the same llik-forced path as t()/cauchy() (rx_pred_ is the
    # log-density), so it needs the same rx_pred_f_/rx_r_ columns -- routed
    # through doCensNormal1() rather than doCensT1().
    .dn <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        add.sd <- c(0, 0.7)
        eta.ka ~ 0.2
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl)
        v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd) + dnorm()
      })
    }
    .manDn <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        add.sd <- c(0, 0.7)
        eta.ka ~ 0.2
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl)
        v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        sdv <- abs(add.sd)
        zz <- (DV - cp) / sdv
        lcens <- log(pnorm((LOQ - cp) / sdv))
        lunc <- -0.5 * log(2 * pi) - log(sdv) - 0.5 * zz * zz
        ll(err) ~ CENSF * lcens + (1 - CENSF) * lunc
      })
    }
    .d <- .censTheo()
    .ctl <- foceiControl(print = 0L, covMethod = "", calcTables = FALSE,
                         maxOuterIterations = 0L)
    expect_equal(as.numeric(.nlmixr(.dn, .censDat(.d, "m3"), est = "focei", .ctl)$objf),
                 as.numeric(.nlmixr(.manDn, .manDat(.d, "m3"), est = "focei", .ctl)$objf),
                 tolerance = 1e-6)
  })
  ## ------------------------------------------------------- 3. eta gradient
  # LAST in the file on purpose: .vaeInnerSetup()'s direct likInner() calls set
  # censEst.h's globalCensFlag, which only a completed FIT resets -- running
  # this before a fit-based test would leak M2/M3/M4 into that fit's
  # $censInformation.
  test_that("the censored inner eta-gradient matches central differences (#992)", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    .addM <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        add.sd <- c(0, 0.7)
        eta.ka ~ 0.2
        eta.cl ~ 0.1
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd) + cauchy()
      })
    }
    # a prop() error makes R itself eta-dependent, so this one only agrees
    # with central differences if d(R)/d(eta) reaches dCensT1 -- FOCE's inner
    # model has no such column in the FOCEi position, so it is appended
    .propM <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        prop.sd <- c(0, 0.3)
        eta.ka ~ 0.2
        eta.cl ~ 0.1
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ prop(prop.sd) + cauchy()
      })
    }
    .d <- .censTheo()
    .N <- length(unique(.d$ID))
    .testSeed(11)
    .etaMat <- matrix(stats::rnorm(.N * 2, 0, 0.1), .N, 2)
    .fdLp <- function(eta, id, h = 1e-5) {
      vapply(seq_along(eta), function(j) {
        .ep <- eta
        .ep[j] <- .ep[j] + h
        .em <- eta
        .em[j] <- .em[j] - h
        (likInner(.ep, id) - likInner(.em, id)) / (2 * h)
      }, numeric(1))
    }
    for (.err in c("add", "prop")) {
      .ui <- rxode2::assertRxUi(if (.err == "add") .addM else .propM)
      for (.mode in c("m3", "m4", "m2")) {
        .dat <- .censDat(.d, .mode, limit = if (.mode == "m2") 0.25 else 0.5)
        suppressWarnings(.vaeInnerSetup(.ui, .dat, .etaMat, vaeControl()))
        for (.id in seq_len(.N)) {
          expect_equal(as.numeric(foceiInnerLp(.etaMat[.id, ], .id)),
                       .fdLp(.etaMat[.id, ], .id), tolerance = 1e-3,
                       info = paste(.err, .mode, "id", .id))
        }
        .vaeInnerFree()
      }
    }
  })

})
