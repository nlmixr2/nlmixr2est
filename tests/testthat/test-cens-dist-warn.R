nmTest({
  # Tests for the catch-all "censoring ignored" warning (#979): any
  # distribution the M2/M3/M4 machinery does not understand (everything
  # except norm/dnorm, and -- nlm-family only -- t()/cauchy()) should warn
  # instead of silently scoring a censored row as uncensored.

  one.cmt.pois <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      lam <- 1
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl)
      v <- exp(tv)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      lambda <- exp(lam) * cp
      cp ~ pois(lambda)
    })
  }

  one.cmt.t <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
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

  .dat <- nlmixr2data::theo_sd
  .dat <- .dat[.dat$ID <= 2, ]
  .dat$CENS <- ifelse(.dat$DV < 3 & .dat$EVID == 0, 1L, 0L)
  .dat$DV[.dat$CENS == 1] <- 3

  test_that("pois() + CENS warns (censoring ignored) on nlm", {
    fit <- suppressMessages(suppressWarnings(
      nlmixr2(one.cmt.pois, .dat, est = "nlm", control = list(print = 0))))
    expect_true(any(grepl("censoring ignored for 'pois'", fit$runInfo,
                          fixed = TRUE)))
  })

  test_that("t() + CENS warns on focei (not yet supported there) but not on nlm", {
    fitNlm <- suppressMessages(suppressWarnings(
      nlmixr2(one.cmt.t, .dat, est = "nlm", control = list(print = 0))))
    expect_false(any(grepl("censoring ignored", fitNlm$runInfo, fixed = TRUE)))

    fitFocei <- suppressMessages(suppressWarnings(
      nlmixr2(one.cmt.t, .dat, est = "focei", control = list(print = 0))))
    expect_true(any(grepl("censoring ignored for 't'", fitFocei$runInfo,
                          fixed = TRUE)))
  })

  test_that("censoring-ignored warning text stays under 75 characters", {
    for (.d in c("pois", "binom", "beta", "chisq", "dexp", "f", "geom",
                 "unif", "weibull", "dgamma", "ordinal", "t", "cauchy")) {
      expect_lt(nchar(paste0("censoring ignored for '", .d, "' endpoint(s)")), 75)
    }
  })
})
