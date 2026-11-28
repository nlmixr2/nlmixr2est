## emvi recovers proportional, combined and eta-dependent residual error, which
## the SGA optimizer overestimated about 2x (#1462).
nmTest({
  .viResidData <- function(errSd, nsub = 40L) {
    .testSeed(42)
    times <- c(0.25, 0.5, 1, 2, 4, 6, 8, 12, 16, 24)
    sim <- rxode2::rxode2({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka * depot - cl / v * center
      cp <- center / v
    })
    ev <- rxode2::et(amt = 320, cmt = "depot") |>
      rxode2::et(times) |>
      rxode2::et(id = seq_len(nsub))
    s <- rxode2::rxSolve(
      sim,
      c(tka = log(1.5), tcl = log(2.8), tv = log(32)),
      ev,
      omega = lotri::lotri(eta.ka + eta.cl + eta.v ~ c(0.1, 0, 0.1, 0, 0, 0.1)),
      returnType = "data.frame"
    )
    etaRes <- rnorm(nsub, 0, sqrt(0.2))[s$id]
    s$DV <- s$cp + errSd(s, etaRes) * rnorm(nrow(s))
    d <- rbind(
      data.frame(ID = seq_len(nsub), TIME = 0, DV = NA, AMT = 320, EVID = 1, CMT = 1),
      data.frame(ID = s$id, TIME = s$time, DV = s$DV, AMT = 0, EVID = 0, CMT = 2)
    )
    d[order(d$ID, d$TIME, -d$EVID), ]
  }

  .viResidFits <- function(mod, d, control = emviControl(print = 0L)) {
    list(
      focei = suppressMessages(suppressWarnings(
        nlmixr2(mod, d, est = "focei", control = foceiControl(print = 0L))
      )),
      emvi = suppressMessages(suppressWarnings(
        nlmixr2(mod, d, est = if (isFALSE(control$pointEstimate)) "fbvi" else "emvi", control = control)
      ))
    )
  }

  test_that("emvi recovers a proportional residual error (#1462)", {
    skip_on_cran()
    d <- .viResidData(function(s, e) 0.1 * s$cp)
    mod <- function() {
      ini({
        tka <- log(1.2)
        tcl <- log(3)
        tv <- log(30)
        eta.ka ~ 0.1
        eta.cl ~ 0.1
        eta.v ~ 0.1
        prop.sd <- 0.2
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        d / dt(depot) <- -ka * depot
        d / dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ prop(prop.sd)
      })
    }
    f <- .viResidFits(mod, d)
    expect_equal(f$emvi$theta[["prop.sd"]], f$focei$theta[["prop.sd"]], tolerance = 0.05)
    expect_equal(unname(diag(f$emvi$omega)), unname(diag(f$focei$omega)), tolerance = 0.15)
    fMf <- suppressMessages(suppressWarnings(
      nlmixr2(mod, d, est = "emvi", control = emviControl(viFamily = "meanField", print = 0L))
    ))
    expect_equal(fMf$theta[["prop.sd"]], f$focei$theta[["prop.sd"]], tolerance = 0.05)
  })

  test_that("fbvi recovers a proportional residual error and FOCEi's SEs (#1462)", {
    skip_on_cran()
    d <- .viResidData(function(s, e) 0.1 * s$cp)
    mod <- function() {
      ini({
        tka <- log(1.2)
        tcl <- log(3)
        tv <- log(30)
        eta.ka ~ 0.1
        eta.cl ~ 0.1
        eta.v ~ 0.1
        prop.sd <- 0.2
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        d / dt(depot) <- -ka * depot
        d / dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ prop(prop.sd)
      })
    }
    f <- .viResidFits(mod, d, fbviControl(print = 0L))
    expect_identical(f$emvi$covMethod, "vi")
    expect_equal(f$emvi$theta[["prop.sd"]], f$focei$theta[["prop.sd"]], tolerance = 0.05)
    expect_equal(unname(diag(f$emvi$omega)), unname(diag(f$focei$omega)), tolerance = 0.15)
    ## mu-referenced SEs integrate the etas out, as FOCEi's do
    .mu <- c("tka", "tcl", "tv")
    expect_equal(f$emvi$parFixedDf[.mu, "SE"], f$focei$parFixedDf[.mu, "SE"], tolerance = 0.2)
  })

  test_that("emvi recovers a combined add + prop residual error (#1462)", {
    skip_on_cran()
    d <- .viResidData(function(s, e) sqrt(0.3^2 + (0.1 * s$cp)^2))
    mod <- function() {
      ini({
        tka <- log(1.2)
        tcl <- log(3)
        tv <- log(30)
        eta.ka ~ 0.1
        eta.cl ~ 0.1
        eta.v ~ 0.1
        add.sd <- 0.5
        prop.sd <- 0.2
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        d / dt(depot) <- -ka * depot
        d / dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd) + prop(prop.sd)
      })
    }
    f <- .viResidFits(mod, d)
    expect_equal(f$emvi$theta[c("add.sd", "prop.sd")], f$focei$theta[c("add.sd", "prop.sd")], tolerance = 0.05)
  })

  test_that("emvi recovers an eta on the residual sd (#1462)", {
    skip_on_cran()
    d <- .viResidData(function(s, e) 0.5 * exp(e))
    mod <- function() {
      ini({
        tka <- log(1.2)
        tcl <- log(3)
        tv <- log(30)
        eta.ka ~ 0.1
        eta.cl ~ 0.1
        eta.v ~ 0.1
        add.sd <- 0.7
        eta.add.sd ~ 0.1
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        d / dt(depot) <- -ka * depot
        d / dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        sd2 <- add.sd * exp(eta.add.sd)
        cp ~ add(sd2)
      })
    }
    f <- .viResidFits(mod, d)
    expect_equal(f$emvi$theta[["add.sd"]], f$focei$theta[["add.sd"]], tolerance = 0.15)
    expect_equal(f$emvi$omega["eta.add.sd", "eta.add.sd"], f$focei$omega["eta.add.sd", "eta.add.sd"], tolerance = 0.2)
  })
  test_that("fbvi with a mu-referenced covariate theta stays finite (#1462)", {
    skip_on_cran()
    d <- .viResidData(function(s, e) 0.1 * s$cp)
    d$WT <- 50 + 10 * (d$ID %% 5)
    mod <- function() {
      ini({
        tka <- log(1.2)
        tcl <- log(3)
        tv <- log(30)
        wtCl <- 0.1
        eta.ka ~ 0.1
        eta.cl ~ 0.1
        eta.v ~ 0.1
        prop.sd <- 0.3
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + wtCl * log(WT / 70) + eta.cl)
        v <- exp(tv + eta.v)
        d / dt(depot) <- -ka * depot
        d / dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ prop(prop.sd)
      })
    }
    f <- .viResidFits(mod, d, fbviControl(print = 0L))
    expect_true(all(is.finite(f$emvi$theta)))
    expect_equal(f$emvi$theta[["prop.sd"]], f$focei$theta[["prop.sd"]], tolerance = 0.1)
    expect_true(all(is.finite(f$emvi$parFixedDf$SE[!is.na(f$emvi$parFixedDf$SE)])))
  })
})
