## Phase 1/3 of plans/vae-els-residual.md: residual-error parameters are
## estimated by two-stage block coordinate descent --
##   stage 1  non-mu-referenced structural thetas, residuals held (driven by dv-f)
##   stage 2  residuals alone, against sum[(y-f)^2/r + log r] over the CACHED
##            (y, f) pairs, so no ODE re-solve and f is fixed
##
## These tests assert the parameters MOVE, not merely that they are finite.  The
## behavior being fixed is a silent freeze: `pow`, `lnorm` and friends were
## classified "other" and left at their ini() value, which any finiteness-only
## check would have passed.

nmTest({
  .powMod <- function() {
    ini({ lka <- 0.45; lcl <- 1; lv <- 3.45; eta.ka ~ 0.6; eta.cl ~ 0.3
      prop.err <- 0.3; pw <- 0.8 })
    model({ ka <- exp(lka + eta.ka); cl <- exp(lcl + eta.cl); v <- exp(lv)
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ pow(prop.err, pw) })
  }
  .lnMod <- function() {
    ini({ lka <- 0.45; lcl <- 1; lv <- 3.45; eta.ka ~ 0.6; eta.cl ~ 0.3
      add.err <- 0.5 })
    model({ ka <- exp(lka + eta.ka); cl <- exp(lcl + eta.cl); v <- exp(lv)
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ lnorm(add.err) })
  }
  .combMod <- function() {
    ini({ lka <- 0.45; lcl <- 1; lv <- 3.45; eta.ka ~ 0.6; eta.cl ~ 0.3
      add.err <- 0.7; prop.err <- 0.1 })
    model({ ka <- exp(lka + eta.ka); cl <- exp(lcl + eta.cl); v <- exp(lv)
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.err) + prop(prop.err) })
  }
  .fit <- function(mod, ro) {
    suppressMessages(suppressWarnings(nlmixr2(
      mod, nlmixr2data::theo_sd, est = "vae",
      control = vaeControl(print = 0L, calcTables = FALSE, residOptimize = ro,
                           itersBurnIn = 40L, iters = 80L, klWarmup = 30L,
                           gammaIter = 60L))))
  }

  test_that("residOptimize defaults to twoStage", {
    expect_equal(vaeControl()$residOptimize, "twoStage")
    expect_equal(vaeControl(residOptimize = "moment")$residOptimize, "moment")
    expect_error(vaeControl(residOptimize = "els"))
  })

  test_that("pow() is estimated, not frozen at its ini() value", {
    skip_on_cran()
    m <- .fit(.powMod(), "moment")
    o <- .fit(.powMod(), "twoStage")
    ## the moment estimator has no closed form here and silently holds both
    expect_equal(m$theta[["prop.err"]], 0.3, tolerance = 1e-8)
    expect_equal(m$theta[["pw"]], 0.8, tolerance = 1e-8)
    ## two-stage moves them and improves the objective
    expect_gt(abs(o$theta[["prop.err"]] - 0.3), 1e-3)
    expect_gt(abs(o$theta[["pw"]] - 0.8), 1e-3)
    expect_lt(o$objf, m$objf)
  })

  test_that("lnorm() is estimated on the log scale, not frozen", {
    skip_on_cran()
    m <- .fit(.lnMod(), "moment")
    o <- .fit(.lnMod(), "twoStage")
    expect_equal(m$theta[["add.err"]], 0.5, tolerance = 1e-8)
    expect_gt(abs(o$theta[["add.err"]] - 0.5), 1e-3)
    ## the frozen value is badly wrong for a log-scale residual, so the
    ## improvement here is large rather than marginal
    expect_lt(o$objf, 0.5 * m$objf)
  })

  test_that("two-stage beats the moment estimator on a combined model", {
    skip_on_cran()
    ## add and prop are near-collinear; a single JOINT solve against the full
    ## outer objective diverges (that is why stage 2 is a separate ELS solve at
    ## fixed f).  Two-stage must beat the closed-form estimator here.
    m <- .fit(.combMod(), "moment")
    o <- .fit(.combMod(), "twoStage")
    expect_true(all(is.finite(c(o$theta[["add.err"]], o$theta[["prop.err"]]))))
    expect_gt(o$theta[["add.err"]], 0)
    expect_gt(o$theta[["prop.err"]], 0)
    expect_lt(o$objf, m$objf)
  })
})
