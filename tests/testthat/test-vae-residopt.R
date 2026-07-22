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

  .bcMod <- function() {
    ini({ lka <- 0.45; lcl <- 1; lv <- 3.45; eta.ka ~ 0.6; eta.cl ~ 0.3
      add.err <- 0.7; lam <- 1 })
    model({ ka <- exp(lka + eta.ka); cl <- exp(lcl + eta.cl); v <- exp(lv)
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.err) + boxCox(lam) })
  }
  .yjMod <- function() {
    ini({ lka <- 0.45; lcl <- 1; lv <- 3.45; eta.ka ~ 0.6; eta.cl ~ 0.3
      add.err <- 0.7; lam <- 1 })
    model({ ka <- exp(lka + eta.ka); cl <- exp(lcl + eta.cl); v <- exp(lv)
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.err) + yeoJohnson(lam) })
  }

  test_that("a boxCox lambda is bounded to (-2, 2)", {
    ## unbounded in the ini() block, and only meaningful on a narrow interval;
    ## SAEM likewise maps lambda through a bounded transform
    p <- .vaeDataPrep(rxode2::assertRxUi(.bcMod()), nlmixr2data::theo_sd, vaeControl())
    i <- match("lam", p$regressNames)
    expect_false(is.na(i))
    expect_equal(p$regressLower[i], -2)
    expect_equal(p$regressUpper[i], 2)
  })

  test_that("a boxCox lambda is estimated and improves the fit", {
    skip_on_cran()
    m <- .fit(.bcMod(), "moment")
    o <- .fit(.bcMod(), "twoStage")
    expect_equal(m$theta[["lam"]], 1, tolerance = 1e-8)   # frozen without the optimizer
    expect_gt(abs(o$theta[["lam"]] - 1), 1e-3)
    expect_gte(o$theta[["lam"]], -2)
    expect_lte(o$theta[["lam"]], 2)
    expect_lt(o$objf, m$objf)
    ## NON-DEGENERACY.  This model previously converged to add.err = 0: the
    ## likelihood floors a zero variance (r == 0 -> r = 1) to stay finite, which
    ## makes a collapsed residual look attractive instead of forbidden.  An
    ## objective-only check does NOT catch it -- the degenerate fit still beat the
    ## moment estimator (68.7 against 181.6) -- so assert the scale itself.
    expect_gt(o$theta[["add.err"]], 0.01)
  })

  test_that("a residual scale parameter cannot reach zero", {
    ## the bound that makes the above hold: absolute scales are floored relative
    ## to the spread of the data, so bobyqa can never propose a zero variance
    p <- .vaeDataPrep(rxode2::assertRxUi(.bcMod()), nlmixr2data::theo_sd, vaeControl())
    lo <- p$regressLower[match("add.err", p$regressNames)]
    expect_gt(lo, 0)
    expect_lt(lo, 0.01)   # far below any plausible estimate
  })

  test_that("a yeoJohnson lambda is estimated and improves the fit", {
    skip_on_cran()
    ## Only dv is transformed in the objective: f leaves the solve ALREADY on the
    ## transformed scale.  Transforming both double-transforms and misfits badly
    ## -- it regressed this model from 131.8 to 383.6 before that was corrected,
    ## which is what makes it worth a test rather than a comment.
    m <- .fit(.yjMod(), "moment")
    o <- .fit(.yjMod(), "twoStage")
    expect_equal(m$theta[["lam"]], 1, tolerance = 1e-8)
    expect_gt(abs(o$theta[["lam"]] - 1), 1e-3)
    expect_gte(o$theta[["lam"]], -2)
    expect_lte(o$theta[["lam"]], 2)
    expect_lt(o$objf, m$objf)
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
