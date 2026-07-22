## Phase 0 of the ELS residual step (plans/vae-els-residual.md): the ELBO step
## carries the per-observation residual VARIANCE r alongside the prediction f.
## grabRFmatFromInner already computed both -- only f was kept -- so this costs no
## extra solve.  An ELS residual objective is sum[(y-f)^2/r + log r], so it needs
## r per observation; these tests pin that r is what it claims to be.
##
## rx_r_ is the VARIANCE (sigma^2), not the SD: nlmixr2's add()/prop() are
## SD-scale parameters that get squared into it.

nmTest({
  .addMod <- function() {
    ini({ lka <- 0.45; lcl <- 1; lv <- 3.45; eta.ka ~ 0.6; add.err <- 0.7 })
    model({ ka <- exp(lka + eta.ka); cl <- exp(lcl); v <- exp(lv)
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.err) })
  }
  .propMod <- function() {
    ini({ lka <- 0.45; lcl <- 1; lv <- 3.45; eta.ka ~ 0.6; prop.err <- 0.15 })
    model({ ka <- exp(lka + eta.ka); cl <- exp(lcl); v <- exp(lv)
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ prop(prop.err) })
  }
  .combMod <- function() {
    ini({ lka <- 0.45; lcl <- 1; lv <- 3.45; eta.ka ~ 0.6
      add.err <- 0.7; prop.err <- 0.15 })
    model({ ka <- exp(lka + eta.ka); cl <- exp(lcl); v <- exp(lv)
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.err) + prop(prop.err) })
  }

  ## one ELBO step at fixed encoder params; returns the per-subject (f, r)
  .step <- function(mod) {
    ui <- rxode2::assertRxUi(mod)
    ctl <- vaeControl(print = 0L, calcTables = FALSE)
    prep <- .vaeDataPrep(ui, nlmixr2data::theo_sd, ctl)
    N <- prep$N; zDim <- prep$zDim
    innerEnv <- .vaeInnerSetup(ui, nlmixr2data::theo_sd, matrix(0, N, zDim), ctl)
    on.exit(.vaeInnerFree(), add = TRUE)
    .testSeed(3)
    params <- .vaeEncoderInitParams(zDim, 8L, ncol(prep$covIn), prep$zPop, rep(0.1, zDim))
    eps <- matrix(stats::rnorm(N * zDim), N, zDim)
    st <- .vaeElboStepInner(params, prep, innerEnv, prep$zPop, prep$omega, prep$a,
                            1, eps, ctl)
    list(st = st, prep = prep)
  }

  test_that("the ELBO step returns r alongside f, one per observation", {
    skip_on_cran()
    r <- .step(.addMod())
    expect_true(!is.null(r$st$rvar))
    expect_length(r$st$rvar, r$prep$N)
    ## same shape as the predictions, subject by subject
    for (i in seq_len(r$prep$N)) {
      expect_length(r$st$rvar[[i]], length(r$st$preds[[i]]))
      expect_true(all(is.finite(r$st$rvar[[i]])))
      expect_true(all(r$st$rvar[[i]] > 0))   # a variance
    }
  })

  ## The residual parameter reaches the solve through the inner problem's
  ## parameter scaling (foceiSetupTheta_/updateTheta), whose round trip is not
  ## bit-exact: the recovered SD sits ~2e-5 ABOVE the ini value (measured at both
  ## 0.7 -> 0.7000153 and 2.5 -> 2.500021, i.e. a roughly constant absolute
  ## offset, not a relative one).  That is far below any estimation tolerance and
  ## irrelevant to an ELS objective, so these check the FORM of r with a
  ## tolerance rather than exact equality (the offset is on the SD, so it roughly
  ## doubles in the variance).
  test_that("r is the additive VARIANCE (add.err^2), not the SD", {
    skip_on_cran()
    r <- .step(.addMod())
    rr <- unlist(r$st$rvar)
    ## constant across observations -- an additive variance does not depend on f
    expect_lt(diff(range(rr)), 1e-8)
    expect_equal(mean(rr), 0.7^2, tolerance = 1e-4)
    ## and it is the SQUARE: the SD would be 0.7, an order of magnitude off
    expect_gt(abs(mean(rr) - 0.7), 0.2)
  })

  test_that("r is the proportional variance (f*prop.err)^2 where f > 0", {
    skip_on_cran()
    r <- .step(.propMod())
    f <- unlist(r$st$preds); rr <- unlist(r$st$rvar)
    ok <- f > 0
    expect_gt(sum(ok), 0)
    expect_equal(rr[ok], (f[ok] * 0.15)^2, tolerance = 1e-3)
    ## f-dependent, unlike the additive case
    expect_gt(diff(range(rr)), 1e-8)
  })

  test_that("r is EXACTLY zero where f is zero (the throttle is not in this column)", {
    skip_on_cran()
    ## theo_sd has predose observations where the prediction is exactly 0, and a
    ## proportional variance is then literally 0.  The f-throttle that keeps the
    ## likelihood finite (handleF) is applied when the likelihood is EVALUATED,
    ## not baked into the stored rx_r_ column -- so the raw r carried here does
    ## contain zeros.
    ##
    ## This is a constraint on the ELS step (plans/vae-els-residual.md): an
    ## objective of the form sum[(y-f)^2/r + log r] must apply the throttle
    ## itself rather than consume this column raw, or it divides by zero and
    ## takes log(0) at every predose record of a proportional model.  The rule is
    ## `if (r == 0) r = 1`, matching handleF's treatment of a zero prediction.
    r <- .step(.propMod())
    f <- unlist(r$st$preds); rr <- unlist(r$st$rvar)
    zero <- f == 0
    expect_gt(sum(zero), 0)
    expect_true(all(rr[zero] == 0))
  })

  test_that("r is the combined variance add.err^2 + (f*prop.err)^2", {
    skip_on_cran()
    ## vaeControl() defaults to combined2, i.e. variances add
    r <- .step(.combMod())
    f <- unlist(r$st$preds); rr <- unlist(r$st$rvar)
    ok <- f > 0
    expect_equal(rr[ok], 0.7^2 + (f[ok] * 0.15)^2, tolerance = 1e-3)
    ## combined1 would be (add + f*prop)^2, which differs by the cross term
    expect_false(isTRUE(all.equal(rr[ok], (0.7 + f[ok] * 0.15)^2, tolerance = 1e-3)))
  })
})
