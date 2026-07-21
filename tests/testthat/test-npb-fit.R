## est="npb" end-to-end: the truncated stick-breaking blocked-Gibbs sampler
## returns a nlmixr2FitData with the posterior mixing distribution, per-subject
## posterior-mean etas, and posterior draws of the population mean (credible
## intervals).  Real fit -> weekly slow batch.

nmTest({
  .npbMod <- function() {
    ini({ tka <- log(1.5); tv <- log(32); tke <- log(0.08)
      eta.ka ~ 0.3; eta.v ~ 0.1; eta.ke ~ 0.1; add.sd <- 0.7 })
    model({ ka <- exp(tka + eta.ka); v <- exp(tv + eta.v); ke <- exp(tke + eta.ke)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - ke * center
      cp <- center / v
      cp ~ add(add.sd) })
  }
  .ctl <- function() npbControl(points = 20L, burnin = 100L, nsamp = 100L, seed = 42L)

  test_that("est='npb' returns a valid fit with a posterior distribution", {
    f <- nlmixr2(.npbMod, nlmixr2data::theo_sd, est = "npb", control = .ctl())
    expect_s3_class(f, "nlmixr2FitData")
    expect_true(inherits(f, "nlmixr2.npb"))
    expect_true(is.finite(as.numeric(f$objf)))
    expect_true(all(c("tka", "tv", "tke", "add.sd") %in% rownames(f$parFixed)))
    N <- length(unique(nlmixr2data::theo_sd$ID))
    expect_equal(ncol(f$env$npbSupport), 3L)
    expect_equal(sum(f$env$npbWeights), 1, tolerance = 1e-6)
    expect_equal(nrow(f$env$npbPosteriorEta), N)
    # posterior draws of the population mean -> credible intervals
    md <- f$env$npbMeanDraws
    expect_equal(dim(md), c(100L, 3L))
    ci <- apply(md, 2, quantile, c(0.025, 0.975))
    expect_true(all(ci[1, ] < ci[2, ]))          # non-degenerate intervals
    expect_equal(f$env$npbK, 20L)
    # eta-space outputs carry the eta names (columns for the matrices, rows for
    # the per-eta R-hat vector)
    .en <- c("eta.ka", "eta.v", "eta.ke")
    expect_equal(colnames(f$env$npbSupport), .en)
    expect_equal(colnames(f$env$npbPosteriorEta), .en)
    expect_equal(colnames(f$env$npbMeanDraws), .en)
    expect_equal(rownames(f$env$npbRhat), .en)
    # npb records + saves its iteration history through the shared scale.h printer
    # (one Scaled + one Back-Transformed row per chain-0 sweep = burnin + nsamp)
    ph <- f$env$parHistData
    expect_s3_class(ph, "data.frame")
    expect_true(all(c("iter", "type", "objf") %in% names(ph)))
    expect_equal(sum(ph$type == "Scaled"), 100L + 100L)
  })

  test_that("est='npb' is reproducible for a fixed seed", {
    f1 <- nlmixr2(.npbMod, nlmixr2data::theo_sd, est = "npb", control = .ctl())
    f2 <- nlmixr2(.npbMod, nlmixr2data::theo_sd, est = "npb", control = .ctl())
    expect_equal(as.numeric(f1$objf), as.numeric(f2$objf), tolerance = 1e-8)
  })

  test_that("est='npb' multi-chain reports Gelman-Rubin R-hat", {
    f <- nlmixr2(.npbMod, nlmixr2data::theo_sd, est = "npb",
                 control = npbControl(points = 20L, burnin = 100L, nsamp = 150L,
                                      nchains = 3L, seed = 42L))
    expect_equal(f$env$npbNchains, 3L)
    # pooled draws across chains, and one R-hat per eta near 1 at convergence
    expect_equal(nrow(f$env$npbMeanDraws), 3L * 150L)
    .rhat <- f$env$npbRhat
    expect_equal(length(.rhat), 3L)
    expect_true(all(is.finite(.rhat)) && all(.rhat >= 1 - 1e-8) && all(.rhat < 1.3))
  })

  test_that("mu-referenced sugar est='mnpb' fits", {
    f <- nlmixr2(.npbMod, nlmixr2data::theo_sd, est = "mnpb", control = .ctl())
    expect_s3_class(f, "nlmixr2FitData")
    expect_true(is.finite(as.numeric(f$objf)))
  })

  test_that("est='npb' does the residual/regressor optimization (like npag)", {
    # tke has no eta and is not mu-referenced -> a "regressor" the residual step
    # re-solves the ODE for; add.sd is a residual-error param.  Both are held under
    # residOptimize="none" and estimated under "alternate".
    .rMod <- function() {
      ini({ tka <- log(1.5); tv <- log(32); tke <- log(0.03)   # truth ~ log(0.08)
        eta.ka ~ 0.3; eta.v ~ 0.1; add.sd <- 0.4 })
      model({ ka <- exp(tka + eta.ka); v <- exp(tv + eta.v); ke <- exp(tke)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - ke * center
        cp <- center / v
        cp ~ add(add.sd) })
    }
    .fit <- function(mode)
      nlmixr2(.rMod, nlmixr2data::theo_sd, est = "npb",
              control = npbControl(points = 30L, burnin = 40L, nsamp = 30L,
                                   residOptimize = mode, seed = 1L))
    fNone <- .fit("none")
    fAlt  <- .fit("alternate")
    # "none" holds the regressor and the residual param at their initial values
    expect_equal(as.numeric(fNone$theta[["tke"]]), log(0.03), tolerance = 1e-6)
    expect_equal(as.numeric(fNone$theta[["add.sd"]]), 0.4, tolerance = 1e-6)
    # "alternate" recovers the structural regressor from its poor start (log 0.03 ->
    # near the log(0.08) truth) and moves the residual param off its start
    expect_true(as.numeric(fAlt$theta[["tke"]]) > log(0.03) + 0.5)
    expect_lt(abs(as.numeric(fAlt$theta[["tke"]]) - log(0.08)), 0.4)
    expect_false(isTRUE(all.equal(as.numeric(fAlt$theta[["add.sd"]]), 0.4)))
  })
})
