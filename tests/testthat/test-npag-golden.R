## Golden comparison against Pmetrics NPAG (the Rust reference engine).  We do NOT
## run Rust/Pmetrics at test time (CRAN constraint); the reference population means
## are frozen in fixtures/npag/pmetrics_reference.csv (see that file for provenance).
## This validates that native npag recovers the same nonparametric mixing
## distribution as the oracle on Theophylline PK.  The -2LL uses a different constant
## convention than Pmetrics, so the comparison is on the population means, not the
## objective.  Real fit -> weekly slow batch.

nmTest({
  test_that("est='npag' matches the Pmetrics NPAG population means on Theophylline", {
    skip_if_not_installed("rxode2")
    .ref <- utils::read.csv(
      testthat::test_path("fixtures", "npag", "pmetrics_reference.csv"),
      comment.char = "#")
    .ref <- .ref[.ref$model == "theo", ]
    .m <- setNames(.ref$pmetrics_mean, .ref$param)   # ka, v, ke

    one.cmt <- function() {
      ini({ tka <- log(1.5); tv <- log(32); tke <- log(0.08)
        eta.ka ~ 0.3; eta.v ~ 0.1; eta.ke ~ 0.1; add.sd <- 0.4 })
      model({ ka <- exp(tka + eta.ka); v <- exp(tv + eta.v); ke <- exp(tke + eta.ke)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - ke * center
        cp <- center / v; cp ~ add(add.sd) })
    }
    f <- nlmixr2(one.cmt, nlmixr2data::theo_sd, est = "npag",
                 control = npagControl(points = 512L, cycles = 20L,
                                       gammaOptimize = TRUE, seed = 1L,
                                       calcTables = FALSE))
    .sp <- f$env$npagSupport; .wt <- f$env$npagWeights
    .th <- f$theta
    # population mean of each parameter = weighted mean over the support of
    # exp(theta + eta) (npag keeps the theta at reference; the location is in the
    # support distribution).  Columns of npagSupport are ordered ka, v, ke.
    .pm <- c(
      ka = sum(.wt * exp(as.numeric(.th[["tka"]]) + .sp[, 1])),
      v  = sum(.wt * exp(as.numeric(.th[["tv"]])  + .sp[, 2])),
      ke = sum(.wt * exp(as.numeric(.th[["tke"]]) + .sp[, 3])))

    # V and Ke are well-identified -> tight agreement with the oracle (~<10%).
    expect_lt(abs(.pm[["v"]]  - .m[["v"]])  / .m[["v"]],  0.10)
    expect_lt(abs(.pm[["ke"]] - .m[["ke"]]) / .m[["ke"]], 0.10)
    # Ka (absorption) is weakly identified in oral PK and grid-dependent, so allow a
    # looser band; it should still be the right order of magnitude.
    expect_lt(abs(.pm[["ka"]] - .m[["ka"]]) / .m[["ka"]], 0.40)
    # nonparametric support of comparable size (oracle nspp = 12)
    expect_true(nrow(.sp) >= 5L && nrow(.sp) <= 25L)
  })

  test_that("est='npag' matches the Pmetrics NPAG population means on Warfarin PK/PD", {
    skip_if_not_installed("rxode2")
    skip_if_not_installed("nlmixr2data")
    .ref <- utils::read.csv(
      testthat::test_path("fixtures", "npag", "pmetrics_reference.csv"),
      comment.char = "#")
    .ref <- .ref[.ref$model == "warfarin", ]
    .m <- setNames(.ref$pmetrics_mean, .ref$param)

    # transit-absorption PK + Emax turnover PD (2 endpoints: cp, pca).  Bounded ini
    # estimates (matching the Pmetrics ab() ranges) + gridBounds="ini" keep the 8-eta
    # grid in range -- an unbounded box collapses this high-dimensional support.
    warf <- function() {
      ini({ tktr <- log(c(0.5, 1.0, 3.0)); tka <- log(c(0.3, 1.0, 2.0))
        tv <- log(c(4, 8, 15)); tcl <- log(c(0.05, 0.15, 0.4))
        temax <- logit(c(0.8, 0.9, 0.99)); tec50 <- log(c(0.3, 1.0, 3.0))
        tkout <- log(c(0.02, 0.05, 0.15)); te0 <- log(c(70, 100, 130))
        eta.ktr ~ 0.5; eta.ka ~ 0.5; eta.v ~ 0.5; eta.cl ~ 0.5
        eta.emax ~ 0.5; eta.ec50 ~ 0.5; eta.kout ~ 0.5; eta.e0 ~ 0.5
        prop.sd <- 0.1; add.sd <- 2.0 })
      model({ ktr <- exp(tktr + eta.ktr); ka <- exp(tka + eta.ka)
        v <- exp(tv + eta.v); cl <- exp(tcl + eta.cl)
        emax <- expit(temax + eta.emax); ec50 <- exp(tec50 + eta.ec50)
        kout <- exp(tkout + eta.kout); e0 <- exp(te0 + eta.e0)
        d/dt(depot) <- -ktr * depot
        d/dt(gut) <- ktr * depot - ka * gut
        d/dt(center) <- ka * gut - cl / v * center
        DCP <- center / v
        d/dt(effect) <- -e0 * kout * (emax * DCP / (ec50 + DCP)) - kout * effect
        cp <- center / v; pca <- effect + e0
        cp ~ prop(prop.sd); pca ~ add(add.sd) })
    }
    f <- nlmixr2(warf, nlmixr2data::warfarin, est = "npag",
                 control = npagControl(points = 400L, cycles = 8L, seed = 1L,
                                       gammaOptimize = FALSE, gridBounds = "ini",
                                       calcTables = FALSE))
    .sp <- f$env$npagSupport; .wt <- f$env$npagWeights; .th <- f$theta
    .tn <- c(ktr="tktr", ka="tka", v="tv", cl="tcl", emax="temax",
             ec50="tec50", kout="tkout", e0="te0")
    .pm <- vapply(seq_along(.tn), function(j) {
      .tr <- if (names(.tn)[j] == "emax") plogis else exp   # emax is on the logit scale
      sum(.wt * .tr(as.numeric(.th[[.tn[j]]]) + .sp[, j]))
    }, numeric(1))
    names(.pm) <- names(.tn)
    # a real 8-parameter nonparametric PK/PD fit is not degenerate ...
    expect_true(nrow(.sp) >= 8L)
    # ... and its population means track the Pmetrics oracle (within 20%, looser than
    # theo -- more parameters, an Emax PD, and a different grid/engine).
    for (.p in names(.m)) {
      expect_lt(abs(.pm[[.p]] - .m[[.p]]) / .m[[.p]], 0.20)
    }
  })
})
