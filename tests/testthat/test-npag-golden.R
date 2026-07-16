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
})
