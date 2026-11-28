# A posthoc (maxOuterIterations = 0) fast fit evaluates the analytic gradient once before
# the final objective.  That diagnostic pass must not spend a subject's one eta nudge
# cascade: warfarin q = 3, subject 27 (no absorption-phase data, two eta.ka basins) needs
# it to leave the spurious basin its first trust solve stops in.  Before the fix the fast
# and adjoint posthoc objectives were 2.1 above fast = FALSE.
nmTest({
  test_that("posthoc fast fit keeps the eta nudge (objective as fast = FALSE)", {
    skip_if_not_installed("nlmixr2data")
    .m <- function() {
      ini({
        tka <- 0.4; tcl <- -2.0; tv <- 2.0; tq <- -1.5; tvp <- 2.5; wtcl <- 0.75; wtv <- 1
        wtq <- 0.75; wtvp <- 1; eta.ka ~ 0.3; eta.cl ~ 0.1; eta.v ~ 0.05; add.sd <- 0.3; prop.sd <- 0.1
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl) * (WT / 70)^wtcl
        v <- exp(tv + eta.v) * (WT / 70)^wtv; q <- exp(tq) * (WT / 70)^wtq
        vp <- exp(tvp) * (WT / 70)^wtvp
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center - q / v * center + q / vp * periph
        d/dt(periph) <- q / v * center - q / vp * periph
        cp <- center / v
        cp ~ add(add.sd) + prop(prop.sd)
      })
    }
    d <- nlmixr2data::warfarin
    d <- d[d$dvid == "cp" & d$id == 27, c("id", "time", "amt", "dv", "evid", "wt")]
    names(d) <- c("ID", "TIME", "AMT", "DV", "EVID", "WT")
    .fit <- function(...) {
      suppressMessages(suppressWarnings(nlmixr2(
        .m,
        d,
        "focei",
        foceiControl(print = 0L, covMethod = "", calcTables = FALSE, maxOuterIterations = 0L, ...)
      )))
    }
    .plain <- .fit()
    .fast <- .fit(fast = TRUE)
    # the global mode (slow absorption); the spurious one is eta.ka ~ -0.8, -2LL ~ 51.59
    expect_lt(.plain$etaObf$eta.ka, -2)
    expect_lt(.fast$etaObf$eta.ka, -2)
    expect_equal(.fast$objf, .plain$objf, tolerance = 1e-6)
    expect_equal(.fast$etaObf$eta.ka, .plain$etaObf$eta.ka, tolerance = 1e-6)
  })
})
