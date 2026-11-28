# FOCE "nonmem" eta-hat Newton (foceEbeNewton): the frozen-variance objective of a
# warfarin subject has a flat, non-convex eta.ka valley between the inner optimizer's
# stopping point and the mode.  Judged by |S| alone (and with the plain Newton step on an
# indefinite Hf) the Newton stalled and the subject fell back to the per-subject finite
# difference (fit$env$nOuterFdInd), which took minutes per gradient.
test_that("FOCE nonmem Newton reaches the frozen-variance mode in a non-convex valley", {
  skip_on_cran()
  d <- nlmixr2data::warfarin
  d <- d[d$dvid == "cp", c("id", "time", "amt", "dv", "evid", "wt")]
  names(d) <- c("ID", "TIME", "AMT", "DV", "EVID", "WT")
  d <- d[d$ID %in% unique(d$ID)[1:8], ]
  m <- function() {
    ini({
      tka <- 0.8512063309; tcl <- -2.069603774; tv <- 2.002401668; tq <- -2.912791681
      tvp <- 2.287699806; wtcl <- 0.7662071484; wtv <- 0.7621879747; wtq <- 1.079585503
      wtvp <- 0.8289089904
      eta.ka ~ 1.499182933; eta.cl ~ 0.2363133025; eta.v ~ 0.03219586673
      add.sd <- 0.2709709107; prop.sd <- 0.09466092473
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
  ctl <- foceiControl(
    print = 0L,
    covMethod = "",
    calcTables = FALSE,
    interaction = FALSE,
    maxOuterIterations = 1L,
    fast = TRUE,
    rxControl = rxode2::rxControl(method = "liblsoda", cores = 1L)
  )
  f <- suppressMessages(suppressWarnings(nlmixr2(m, d, "focei", ctl)))
  expect_equal(f$env$nOuterFdInd, 0L)
  expect_equal(f$env$nNewtonFail[["maxit"]], 0L)
})
