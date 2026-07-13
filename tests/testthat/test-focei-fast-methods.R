# The "*f" convenience methods: each is its base method with fast=TRUE default.

nmTest({
  .fastm_one_cmt <- function() {
    ini({
      tka <- log(1.5); tcl <- log(2.7); tv <- log(31.5)
      eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
      d/dt(depot)  <- -ka * depot
      d/dt(center) <-  ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }

  test_that("all nine *f methods are registered and default to fast=TRUE", {
    .fm <- c("focef", "focepf", "foceif", "mfocef", "mfocepf", "mfoceif",
             "ifocef", "ifocepf", "ifoceif")
    expect_true(all(.fm %in% nlmixr2AllEst()))
    # each *f control forces fast=TRUE, even from an empty control
    for (m in .fm) {
      .ctl <- getValidNlmixrCtl(structure(list(NULL), class = c(m, "getValidNlmixrControl")))
      expect_true(isTRUE(.ctl$fast), info = m)
    }
    # a derivative-free outerOpt still downgrades fast (runs in the base constructor)
    expect_warning(
      .ctl <- getValidNlmixrCtl(structure(list(foceiControl(outerOpt = "bobyqa")),
                                          class = c("foceif", "getValidNlmixrControl"))),
      "derivative-free")
    expect_false(isTRUE(.ctl$fast))
  })

  test_that("mu-referenced families use the analytic gradient (fast matches FD)", {
    skip_on_cran()
    skip_on_ci()
    skip_if_not_installed("nlmixr2data")
    d <- nlmixr2data::theo_sd
    # mu-referenced covariate coefficient cl.wt (regression-updated); the analytic
    # gradient must exclude it from the outer optimizer's parameter set (A1b).
    mc <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; cl.wt <- 0.01; eta.ka ~ 0.3; eta.cl ~ 0.1; add.sd <- 0.7 })
      model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl + cl.wt * (WT - 70)); v <- exp(tv)
              d/dt(depot) <- -ka * depot; d/dt(center) <- ka * depot - cl / v * center
              cp <- center / v; cp ~ add(add.sd) })
    }
    for (est in c("mfocei", "ifocei", "mfoce")) {
      f0 <- suppressMessages(suppressWarnings(nlmixr2(mc, d, est, foceiControl(print = 0L, covMethod = "", fast = FALSE))))
      fF <- suppressMessages(suppressWarnings(nlmixr2(mc, d, est, foceiControl(print = 0L, covMethod = "", fast = TRUE))))
      expect_equal(fF$objf, f0$objf, tolerance = 0.05, info = est)
      # analytic gradient actually consumed on the mu-profiled parameter set
      .gt <- fF$parHistData$type
      expect_gt(sum(.gt == "Analytic Gradient"), 0)
      expect_equal(sum(.gt %in% c("Gill83 Gradient", "Mixed Gradient",
                                  "Forward Difference", "Central Difference")), 0)
      expect_match(fF$extra, "grad: analytic", info = est)
      expect_match(fF$extra, if (grepl("^irls", est)) "mu: irls" else "mu: lin", info = est)
    }
  })

  test_that("est='foceif' equals est='focei' + fast=TRUE", {
    skip_on_cran()
    skip_on_ci()
    skip_if_not_installed("nlmixr2data")
    d <- nlmixr2data::theo_sd
    fF   <- suppressMessages(nlmixr2(.fastm_one_cmt, d, "foceif", foceiControl(print = 0L, covMethod = "")))
    fRef <- suppressMessages(nlmixr2(.fastm_one_cmt, d, "focei", foceiControl(print = 0L, covMethod = "", fast = TRUE)))
    expect_equal(fF$objf, fRef$objf, tolerance = 0.02)
    expect_equal(unname(fixef(fF)), unname(fixef(fRef)), tolerance = 1e-2)
    # both consumed the analytic gradient and default to lbfgsb3c
    expect_match(fF$extra, "outer: lbfgsb3c; grad: analytic")
    expect_match(fRef$extra, "outer: lbfgsb3c; grad: analytic")
  })
})
