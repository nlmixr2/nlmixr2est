# numericGrad()'s outer finite-difference fallbacks: the small-gradient central
# confirmation (gradCalcCentralSmall) and the gradTrim clamp.  These only run
# under a gradient-based outerOpt, so every fit below pins one.
nmTest({
  test_that("gradTrim/gradCalcCentralSmall round-trip and are validated", {
    ctl <- foceiControl()
    expect_equal(ctl$gradTrim, Inf)
    expect_equal(ctl$gradCalcCentralSmall, 1e-4)
    expect_equal(ctl$gradCalcCentralLarge, 1e4)
    ctl2 <- foceiControl(gradTrim = 10, gradCalcCentralSmall = 1e-2)
    expect_equal(ctl2$gradTrim, 10)
    expect_equal(ctl2$gradCalcCentralSmall, 1e-2)
    expect_error(foceiControl(gradCalcCentralSmall = -1))
    expect_error(foceiControl(gradCalcCentralLarge = Inf))
  })

  test_that("the small-gradient central confirmation keeps the gradient scale", {
    skip_on_cran()
    ode.cmt <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45; add.sd <- 0.7
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        d/dt(depot) = -ka * depot
        d/dt(center) = ka * depot - cl / v * center
        cp = center / v
        cp ~ add(add.sd)
      })
    }
    run <- function(...) {
      suppressWarnings(suppressMessages(
        nlmixr2(ode.cmt, theo_sd, "focei",
                foceiControl(covMethod = "", calcTables = FALSE, print = 1L,
                             outerOpt = "nlminb", maxOuterIterations = 8L, ...))))
    }
    gradRows <- function(fit, lvl) {
      ph <- fit$parHistData
      ph[as.character(ph$type) == lvl, setdiff(names(ph), c("iter", "type", "objf")),
         drop = FALSE]
    }

    base <- run()
    # a plain forward difference is what the default takes
    fwd <- gradRows(base, "Forward Difference")
    expect_true(nrow(fwd) > 0)

    # raising gradCalcCentralSmall above every |forward difference| sends EVERY
    # component through the central-difference confirmation instead -- the branch
    # reports itself as a mixed gradient, which is the evidence it ran
    forced <- run(gradCalcCentralSmall = 1e10)
    mixed <- gradRows(forced, "Mixed Gradient")
    expect_true(nrow(mixed) > 0)
    expect_equal(nrow(gradRows(forced, "Forward Difference")), 0L)

    # the confirmation must return a derivative, not the objective it differences
    # against.  Overwriting the saved objective at theta+delta with the forward
    # gradient turned this into roughly -objf/(2*h): ~1e6 here, and every
    # component negative.
    expect_true(all(is.finite(unlist(mixed))))
    expect_lt(max(abs(unlist(mixed))), 1e3)
    expect_lt(max(abs(unlist(mixed))), 50 * max(abs(unlist(fwd))))
    expect_true(any(unlist(mixed) > 0))
  })

  # NOTE: this does not reproduce the `g < gradTrim` sign slip itself.  That needs
  # a component whose forward difference clears gradTrim while its recomputed
  # central difference lands back inside the band, which is not arrangeable from
  # a control setting.  It does pin the clamp bounds, so a lost or re-broken
  # lower clamp shows up here.
  test_that("a finite gradTrim bounds the outer gradient symmetrically", {
    skip_on_cran()
    one.cmt <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45; add.sd <- 0.7
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    trim <- 2
    fit <- suppressWarnings(suppressMessages(
      nlmixr2(one.cmt, theo_sd, "focei",
              foceiControl(covMethod = "", calcTables = FALSE, print = 1L,
                           outerOpt = "nlminb", maxOuterIterations = 8L,
                           gradTrim = trim))))
    ph <- fit$parHistData
    # gradTrim only applies to numericGrad(); the Gill83 iteration is a separate
    # path and is not clamped
    gr <- ph[as.character(ph$type) %in%
               c("Mixed Gradient", "Forward Difference", "Central Difference"),
             setdiff(names(ph), c("iter", "type", "objf")), drop = FALSE]
    expect_true(nrow(gr) > 0)
    expect_true(all(abs(unlist(gr)) <= trim + 1e-8))
    # clamping must not collapse the gradient onto one bound: both signs and
    # untouched interior values have to survive
    expect_true(any(unlist(gr) > 0))
    expect_true(any(abs(unlist(gr)) < trim - 1e-8))
  })
})
