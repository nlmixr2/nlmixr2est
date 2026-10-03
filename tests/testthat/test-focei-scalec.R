nmTest({
  # A fixed or mu-profiled theta is not an optimizer parameter, so the
  # parameters after it have a lower optimizer index than parameter index.
  .mod <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }

  # The scaleC the optimizer used for each parameter, from the first
  # evaluation that moved it: unscaled - init = (scaled - scaled init) * scaleC.
  .usedScaleC <- function(fit, pars) {
    .ph <- fit$parHistData
    .x <- as.matrix(.ph[.ph$type == "Scaled", pars])
    .u <- as.matrix(.ph[.ph$type == "Unscaled", pars])
    .dx <- sweep(.x, 2, .x[1, ])
    .at <- cbind(apply(.dx != 0, 2, which.max), seq_along(pars))
    (.u[.at] - .u[1, ]) / .dx[.at]
  }

  .ctl <- function(...) {
    foceiControl(print = 0, maxOuterIterations = 2L, covMethod = "", calcTables = FALSE, outerOpt = "lbfgsb3c", ...)
  }

  # add.sd gets 0.5 * |init|; an omega parameter, chol(omega^-1) with sqrt
  # diagonals, gets 1/|init| = omega^(1/4)
  .omegaC <- c(o1 = 0.6^0.25, o2 = 0.3^0.25, o3 = 0.1^0.25)

  test_that("mfocei scales each parameter by its own scaleC", {
    skip_on_cran()
    # the band misses add.sd's 0.35 and two of the omegas' values, and guards
    # none of them: scaleCband applies to linear thetas, in R
    .f <- suppressMessages(suppressWarnings(nlmixr(.mod, theo_sd, "mfocei", control = .ctl(scaleCband = c(0.8, 10)))))
    .want <- c(add.sd = 0.35, .omegaC)
    expect_equal(.usedScaleC(.f, names(.want)), .want, tolerance = 1e-6)
    expect_equal(.f$scaleInfo$scaleC, c(NA, NA, NA, unname(.want)), tolerance = 1e-6)
  })

  # tcl lies inside expit(., 1, 100) and prop.sd is a residual error; R gives
  # them 21.3 (inside expit's own band) and 0.5 * 0.1, both outside scaleCband
  .modR <- function() {
    ini({
      tka <- 0.45
      tcl <- 3
      tv <- 3.45
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      prop.sd <- 0.1
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- expit(tcl + eta.cl, 1, 100) / 30
      v <- exp(tv + eta.v)
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ prop(prop.sd)
    })
  }

  test_that("FOCEi scales each theta by the scaleC R gives it", {
    skip_on_cran()
    .want <- setNames(rxode2::rxode2(.modR)$scaleCtheta, c("tka", "tcl", "tv", "prop.sd"))
    expect_equal(unname(.want), c(1, 21.309126, 1, 0.05), tolerance = 1e-6)
    # FOCEi guarded tcl and prop.sd to scaleCband a second time, in C++, and
    # scaled them by |init|: 3 and 0.1
    .f <- suppressMessages(suppressWarnings(nlmixr(.modR, theo_sd, "focei", control = .ctl())))
    expect_equal(.usedScaleC(.f, names(.want)), .want, tolerance = 1e-6)
    expect_equal(.f$scaleInfo$scaleC[1:4], unname(.want), tolerance = 1e-6)
  })

  test_that("FOCEi scales a theta by the scaleC the user gives it", {
    skip_on_cran()
    # add.sd's 0.02 was replaced by |init| = 0.7
    .f <- suppressMessages(suppressWarnings(nlmixr(.mod, theo_sd, "focei", control = .ctl(scaleC = c(1, 1, 1, 0.02)))))
    .want <- c(tka = 1, tcl = 1, tv = 1, add.sd = 0.02)
    expect_equal(.usedScaleC(.f, names(.want)), .want, tolerance = 1e-6)
  })

  # beta multiplies a covariate that is 0 throughout, so no step finds a slope
  # for it (a population model: no etas to re-optimize between the legs)
  .modZ <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      beta <- 0.5
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl + beta * ZERO)
      v <- exp(tv)
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }
  .dZ <- theo_sd
  .dZ$ZERO <- 0

  test_that("a zero gradient is searched again with scaleC0, then 1/scaleC0", {
    skip_on_cran()
    # both searches run and the second one's scale stays; scaleCband turned
    # both into |init| = 0.5
    .f <- suppressMessages(suppressWarnings(nlmixr(.modZ, .dZ, "focei", control = .ctl(scaleC0 = 1000))))
    expect_equal(.f$scaleInfo$scaleC[4], 1e-3)
  })

  test_that("a fixed theta does not shift the scaleC of the parameters after it", {
    skip_on_cran()
    .m <- .mod |> rxode2::ini(tka = fix(0.45))
    .f <- suppressMessages(suppressWarnings(nlmixr(.m, theo_sd, "focei", control = .ctl(literalFix = FALSE))))
    .want <- c(tcl = 1, tv = 1, add.sd = 0.35, .omegaC)
    expect_equal(.usedScaleC(.f, names(.want)), .want, tolerance = 1e-6)
  })

  test_that("a fixed theta does not shift the bound codes of the parameters after it", {
    skip_on_cran()
    # tcl's lower bound is above its estimate, so the fit ends on it and the
    # covariance step must report the boundary.  tcl's bound code (lower and
    # upper) is its own, not that of tv after it (none), which would leave only
    # its upper bound checked.  add.sd comes before tcl so that its own
    # lower-bound check does not carry over to tcl.
    .m <- function() {
      ini({
        tka <- fix(0.45)
        add.sd <- 0.7
        tcl <- c(1.1, 1.2, 5)
        tv <- 3.45
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        d / dt(depot) <- -ka * depot
        d / dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }
    .f <- suppressMessages(suppressWarnings(nlmixr(
      .m,
      theo_sd,
      "focei",
      control = foceiControl(print = 0, calcTables = FALSE, literalFix = FALSE)
    )))
    expect_equal(fixef(.f)[["tcl"]], 1.1, tolerance = 1e-4)
    expect_match(.f$covMethod, "^Boundary issue")
    expect_match(.f$covMethod, "\"tcl\"")
  })

  .agqObjf <- function(scaleC) {
    .ctl <- agqControl(print = 0, maxOuterIterations = 0L, covMethod = "", calcTables = FALSE, scaleC = scaleC)
    suppressMessages(suppressWarnings(nlmixr(.mod, theo_sd, "agq", control = .ctl)))$objf
  }

  test_that("a scaleC longer than the parameter vector is only read up to it", {
    skip_on_cran()
    # the values past the parameters were copied over the AGQ nodes that follow
    # the scaleC block (and a long enough vector past the buffer, a crash)
    expect_identical(.agqObjf(rep(2, 10)), .agqObjf(rep(2, 7)))
  })

  test_that("ui$scaleCtheta has one value per estimated theta", {
    .ui <- rxode2::rxUiDecompress(rxode2::rxode2(.mod))
    # it was followed by an NA for each eta row of iniDf
    expect_identical(.ui$scaleCtheta, c(1, 1, 1, 0.35))
    # and a longer foceiControl(scaleC=) was returned whole
    assign("control", foceiControl(scaleC = rep(2, 10)), envir = .ui)
    expect_warning(.sc <- .ui$scaleCtheta, "more options than estimated")
    expect_identical(.sc, rep(2, 4))
  })

  # With literalFix = TRUE a fixed tka leaves the model, so the optimizer moves
  # the same parameters as with literalFix = FALSE, where tka keeps its row of
  # $scaleInfo (and of the parameter vectors behind it)
  .scaleInfoFixed <- function(literalFix, ...) {
    .m <- .mod |> rxode2::ini(tka = fix(0.45))
    .ctlS <- foceiControl(print = 0, calcTables = FALSE, literalFix = literalFix, ...)
    suppressMessages(suppressWarnings(nlmixr(.m, theo_sd, "focei", control = .ctlS)))$scaleInfo
  }

  test_that("$scaleInfo reports each parameter's own initial gradient search", {
    skip_on_cran()
    .cols <- c("Initial Gradient", "Forward aEps", "Forward rEps", "Central aEps", "Central rEps")
    .a <- .scaleInfoFixed(FALSE, maxOuterIterations = 1L, covMethod = "", outerOpt = "nlminb")
    .b <- .scaleInfoFixed(TRUE, maxOuterIterations = 1L, covMethod = "", outerOpt = "nlminb")
    # the fixed tka has no search; every other row is its own parameter's
    expect_equal(as.character(.a[["Initial Gradient"]][1]), "Not Assessed")
    expect_true(all(is.na(unlist(.a[1, .cols[-1]]))))
    expect_true(all(as.character(.b[["Initial Gradient"]]) != "Not Assessed"))
    expect_equal(as.character(.a[["Initial Gradient"]][-1]), as.character(.b[["Initial Gradient"]]))
    expect_equal(.a[-1, .cols[-1]], .b[, .cols[-1]], ignore_attr = TRUE, tolerance = 1e-10)
  })

  test_that("$scaleInfo reports each parameter's own covariance step search", {
    skip_on_cran()
    .cols <- c("Covariance Gradient", "Covariance aEps", "Covariance rEps")
    # the covariance step at the initial estimates differences the thetas only
    .a <- .scaleInfoFixed(FALSE, maxOuterIterations = 0L)
    .b <- .scaleInfoFixed(TRUE, maxOuterIterations = 0L)
    # the fixed tka has no search; every other theta's row is its own
    expect_equal(as.character(.a[["Covariance Gradient"]][1]), "Not Assessed")
    expect_true(all(is.na(unlist(.a[1, .cols[-1]]))))
    expect_true(all(as.character(.b[["Covariance Gradient"]][1:3]) != "Not Assessed"))
    expect_equal(as.character(.a[["Covariance Gradient"]][-1]), as.character(.b[["Covariance Gradient"]]))
    expect_equal(.a[-1, .cols[-1]], .b[, .cols[-1]], ignore_attr = TRUE, tolerance = 1e-10)
  })

  test_that("the first omega parameter is scaled by its diagXform", {
    skip_on_cran()
    .om <- c("o1", "o2", "o3")
    # "log": each diagonal of chol(omega^-1) is exp(x), scaled by 1/2
    .f <- suppressMessages(suppressWarnings(nlmixr(.mod, theo_sd, "focei", control = .ctl(diagXform = "log"))))
    expect_equal(.usedScaleC(.f, .om), c(o1 = 0.5, o2 = 0.5, o3 = 0.5), tolerance = 1e-6)
    expect_equal(.f$scaleInfo$scaleC[5:7], rep(0.5, 3))
    # "identity": the diagonal 1/sqrt(omega) itself, scaled by 1/(2|init|)
    .f <- suppressMessages(suppressWarnings(nlmixr(.mod, theo_sd, "focei", control = .ctl(diagXform = "identity"))))
    .want <- setNames(sqrt(c(0.6, 0.3, 0.1)) / 2, .om)
    expect_equal(.usedScaleC(.f, .om), .want, tolerance = 1e-6)
    expect_equal(.f$scaleInfo$scaleC[5:7], unname(.want), tolerance = 1e-12)
  })

  test_that("a zero gradient away from the scale's anchor keeps its scale", {
    skip_on_cran()
    # An outer optimizer whose first gradient is not at the starting values:
    # beta has moved by 1 on the optimizer's scale.  A new scaleC there would
    # move beta under the optimizer (and the search would difference about the
    # objective at the old point), so the same point must still be beta = 2.5.
    .opt <- function(par, fn, gr, lower, upper, control, ...) {
      .p <- par
      .p[4] <- .p[4] + 1
      fn(.p)
      gr(.p)
      .v <- fn(.p)
      list(x = .p, value = .v, convergence = 0L, message = "")
    }
    .f <- suppressMessages(suppressWarnings(nlmixr(
      .modZ,
      .dZ,
      "focei",
      control = foceiControl(
        print = 0,
        maxOuterIterations = 1L,
        covMethod = "",
        calcTables = FALSE,
        outerOpt = .opt
      )
    )))
    .ph <- .f$parHistData
    .beta <- .ph$beta[.ph$type == "Unscaled"]
    expect_length(.beta, 2L)
    # beta = 0.5 + 1 * its scaleC, 1/0.5
    expect_equal(.beta, c(2.5, 2.5))
  })
})
