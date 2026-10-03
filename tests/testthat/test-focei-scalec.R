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
    # the band misses add.sd's 0.35 and two of the omegas' values: only the
    # theta is guarded (to |init|), the omegas keep theirs
    .f <- suppressMessages(suppressWarnings(nlmixr(.mod, theo_sd, "mfocei", control = .ctl(scaleCband = c(0.8, 10)))))
    .want <- c(add.sd = 0.7, .omegaC)
    expect_equal(.usedScaleC(.f, names(.want)), .want, tolerance = 1e-6)
    expect_equal(.f$scaleInfo$scaleC, c(NA, NA, NA, unname(.want)), tolerance = 1e-6)
  })

  test_that("a fixed theta does not shift the scaleC of the parameters after it", {
    skip_on_cran()
    .m <- .mod |> rxode2::ini(tka = fix(0.45))
    .f <- suppressMessages(suppressWarnings(nlmixr(.m, theo_sd, "focei", control = .ctl(literalFix = FALSE))))
    .want <- c(tcl = 1, tv = 1, add.sd = 0.35, .omegaC)
    expect_equal(.usedScaleC(.f, names(.want)), .want, tolerance = 1e-6)
  })

  test_that("a scaleC longer than the parameter vector is only read up to it", {
    skip_on_cran()
    # the extra values were copied past the end of the buffer (a crash)
    .f <- suppressMessages(suppressWarnings(nlmixr(.mod, theo_sd, "focei", control = .ctl(scaleC = rep(2, 2000)))))
    expect_equal(.f$scaleInfo$scaleC, rep(2, 7))
  })
})
