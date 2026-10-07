nmTest({
  # lag() of a calculated variable (c0): the variable is output ahead of the
  # prediction, and has no symbolic sensitivity

  .lagMod <- function() {
    ini({
      tka <- 0.45
      tcl <- 1.0
      tv <- 3.45
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl)
      v <- exp(tv)
      d/dt(depot) <- -ka * depot
      d/dt(central) <- ka * depot - cl / v * central
      c0 <- central / v
      cp <- 0.5 * c0 + 0.5 * lag(c0)
      cp ~ add(add.sd)
    })
  }

  .lagDat <- local({
    .m <- rxode2::rxode2(
      "ka=exp(tka)\ncl=exp(tcl)\nv=exp(tv)\nd/dt(depot)=-ka*depot\nd/dt(central)=ka*depot-cl/v*central\nc0=central/v\ncp=0.5*c0+0.5*lag(c0)\n"
    )
    .ev <- rxode2::et(amt = 320, cmt = "depot") |> rxode2::et(seq(0.5, 24, by = 1.5))
    .s <- rxode2::rxSolve(.m, .ev, params = c(tka = 0.6, tcl = 1.1, tv = 3.6), returnType = "data.frame")
    rxode2::rxWithSeed(1234, {
      .d <- rbind(
        data.frame(ID = 1:4, TIME = 0, DV = NA, AMT = 320, EVID = 1),
        data.frame(
          ID = rep(1:4, each = nrow(.s)),
          TIME = rep(.s$time, 4),
          DV = rep(.s$cp, 4) + stats::rnorm(4 * nrow(.s), 0, 0.3),
          AMT = 0,
          EVID = 0
        )
      )
      .d[order(.d$ID, .d$TIME, -.d$EVID), ]
    })
  })

  .lagEtaMod <- function() {
    ini({
      tka <- 0.45
      tcl <- 1.0
      tv <- 3.45
      add.sd <- 0.7
      eta.f ~ 0.1
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl)
      v <- exp(tv)
      d/dt(depot) <- -ka * depot
      d/dt(central) <- ka * depot - cl / v * central
      c0 <- central / v
      cp <- (0.5 * c0 + 0.5 * lag(c0)) * exp(eta.f)
      cp ~ add(add.sd)
    })
  }

  # the analytic nlm gradient of a model at its initial estimates, and its
  # central difference
  .lagGradFd <- function(mod) {
    .x <- suppressMessages(nlmObjectiveSetup(
      mod,
      .lagDat,
      control = nlmControl(print = 0L),
      gradient = TRUE,
      scale = "natural"
    ))
    on.exit(.nlmFreeEnv())
    list(
      grad = nlmLikEvalC_(.x)$grad,
      fd = vapply(
        seq_along(.x),
        function(i) {
          .e <- replace(numeric(length(.x)), i, 1e-5)
          (nlmSolveR(.x + .e) - nlmSolveR(.x - .e)) / 2e-5
        },
        numeric(1)
      )
    )
  }

  test_that("the nlm family carries theta sensitivities through a lagged calculated variable (issue 1140)", {
    skip_on_cran()
    # c0 is a bare symbol to symengine, so its sensitivity is its own lhs and
    # d(lag(c0))/d(theta) = lag(d(c0)/d(theta)); no theta is finite-differenced
    .s <- suppressMessages(rxode2::rxode2(.lagMod)$nlmEnv)
    expect_equal(.s$.eventTheta, rep(0L, 4))
    expect_true(any(grepl("lag(rx_lsens_1_1_)", .s$..nlmS, fixed = TRUE)))
    .r <- .lagGradFd(.lagMod)
    expect_equal(.r$grad, .r$fd, tolerance = 5e-3)
    # the gradient methods reach the optimum least squares (nls) finds
    .nls <- .nlmixr(.lagMod, .lagDat, est = "nls", control = nlsControl(print = 0L, solveType = "fun"))
    .nlm <- .nlmixr(.lagMod, .lagDat, est = "nlm", control = nlmControl(print = 0L))
    expect_equal(unname(.nlm$theta[1:3]), unname(.nls$theta[1:3]), tolerance = 1e-3)
    .n1qn1 <- .nlmixr(.lagMod, .lagDat, est = "n1qn1", control = n1qn1Control(print = 0L))
    expect_equal(unname(.n1qn1$theta[1:3]), unname(.nls$theta[1:3]), tolerance = 2e-2)
    # nls fits it with its own gradient too
    .nlsGrad <- .nlmixr(.lagMod, .lagDat, est = "nls", control = nlsControl(print = 0L))
    expect_equal(.nlsGrad$theta, .nls$theta, tolerance = 1e-4)
  })

  test_that("lagged sensitivities chain through diff() and into the ODEs (issue 1140)", {
    skip_on_cran()
    # c1 lags c0, and diff(c1) needs the sensitivity of c1
    .diffMod <- .lagMod |>
      rxode2::model(c1 <- 2 * c0 + lag(c0), append = c0) |>
      rxode2::model(cp <- 0.5 * c0 + 0.5 * lag(c0) + 0.1 * diff(c1))
    expect_equal(suppressMessages(rxode2::rxode2(.diffMod)$nlmEnv$.eventTheta), rep(0L, 4))
    .r <- .lagGradFd(.diffMod)
    expect_equal(.r$grad, .r$fd, tolerance = 5e-3)
    # an ODE that uses c0 gets c0's definition for its sensitivities
    .odeMod <- .lagMod |>
      rxode2::model(d / dt(eff) <- c0 - eff, append = c0) |>
      rxode2::model(cp <- eff + 0.5 * lag(c0))
    expect_equal(suppressMessages(rxode2::rxode2(.odeMod)$nlmEnv$.eventTheta), rep(0L, 4))
    .r <- .lagGradFd(.odeMod)
    expect_equal(.r$grad, .r$fd, tolerance = 5e-3)
    # an ODE with lag(c0) has no sensitivity ODE: the thetas are finite-differenced
    .lagOdeMod <- .lagMod |>
      rxode2::model(d / dt(eff) <- lag(c0) - eff, append = c0) |>
      rxode2::model(cp <- eff + 0.5 * c0)
    expect_equal(suppressMessages(rxode2::rxode2(.lagOdeMod)$nlmEnv$.eventTheta), rep(1L, 4))
  })

  test_that("lag(v, 1) and a redefined lagged variable (issue 1140)", {
    skip_on_cran()
    # symengine keeps lag(c0, 1) as lag(c0, 1.0); it is still followed analytically
    .lag1Mod <- .lagMod |>
      rxode2::model(cp <- 0.5 * c0 + 0.5 * lag(c0, 1))
    expect_equal(suppressMessages(rxode2::rxode2(.lag1Mod)$nlmEnv$.eventTheta), rep(0L, 4))
    .r <- .lagGradFd(.lag1Mod)
    expect_equal(.r$grad, .r$fd, tolerance = 5e-3)
    # c1 reads the first c0, but symengine inlines it with the last one, so a
    # lagged variable defined twice is finite-differenced
    .redefMod <- .lagMod |>
      rxode2::model(c1 <- 2 * c0 + lag(c0), append = c0) |>
      rxode2::model(c0 <- c0 * exp(tka), append = c1) |>
      rxode2::model(cp <- 0.3 * c1 + 0.5 * lag(c0) + c0)
    expect_equal(suppressMessages(rxode2::rxode2(.redefMod)$nlmEnv$.eventTheta), rep(1L, 4))
  })

  test_that("the lagged definitions are matched by name (issue 1140)", {
    .s <- new.env(parent = emptyenv())
    .s$..laggedVars <- "c.0"
    .s$..lhs <- c("cx0=1", "c.0=2*central")
    expect_identical(.nlmFamilyLagDefs(.s), "c.0=2*central")
    .s <- suppressMessages(rxode2::rxode2(.lagMod)$nlmEnv)
    expect_identical(.nlmFamilyLagDefs(.s), "c0=exp(-THETA[3])*central")
  })

  test_that("the predictions of a lagged model depend on the thetas it uses (issue 1140)", {
    # used when a lagged variable is finite-differenced: the build counts the
    # thetas the predictions, ODEs and calculated variables use: THETA[1-3]
    # enter the ODEs and c0, THETA[4] the error model
    .s <- suppressMessages(rxode2::rxode2(.lagMod)$nlmEnv)
    expect_identical(.nlmFamilyThetaUsed(.s), rep(TRUE, 4))
    .s$..maxTheta <- 5L
    expect_identical(.nlmFamilyThetaUsed(.s), c(rep(TRUE, 4), FALSE))
  })

  test_that("the table step stops when the solve has no prediction column (issue 1140)", {
    # columns are found by name, never guessed by position
    .df <- list(ID = 1L, time = 0, c0 = 1, rx_r_ = 1, rxLambda = 1, rxYj = 2, rxLow = 0, rxHi = 1)
    .ires <- function(df) {
      .Call(`_nlmixr2est_iresCalc`, df, 1, 0L, NULL, NULL, character(0), character(0), character(0), NULL, list())
    }
    expect_error(.ires(.df), "'rx_pred_' not found in the solved data.frame")
    expect_error(.ires(.df[names(.df) != "time"]), "'time' not found in the solved data.frame")
    .df$rx_pred_ <- 1
    expect_error(.ires(.df[names(.df) != "rx_r_"]), "'rx_r_' not found in the solved data.frame")
  })

  test_that("the fit table predicts with lag() of a calculated variable (issue 1140)", {
    skip_on_cran()
    # c0 is output ahead of the prediction: the table finds rx_pred_ and rx_r_
    # by name
    .fit <- .nlmixr(.lagMod, .lagDat, est = "nlm", control = nlmControl(print = 0L, solveType = "fun"))
    expect_equal(.fit$IPRED, .fit$cp, tolerance = 1e-12)
    expect_equal(.fit$IWRES, (.fit$DV - .fit$cp) / .fit$theta[["add.sd"]], tolerance = 1e-12)
    # a mixed-effects fit (the CWRES table); the record before the first
    # observation is the dose, where c0 is 0
    .fit <- .nlmixr(
      .lagEtaMod,
      .lagDat,
      est = "focei",
      control = foceiControl(print = 0L, maxOuterIterations = 0L, covMethod = "")
    )
    .d <- as.data.frame(.fit)
    .prev <- ave(.d$c0, .d$ID, FUN = function(x) c(0, x[-length(x)]))
    .ipred <- (0.5 * .d$c0 + 0.5 * .prev) * exp(.d$eta.f)
    expect_equal(.d$IPRED, .ipred, tolerance = 1e-12)
    expect_equal(.d$IWRES, (.d$DV - .ipred) / .fit$theta[["add.sd"]], tolerance = 1e-12)
    # c0 does not depend on eta.f
    expect_equal(.d$PRED, 0.5 * .d$c0 + 0.5 * .prev, tolerance = 1e-12)
  })
})
