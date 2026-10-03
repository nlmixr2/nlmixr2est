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

  test_that("the nlm family finite-differences the gradient of a lagged calculated variable (issue 1140)", {
    skip_on_cran()
    # c0 is a bare symbol to symengine, so the analytic gradient of every
    # structural theta was exactly 0 (against 9.6, -190, -215)
    .x <- suppressMessages(nlmObjectiveSetup(
      .lagMod,
      .lagDat,
      control = nlmControl(print = 0L),
      gradient = TRUE,
      scale = "natural"
    ))
    .g <- nlmLikEvalC_(.x)$grad
    .fd <- vapply(
      seq_along(.x),
      function(i) {
        .e <- replace(numeric(length(.x)), i, 1e-5)
        (nlmSolveR(.x + .e) - nlmSolveR(.x - .e)) / 2e-5
      },
      numeric(1)
    )
    .nlmFreeEnv()
    expect_equal(.g, .fd, tolerance = 1e-2)
    # so the gradient methods stayed at the initial structural estimates
    # (0.45, 1, 3.45); least squares (nls) has the same optimum
    .nls <- .nlmixr(.lagMod, .lagDat, est = "nls", control = nlsControl(print = 0L, solveType = "fun"))
    .nlm <- .nlmixr(.lagMod, .lagDat, est = "nlm", control = nlmControl(print = 0L))
    expect_equal(unname(.nlm$theta[1:3]), unname(.nls$theta[1:3]), tolerance = 1e-3)
    .n1qn1 <- .nlmixr(.lagMod, .lagDat, est = "n1qn1", control = n1qn1Control(print = 0L))
    expect_equal(unname(.n1qn1$theta[1:3]), unname(.nls$theta[1:3]), tolerance = 2e-2)
    # nls with its gradient stopped with "none of the predictions depend on 'THETA'"
    .nlsGrad <- .nlmixr(.lagMod, .lagDat, est = "nls", control = nlsControl(print = 0L))
    expect_equal(.nlsGrad$theta, .nls$theta, tolerance = 1e-4)
  })

  test_that("the fit table predicts with lag() of a calculated variable (issue 1140)", {
    skip_on_cran()
    # The tables took the prediction to be the column after time, which is c0
    # here: IPRED was c0, and the residual variance the column after it.
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
