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
