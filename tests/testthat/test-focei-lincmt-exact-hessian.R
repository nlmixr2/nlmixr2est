nmTest({
  # #1498: with second-order linCmt() reads (nlmixr2/rxode2lincmt#17) a
  # linCmt() model gets the exact inner Hessian and the analytic 'fast'
  # gradient instead of the Shi21 finite difference.  Ground truth is the same
  # model written as ODEs.  Several subjects on two cores: the inner model and
  # the second-order model are solved concurrently for different subjects,
  # which is what exposed the shared rx->ndiff (odeSwapSolveInd).  (An
  # observation at the dose time is checked in test-focei-lincmt-detect-1103.R:
  # this oral model's log-density is undefined there.)

  .exactHessData <- function() {
    set.seed(1498)
    do.call(
      rbind,
      lapply(1:6, function(i) {
        d <- data.frame(
          ID = i,
          TIME = c(0, 0.5, 1, 2, 4, 8, 12, 24),
          EVID = c(1, rep(0, 7)),
          AMT = c(100, rep(0, 7)),
          CMT = c(1, rep(2, 7))
        )
        k <- 2 * exp(rnorm(1, 0, 0.3)) / (20 * exp(rnorm(1, 0, 0.2)))
        v <- 20 * exp(rnorm(1, 0, 0.2))
        cp <- 100 / v * 1.2 / (1.2 - k) * (exp(-k * d$TIME) - exp(-1.2 * d$TIME))
        d$DV <- ifelse(d$EVID == 0, cp * (1 + 0.15 * rnorm(8)), 0)
        d
      })
    )
  }

  .lagged <- function(f) {
    # the same model with a lag on the depot: the second-order reads are NA
    # there (nlmixr2/rxode2lincmt#19), so the inner Hessian falls back to the
    # finite difference for those subjects
    rxode2::model(rxode2::rxode2(f), alag(depot) <- 0.25, append = TRUE)
  }

  .llOde <- function() {
    ini({
      tka <- log(1.2); tcl <- log(2); tv <- log(20)
      eta.cl ~ 0.09; eta.v ~ 0.04
      lsd <- log(0.15)
    })
    model({
      ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v); sd <- exp(lsd)
      d/dt(depot) <- -ka * depot
      d/dt(central) <- ka * depot - cl / v * central
      cp <- central / v
      ll(err) ~ -0.5 * log(2 * pi) - log(sd * cp) - 0.5 * ((DV - cp) / (sd * cp))^2
    })
  }

  # the linCmt() twin of an ODE model
  .linOf <- function(f) {
    suppressMessages(rxode2::model(f, -d/dt(depot), -d/dt(central), cp <- linCmt()))
  }

  test_that("ll() + linCmt(): exact inner Hessian and fast fit match the ODE twin (#1498)", {
    skip_on_cran()
    skip_if_not(.foceiLinCmtHasSecondOrder())
    .llLin <- .linOf(.llOde)
    d <- .exactHessData()
    .ctl <- function(cores, maxOuter = NULL) {
      args <- list(
        fast = TRUE,
        linCmtSensCarry = "none",
        print = 0L,
        covMethod = "",
        rxControl = rxode2::rxControl(cores = cores)
      )
      if (!is.null(maxOuter)) {
        args$maxOuterIterations <- maxOuter
      }
      do.call(foceiControl, args)
    }
    ui <- rxode2::.copyUi(suppressMessages(nlmixr2(.llLin)))
    assign("control", .ctl(2L), envir = ui)
    s <- suppressWarnings(suppressMessages(ui$foceEnv))
    expect_false(is.null(s$..HdEta2))
    # The linCmt() fits (two cores) run before the ODE twins (one core): a
    # multi-threaded fit after a fast ll() ODE fit in the same process
    # crashes on main too, unrelated to linCmt().
    fL0 <- suppressMessages(suppressWarnings(nlmixr2(.llLin, d, "focei", .ctl(2L, 0L))))
    fL <- suppressMessages(suppressWarnings(nlmixr2(.llLin, d, "focei", .ctl(2L))))
    fLlag <- suppressMessages(suppressWarnings(nlmixr2(.lagged(.llLin), d, "focei", .ctl(2L, 0L))))
    fO0 <- suppressMessages(suppressWarnings(nlmixr2(.llOde, d, "focei", .ctl(1L, 0L))))
    fO <- suppressMessages(suppressWarnings(nlmixr2(.llOde, d, "focei", .ctl(1L))))
    fOlag <- suppressMessages(suppressWarnings(nlmixr2(.lagged(.llOde), d, "focei", .ctl(1L, 0L))))
    # objective and EBEs at the initial estimates
    expect_true(fL0$control$fast)
    expect_equal(fL0$objf, fO0$objf, tolerance = 1e-5)
    expect_equal(as.matrix(fL0$eta[, -1]), as.matrix(fO0$eta[, -1]), tolerance = 1e-4)
    # the converged fit
    expect_equal(fL$objf, fO$objf, tolerance = 1e-3)
    expect_equal(unname(fL$theta), unname(fO$theta), tolerance = 1e-2)
    # a lag: the finite-difference fallback, reported, at its own accuracy
    expect_true(any(grepl("exact eta Hessian not finite", fLlag$runInfo, fixed = TRUE)))
    expect_lt(abs(fLlag$objf - fOlag$objf), 0.01)
    # the exact Hessian is the objective's Hessian: the objective matches the
    # finite-difference one (fast = FALSE) at its accuracy.  linCmt() and the ODE
    # twin share the Hessian assembly, so their agreement does not show this.
    .fd <- foceiControl(fast = FALSE, print = 0L, covMethod = "", maxOuterIterations = 0L,
                        rxControl = rxode2::rxControl(cores = 1L))
    fOfd0 <- suppressMessages(suppressWarnings(nlmixr2(.llOde, d, "focei", .fd)))
    expect_lt(abs(fO0$objf - fOfd0$objf), 0.01)
    expect_lt(abs(fL0$objf - fOfd0$objf), 0.01)
  })

  test_that("ll() + linCmt() with a time-varying covariate: the eta-space carry's exact Hessian (#1498)", {
    skip_on_cran()
    skip_if_not(.foceiLinCmtHasSecondOrder())
    .tvOde <- suppressMessages(rxode2::model(.llOde, cl <- exp(tcl + eta.cl) * (wt / 70)^0.75))
    .tvLin <- .linOf(.tvOde)
    set.seed(7)
    d <- do.call(
      rbind,
      lapply(1:6, function(i) {
        tt <- c(0, 0.5, 1, 2, 4, 6, 8, 12, 12, 14, 16, 20, 24, 36)
        ev <- c(1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0)
        w0 <- runif(1, 50, 90)
        data.frame(
          ID = i,
          TIME = tt,
          EVID = ev,
          AMT = ifelse(ev == 1, 100, 0),
          CMT = ifelse(ev == 1, 1, 2),
          wt = ifelse(tt < 8, w0, w0 * 1.3),
          DV = ifelse(ev == 0, 3 * exp(-0.1 * tt) * exp(rnorm(length(tt), 0, 0.2)), 0)
        )
      })
    )
    for (interp in c("locf", "nocb")) {
      # the ODE twin needs tight tolerances across the covariate step to be
      # a reference at this precision
      .ctl <- function(cores) {
        foceiControl(
          fast = TRUE,
          print = 0L,
          covMethod = "",
          maxOuterIterations = 0L,
          rxControl = rxode2::rxControl(
            covsInterpolation = interp,
            cores = cores,
            atol = 1e-10,
            rtol = 1e-10,
            atolSens = 1e-10,
            rtolSens = 1e-10
          )
        )
      }
      # one core: a multi-threaded fit after the ODE fits above would hit #1502
      fL <- suppressMessages(suppressWarnings(nlmixr2(.tvLin, d, "focei", .ctl(1L))))
      fO <- suppressMessages(suppressWarnings(nlmixr2(.tvOde, d, "focei", .ctl(1L))))
      expect_false(any(grepl("exact eta Hessian not finite", fL$runInfo, fixed = TRUE)))
      expect_equal(fL$objf, fO$objf, tolerance = 1e-6)
    }
    # an eta and the time-varying covariate on ka: the carry's ka row is 2 * ncmt, not
    # the last of the seven slots (a 1-cmt ka derivative read zero there)
    .kaOde <- suppressMessages(rxode2::ini(
      rxode2::model(.llOde, ka <- exp(tka + eta.ka) * (wt / 70)^0.5, v <- exp(tv)),
      eta.ka ~ 0.09
    ))
    .kaLin <- .linOf(.kaOde)
    fL <- suppressMessages(suppressWarnings(nlmixr2(.kaLin, d, "focei", .ctl(1L))))
    fO <- suppressMessages(suppressWarnings(nlmixr2(.kaOde, d, "focei", .ctl(1L))))
    expect_false(any(grepl("exact eta Hessian not finite", fL$runInfo, fixed = TRUE)))
    expect_equal(fL$objf, fO$objf, tolerance = 1e-6)
  })
})
