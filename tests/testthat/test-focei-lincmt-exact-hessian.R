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
    do.call(rbind, lapply(1:6, function(i) {
      d <- data.frame(ID = i, TIME = c(0, 0.5, 1, 2, 4, 8, 12, 24),
                      EVID = c(1, rep(0, 7)), AMT = c(100, rep(0, 7)), CMT = c(1, rep(2, 7)))
      k <- 2 * exp(rnorm(1, 0, 0.3)) / (20 * exp(rnorm(1, 0, 0.2)))
      v <- 20 * exp(rnorm(1, 0, 0.2))
      cp <- 100 / v * 1.2 / (1.2 - k) * (exp(-k * d$TIME) - exp(-1.2 * d$TIME))
      d$DV <- ifelse(d$EVID == 0, cp * (1 + 0.15 * rnorm(8)), 0)
      d
    }))
  }

  .llLin <- function() {
    ini({
      tka <- log(1.2); tcl <- log(2); tv <- log(20)
      eta.cl ~ 0.09; eta.v ~ 0.04
      lsd <- log(0.15)
    })
    model({
      ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v); sd <- exp(lsd)
      cp <- linCmt()
      ll(err) ~ -0.5 * log(2 * pi) - log(sd * cp) - 0.5 * ((DV - cp) / (sd * cp))^2
    })
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

  test_that("ll() + linCmt(): exact inner Hessian and fast fit match the ODE twin (#1498)", {
    skip_on_cran()
    skip_if_not(.foceiLinCmtHasSecondOrder())
    d <- .exactHessData()
    .ctl <- function(cores, maxOuter = NULL) {
      args <- list(fast = TRUE, linCmtSensCarry = "none", print = 0L, covMethod = "",
                   rxControl = rxode2::rxControl(cores = cores))
      if (!is.null(maxOuter)) args$maxOuterIterations <- maxOuter
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
    fO0 <- suppressMessages(suppressWarnings(nlmixr2(.llOde, d, "focei", .ctl(1L, 0L))))
    fO <- suppressMessages(suppressWarnings(nlmixr2(.llOde, d, "focei", .ctl(1L))))
    # objective and EBEs at the initial estimates
    expect_true(fL0$control$fast)
    expect_equal(fL0$objf, fO0$objf, tolerance = 1e-5)
    expect_equal(as.matrix(fL0$eta[, -1]), as.matrix(fO0$eta[, -1]), tolerance = 1e-4)
    # the converged fit
    expect_equal(fL$objf, fO$objf, tolerance = 1e-3)
    expect_equal(unname(fL$theta), unname(fO$theta), tolerance = 1e-2)
  })
})
