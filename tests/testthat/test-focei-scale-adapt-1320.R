nmTest({
  # -2LL of independent normals: the marginal SE of x[i] is se[i]
  .se <- c(0.05, 2, 30)
  .mu <- c(1, -3, 20)
  .fn <- function(x) sum(((x - .mu) / .se)^2)

  test_that("scaleAdapt validates and round-trips", {
    expect_equal(foceiControl()$scaleAdapt, 0L)
    expect_equal(foceiControl(scaleAdapt = 3)$scaleAdapt, 3L)
    expect_error(foceiControl(scaleAdapt = -1))
    expect_error(foceiControl(scaleAdapt = 1.5))
  })

  test_that(".bobyqaCurvature gives each coordinate's marginal SE", {
    .x <- c(0, 0, 0)
    .c <- .bobyqaCurvature(.fn, .x, .fn(.x), 0.1, rep(-Inf, 3), rep(Inf, 3))
    expect_equal(.c$scale, .se, tolerance = 1e-8)
    expect_equal(.c$nEval, 6L)
    # a quadratic's diagonal Newton decrement is the whole decrease here
    expect_equal(.c$dec, .fn(.x), tolerance = 1e-8)
    # one-sided at a bound, no curvature on a flat coordinate
    .flat <- function(x) (x[1] - 1)^2
    .c <- .bobyqaCurvature(.flat, c(0, 0), 1, 0.1, c(-0.05, -Inf), c(Inf, Inf))
    expect_equal(.c$scale, c(sqrt(2 / 2), NA))
  })

  test_that(".bobyqaAdapt restarts from a stopped search to the minimum", {
    .ctl <- list(rhobeg = 0.2, rhoend = 1e-6, npt = 7L, maxfun = 2000L, iprint = 0L)
    .x <- c(0, 0, 0)
    .stopped <- list(par = .x, fval = .fn(.x), feval = 1L, ierr = 0L, msg = "stopped")
    .r <- .bobyqaAdapt(.fn, -Inf, Inf, .ctl, .stopped, tol = 0.01, maxLeg = 3L)
    expect_gte(.r$nScaleAdapt, 1L)
    expect_equal(.r$par, .mu, tolerance = 1e-4)
    expect_lt(.r$fval, 1e-6)
    # at the minimum the decrement stops it before any restart
    .done <- list(par = .mu, fval = 0, feval = 1L, ierr = 0L, msg = "done")
    .r <- .bobyqaAdapt(.fn, -Inf, Inf, .ctl, .done, tol = 0.01, maxLeg = 3L)
    expect_equal(.r$nScaleAdapt, 0L)
    expect_equal(.r$par, .mu)
    # maxLeg = 0 never restarts
    expect_equal(.bobyqaAdapt(.fn, -Inf, Inf, .ctl, .stopped, tol = 0.01, maxLeg = 0L)$nScaleAdapt, 0L)
  })

  test_that("scaleAdapt takes the #1320 Emax fit past the default scale's stop", {
    skip_on_cran()
    rxode2::setRxThreads(1L)
    on.exit(rxode2::setRxThreads(NULL), add = TRUE)
    set.seed(20261108)
    .arms <- expand.grid(ARM = 1:3, ID = 1:8)
    .arms$DOSE <- c(0, 50, 200)[.arms$ARM]
    .arms$N <- sample(40:200, nrow(.arms), replace = TRUE)
    .arms$OCC <- .arms$ARM
    .d <- merge(.arms, data.frame(TIME = c(0, 2, 4, 8, 12)))
    .d <- .d[order(.d$ID, .d$OCC, .d$TIME), ]
    .e0 <- 10 + stats::rnorm(8, 0, 1)[.d$ID]
    .emaxArm <- 8 * exp(stats::rnorm(nrow(.arms), 0, 0.1))[match(paste(.d$ID, .d$OCC), paste(.arms$ID, .arms$OCC))]
    .mean <- .e0 - .emaxArm * .d$DOSE / (40 + .d$DOSE) * .d$TIME / (3 + .d$TIME)
    .d$W <- .d$N / 9
    .d$DV <- (.mean + stats::rnorm(nrow(.d), 0, 1 / sqrt(.d$W))) * sqrt(.d$W)
    .mod <- function() {
      ini({
        tve0 <- 9
        tvemax <- 6
        led50 <- log(30)
        lt50 <- log(2)
        eta.e0 ~ 1
        iov.emax ~ 0.05 | OCC
        sdw <- 1
      })
      model({
        E0 <- tve0 + eta.e0
        EMAX <- tvemax * exp(iov.emax)
        ONSET <- 0
        if (TIME > 0) ONSET <- TIME / (exp(lt50) + TIME)
        EFF <- EMAX * DOSE / (exp(led50) + DOSE) * ONSET
        pred <- (E0 - EFF) * sqrt(W)
        pred ~ add(sdw)
      })
    }
    .f <- suppressMessages(suppressWarnings(nlmixr2(
      .mod,
      .d,
      est = "focei",
      control = foceiControl(print = 0, covMethod = "", calcTables = FALSE, scaleAdapt = 3L)
    )))
    # the default scale's first search is unchanged; a rescaled restart ran
    expect_equal(.f$scaleInfo$scaleC[1:2], c(1 / 9, 1 / 6))
    expect_gte(.f$env$optReturn$nScaleAdapt, 1L)
    # and reached the optimum FeRx and |init| scaling agree on (161.24)
    expect_lt(.f$objf, 165)
    expect_gt(.f$theta[["tvemax"]], 7.5)
  })
})
