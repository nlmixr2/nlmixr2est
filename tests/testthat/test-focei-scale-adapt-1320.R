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
    # an asymmetric step near a bound; no curvature on a flat coordinate
    .flat <- function(x) (x[1] - 1)^2
    .c <- .bobyqaCurvature(.flat, c(0, 0), 1, 0.1, c(-0.05, -Inf), c(Inf, Inf))
    expect_equal(.c$scale, c(1, NA))
    # at a bound both probes go to the open side
    .c <- .bobyqaCurvature(.flat, 0, 1, 0.1, 0, Inf)
    expect_equal(.c$scale, 1)
    expect_equal(.c$dec, 1)
    # a Newton step that leaves the box predicts only the decrease inside it
    .wall <- function(x) (x + 5)^2
    .c <- .bobyqaCurvature(.wall, 1e-5, .wall(1e-5), 0.1, 0, Inf)
    expect_lt(.c$dec, 1e-3)
    # no positive curvature: the decrease a probe saw still counts
    .bump <- function(x) -10 * exp(-((x - 50) / 20)^2)
    .c <- .bobyqaCurvature(.bump, 0, .bump(0), 0.1, -Inf, Inf)
    expect_true(is.na(.c$scale))
    expect_gt(.c$dec, 0)
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
    # a stop on a concave flank still restarts
    .flank <- function(x) x[1]^2 - 10 * exp(-((x[2] - 50) / 20)^2)
    .st <- list(par = c(0, 0), fval = .flank(c(0, 0)), feval = 1L, ierr = 0L, msg = "stopped")
    .r <- .bobyqaAdapt(.flank, -Inf, Inf, .ctl, .st, tol = 1e-4, maxLeg = 3L)
    expect_gte(.r$nScaleAdapt, 1L)
    expect_lt(.r$fval, .st$fval)
    # maxLeg = 0 never restarts
    expect_equal(.bobyqaAdapt(.fn, -Inf, Inf, .ctl, .stopped, tol = 0.01, maxLeg = 0L)$nScaleAdapt, 0L)
  })

  test_that(".bobyqaAdapt keeps bobyqa's box and start with finite bounds", {
    .ctl <- list(rhobeg = 0.2, rhoend = 1e-6, npt = 5L, maxfun = 2000L, iprint = 0L)
    .fn2 <- function(x) ((x[1] - 1) / 500)^2 + ((x[2] - 20) / 30)^2
    .st <- list(par = c(0, 0), fval = .fn2(c(0, 0)), feval = 1L, ierr = 0L, msg = "stopped")
    .seen <- NULL
    .f <- function(x) {
      .seen <<- rbind(.seen, x)
      .fn2(x)
    }
    expect_no_warning(.r <- .bobyqaAdapt(.f, c(-5, 0), c(5, Inf), .ctl, .st, tol = 0.01, maxLeg = 3L))
    expect_true(all(.seen[, 1] >= -5 & .seen[, 1] <= 5 & .seen[, 2] >= 0))
    expect_lt(.r$fval, .st$fval)
    expect_equal(.r$par[2], 20, tolerance = 1e-3)
  })

  test_that(".bobyqaAdapt stays inside the evaluation budget", {
    .n <- 0L
    .f <- function(x) {
      .n <<- .n + 1L
      .fn(x)
    }
    .x <- c(0, 0, 0)
    .st <- list(par = .x, fval = .fn(.x), feval = 1L, ierr = 0L, msg = "stopped")
    .ctl <- list(rhobeg = 0.2, rhoend = 1e-6, npt = 7L, maxfun = 60L, iprint = 0L)
    expect_no_warning(.r <- .bobyqaAdapt(.f, -Inf, Inf, .ctl, .st, tol = 0.01, maxLeg = 3L))
    expect_lte(.r$feval, 60L)
    expect_lte(.n, 60L)
    # too little budget left for a probe and a restart: nothing runs
    .n <- 0L
    .ctl$maxfun <- 2L * 3L + 7L + 3L
    .r <- .bobyqaAdapt(.f, -Inf, Inf, .ctl, .st, tol = 0.01, maxLeg = 3L)
    expect_equal(.n, 0L)
    expect_equal(.r$nScaleAdapt, 0L)
  })

  test_that(".bobyqaAdapt keeps its probes out of the history and the inner state", {
    .log <- character(0)
    .rec <- function(on) {
      .log <<- c(.log, paste0("record", on))
      TRUE
    }
    .hold <- function(on) {
      .log <<- c(.log, paste0("hold", on))
      TRUE
    }
    .done <- list(par = .mu, fval = 0, feval = 1L, ierr = 0L, msg = "done")
    .ctl <- list(rhobeg = 0.2, rhoend = 1e-6, npt = 7L, maxfun = 2000L, iprint = 0L)
    .bobyqaAdapt(.fn, -Inf, Inf, .ctl, .done, tol = 0.01, maxLeg = 3L, record = .rec, hold = .hold)
    expect_equal(.log, c("holdTRUE", "recordFALSE", "holdFALSE", "recordTRUE", "holdFALSE"))
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
