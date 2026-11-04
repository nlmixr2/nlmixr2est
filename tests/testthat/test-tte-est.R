nmTest({
  skip_if_not(.nlmixr2estRxode2HasTte(), "rxode2 without tte() endpoints")

  # Monolix-style rows of a constant-hazard endpoint with a treatment effect
  .tteGlmData <- function() {
    withr::with_seed(1099, {
      .n <- 60L
      .trt <- rep(0:1, each = .n)
      .h <- exp(-2 + -0.4 * .trt)
      .time <- pmin(stats::rexp(length(.trt), rate = .h), 5)
      .event <- as.integer(.time < 5)
    })
    .id <- seq_along(.trt)
    .d <- rbind(
      data.frame(id = .id, time = 0, dv = 0, trt = .trt),
      data.frame(id = .id, time = .time, dv = .event, trt = .trt)
    )
    list(
      data = .d[order(.d$id, .d$time), ],
      glm = data.frame(time = .time, event = .event, trt = .trt)
    )
  }

  test_that("tte() endpoint fit matches the Poisson GLM", {
    .d <- .tteGlmData()
    .mod <- function() {
      ini({
        log_h0 <- log(0.2)
        log_hr <- -0.1
      })
      model({
        haz <- exp(log_h0 + log_hr * trt)
        ev ~ tte(haz)
      })
    }
    .fit <- suppressMessages(.nlmixr(.mod, .d$data, est = "bobyqa", bobyqaControl(print = 0)))
    .glm <- glm(event ~ trt, offset = log(time), family = poisson(), data = .d$glm)
    expect_equal(unname(.fit$theta[c("log_h0", "log_hr")]), unname(coef(.glm)), tolerance = 0.01)
    expect_equal(
      unname(.fit$parFixedDf[c("log_h0", "log_hr"), "SE"]),
      unname(summary(.glm)$coefficients[, "Std. Error"]),
      tolerance = 0.1
    )
  })

  test_that("tte() mixed-effects fits match the exact marginal likelihood", {
    skip_on_cran()
    .mod <- function() {
      ini({
        llam <- log(0.05)
        g <- 0.05
        eta.lam ~ 0.3
      })
      model({
        lam <- exp(llam + eta.lam)
        haz <- lam * exp(g * time)
        ev ~ tte(haz)
      })
    }
    .s <- withr::with_seed(11, rxode2::rxSolve(.mod(), rxode2::et(c(0, 30)) |> rxode2::et(id = 1:300)))
    .d <- as.data.frame(.s)[, c("id", "time", "sim")]
    names(.d)[3] <- "dv"
    # exact marginal maximum likelihood (Gauss-Hermite over the eta)
    .ev <- split(.d$time[.d$dv == 1], factor(.d$id[.d$dv == 1], levels = unique(.d$id)))
    .end <- tapply(.d$time, .d$id, max)
    .J <- matrix(0, 40, 40)
    for (.i in 1:39) .J[.i, .i + 1] <- .J[.i + 1, .i] <- sqrt(.i)
    .e <- eigen(.J, symmetric = TRUE)
    .nll <- function(p) {
      -sum(vapply(seq_along(.end), function(i) {
        .t <- .ev[[names(.end)[i]]]
        .lam <- exp(p[1] + .e$values * exp(p[3] / 2))
        .l <- vapply(.lam, function(L) sum(log(L) + p[2] * .t) - L / p[2] * (exp(p[2] * .end[i]) - 1), numeric(1))
        max(.l) + log(sum(.e$vectors[1, ]^2 * exp(.l - max(.l))))
      }, numeric(1)))
    }
    .o <- stats::optim(c(log(0.05), 0.05, log(0.3)), .nll, method = "BFGS")$par
    .fit <- suppressMessages(.nlmixr(.mod, .d, est = "focei", foceiControl(print = 0)))
    expect_equal(unname(.fit$theta[c("llam", "g")]), .o[1:2], tolerance = 0.02)
    expect_equal(unname(.fit$omega[1, 1]), exp(.o[3]), tolerance = 0.05)
    .fit <- suppressMessages(.nlmixr(.mod, .d, est = "saem", saemControl(print = 0)))
    expect_equal(unname(.fit$theta[c("llam", "g")]), .o[1:2], tolerance = 0.05)
    expect_equal(unname(.fit$omega[1, 1]), exp(.o[3]), tolerance = 0.1)
  })

  test_that("tte() endpoint data are checked", {
    .mod <- function() {
      ini({
        lam <- 0.1
      })
      model({
        haz <- lam
        ev ~ tte(haz)
      })
    }
    .d <- data.frame(id = 1, time = c(0, 2, 5), dv = c(0, 2, 0))
    expect_error(.nlmixr(.mod, .d, est = "bobyqa"), "must be 0 or 1")
    .d <- data.frame(id = rep(1:2, each = 2), time = c(1, 5, 0, 5), dv = c(1, 0, 0, 0))
    # the fit keeps its warnings in $runInfo
    .fit <- suppressMessages(.nlmixr(.mod, .d, est = "bobyqa", bobyqaControl(print = 0)))
    expect_true(any(grepl("start of observation", .fit$runInfo, fixed = TRUE)))
  })

  test_that("fo cannot estimate a tte() endpoint", {
    .mod <- function() {
      ini({
        llam <- log(0.1)
        eta.lam ~ 0.1
      })
      model({
        haz <- exp(llam + eta.lam)
        ev ~ tte(haz)
      })
    }
    .d <- data.frame(id = rep(1:2, each = 2), time = c(0, 5, 0, 5), dv = c(0, 1, 0, 0))
    expect_error(suppressMessages(.nlmixr(.mod, .d, est = "fo")))
  })
})
