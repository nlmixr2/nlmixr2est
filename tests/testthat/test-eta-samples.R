# etaSamples() draws each subject's etas from p(eta | y, theta, Omega, sigma)
# at a fit's estimates (#1410).

test_that(".etaSamplesDf is long by draw then ID", {
  .d <- array(seq_len(2 * 3 * 4), c(2, 3, 4))
  .s <- .etaSamplesDf(.d, factor(c("a", "b")), c("e1", "e2", "e3"))
  expect_equal(names(.s), c("ID", "draw", "e1", "e2", "e3"))
  expect_equal(as.character(.s$ID), rep(c("a", "b"), 4))
  expect_equal(.s$draw, rep(1:4, each = 2))
  expect_equal(unlist(.s[4, c("e1", "e2", "e3")], use.names = FALSE), .d[2, , 2])
})

test_that(".etaSamplesChol factors a singular covariance", {
  .c <- matrix(c(2, 0, 1, 0, 0, 0, 1, 0, 1), 3)
  .r <- .etaSamplesChol(.c)
  expect_equal(crossprod(.r), .c)
  expect_equal(crossprod(.etaSamplesChol(diag(c(4, 1)))), diag(c(4, 1)))
})

# A linear-Gaussian model, y = tv + eta.a + eta.b * t + eps, whose conditional
# eta distribution is normal with a closed form
.etaSamplesLin <- local({
  .cache <- NULL
  function() {
    if (!is.null(.cache)) {
      return(.cache)
    }
    .testSeed(1410)
    .n <- 30L
    .t <- c(0.5, 1, 2)
    .d <- data.frame(ID = rep(seq_len(.n), each = length(.t)), TIME = rep(.t, .n))
    .ea <- stats::rnorm(.n, 0, 1)
    .eb <- stats::rnorm(.n, 0, 0.5)
    .d$DV <- 2 + .ea[.d$ID] + .eb[.d$ID] * .d$TIME + stats::rnorm(nrow(.d), 0, 1.5)
    # nolint start: object_usage_linter.
    .mod <- function() {
      ini({
        tv <- 2
        eta.a ~ 1
        eta.b ~ 0.25
        add.sd <- 1.5
      })
      model({
        v <- tv + eta.a + eta.b * t
        v ~ add(add.sd)
      })
    }
    # nolint end
    .fit <- .nlmixr(.mod, .d, "focei", control = foceiControl(print = 0L, calcTables = FALSE))
    .om <- as.matrix(.fit$omega)
    .s2 <- .fit$theta[["add.sd"]]^2
    .tv <- .fit$theta[["tv"]]
    .x <- cbind(1, .t)
    .prec <- solve(.om) + crossprod(.x) / .s2
    .cov <- solve(.prec)
    .mean <- t(vapply(
      seq_len(.n),
      function(i) {
        as.numeric(.cov %*% crossprod(.x, .d$DV[.d$ID == i] - .tv) / .s2)
      },
      numeric(2)
    ))
    .cache <<- list(fit = .fit, mean = .mean, cov = .cov)
    .cache
  }
})

.etaSamplesMoments <- function(smp) {
  .ids <- unique(smp$ID)
  list(
    mean = t(vapply(.ids, function(i) colMeans(smp[smp$ID == i, c("eta.a", "eta.b")]), numeric(2))),
    cov = lapply(.ids, function(i) stats::cov(smp[smp$ID == i, c("eta.a", "eta.b")]))
  )
}

test_that("laplace draws match the linear-Gaussian conditional distribution", {
  skip_on_cran()
  .l <- .etaSamplesLin()
  .s <- etaSamples(.l$fit, n = 4000L, method = "laplace", seed = 3L)
  expect_equal(attr(.s, "method"), "laplace")
  expect_equal(nrow(.s), 4000L * 30L)
  .m <- .etaSamplesMoments(.s)
  .se <- sqrt(diag(.l$cov) / 4000)
  expect_lt(max(abs(sweep(.m$mean - .l$mean, 2, .se, "/"))), 5)
  for (.c in .m$cov) {
    expect_equal(.c, .l$cov, tolerance = 0.1, ignore_attr = TRUE)
  }
  # the mode the draws are centred on is the EBE
  expect_equal(.l$mean, as.matrix(.l$fit$eta[, c("eta.a", "eta.b")]), tolerance = 1e-3, ignore_attr = TRUE)
  # the same seed gives the same draws, whatever the thread count
  .old <- rxode2::getRxThreads()
  on.exit(rxode2::setRxThreads(.old), add = TRUE)
  rxode2::setRxThreads(1L)
  .a <- etaSamples(.l$fit, n = 50L, method = "laplace", seed = 3L)
  rxode2::setRxThreads(2L)
  expect_identical(etaSamples(.l$fit, n = 50L, method = "laplace", seed = 3L), .a)
  expect_identical(etaSamples(.l$fit, n = 50L, seed = 3L), .a)
  expect_false(identical(etaSamples(.l$fit, n = 50L, method = "laplace", seed = 4L), .a))
})

test_that("imp draws match the linear-Gaussian conditional distribution", {
  skip_on_cran()
  .l <- .etaSamplesLin()
  .s <- etaSamples(.l$fit, n = 1000L, method = "imp", seed = 3L)
  expect_equal(attr(.s, "method"), "imp")
  .m <- .etaSamplesMoments(.s)
  .sd <- sqrt(diag(.l$cov))
  expect_lt(max(abs(sweep(.m$mean - .l$mean, 2, .sd, "/"))), 0.2)
  .v <- t(vapply(.m$cov, diag, numeric(2)))
  expect_lt(max(abs(log(sweep(.v, 2, diag(.l$cov), "/")))), 0.35)
})

test_that("draws restore the spread a shrunken EBE loses", {
  skip_on_cran()
  .testSeed(14101)
  .n <- 40L
  .d <- data.frame(ID = seq_len(.n), TIME = 1)
  .d$DV <- 2 + stats::rnorm(.n, 0, 1) + stats::rnorm(.n, 0, 2)
  # nolint start: object_usage_linter.
  .mod <- function() {
    ini({
      tv <- 2
      eta.v ~ fix(1)
      add.sd <- fix(2)
    })
    model({
      v <- tv + eta.v
      v ~ add(add.sd)
    })
  }
  # nolint end
  .fit <- .nlmixr(.mod, .d, "focei", control = foceiControl(print = 0L, calcTables = FALSE))
  .s <- etaSamples(.fit, n = 500L)
  # the EBEs keep 1 / (1 + 4) of Omega's variance; the draws keep all of it
  expect_lt(stats::var(.fit$eta$eta.v), 0.5)
  expect_equal(stats::var(.s$eta.v), 1, tolerance = 0.15)
  # within 5 Monte-Carlo SEs of the EBE (conditional variance 0.8)
  expect_lt(max(abs(as.numeric(tapply(.s$eta.v, .s$ID, mean)) - .fit$eta$eta.v)), 5 * sqrt(0.8 / 500))
})

test_that("a SAEM fit gets its modes from a posthoc step and simulates per draw", {
  skip_on_cran()
  # nolint start: object_usage_linter.
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
      linCmt() ~ add(add.sd)
    })
  }
  # nolint end
  .fit <- .nlmixr(.mod, nlmixr2data::theo_sd, "saem", control = saemControl(print = 0L, calcTables = FALSE))
  .s <- etaSamples(.fit, n = 20L)
  expect_equal(names(.s), c("ID", "draw", "eta.ka", "eta.cl", "eta.v"))
  expect_equal(nrow(.s), 20L * 12L)
  expect_true(all(is.finite(as.matrix(.s[, -(1:2)]))))
  .sim <- suppressMessages(suppressWarnings(
    rxode2::rxSolve(.fit, params = .s, nStud = 20L, simVariability = FALSE)
  ))
  .cl <- tapply(.sim$cl, .sim$sim.id, function(x) x[1])
  expect_equal(as.numeric(.cl), exp(.fit$theta[["tcl"]] + .s$eta.cl[.s$ID == .s$ID[1]]), tolerance = 1e-6)
})

test_that("etaSamples refuses fits without a conditional eta distribution", {
  skip_on_cran()
  expect_error(etaSamples(1), "nlmixr2 fit")
  # nolint start: object_usage_linter.
  .mod <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl)
      v <- exp(tv)
      linCmt() ~ add(add.sd)
    })
  }
  # nolint end
  .fit <- .nlmixr(.mod, nlmixr2data::theo_sd, "nlm", control = nlmControl(print = 0L, calcTables = FALSE))
  expect_error(etaSamples(.fit), "no random effects")
})
