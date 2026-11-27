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
  expect_false(identical(etaSamples(.l$fit, n = 50L, method = "laplace", seed = 4L), .a))
})

test_that("mcmc draws match the linear-Gaussian conditional distribution", {
  skip_on_cran()
  .l <- .etaSamplesLin()
  .s <- etaSamples(.l$fit, n = 3000L, method = "mcmc", seed = 3L, thin = 2L)
  expect_equal(attr(.s, "method"), "mcmc")
  .acc <- attr(.s, "acceptance")
  expect_length(.acc, 30L)
  expect_true(all(.acc > 0.15 & .acc < 0.6))
  .m <- .etaSamplesMoments(.s)
  .sd <- sqrt(diag(.l$cov))
  # correlated draws: a looser bound than the independent laplace draws
  expect_lt(max(abs(sweep(.m$mean - .l$mean, 2, .sd, "/"))), 0.25)
  .v <- t(vapply(.m$cov, diag, numeric(2)))
  expect_lt(max(abs(log(sweep(.v, 2, diag(.l$cov), "/")))), 0.35)
  expect_lt(abs(mean(log(sweep(.v, 2, diag(.l$cov), "/")))), 0.05)
  # the same seed gives the same chain, whatever the thread count
  .old <- rxode2::getRxThreads()
  on.exit(rxode2::setRxThreads(.old), add = TRUE)
  rxode2::setRxThreads(1L)
  .a <- etaSamples(.l$fit, n = 30L, method = "mcmc", seed = 3L, nBurn = 20L)
  rxode2::setRxThreads(2L)
  expect_identical(etaSamples(.l$fit, n = 30L, method = "mcmc", seed = 3L, nBurn = 20L), .a)
  expect_identical(etaSamples(.l$fit, n = 30L, seed = 3L, nBurn = 20L), .a)
})

# The theophylline one-compartment fit, cached per estimation method
.etaSamplesTheo <- local({
  .cache <- list()
  function(est) {
    if (!is.null(.cache[[est]])) {
      return(.cache[[est]])
    }
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
    .ctl <- switch(
      est,
      saem = saemControl(print = 0L, calcTables = FALSE),
      emvi = emviControl(print = 0L, calcTables = FALSE)
    )
    .cache[[est]] <<- .nlmixr(.mod, nlmixr2data::theo_sd, est, control = .ctl)
    .cache[[est]]
  }
})

test_that("mcmc samples a SAEM fit's conditional distribution", {
  skip_on_cran()
  .fit <- .etaSamplesTheo("saem")
  .s <- etaSamples(.fit, n = 400L, method = "mcmc")
  expect_equal(attr(.s, "method"), "mcmc")
  expect_equal(nrow(.s), 400L * 12L)
  expect_true(all(is.finite(as.matrix(.s[, -(1:2)]))))
  expect_true(all(attr(.s, "acceptance") > 0.1))
  # the chains centre near the modes (the posterior is close to normal here)
  .m <- as.matrix(aggregate(cbind(eta.ka, eta.cl, eta.v) ~ ID, .s, mean)[, -1])
  .sd <- as.matrix(aggregate(cbind(eta.ka, eta.cl, eta.v) ~ ID, .s, sd)[, -1])
  expect_lt(max(abs(.m - as.matrix(.fit$eta[, -1])) / .sd), 1)
})

test_that("auto reuses a SAEM fit's own chains", {
  skip_on_cran()
  .fit <- .etaSamplesTheo("saem")
  .s <- etaSamples(.fit, n = 300L)
  expect_equal(attr(.s, "method"), "saem")
  expect_equal(names(.s), c("ID", "draw", "eta.ka", "eta.cl", "eta.v"))
  expect_equal(nrow(.s), 300L * 12L)
  expect_identical(etaSamples(.fit, n = 300L, method = "fit"), .s)
  # the last draw is the last chain's state at the last iteration
  .d <- dim(.fit$phiM)
  expect_equal(
    unlist(.s[.s$draw == 300L, c("eta.ka", "eta.cl", "eta.v")], use.names = FALSE),
    as.vector(.fit$phiM[, .d[2], .d[3], ] - .fit$saem$mprior_phi)
  )
  # chain draws centre near the EBEs
  .m <- as.matrix(aggregate(cbind(eta.ka, eta.cl, eta.v) ~ ID, .s, mean)[, -1])
  .sd <- as.matrix(aggregate(cbind(eta.ka, eta.cl, eta.v) ~ ID, .s, sd)[, -1])
  expect_lt(max(abs(.m - as.matrix(.fit$eta[, -1])) / .sd), 1)
  # more draws than the chains hold: auto samples, "fit" refuses
  .n <- .d[2] * 150L + 1L
  expect_equal(attr(etaSamples(.fit, n = .n, nBurn = 10L), "method"), "mcmc")
  expect_error(etaSamples(.fit, n = .n, method = "fit"), "holds no")
})

test_that("SAEM chain reuse follows a plateau-stopped EM phase", {
  skip_on_cran()
  .fit <- .etaSamplesTheo("saem")
  .niter <- .fit$env$saemControl$mcmc$niter
  # the chains as plateauStop leaves them: EM stopped after 100 iterations
  .stop <- .niter[1] + 100L
  .mock <- list(
    phiM = .fit$phiM[,, seq_len(.stop), , drop = FALSE],
    saem = .fit$saem,
    env = list(saemControl = .fit$env$saemControl),
    eta = .fit$eta,
    ui = .fit$ui
  )
  .d <- .etaSamplesSaem(.mock, 150L, 3L)
  expect_equal(dim(.d), c(12L, 3L, 150L))
  # the last draw is the last chain at the stopping iteration
  expect_equal(.d[,, 150L], .fit$phiM[, 3L, .stop, ] - .fit$saem$mprior_phi, ignore_attr = TRUE)
  # 50 iterations x 3 chains hold no more than 150 draws
  expect_null(.etaSamplesSaem(.mock, 151L, 3L))
  # stopped inside burn-in: no EM draws
  .mock$phiM <- .fit$phiM[,, seq_len(.niter[1]), , drop = FALSE]
  expect_null(.etaSamplesSaem(.mock, 1L, 3L))
})

test_that("SAEM chain draws follow the etas' declared order", {
  skip_on_cran()
  # nolint start: object_usage_linter.
  .mod <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      eta.v ~ 0.1
      eta.cl ~ 0.3
      eta.ka ~ 0.6
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
  .s <- etaSamples(.fit, n = 300L)
  expect_equal(attr(.s, "method"), "saem")
  expect_equal(names(.s), c("ID", "draw", "eta.v", "eta.cl", "eta.ka"))
  .m <- as.matrix(aggregate(cbind(eta.v, eta.cl, eta.ka) ~ ID, .s, mean)[, -1])
  .sd <- as.matrix(aggregate(cbind(eta.v, eta.cl, eta.ka) ~ ID, .s, sd)[, -1])
  expect_lt(max(abs(.m - as.matrix(.fit$eta[, c("eta.v", "eta.cl", "eta.ka")])) / .sd), 1)
  # eta.ka is the widest conditional distribution, eta.v the narrowest
  expect_true(all(.sd[, "eta.ka"] > .sd[, "eta.v"]))
})

test_that("auto reuses an emvi fit's variational posterior", {
  skip_on_cran()
  .fit <- .etaSamplesTheo("emvi")
  .s <- etaSamples(.fit, n = 2000L, seed = 5L)
  expect_equal(attr(.s, "method"), "vi")
  expect_equal(nrow(.s), 2000L * 12L)
  expect_identical(etaSamples(.fit, n = 2000L, seed = 5L), .s)
  .v <- .fit$env$viState
  .m <- as.matrix(aggregate(cbind(eta.ka, eta.cl, eta.v) ~ ID, .s, mean)[, -1])
  .sd <- as.matrix(aggregate(cbind(eta.ka, eta.cl, eta.v) ~ ID, .s, sd)[, -1])
  expect_lt(max(abs(.m - .v$mu) / .sd), 5 / sqrt(2000))
  # full-rank family: covariance L L' from the packed lower factor
  expect_equal(.v$family, "fullRank")
  .l <- matrix(0, 3, 3)
  .l[upper.tri(.l, diag = TRUE)] <- .v$Lpack[1, ]
  .l <- t(.l)
  expect_equal(
    stats::cov(.s[.s$ID == .s$ID[1], c("eta.ka", "eta.cl", "eta.v")]),
    tcrossprod(.l),
    tolerance = 0.1,
    ignore_attr = TRUE
  )
})

test_that("a fit without its own draws refuses method = 'fit'", {
  skip_on_cran()
  .l <- .etaSamplesLin()
  expect_error(etaSamples(.l$fit, n = 10L, method = "fit"), "holds no")
  expect_equal(attr(etaSamples(.l$fit, n = 10L, nBurn = 10L), "method"), "mcmc")
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
  # near the EBE: 0.3 is ~4 Monte-Carlo SEs of 500 autocorrelated draws
  # (conditional variance 0.8)
  expect_lt(max(abs(as.numeric(tapply(.s$eta.v, .s$ID, mean)) - .fit$eta$eta.v)), 0.3)
})

test_that("a SAEM fit gets its modes from a posthoc step and simulates per draw", {
  skip_on_cran()
  .fit <- .etaSamplesTheo("saem")
  .s <- etaSamples(.fit, n = 20L, method = "laplace")
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
