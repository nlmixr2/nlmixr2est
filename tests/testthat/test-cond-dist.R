# addCondDist() samples each subject's conditional eta distribution at a fit's
# final estimates, for any estimation method (#1378).

test_that(".weightedDraws gives weighted moments and systematic draws", {
  .s <- list(cbind(c(0, 1, 2, 3), c(1, 1, 1, 1)), cbind(c(5, 5, 5, 5), c(-1, 1, -1, 1)))
  .w <- list(c(0.1, 0.2, 0.3, 0.4), c(1, 1, 1, 1))
  .r <- .weightedDraws(.s, .w, 10L)
  .r2 <- .weightedDraws(.s, .w, 8L)
  expect_equal(.r$mean[1, ], c(sum(0:3 * .w[[1]]), 1))
  expect_equal(.r$sd[1, 1], sqrt(sum(.w[[1]] * (0:3 - 2)^2)))
  expect_equal(.r$sd[1, 2], 0)
  expect_equal(.r$mean[2, ], c(5, 0))
  expect_equal(.r$ess, c(1 / sum(.w[[1]]^2), 4))
  expect_equal(dim(.r$draws), c(2L, 2L, 10L))
  # draws follow the weights: 1, 2, 3, 4 of the 10 draws at 0, 1, 2, 3
  expect_equal(as.vector(table(.r$draws[1, 1, ])), 1:4)
  expect_equal(mean(.r2$draws[2, 2, ]), 0)
})

test_that(".condDistSummary is long by ID then ETA", {
  .m <- matrix(1:4, 2)
  .sd <- matrix(5:8, 2)
  .s <- .condDistSummary(.m, .sd, c(10, 20), c("eta.a", "eta.b"))
  expect_equal(names(.s), c("ID", "ETA", "MEAN", "SD"))
  expect_equal(.s$ID, c(10, 10, 20, 20))
  expect_equal(as.character(.s$ETA), c("eta.a", "eta.b", "eta.a", "eta.b"))
  expect_equal(.s$MEAN, c(1, 3, 2, 4))
  expect_equal(.s$SD, c(5, 7, 6, 8))
  expect_null(attr(.s, "samples"))
  .d <- array(0, c(2, 2, 3))
  .s <- .condDistSummary(.m, .sd, c(10, 20), c("eta.a", "eta.b"), .d)
  expect_equal(dimnames(attr(.s, "samples")), list(c("10", "20"), c("eta.a", "eta.b"), NULL))
})

test_that("addCondDist samples a FOCEi fit's conditional eta distribution", {
  skip_on_cran()

  # nolint start: object_usage_linter.
  .mod <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      add.sd <- 0.7
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }
  # nolint end
  .dat <- nlmixr2data::theo_sd
  fit <- .nlmixr(.mod, .dat, est = "focei", control = foceiControl(print = 0))
  expect_null(fit$condDist)
  .theta <- fit$theta
  .eta <- fit$eta
  fit <- addCondDist(fit, keepSamples = TRUE)
  .cd <- fit$condDist
  .etaNames <- c("eta.ka", "eta.cl", "eta.v")
  .n <- length(unique(.dat$ID))
  expect_s3_class(.cd, "data.frame")
  expect_equal(names(.cd), c("ID", "ETA", "MEAN", "SD"))
  expect_equal(nrow(.cd), .n * 3L)
  expect_equal(attr(.cd, "method"), "imp")
  expect_true(all(is.finite(.cd$SD) & .cd$SD > 0))
  expect_true(all(attr(.cd, "ess") > 100))
  expect_equal(dim(attr(.cd, "samples")), c(.n, 3L, 300L))
  # the fit itself is unchanged
  expect_equal(fit$theta, .theta)
  expect_equal(fit$eta, .eta)
  # the conditional mean is near the conditional mode the fit reports
  for (.e in .etaNames) {
    .z <- (.cd$MEAN[.cd$ETA == .e] - .eta[[.e]]) / .cd$SD[.cd$ETA == .e]
    expect_lt(mean(abs(.z)), 0.5)
  }
  # the draws carry the same moments
  .draws <- attr(.cd, "samples")
  expect_equal(as.vector(t(apply(.draws, c(1, 2), mean))), .cd$MEAN, tolerance = 0.05)

  # reproducible under a seed; another seed agrees within Monte-Carlo error
  expect_equal(addCondDist(fit, keepSamples = TRUE)$condDist, .cd)
  .cd2 <- addCondDist(fit, seed = 7L)$condDist
  expect_null(attr(.cd2, "samples"))
  expect_lt(max(abs(.cd2$MEAN - .cd$MEAN) / .cd$SD), 0.5)
  expect_equal(.cd2$SD, .cd$SD, tolerance = 0.25)
})

test_that("addCondDist refuses a fit with no random effects", {
  expect_error(addCondDist(list()), "nlmixr2 fit")
})
