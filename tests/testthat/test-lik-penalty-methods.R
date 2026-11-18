## Likelihood penalties (inst/include/nlmixr2estLikPenalty.h) across every
## estimation method.  Three probes with the built-in penalties:
##   - a constant subject log-likelihood c (etaL2 with lambda 0) moves the
##     reported objective by -2*N*c when counted in it,
##   - a constant population penalty P (thetaL2 offset) moves it by P,
##   - a strong ridge on tka pulls the estimate toward 0 for a method that
##     estimates tka.

.likPenEta <- function() {
  ini({
    tka <- 0.45; tcl <- 1; tv <- 3.45
    eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
    linCmt() ~ add(add.sd)
  })
}

.likPenPop <- function() {
  ini({
    tka <- 0.45; tcl <- 1; tv <- 3.45
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka); cl <- exp(tcl); v <- exp(tv)
    linCmt() ~ add(add.sd)
  })
}

.likPenProbe <- function(est, mod, ridge = TRUE, tol = 0.05) {
  .cf <- paste0(est, "Control")
  .ctl <- if (exists(.cf, mode = "function")) {
    do.call(.cf, list(print = 0L))
  } else {
    getValidNlmixrControl(list(print = 0L), est)
  }
  .fit <- function(...) .nlmixr(mod, theo_sd, est = est, control = .ctl, ...)
  .nsub <- length(unique(theo_sd$ID))
  f0 <- .fit()
  f1 <- .fit(likPenalty = list(nlmixr2LikPenaltySpec("etaL2", par = c(0, 1), inObjf = TRUE)))
  expect_equal(f1$objf - f0$objf, -2 * .nsub, tolerance = tol, scale = 1, label = paste(est, "subject"))
  f2 <- .fit(likPenalty = list(nlmixr2LikPenaltySpec("thetaL2", thetas = "tka", par = c(0, 0, 100), inObjf = TRUE)))
  expect_equal(f2$objf - f0$objf, 100, tolerance = tol, scale = 1, label = paste(est, "population"))
  if (ridge) {
    f3 <- .fit(likPenalty = list(nlmixr2LikPenaltySpec("thetaL2", thetas = "tka", par = 1000)))
    expect_lt(abs(fixef(f3)[["tka"]]), 0.1, label = paste(est, "ridge"))
    expect_true(is.data.frame(f3$likPenalty))
  }
  expect_equal(.likPenaltyInfo()[["n"]], 0L)
}

test_that("likelihood penalties reach the FOCEi family", {
  skip_on_cran()
  for (.est in c("foce", "fo", "focep", "laplace", "agq")) {
    .likPenProbe(.est, .likPenEta)
  }
  ## posthoc estimates no thetas
  .likPenProbe("posthoc", .likPenEta, ridge = FALSE)
})

test_that("likelihood penalties reach the EM and variational methods", {
  skip_on_cran()
  for (.est in c("imp", "impmap", "qrpem", "saem", "emvi", "fbvi")) {
    .likPenProbe(.est, .likPenEta, tol = 0.5)
  }
  ## vae's objective is noisier between fits
  .likPenProbe("vae", .likPenEta, tol = 1)
})

test_that("likelihood penalties reach the nonparametric methods", {
  skip_on_cran()
  for (.est in c("npag", "npb")) {
    .likPenProbe(.est, .likPenEta, tol = 0.5)
  }
})

test_that("likelihood penalties reach the nlm family", {
  skip_on_cran()
  for (.est in c("nlm", "nlminb", "optim", "bobyqa", "newuoa", "n1qn1", "lbfgsb3c", "trust")) {
    .likPenProbe(.est, .likPenPop, tol = 0.5)
  }
})

test_that("methods without likelihood penalty support refuse them", {
  skip_on_cran()
  .pen <- list(nlmixr2LikPenaltySpec("thetaL2", thetas = "tka", par = 1))
  expect_error(
    .nlmixr(.likPenEta, theo_sd, est = "nlme", control = list(print = 0L), likPenalty = .pen),
    "likelihood penalties"
  )
  expect_equal(.likPenaltyInfo()[["n"]], 0L)
  ## simulation ignores them
  expect_error(.nlmixr(.likPenEta, theo_sd, est = "rxSolve", likPenalty = .pen), NA)
})
