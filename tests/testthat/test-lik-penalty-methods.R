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
  expect_lt(abs(f1$objf - f0$objf + 2 * .nsub), tol, label = paste(est, "subject"))
  f2 <- .fit(likPenalty = list(nlmixr2LikPenaltySpec("thetaL2", thetas = "tka", par = c(0, 0, 100), inObjf = TRUE)))
  expect_lt(abs(f2$objf - f0$objf - 100), tol, label = paste(est, "population"))
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

test_that("saem penalizes a theta without a random effect", {
  skip_on_cran()
  .m <- function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45
      eta.cl ~ 0.3; eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }
  .f <- .nlmixr(
    .m,
    theo_sd,
    est = "saem",
    control = saemControl(print = 0L),
    likPenalty = list(nlmixr2LikPenaltySpec("thetaL2", thetas = "tka", par = 1000))
  )
  expect_lt(abs(fixef(.f)[["tka"]]), 0.1)
})

test_that("saem penalties: bounded thetas, ll() models, sampled phi0 and mode etas", {
  skip_on_cran()
  .ctl <- saemControl(print = 0L, nBurn = 100, nEm = 100, seed = 1)
  .ridge <- list(nlmixr2LikPenaltySpec("thetaL2", thetas = "tka", par = 1000))
  ## a bounded theta is estimated as rxBoundedTr.tka; the penalty still finds it
  ## and sees its natural value
  .mB <- function() {
    ini({
      tka <- c(-5, 0.45, 5); tcl <- 1; tv <- 3.45
      eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }
  .f <- .nlmixr(.mB, theo_sd, est = "saem", control = .ctl, likPenalty = .ridge)
  expect_lt(abs(fixef(.f)[["tka"]]), 0.05)
  expect_equal(.f$likPenalty$population, 1000 * fixef(.f)[["tka"]]^2, tolerance = 1e-6)
  ## a compiled spec carried in the model (an external pointer) survives the
  ## bounded transform's model rebuild
  .ui <- nlmixr2LikPenaltyAdd(
    .mB,
    nlmixr2LikPenaltySpec(.Call(`_nlmixr2est_likPenaltyBuiltin`, "thetaL2"), thetas = "tka", par = 1000)
  )
  .f <- .nlmixr(.ui, theo_sd, est = "saem", control = .ctl)
  expect_lt(abs(fixef(.f)[["tka"]]), 0.05)
  ## general likelihood: the pooled phi1 refinement keeps the penalty
  .mL <- function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45
      eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
      lsd <- -0.3
    })
    model({
      ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
      cp <- linCmt(); sd <- exp(lsd)
      ll(err) ~ -log(sd) - 0.5 * log(2 * pi) - 0.5 * ((DV - cp) / sd)^2
    })
  }
  .f <- .nlmixr(.mL, theo_sd, est = "saem", control = .ctl, likPenalty = .ridge)
  expect_lt(abs(fixef(.f)[["tka"]]), 0.05)
  ## nonMuTheta = "eta": the penalty shapes the phi0 sampling target
  .m0 <- function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45
      eta.cl ~ 0.3; eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }
  .ctl0 <- saemControl(print = 0L, nBurn = 100, nEm = 100, seed = 1, nonMuTheta = "eta")
  .f0 <- .nlmixr(.m0, theo_sd, est = "saem", control = .ctl0)
  .f1 <- .nlmixr(.m0, theo_sd, est = "saem", control = .ctl0, likPenalty = .ridge)
  expect_lt(abs(fixef(.f1)[["tka"]]), 0.25 * abs(fixef(.f0)[["tka"]]))
  ## the reported (mode) etas carry the subject penalty
  .f <- .nlmixr(
    .mB,
    theo_sd,
    est = "saem",
    control = .ctl,
    likPenalty = list(nlmixr2LikPenaltySpec("etaL2", par = c(1000, 0)))
  )
  expect_lt(max(abs(as.matrix(.f$eta[, -1]))), 0.05)
})

test_that("saem orders the subject penalty etas as ini() does", {
  .m <- function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45
      eta.v ~ 0.1; eta.cl ~ 0.3; eta.ka ~ 0.6
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }
  .o <- rxode2::rxode2(.m)$saemOmegaTrans
  expect_equal(unname(.o), c(3, 2, 1))
  .e1 <- matrix(c(1, 2, 3), 1, 3)
  expect_equal(.saemLikPenaltyEta(.e1, as.integer(.o) - 1L), matrix(c(3, 2, 1), 1, 3))
})
