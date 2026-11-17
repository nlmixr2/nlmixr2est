## External likelihood-contribution API: a contributor registers begin/obs/end
## hooks that nlmixr2est cycles in series inside likInner0's per-observation
## loop.  Validated with a test contributor (src/inner.cpp) that (a) counts the
## begin/obs/end calls and (b) adds a constant c to each observation's log-
## likelihood -- which must shift the objective by exactly -2*c*nObs, since a
## constant LL addition does not move the inner optimum.

test_that("likelihood-contribution hooks fire per observation and fold into the objective", {
  skip_on_cran()

  one.compartment <- function() {
    ini({
      tka <- fix(0.45); tcl <- fix(1); tv <- fix(3.45)
      eta.ka ~ fix(0.6); eta.cl ~ fix(0.3); eta.v ~ fix(0.1)
      add.sd <- fix(0.7)
    })
    model({
      ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }

  .old <- rxode2::getRxThreads()
  on.exit(rxode2::setRxThreads(.old), add = TRUE)
  rxode2::setRxThreads(1L) # test contributor uses global accumulators

  .nObs <- sum(theo_sd$EVID == 0)
  .nsub <- length(unique(theo_sd$ID))

  f0 <- .nlmixr(one.compartment, theo_sd, est = "focei", control = foceiControl(print = 0L))

  cc <- 0.01
  .Call("_nlmixr2est_registerTestContrib", PACKAGE = "nlmixr2est")
  on.exit(.Call("_nlmixr2est_removeTestContrib", PACKAGE = "nlmixr2est"), add = TRUE)
  .Call("_nlmixr2est_setTestContribAddLL", cc, PACKAGE = "nlmixr2est")
  f1 <- .nlmixr(one.compartment, theo_sd, est = "focei", control = foceiControl(print = 0L))
  res <- .Call("_nlmixr2est_getTestContrib", PACKAGE = "nlmixr2est")
  names(res) <- c("nObs", "sumDLLdf", "sumErr", "sumF", "nBegin", "nEnd")

  ## a constant c added per obs shifts the objective by exactly -2*c*nObs
  expect_equal(f1$objf - f0$objf, -2 * cc * .nObs, tolerance = 1e-3)

  ## the obs hook fires once per observation within each subject bracket
  expect_gt(res[["nObs"]], 0)
  expect_equal(res[["nObs"]] / res[["nBegin"]], .nObs / .nsub) # obs per subject
  ## begin/end are balanced
  expect_equal(res[["nBegin"]], res[["nEnd"]])

  ## after removing the contributor, a fit does not invoke it
  .Call("_nlmixr2est_removeTestContrib", PACKAGE = "nlmixr2est")
  .Call("_nlmixr2est_registerTestContrib", PACKAGE = "nlmixr2est") # resets counters
  .Call("_nlmixr2est_removeTestContrib", PACKAGE = "nlmixr2est")
  f2 <- .nlmixr(one.compartment, theo_sd, est = "focei", control = foceiControl(print = 0L))
  res2 <- .Call("_nlmixr2est_getTestContrib", PACKAGE = "nlmixr2est")
  expect_equal(res2[[1]], 0) # nObs: not invoked
})

## Likelihood penalty ABI (inst/include/nlmixr2estLikPenalty.h): per-subject and
## population terms, activated per fit from a model spec or `likPenalty=`.

test_that("likelihood penalty entry points are appended to the pointer table", {
  .p <- .nlmixr2estLikContribPtrs()
  expect_length(.p, 7L)
  expect_equal(names(.p)[6:7], c("registerLikPenalty", "removeLikPenalty"))
  expect_true(all(vapply(.p, inherits, logical(1), "externalptr")))
})

test_that("likelihood penalty specs validate and attach to a model", {
  expect_error(nlmixr2LikPenaltySpec(1))
  expect_error(nlmixr2LikPenaltySpec("thetaL2", data = 1))
  .s <- nlmixr2LikPenaltySpec("thetaL2", thetas = "tka", par = 10)
  expect_s3_class(.s, "nlmixr2LikPenaltySpec")
  expect_output(print(.s), "thetaL2")
  mod <- function() {
    ini({
      tka <- 0.45
      eta.ka ~ 0.6
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      ka ~ add(add.sd)
    })
  }
  .ui <- nlmixr2LikPenaltyAdd(mod, .s)
  expect_length(.ui$nlmixr2LikPenalty, 1L)
  .ui <- nlmixr2LikPenaltyAdd(.ui, .s)
  expect_length(.ui$nlmixr2LikPenalty, 2L)
  .ui <- nlmixr2LikPenaltyAdd(.ui, .s, replace = TRUE)
  expect_length(.ui$nlmixr2LikPenalty, 1L)
  expect_error(nlmixr2LikPenaltyAdd(mod, list(1)))
  expect_error(.nlmixr2LikPenaltyBuiltin("noSuchPenalty"))
})

test_that("likelihood penalties shift and shape a focei fit, then turn off", {
  skip_on_cran()
  one.cmt <- function() {
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
  .ctl <- foceiControl(print = 0L, calcTables = FALSE)
  .nsub <- length(unique(theo_sd$ID))
  f0 <- .nlmixr(one.cmt, theo_sd, est = "focei", control = .ctl)

  ## a constant subject log-likelihood c moves the objective by -2*N*c ...
  f1 <- .nlmixr(
    one.cmt,
    theo_sd,
    est = "focei",
    control = .ctl,
    likPenalty = list(nlmixr2LikPenaltySpec("etaL2", par = c(0, 1), inObjf = TRUE))
  )
  expect_equal(f1$objf - f0$objf, -2 * .nsub, tolerance = 1e-4)
  expect_equal(f1$likPenalty$subjectLL, .nsub)
  ## ... and not at all when it is kept out of the reported objective
  f1b <- .nlmixr(
    one.cmt,
    theo_sd,
    est = "focei",
    control = .ctl,
    likPenalty = list(nlmixr2LikPenaltySpec("etaL2", par = c(0, 1)))
  )
  expect_equal(f1b$objf, f0$objf, tolerance = 1e-4)

  ## a constant population penalty P moves it by P
  f2 <- .nlmixr(
    one.cmt,
    theo_sd,
    est = "focei",
    control = .ctl,
    likPenalty = list(nlmixr2LikPenaltySpec("thetaL2", thetas = "tka", par = c(0, 0, 100), inObjf = TRUE))
  )
  expect_equal(f2$objf - f0$objf, 100, tolerance = 1e-4)

  ## a ridge on tka pulls it to 0; the reported objective stays unpenalized
  .ui <- nlmixr2LikPenaltyAdd(one.cmt, nlmixr2LikPenaltySpec("thetaL2", thetas = "tka", par = 1000))
  f3 <- .nlmixr(.ui, theo_sd, est = "focei", control = .ctl)
  expect_lt(abs(fixef(f3)[["tka"]]), 0.05)
  ## unpenalized objective at the penalized estimates: worse than the optimum
  expect_gt(f3$objf, f0$objf)
  expect_false(f3$likPenalty$inObjf)
  expect_true(any(grepl("likelihood penalty excluded", f3$runInfo)))

  ## an eta shrinkage moves the etas toward 0
  f4 <- .nlmixr(
    one.cmt,
    theo_sd,
    est = "focei",
    control = .ctl,
    likPenalty = list(nlmixr2LikPenaltySpec("etaL2", par = 50))
  )
  expect_lt(sum(f4$eta$eta.ka^2), sum(f0$eta$eta.ka^2))

  ## inactive between fits, and after a refused method
  expect_equal(.likPenaltyInfo()[["n"]], 0L)
  expect_error(
    suppressWarnings(suppressMessages(nlmixr2(
      .ui,
      theo_sd,
      est = "nls",
      control = nlsControl(print = 0L)
    ))),
    "likelihood penalties"
  )
  expect_equal(.likPenaltyInfo()[["n"]], 0L)
  expect_null(nlmixr2global$likPenaltyActive)
})

test_that("a population penalty keeps the analytic outer gradient exact", {
  skip_on_cran()
  d <- data.frame(
    ID = rep(1:3, each = 3),
    TIME = rep(c(0.5, 1, 2), 3),
    DV = c(1.4, 1.5, 1.3, 1.7, 1.2, 1.6, 1.1, 1.9, 1.5),
    EVID = 0L
  )
  mod <- function() {
    ini({
      level <- 0.3
      etaLevel ~ 0.1
      error <- fix(0.2)
    })
    model({
      prediction <- exp(level + etaLevel)
      prediction ~ prop(error)
    })
  }
  .pen <- list(nlmixr2LikPenaltySpec("thetaL2", thetas = "level", par = c(20, -0.5), inObjf = TRUE))
  .run <- function(fast) {
    .nlmixr(
      mod,
      d,
      est = "focei",
      likPenalty = .pen,
      control = foceiControl(print = 0L, fast = fast, calcTables = FALSE)
    )
  }
  .t <- .run(TRUE)
  .f <- .run(FALSE)
  .b <- .nlmixr(mod, d, est = "focei", control = foceiControl(print = 0L, calcTables = FALSE))
  expect_gt(as.integer(.t$env$nAnalyticGradDirect), 0L)
  ## same optimum as the finite-difference gradient; the analytic arm reaches a
  ## penalized objective at least as low
  expect_equal(fixef(.t)[["level"]], fixef(.f)[["level"]], tolerance = 5e-3)
  expect_lt(.t$objf, .f$objf + 1e-4)
  ## the penalty really moved the estimate toward its center
  expect_lt(fixef(.f)[["level"]], fixef(.b)[["level"]] - 0.05)
})

test_that("nlm objective setup carries the model's likelihood penalty", {
  skip_on_cran()
  pop <- function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka); cl <- exp(tcl); v <- exp(tv)
      linCmt() ~ add(add.sd)
    })
  }
  .p <- suppressMessages(nlmObjectiveSetup(pop, theo_sd))
  .v0 <- nlmSolveR(.p)
  .nlmFreeEnv()
  .p <- suppressMessages(nlmObjectiveSetup(
    pop,
    theo_sd,
    likPenalty = list(
      nlmixr2LikPenaltySpec("thetaL2", thetas = "tka", par = c(0, 0, 10))
    )
  ))
  .v1 <- nlmSolveR(.p)
  .nlmFreeEnv()
  ## nlm's objective is -LL, so a -2LL penalty of 10 adds 5
  expect_equal(.v1 - .v0, 5, tolerance = 1e-8)
  expect_equal(.likPenaltyInfo()[["n"]], 0L)
})
