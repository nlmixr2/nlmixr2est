nmTest({
  test_that("nlmixr interface for solving works", {
    one.cmt <- function() {
      ini({
        ## You may label each parameter with a comment
        tka <- 0.45 # Ka
        tcl <- log(c(0, 2.7, 100)) # Log Cl
        ## This works with interactive models
        ## You may also label the preceding line with label("label text")
        tv <- 3.45; label("log V")
        ## the label("Label name") works with all models
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

    f <- .nlmixr(one.cmt, nlmixr2data::theo_sd, "rxSolve")

    expect_s3_class(f, "rxSolve")

    f2 <- .nlmixr(one.cmt, nlmixr2data::theo_sd, "rxSolve", rxControl(returnType = "data.frame"))

    expect_s3_class(f, "data.frame")

    one.cmt <- function() {
      ini({
        ## You may label each parameter with a comment
        tka <- 0.45 # Ka
        tcl <- log(c(0, 2.7, 100)) # Log Cl
        ## This works with interactive models
        ## You may also label the preceding line with label("label text")
        tv <- 3.45; label("log V")
        ## the label("Label name") works with all models
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

    expect_error(.nlmixr(one.cmt, nlmixr2data::theo_sd, "matt"))
  })

  test_that("rxSolve will warn when necessary", {
    library(rxode2)
    eventTable <- et(amt = 320, evid = 1, cmt = 20, time = 0) |> # nolint: object_name_linter.
      et(2, 4)

    # Now define the nlmixr2/rxode2 model used for both estimation and simulation
    mod <- function() {
      # Parameters
      ini({
        tka <- 0.45; label("Ka (first order absorption)")
        trate <- 0.4; label("Zero order rate")
        tcl <- 1; label("Cl")
        tv <- 3.45; label("V")
        fDepot <- logit(0.5); label("amount of dose in first order absorption")
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      # and a model block with the error specification and model specification
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center

        f(depot) <- expit(fDepot)
        f(center) <- 1 - expit(fDepot)
        rate(center) <- exp(trate)

        cp <- center / v
        cp ~ add(add.sd)
      })
    }

    mod <- mod()

    expect_warning(
      suppressMessages(
        nlmixr2(mod, eventTable, "rxSolve", control = rxControl(addDosing = TRUE))
      )
    )
  })
})

nmTest({
  # nlmixr2/rxode2 issue 1289: `$simInfo` re-runs the pre-process hooks and re-derives
  # the simulation model, which dominated the cost of solving a fit and grew
  # process memory that was never returned -- and a plain `rxSolve(fit, ...)`
  # never reads any of it.

  .simInfoEv <- function() {
    rxode2::et(rxode2::et(amt = 320, cmt = "depot"), seq(0, 24, by = 4))
  }

  test_that("solving a fit skips simulation information it cannot use", {
    skip_if(is.null(one.compartment.fit.saem))
    local_mocked_bindings(
      .simInfo = function(object) stop("simulation information was derived")
    )
    expect_no_error(
      suppressMessages(rxode2::rxSolve(one.compartment.fit.saem, .simInfoEv()))
    )
  })

  test_that("solving a fit with uncertainty still uses its simulation information", {
    skip_if(is.null(one.compartment.fit.saem))
    .acc <- new.env(parent = emptyenv())
    .acc$msg <- character(0)
    withCallingHandlers(
      .sim <- rxode2::rxSolve(one.compartment.fit.saem, .simInfoEv(), nStud = 2),
      message = function(m) {
        .acc$msg <- c(.acc$msg, conditionMessage(m))
        invokeRestart("muffleMessage")
      }
    )
    expect_true(any(grepl("population uncertainty", .acc$msg, fixed = TRUE)))
    expect_true(any(grepl("from the number of observations", .acc$msg, fixed = TRUE)))
    expect_true(any(grepl("from the number of subjects", .acc$msg, fixed = TRUE)))
    expect_equal(length(unique(.sim$sim.id)), 2L)
  })
})

test_that(".predictZeroNestedEtas() zeroes every eta of a nested omega", {
  .d <- data.frame(ID = 1, TIME = 0:1)
  .nested <- list(
    omega = lotri::lotri(lotri::lotri(eta.v ~ 0.1) | id(nu = 100), lotri::lotri(iov.cl ~ 0.1) | occ(nu = 100))
  )
  expect_identical(.predictZeroNestedEtas(.nested, .d), data.frame(ID = 1, TIME = 0:1, eta.v = 0, iov.cl = 0))
  # a matrix omega's etas rxode2 zeroes itself
  .flat <- list(omega = lotri::lotri(eta.v ~ 0.1))
  expect_identical(.predictZeroNestedEtas(.flat, .d), .d)
  expect_identical(.predictZeroNestedEtas(list(omega = NULL), .d), .d)
})

# A fit with IOV as far as the eta helpers read it: subject etas and
# occasion estimates for subjects "1" and "2", occasions 1 and 2
.rxsolveIovFitStub <- function() {
  .iov <- data.frame(ID = factor(c("1", "1", "2", "2")), occ = c(1L, 2L, 1L, 2L), iov.cl = c(0.1, 0.2, 0.3, 0.4))
  list(env = list(iov = list(occ = .iov)), eta = data.frame(ID = factor(c("1", "2")), eta.v = c(0.5, 0.6)))
}

test_that(".iovEtaRows() leaves a row with no occasion NA, not unseen", {
  .r <- .iovEtaRows(.rxsolveIovFitStub(), c("1", "1", "2"), list(occ = c(1, NA, 2)))
  expect_identical(.r$etas, data.frame(iov.cl = c(0.1, NA, 0.4)))
  expect_identical(.r$unseen, 0L)
})

test_that(".augPredIovEvents() warns about a whole occasion the fit did not see, not a fractional one", {
  .fit <- .rxsolveIovFitStub()
  .fit$idLvl <- c("1", "2")
  # the events' subject levels are the fit's integer codes
  .ev <- data.frame(ID = factor(c("1", "2", "2"), levels = c("1", "2")), occ = c(1, 1.5, 3))
  .w <- capture_warnings(.r <- .augPredIovEvents(.fit, .ev))
  expect_identical(.w, "1 subject-occasion pair(s) not in the fit get occasion etas of 0")
  expect_identical(.r$iov.cl, c(0.1, 0, 0))
  .ev$occ[3] <- 2
  expect_identical(capture_warnings(.r <- .augPredIovEvents(.fit, .ev)), character(0))
  expect_identical(.r$iov.cl, c(0.1, 0, 0.4))
})

test_that(".subjectEtaRows() joins on the subject, and counts unseen subjects", {
  .r <- .subjectEtaRows(.rxsolveIovFitStub(), c("2", "1", "3", "3"))
  expect_identical(.r$etas, data.frame(eta.v = c(0.6, 0.5, 0, 0)))
  expect_identical(.r$unseen, 1L)
})

test_that(".iovEtaRows() joins on subject and occasion, and counts unseen pairs", {
  .r <- .iovEtaRows(.rxsolveIovFitStub(), c("2", "1", "1", "3", "1"), list(occ = c(1, 2, 3, 1, 3)))
  # subject 3 and subject 1's occasion 3 are not in the fit: eta 0, two pairs
  expect_identical(.r$etas, data.frame(iov.cl = c(0.3, 0.2, 0, 0, 0)))
  expect_identical(.r$unseen, 2L)
})

test_that(".nlmixrIpredEtaData() adds the occasion etas and warns about unseen pairs", {
  .d <- data.frame(id = c(1, 1, 2), TIME = 0:2, OCC = c(1, 2, 2))
  .r <- .nlmixrIpredEtaData(.rxsolveIovFitStub(), .d)
  expect_identical(.r$iov.cl, c(0.1, 0.2, 0.4))
  # the subject etas too, so the solve does not take them from its params
  expect_identical(.r$eta.v, c(0.5, 0.5, 0.6))
  .d$OCC[3] <- 5
  .w <- capture_warnings(.r <- .nlmixrIpredEtaData(.rxsolveIovFitStub(), .d))
  expect_identical(.w, "1 subject-occasion pair(s) not in the fit get occasion etas of 0")
  expect_identical(.r$iov.cl, c(0.1, 0.2, 0))
  .e <- expect_error(.nlmixrIpredEtaData(.rxsolveIovFitStub(), .d[, c("id", "TIME")]))
  expect_identical(conditionMessage(.e), "the data needs occ for the fit's individual etas")
})

test_that(".nlmixrIpredEtaData() warns only for subjects and pairs it predicts", {
  # subject 3 has a dose and no observation: it is not predicted
  .d <- data.frame(ID = c(1, 2, 3), TIME = 0, EVID = c(0, 0, 1), AMT = c(0, 0, 1), DV = c(1, 1, NA), OCC = 1)
  expect_identical(capture_warnings(.r <- .nlmixrIpredEtaData(.rxsolveIovFitStub(), .d)), character(0))
  expect_identical(.r$iov.cl, c(0.1, 0.3, 0))
  .d$EVID[3] <- 0
  .d$DV[3] <- 1
  expect_identical(
    capture_warnings(.nlmixrIpredEtaData(.rxsolveIovFitStub(), .d)),
    "1 subject-occasion pair(s) not in the fit get occasion etas of 0"
  )
  # a fit without IOV counts subjects
  .fit <- .rxsolveIovFitStub()
  .fit$env$iov <- NULL
  .d$EVID[3] <- 1
  .d$DV[3] <- NA
  expect_identical(capture_warnings(.r <- .nlmixrIpredEtaData(.fit, .d)), character(0))
  expect_identical(.r$eta.v, c(0.5, 0.6, 0))
  .d$EVID[3] <- 0
  .d$DV[3] <- 1
  expect_identical(capture_warnings(.nlmixrIpredEtaData(.fit, .d)), "1 subject(s) not in the fit get etas of 0")
})
