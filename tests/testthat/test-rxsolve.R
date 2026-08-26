nmTest({
  test_that("nlmixr interface for solving works", {

    one.cmt <- function() {
      ini({
        ## You may label each parameter with a comment
        tka <- 0.45 # Ka
        tcl <- log(c(0, 2.7, 100)) # Log Cl
        ## This works with interactive models
        ## You may also label the preceding line with label("label text")
        tv <- 3.45
        label("log V")
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

    f2 <- .nlmixr(one.cmt, nlmixr2data::theo_sd, "rxSolve", rxControl(returnType="data.frame"))

    expect_s3_class(f, "data.frame")

    one.cmt <- function() {
      ini({
        ## You may label each parameter with a comment
        tka <- 0.45 # Ka
        tcl <- log(c(0, 2.7, 100)) # Log Cl
        ## This works with interactive models
        ## You may also label the preceding line with label("label text")
        tv <- 3.45
        label("log V")
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
    eventTable <- et(amt=320, evid=1, cmt=20, time = 0) |> # nolint: object_name_linter.
      et(2, 4)

    # Now define the nlmixr2/rxode2 model used for both estimation and simulation
    mod <- function() {
      # Parameters
      ini({
        tka <- 0.45
        label("Ka (first order absorption)")
        trate <- 0.4
        label("Zero order rate")
        tcl <- 1
        label("Cl")
        tv <- 3.45
        label("V")
        fDepot <- logit(0.5)
        label("amount of dose in first order absorption")
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
        d/dt(depot) = -ka * depot
        d/dt(center) = ka * depot - cl / v * center

        f(depot) <- expit(fDepot)
        f(center) <- 1-expit(fDepot)
        rate(center) <- exp(trate)

        cp = center / v
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
