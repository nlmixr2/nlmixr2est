# The analytic-gradient scope probes read inter-occasion variability from the ui
# they are given (.uiHasIov()), not from .uiIovEnv, which keeps the last fit's
# IOV rewrite.  After an IOV fit they said a model without IOV was out of scope:
# test-vae-iov.R then test-vae-nonmutheta-grad.R in one process failed 5 checks.

nmTest({
  .isNoIov <- function() {
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
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }

  .isIov <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      iov.cl ~ 0.08 | occ
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl + iov.cl)
      v <- exp(tv + eta.v)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }

  # an ll() endpoint: in scope through the log-density route
  .isLl <- function() {
    ini({
      lka <- 0.4
      lcl <- -0.9
      lb <- log(3)
      lsd <- log(1.2)
      eta.b ~ 0.1
    })
    model({
      ka <- exp(lka)
      cl <- exp(lcl)
      d/dt(depot) <- -ka * depot
      d/dt(central) <- ka * depot - cl * central
      mu <- exp(lb + eta.b) * central
      sd <- exp(lsd)
      ll(cp) ~ -0.5 * log(2 * pi) - log(sd) - 0.5 * ((DV - mu) / sd)^2
    })
  }

  # theo_sd with two occasions per subject
  .isData <- nlmixr2data::theo_sd
  .isData$occ <- ifelse(.isData$TIME < 6, 1L, 2L)

  # .uiApplyIov() and a fit both leave their rewrite in .uiIovEnv; put it back
  .isKeepIovEnv <- function(envir = parent.frame()) {
    .old <- as.list(.uiIovEnv, all.names = TRUE)
    withr::defer(
      {
        rm(list = ls(.uiIovEnv, all.names = TRUE), envir = .uiIovEnv)
        list2env(.old, envir = .uiIovEnv)
      },
      envir = envir
    )
  }

  test_that(".uiHasIov() reads the ui, before and after the IOV rewrite", {
    .isKeepIovEnv()
    expect_false(.uiHasIov(rxode2::assertRxUi(.isNoIov)))
    .ui <- rxode2::assertRxUi(.isIov)
    expect_true(.uiHasIov(.ui))
    .rw <- .uiApplyIov(.ui, "focei", .isData, foceiControl())$ui
    expect_true(length(.getSyntheticIovEtaNames(.rw)) > 0L)
    expect_true(.uiHasIov(.rw))
  })

  test_that("after an IOV fit the probes still accept a model without IOV, and refuse one with it", {
    skip_on_cran()
    .isKeepIovEnv()
    suppressMessages(suppressWarnings(nlmixr2(
      .isIov,
      .isData,
      "focei",
      foceiControl(print = 0L, maxOuterIterations = 0L, covMethod = "", calcTables = FALSE)
    )))
    # what the probes used to read: the IOV fit's rewrite, still in place
    expect_identical(.uiIovEnv$iovVars, "iov.cl")
    .ui <- rxode2::assertRxUi(.isNoIov)
    expect_true(.vaeGradInScope(.ui))
    expect_false(is.null(.foceiOuterDirs(.ui, "vae")))
    .ll <- rxode2::assertRxUi(.isLl)
    expect_true(.foceiLLGradInScope(.ll, "vae"))
    expect_true(.vaeGradInScope(.ll))
    # and with no IOV fit behind it, a model with IOV is still out of scope
    .uiIovEnv$iovVars <- NULL
    expect_null(.foceiOuterDirs(rxode2::assertRxUi(.isIov), "vae"))
    expect_false(.vaeGradInScope(rxode2::assertRxUi(.isIov)))
  })
})
