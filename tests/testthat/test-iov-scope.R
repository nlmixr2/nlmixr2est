# The analytic-gradient scope probes (.foceiOuterDirs(), .foceiLLGradInScope(),
# .foceiAnalyticGradSetup()) take inter-occasion variability from the ui they are
# given (.uiHasIov()), before or after the IOV rewrite, never from .uiIovEnv,
# which holds the last fit's rewrite.

# A one-compartment model with the given ini() eta lines and model() parameter
# lines; .isEtaModel("", ...) has no eta
.isEtaModel <- function(etas, pars) {
  eval(parse(
    text = sprintf(
      "function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        %s
        add.sd <- 0.7
      })
      model({
        %s
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }",
      etas,
      pars
    )
  ))
}

.isNoIov <- .isEtaModel(
  "eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1",
  "ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)"
)
# IOV at two levels, occasion and visit
.isIov <- .isEtaModel(
  "eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1; iov.cl ~ 0.08 | occ; iov.v ~ 0.05 | visit",
  "ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl + iov.cl); v <- exp(tv + eta.v + iov.v)"
)
# etas whose names look like the rewrite's: a user's, and SAEM's pseudo-eta
.isRxCl <- .isEtaModel(
  "eta.ka ~ 0.6; rx.cl ~ 0.3; eta.v ~ 0.1",
  "ka <- exp(tka + eta.ka); cl <- exp(tcl + rx.cl); v <- exp(tv + eta.v)"
)
.isRxEta <- .isEtaModel(
  "eta.ka ~ 0.6; rx.eta.tcl ~ 0.3; eta.v ~ 0.1",
  "ka <- exp(tka + eta.ka); cl <- exp(tcl + rx.eta.tcl); v <- exp(tv + eta.v)"
)
# a repeated (same()) block at the id level, which is not IOV
.isSame <- .isEtaModel(
  "eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ same()",
  "ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)"
)
.isNoEta <- .isEtaModel("", "ka <- exp(tka); cl <- exp(tcl); v <- exp(tv)")

# an ll() endpoint, in scope through the log-density route, with or without IOV
.isLlModel <- function(iov) {
  eval(parse(
    text = sprintf(
      "function() {
      ini({
        lka <- 0.4
        lcl <- -0.9
        lb <- log(3)
        lsd <- log(1.2)
        eta.b ~ 0.1
        %s
      })
      model({
        ka <- exp(lka)
        cl <- exp(lcl)
        d/dt(depot) <- -ka * depot
        d/dt(central) <- ka * depot - cl * central
        mu <- exp(lb + eta.b %s) * central
        sd <- exp(lsd)
        ll(cp) ~ -0.5 * log(2 * pi) - log(sd) - 0.5 * ((DV - mu) / sd)^2
      })
    }",
      if (iov) "iov.b ~ 0.05 | occ" else "",
      if (iov) "+ iov.b" else ""
    )
  ))
}

# theo_sd with two occasions and two visits per subject
.isData <- nlmixr2data::theo_sd
.isData$occ <- ifelse(.isData$TIME < 6, 1L, 2L)
.isData$visit <- ifelse(.isData$TIME < 3, 1L, 2L)

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
  .ui <- rxode2::assertRxUi(.isIov)
  expect_identical(.uiIovNames(.ui), c("iov.cl", "iov.v"))
  expect_true(.uiHasIov(.ui))
  # both rewrites put the occasion etas at the id level, and give the magnitude
  # thetas an nlmixr2iov back-transformation
  for (.m in c("theta", "omega")) {
    .rw <- .uiApplyIov(.ui, "focei", .isData, foceiControl(iovMethod = .m))$ui
    expect_identical(
      .getSyntheticIovEtaNames(.rw),
      c("rx.iov.cl.1", "rx.iov.cl.2", "rx.iov.v.1", "rx.iov.v.2"),
      label = .m
    )
    expect_identical(.uiIovNames(.rw), c("iov.cl", "iov.v"), label = .m)
    expect_true(.uiHasIov(.rw), label = .m)
  }
  # no IOV: no occasion level, whatever the etas are called; a same() block at
  # the id level repeats an eta and is not an occasion
  for (.mod in list(.isNoIov, .isRxCl, .isRxEta, .isSame, .isNoEta, .isLlModel(FALSE))) {
    .u <- rxode2::assertRxUi(.mod)
    expect_identical(.uiIovNames(.u), character(0))
    expect_false(.uiHasIov(.u))
  }
  expect_identical(.uiIovNames(rxode2::assertRxUi(.isLlModel(TRUE))), "iov.b")
})

test_that(".uiIovNames() names the variables of a correlated occasion block, not its pair", {
  .ui <- rxode2::assertRxUi(.isEtaModel(
    "eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1; iov.cl + iov.v ~ c(0.08, 0.01, 0.05) | occ",
    "ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl + iov.cl); v <- exp(tv + eta.v + iov.v)"
  ))
  expect_identical(.uiIovNames(.ui), c("iov.cl", "iov.v"))
})

test_that("the after-fit analytic covariance declines IOV and a repeated block by their own reasons", {
  .isKeepIovEnv()
  .msg <- function(ui) testthat::capture_messages(.foceiCovAnalyticCalc(list(finalUi = ui)))
  .iov <- "covType=\"analytic\": inter-occasion variability (IOV) is out of analytic-covariance scope; using the finite-difference covariance instead\n"
  .ui <- suppressMessages(rxode2::assertRxUi(.isIov))
  .rw <- suppressMessages(.uiApplyIov(.ui, "focei", .isData, foceiControl())$ui)
  .same <- suppressMessages(rxode2::assertRxUi(.isSame))
  expect_identical(.msg(.ui), .iov)
  expect_identical(.msg(.rw), .iov)
  expect_identical(
    .msg(.same),
    "covType=\"analytic\": a repeated (SAME) omega block is out of analytic-covariance scope; using the finite-difference covariance instead\n"
  )
})

test_that(".iovOccRows() takes the level from the condition's base", {
  .idf <- data.frame(
    name = c("t", "a", "b", "c", "d", "e", "(c,d)"),
    condition = c(NA, "id", "id:same:a", "occ", "occ:same:c", "visit", "occ"),
    err = NA_character_
  )
  expect_identical(.iovOccRows(.idf), c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE))
  # an error-model row is never one
  .idf$err[4] <- "add"
  expect_false(.iovOccRows(.idf)[4])
})

nmTest({
  # what the three probes give a model: its outer directions, the log-density
  # scope, and whether the analytic gradient sets itself up
  .isProbes <- function(mod, caller = "focei") {
    .ui <- rxode2::assertRxUi(mod)
    .idf <- .ui$iniDf
    .th <- setNames(.idf$est[!is.na(.idf$ntheta)], .idf$name[!is.na(.idf$ntheta)])
    list(
      dirs = .foceiOuterDirs(.ui, caller)$dirs,
      ll = .foceiLLGradInScope(.ui, caller),
      setup = !is.null(.foceiAnalyticGradSetup(.ui, .th, .ui$omega, caller = caller))
    )
  }
  .isEtaDirs <- c("ETA_1_", "ETA_2_", "ETA_3_")

  test_that("the probes judge the model they are given, whatever fit came before", {
    .isKeepIovEnv()
    # an IOV fit leaves its rewrite in .uiIovEnv
    suppressMessages(suppressWarnings(nlmixr2(
      .isIov,
      .isData,
      "focei",
      foceiControl(print = 0L, maxOuterIterations = 0L, covMethod = "", calcTables = FALSE)
    )))
    expect_identical(.uiIovEnv$iovVars, c("iov.cl", "iov.v"))
    for (.state in c("after an IOV fit", "with no IOV fit")) {
      if (.state == "with no IOV fit") {
        .uiIovEnv$iovVars <- NULL
      }
      # no IOV: every probe accepts the model, with its three eta directions
      for (.mod in list(.isNoIov, .isRxCl, .isSame)) {
        expect_identical(.isProbes(.mod), list(dirs = .isEtaDirs, ll = FALSE, setup = TRUE), label = .state)
      }
      expect_identical(.isProbes(.isNoIov, "vae")$dirs, .isEtaDirs, label = .state)
      expect_true(.foceiLLGradInScope(rxode2::assertRxUi(.isLlModel(FALSE)), "vae"), label = .state)
      # IOV: every probe refuses it, the log-density route included
      expect_identical(.isProbes(.isIov), list(dirs = NULL, ll = FALSE, setup = FALSE), label = .state)
      expect_false(.foceiLLGradInScope(rxode2::assertRxUi(.isLlModel(TRUE)), "vae"), label = .state)
      # no eta: nothing to differentiate
      expect_identical(.isProbes(.isNoEta), list(dirs = NULL, ll = FALSE, setup = FALSE), label = .state)
      # the rewritten ui a fit passes them: its occasion etas are at the id level and
      # mu-referenced, so only the IOV gates refuse it
      for (.m in c("theta", "omega")) {
        .rw <- .uiApplyIov(rxode2::assertRxUi(.isIov), "focei", .isData, foceiControl(iovMethod = .m))$ui
        expect_identical(.isProbes(.rw), list(dirs = NULL, ll = FALSE, setup = FALSE), label = paste(.state, .m))
      }
    }
  })
})
