# The tables of a finished fit with inter-occasion variability.  Its thetas,
# etas and data are the user's model's, so addNpde() and addTable() refit it
# with no iterations at its etas; their columns must be the fit's own.

# The fits shared by the tests below, each made once per run
.iovTablesCache <- new.env(parent = emptyenv())

# The cached `.iovFit(...)` named `key`
.iovTablesFit <- function(key, ...) {
  if (!exists(key, envir = .iovTablesCache, inherits = FALSE)) {
    assign(key, .iovFit(...), envir = .iovTablesCache)
  }
  get(key, envir = .iovTablesCache, inherits = FALSE)
}

# The columns _nlmixr2est_npdeCalc() (src/npde.cpp) returns
.iovTablesNpdeCols <- c("EPRED", "ERES", "NPDE", "NPD", "PDE", "PD")

nmTest({
  test_that("addNpde() on a fit with IOV gives the NPDE the fit computes itself", {
    for (.est in c("focei", "saem")) {
      .with <- .iovTablesFit(paste0(.est, "-npde"), .est, table = tableControl(npde = TRUE))
      .fit <- .iovTablesFit(.est, .est)
      # a fit asked for NPDE is finalized as one that is not: npde's vpcSim()
      # in its table step runs a simulation's hooks, which clear the IOV
      # rewrite's state unless they put it back
      expect_identical(setdiff(names(.with), names(.fit)), .iovTablesNpdeCols, label = .est)
      expect_identical(setdiff(names(.fit), names(.with)), character(0), label = .est)
      expect_identical(.with$ui$iniDf$name, .fit$ui$iniDf$name, label = .est)
      expect_identical(rownames(.with$parFixedDf), rownames(.fit$parFixedDf), label = .est)
      expect_identical(names(.with$ranef), c("ID", "eta.v"), label = .est)
      expect_equal(.with$fixef, .fit$fixef, tolerance = 1e-10, label = .est)
      .new <- suppressMessages(addNpde(.fit, updateObject = FALSE))
      # its time is recorded as for a fit without IOV, though the refit's
      # nlmixr2() resets the timing globals
      expect_true("NPDE" %in% names(.fit$time), label = .est)
      expect_identical(setdiff(names(.new), names(.fit)), .iovTablesNpdeCols, label = .est)
      expect_identical(nrow(.new), nrow(.fit), label = .est)
      for (.c in .iovTablesNpdeCols) {
        expect_equal(.new[[.c]], .with[[.c]], tolerance = 1e-10, label = paste(.est, .c))
      }
      expect_true(all(is.finite(.new$NPDE)), label = .est)
    }
  })

  test_that("addNpde() on an IOV fit with a table asking for doses keeps the fit's rows", {
    # the dosing rows change the simulation's draws, as in the fit's own NPDE
    .with <- .iovTablesFit(
      "focei-npde-dosing",
      "focei",
      table = tableControl(npde = TRUE, addDosing = TRUE)
    )
    .fit <- .iovTablesFit("focei", "focei")
    .new <- suppressMessages(addNpde(.fit, updateObject = FALSE, table = tableControl(addDosing = TRUE)))
    expect_identical(nrow(.new), nrow(.fit))
    expect_equal(.new$NPDE, .with$NPDE[.with$EVID == 0], tolerance = 1e-10)
  })

  test_that("addNpde() on an IOV fit works after simulate()", {
    skip_if(
      .rxode2WritesNestingData(),
      "rxode2 5.1.8's rxSolve() turns the caller's occasion column into a factor"
    )
    .fit <- .iovTablesFit("focei", "focei")
    invisible(suppressWarnings(suppressMessages(simulate(.fit))))
    # simulate() leaves the fit's data as it was, so the refit still sees a
    # numeric occasion
    expect_identical(.fit$origData$occ, .iovFitData()$occ)
    .new <- suppressMessages(addNpde(.fit, updateObject = FALSE))
    expect_true(all(is.finite(.new$NPDE)))
  })

  test_that("addTable() on a fit with IOV gives the fit's own tables", {
    # a fit run without tables gets the same ones
    for (.est in c("focei", "saem")) {
      .fit <- .iovTablesFit(.est, .est)
      .num <- names(Filter(is.numeric, as.data.frame(.fit)))
      expect_true(all(c("IPRED", "PRED", "IWRES", "iov.cl", "eta.v") %in% .num), label = .est)
      .noTab <- .iovTablesFit(paste0(.est, "-notab"), .est, control = list(calcTables = FALSE))
      for (.f in list(.fit, .noTab)) {
        .t <- suppressMessages(addTable(.f, updateObject = FALSE))
        expect_identical(names(.t), names(.fit), label = .est)
        expect_identical(class(.t)[1:3], class(.fit)[1:3], label = .est)
        expect_identical(.t$ID, .fit$ID, label = .est)
        for (.c in .num) {
          expect_equal(.t[[.c]], .fit[[.c]], tolerance = 1e-10, label = paste(.est, .c))
        }
        expect_equal(.t$objf, .fit$objf, label = .est)
      }
    }
  })

  test_that("an IOV fit with a subject without observations is finished, and addTable() gives its table", {
    # the fit drops subject 3 and its table reinserts it, through a
    # population predict(), whose hooks cleared the IOV rewrite's state; and
    # addTable()'s refit already has the subject labels and reinserted rows
    for (.est in c("focei", "saem")) {
      .full <- .iovTablesFit(.est, .est)
      .fit <- .iovTablesFit(paste0(.est, "-noobs"), .est, data = .iovFitNoObsData())
      expect_identical(names(.fit), names(.full), label = .est)
      expect_identical(.fit$ui$iniDf$name, .full$ui$iniDf$name, label = .est)
      expect_identical(rownames(.fit$parFixedDf), rownames(.full$parFixedDf), label = .est)
      expect_identical(names(.fit$ranef), c("ID", "eta.v"), label = .est)
      expect_identical(nrow(.fit), 132L, label = .est)
      expect_identical(.fit$ID, .full$ID, label = .est)
      expect_true(all(is.na(.fit$IPRED[.fit$ID == "3"])), label = .est)
      expect_true(all(is.finite(.fit$IPRED[.fit$ID != "3"])), label = .est)
      .num <- names(Filter(is.numeric, as.data.frame(.fit)))
      .noTab <- .iovTablesFit(
        paste0(.est, "-noobs-notab"),
        .est,
        control = list(calcTables = FALSE),
        data = .iovFitNoObsData()
      )
      for (.f in list(.fit, .noTab)) {
        .t <- suppressMessages(addTable(.f, updateObject = FALSE))
        expect_identical(names(.t), names(.fit), label = .est)
        expect_identical(.t$ID, .fit$ID, label = .est)
        for (.c in .num) {
          expect_equal(.t[[.c]], .fit[[.c]], tolerance = 1e-10, label = paste(.est, .c))
        }
      }
    }
  })

  test_that("saem's two-level IOV: addNpde() and addTable() give the fit's columns", {
    .ctl <- list(iovMethod = "twoLevel")
    .fit <- .iovTablesFit("saem-2l", "saem", control = .ctl, model = .iovFitMuMod)
    .with <- .iovTablesFit(
      "saem-2l-npde",
      "saem",
      control = .ctl,
      model = .iovFitMuMod,
      table = tableControl(npde = TRUE)
    )
    expect_identical(.fit$env$iovNative, "twoLevel")
    expect_identical(.with$env$iovNative, "twoLevel")
    expect_identical(setdiff(names(.with), names(.fit)), .iovTablesNpdeCols)
    .new <- suppressMessages(addNpde(.fit, updateObject = FALSE))
    for (.c in .iovTablesNpdeCols) {
      expect_equal(.new[[.c]], .with[[.c]], tolerance = 1e-10, label = .c)
    }
    .t <- suppressMessages(addTable(.fit, updateObject = FALSE))
    expect_identical(names(.t), names(.fit))
    for (.c in names(Filter(is.numeric, as.data.frame(.fit)))) {
      expect_equal(.t[[.c]], .fit[[.c]], tolerance = 1e-10, label = .c)
    }
  })

  test_that("saem's collapsed IOV: addNpde() and addTable() refuse to refit", {
    # its tables are those of the collapsed model it fitted, which the focei
    # refit does not rebuild
    .fit <- .iovTablesFit("saem-coll", "saem", control = list(iovMethod = "collapsed"), model = .iovFitMuMod)
    expect_identical(.fit$env$iovNative, "collapsed")
    .why <- paste0(
      " can't recompute the tables of a saem fit with iovMethod = \"collapsed\"; ",
      "ask for them when fitting (`table = tableControl()`, `calcTables = TRUE`)"
    )
    .e <- expect_error(suppressMessages(addNpde(.fit, updateObject = FALSE)))
    expect_identical(conditionMessage(.e), paste0("addNpde()", .why))
    .e <- expect_error(suppressMessages(addTable(.fit, updateObject = FALSE)))
    expect_identical(conditionMessage(.e), paste0("addTable()", .why))
  })

  test_that("addTable() on an IOV fit started from an etaMat refits at the fit's etas", {
    # the fit's control keeps the starting etaMat; a refit at it is off by ~4
    .em <- .iovTablesFit("focei", "focei")$etaMat
    .em[] <- 0.4
    .fit <- .iovFit("focei", control = list(etaMat = .em))
    expect_false(is.null(.fit$foceiControl$etaMat))
    .t <- suppressMessages(addTable(.fit, updateObject = FALSE))
    expect_equal(.t$IPRED, .fit$IPRED, tolerance = 1e-10)
  })

  test_that("addTable() and addNpde() on an ODE IOV fit agree with the fit to the solver's tolerance", {
    # the refit solves the ODEs again; focei solves them to rtol = 1e-3
    # (sigdig = 3), and the two solves agree to about 2e-6
    .fit <- .iovFit("focei", model = .iovFitOdeMod)
    .with <- .iovFit("focei", model = .iovFitOdeMod, table = tableControl(npde = TRUE))
    .t <- suppressMessages(addTable(.fit, updateObject = FALSE))
    expect_identical(names(.t), names(.fit))
    for (.c in names(Filter(is.numeric, as.data.frame(.fit)))) {
      expect_equal(.t[[.c]], .fit[[.c]], tolerance = 1e-4, label = .c)
    }
    .new <- suppressMessages(addNpde(.fit, updateObject = FALSE))
    for (.c in .iovTablesNpdeCols) {
      expect_equal(.new[[.c]], .with[[.c]], tolerance = 1e-4, label = .c)
    }
  })

  test_that("addTable() on an IOV fit keeps and drops the columns its table asks for", {
    .fit <- .iovTablesFit("focei", "focei")
    expect_false("WT" %in% names(.fit))
    expect_true("IWRES" %in% names(.fit))
    .t <- suppressMessages(addTable(.fit, updateObject = FALSE, table = tableControl(keep = "WT", drop = "IWRES")))
    expect_equal(.t$WT, .iovFitData()$WT[.iovFitData()$EVID == 0])
    expect_false("IWRES" %in% names(.t))
  })
})
