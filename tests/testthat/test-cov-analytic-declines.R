# Every way covMethod="analytic" can decline gives a reason.  Each helper the live
# hook (.foceiCalcRanalytic) depends on is made to fail -- return NULL, noting
# nothing -- one at a time, on each assembly: add/prop FOCEI, add + prop FOCE (the
# frozen-R0 eta = 0 solve), the general (f,R) assembly under FOCEI and FOCE, and
# AGQ.  A decline must record a reason other than the "no reason recorded"
# sentinel, the same in covTried, and every helper must make the covariance decline
# on some assembly.  A failed batched solve is recovered by the per-subject solve on
# the add/prop assembly; on the general (f,R) assembly the per-subject solve does
# not rebuild the sigma slots and errors (pre-existing on origin/main 011375065,
# reported for its own fix), which is pinned here so that a fix updates this test.
# The reason texts themselves are pinned in test-cov-analytic-reason.R.

nmTest({
  .dcModAdd <- function() {
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

  .dcModAddProp <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.5
      prop.sd <- 0.1
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd) + prop(prop.sd)
    })
  }

  .dcModC1 <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.5
      prop.sd <- 0.1
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd) + prop(prop.sd) + combined1()
    })
  }

  .dcFit <- function(model, est, ...) {
    suppressMessages(suppressWarnings(nlmixr2(
      model,
      nlmixr2data::theo_sd,
      est,
      foceiControl(covMethod = "analytic", print = 0L, maxOuterIterations = 0L, calcTables = FALSE, ...),
      table = tableControl(cwres = FALSE, npde = FALSE)
    )))
  }

  .dcNull <- function(...) NULL

  # helpers before the assemblies split, then those of the solves and assemblies
  .dcLeaves <- list(
    ".foceiEtaThetaMap",
    ".foceiAnalyticErrFull",
    ".foceiAnalyticDirections",
    ".foceiOmegaPairs",
    ".omegaVarCovDeriv",
    ".foceiAnalyticAugModelDirs",
    ".nmSubjectObs",
    ".foceiAnalyticSolveAll",
    ".foceiAnalyticRxSolve",
    ".foceiAnalyticSolveAllFD3",
    c(".foceiAnalyticSolveAllFD3", ".foceiAnalyticSolveSubjectFD3"),
    c(".foceiAnalyticSolveAllFD3", ".foceiAnalyticSolveFA"),
    ".foceiAnalyticExpandSigma",
    ".foceiAnalyticSubjectR",
    ".foceiAnalyticSubjectRfoce",
    "foceiRAllFR_",
    "foceiRAllFoceFR_",
    ".foceiAnalyticSgNames",
    ".foceiAnalyticSolve"
  )
  .dcBatch <- ".foceiAnalyticSolveAllFD3"
  .dcBatchError <- "an error (number of items to replace is not a multiple of replacement length)"

  .dcPaths <- list(
    addFocei = list(model = .dcModAdd, est = "focei", leaves = .dcLeaves, recovers = TRUE),
    addPropFoce = list(model = .dcModAddProp, est = "foce", leaves = .dcLeaves, recovers = TRUE),
    generalFocei = list(model = .dcModC1, est = "focei", leaves = .dcLeaves, recovers = FALSE),
    generalFoce = list(model = .dcModC1, est = "foce", leaves = .dcLeaves, recovers = FALSE),
    # AGQ shares the add/prop assembly; only the node solve and assembly are its own
    agq = list(
      model = .dcModAdd,
      est = "focei",
      agq = 3L,
      leaves = list(".foceiAnalyticSolveFA", ".foceiAnalyticSubjectR", ".foceiAnalyticRxSolve"),
      recovers = TRUE
    )
  )

  test_that("every decline of the analytic covariance records a reason", {
    skip_on_cran()
    .declined <- character(0)
    for (.p in names(.dcPaths)) {
      .path <- .dcPaths[[.p]]
      .model <- .path$model
      .args <- if (is.null(.path$agq)) list() else list(nAGQ = .path$agq)
      # the path must give the analytic covariance before any helper fails
      .base <- do.call(.dcFit, c(list(.model, .path$est), .args))
      expect_false(exists(".analyticReason", envir = .base$env, inherits = FALSE), label = .p)
      for (.leaf in .path$leaves) {
        .name <- paste(.leaf, collapse = "+")
        .r <- local({
          .m <- stats::setNames(rep(list(.dcNull), length(.leaf)), .leaf)
          do.call(local_mocked_bindings, .m)
          do.call(.dcFit, c(list(.model, .path$est), .args))
        })
        .label <- paste(.p, .name)
        .gone <- exists(".analyticReason", envir = .r$env, inherits = FALSE)
        if (identical(.name, .dcBatch)) {
          expect_identical(.gone, !.path$recovers, label = .label)
          if (.gone) {
            expect_identical(.r$env$.analyticReason, .dcBatchError, label = .label)
          }
        }
        if (!.gone) {
          next
        }
        .why <- .r$env$.analyticReason
        expect_false(identical(.why, "no reason recorded"), label = .label)
        expect_true(nzchar(.why), label = .label)
        expect_identical(.r$env$covTried$method[1], "analytic", label = .label)
        expect_identical(.r$env$covTried$outcome[1], .why, label = .label)
        expect_false(identical(.r$env$.analyticReasonShort, "no reason recorded"), label = .label)
        expect_lt(nchar(.r$env$.analyticWarnNoFallback), 75L, label = .label)
        .declined <- c(.declined, .name)
      }
    }
    .all <- unique(vapply(
      unlist(lapply(.dcPaths, `[[`, "leaves"), recursive = FALSE),
      paste,
      character(1),
      collapse = "+"
    ))
    expect_setequal(unique(.declined), .all)
  })
})
