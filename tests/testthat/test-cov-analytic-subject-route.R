# The per-subject augmented solve, the fallback when the batched one fails, gives
# the batched route's covariance on each assembly: add/prop FOCEI, add + prop
# FOCE, and the general (f,R) assembly under FOCEI, FOCE and foce+, where the
# per-subject route rebuilds the sigma slots from rsig as the batched one does,
# for every sigma (combined1) or a subset of them (boxCox, whose lambda has no
# sigma direction).  The batched solve is mocked to fail, so the fit takes the
# real fallback, and a counter shows the per-subject solve ran.  Both routes are
# evaluated at one fit's estimates (maxOuterIterations = 0).  They difference the
# third-order tensor differently (Richardson against Shi): the largest relative
# difference measured is 2e-6 in the covariance and 1.1e-6 in an SE, while
# dropping the sigma slots' eta-derivatives moves the covariance by 6e-5 to
# 3e-4, so 1e-5 separates the two.
#
# The SD-scale IOV assembly never batches; its per-subject solve is from a caller
# that drops no sigma direction, so it must give the same covariance with Route A
# on as with it off (FOCEI_NO_RSIG).

nmTest({
  .srBase <- function() {
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

  .srAddProp <- function() {
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

  .srC1 <- function() {
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

  # combined1 with both sigmas fixed: its solve carries Rsig, but no sigma
  # direction is dropped
  .srC1Fixed <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- fix(0.5)
      prop.sd <- fix(0.1)
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

  .srBoxCox <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.5
      lambda <- 0.5
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd) + boxCox(lambda)
    })
  }

  .srFit <- function(model, est, data = nlmixr2data::theo_sd, ...) {
    suppressMessages(suppressWarnings(nlmixr2(
      model,
      data,
      est,
      foceiControl(print = 0L, calcTables = FALSE, ...),
      table = tableControl(cwres = FALSE, npde = FALSE)
    )))
  }

  # the covariance at the fit's estimates, counting the per-subject solves; with
  # batch = FALSE the batched solve fails, as it does when it cannot solve
  .srCov <- function(fit, est, batch, data = nlmixr2data::theo_sd, ...) {
    .real <- .foceiAnalyticSolveSubjectFD3
    .n <- new.env(parent = emptyenv())
    .n$calls <- 0L
    local_mocked_bindings(.foceiAnalyticSolveSubjectFD3 = function(...) {
      .n$calls <- .n$calls + 1L
      .real(...)
    })
    if (!batch) {
      local_mocked_bindings(.foceiAnalyticSolveAllFD3 = function(...) NULL)
    }
    .f <- .srFit(fit$finalUi, est, data, covMethod = "analytic", maxOuterIterations = 0L, ...)
    list(fit = .f, calls = .n$calls)
  }

  .srExpectSame <- function(a, b, label) {
    expect_false(exists(".analyticReason", envir = a$env, inherits = FALSE), label = label)
    expect_false(exists(".analyticReason", envir = b$env, inherits = FALSE), label = label)
    expect_identical(a$covMethod, b$covMethod, label = label)
    expect_identical(dimnames(a$cov), dimnames(b$cov), label = label)
    expect_equal(sqrt(diag(a$cov)), sqrt(diag(b$cov)), tolerance = 1e-5, label = label)
    expect_equal(a$cov, b$cov, tolerance = 1e-5, label = label)
    expect_equal(stats::cov2cor(a$cov), stats::cov2cor(b$cov), tolerance = 1e-5, label = label)
  }

  .srCase <- list(
    addFocei = list(model = .srBase, est = "focei", label = "analytic (full)"),
    addPropFoce = list(model = .srAddProp, est = "foce", label = "|r| (analytic)"),
    generalFocei = list(model = .srC1, est = "focei", label = "analytic (full)"),
    generalFoce = list(model = .srC1, est = "foce", label = "analytic (full)"),
    generalFocePlus = list(model = .srC1, est = "foce", ctl = list(foce = "foce+"), label = "analytic (full)"),
    boxCoxFocei = list(model = .srBoxCox, est = "focei", label = "analytic (full)"),
    fixedSigmaFocei = list(model = .srC1Fixed, est = "focei", label = "analytic (full)")
  )

  test_that("the per-subject solve gives the batched solve's analytic covariance", {
    skip_on_cran()
    .nsub <- length(unique(nlmixr2data::theo_sd$ID))
    for (.label in names(.srCase)) {
      .c <- .srCase[[.label]]
      .ctl <- if (is.null(.c$ctl)) list() else .c$ctl
      .fit <- do.call(.srFit, c(list(.c$model, .c$est, covMethod = ""), .ctl))
      .batched <- do.call(.srCov, c(list(.fit, .c$est, batch = TRUE), .ctl))
      .subject <- do.call(.srCov, c(list(.fit, .c$est, batch = FALSE), .ctl))
      expect_identical(.batched$calls, 0L, label = .label)
      expect_identical(.subject$calls, .nsub, label = .label)
      expect_identical(.batched$fit$covMethod, .c$label, label = .label)
      .srExpectSame(.subject$fit, .batched$fit, .label)
    }
  })

  test_that("setCov(fit, \"analytic\") installs the analytic covariance with every sigma fixed", {
    skip_on_cran()
    # the solve carries Rsig and no sigma direction is dropped: the analytic
    # covariance is computed, as the batched route gives it
    .fit <- .srFit(.srC1Fixed, "focei", covMethod = "")
    suppressMessages(suppressWarnings(setCov(.fit, "analytic")))
    expect_identical(.fit$covMethod, "analytic")
    expect_true(all(is.finite(sqrt(diag(.fit$cov)))))
  })

  test_that("the SD-scale IOV covariance is the same with Route A on and off", {
    skip_on_cran()
    .d <- .covIovData()
    .fit <- .srFit(.covIovModel, "focei", .d, covMethod = "", iovXform = "sd")
    .on <- .srCov(.fit, "focei", batch = TRUE, data = .d, covFull = TRUE, iovXform = "sd")
    .off <- local({
      withr::local_envvar(FOCEI_NO_RSIG = "1")
      .srCov(.fit, "focei", batch = TRUE, data = .d, covFull = TRUE, iovXform = "sd")
    })
    # the IOV rescale never batches: every subject is solved on its own
    expect_identical(.on$calls, 24L)
    expect_identical(.on$fit$covMethod, "analytic (full)")
    expect_true("iov.cl" %in% rownames(.on$fit$cov))
    expect_true(any(grepl("^om\\.", rownames(.on$fit$cov))))
    .srExpectSame(.on$fit, .off$fit, "sd iov")
  })
})
