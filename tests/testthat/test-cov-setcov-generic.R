# setCov() dispatches on the covariance method, so another package can add one;
# setCov() itself owns the install (PD guard, covList stash, SE refresh).

nmTest({
  .oneCmt <- function() {
    ini({
      tka <- 0.45
      tcl <- 1.0
      tv <- 3.45
      add.sd <- 0.7
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }

  .fitOnce <- function() {
    suppressWarnings(nlmixr2(
      .oneCmt,
      nlmixr2data::theo_sd,
      est = "focei",
      control = foceiControl(print = 0, covMethod = "r,s", calcTables = FALSE)
    ))
  }

  .register <- function(name, fn) {
    .ns <- asNamespace("nlmixr2est")
    registerS3method("setCov", name, fn, envir = .ns)
    withr::defer(rm(list = paste0("setCov.", name), envir = .ns[[".__S3MethodsTable__."]]), envir = parent.frame())
  }

  test_that("a registered method's matrix is installed, stashed and swappable", {
    .fit <- .fitOnce()
    .cov0 <- .fit$cov
    .method0 <- .fit$covMethod
    .se0 <- .fit$parFixedDf$SE
    .calls <- 0L
    .seenControl <- NULL
    .register("testCov", function(fit, method, control = list(n = 7L), ...) {
      .calls <<- .calls + 1L
      .seenControl <<- control
      expect_true(inherits(fit, "nlmixr2FitCore"))
      fit$cov * 4
    })
    expect_true("testCov" %in% setCovAllMethods())

    suppressMessages(setCov(.fit, "testCov"))
    expect_equal(.calls, 1L)
    expect_equal(.seenControl, list(n = 7L))
    expect_identical(.fit$covMethod, "testCov")
    expect_equal(.fit$cov, .cov0 * 4)
    expect_equal(.fit$parFixedDf$SE, .se0 * 2)
    expect_true(.method0 %in% names(.fit$env$covList))
    expect_equal(.fit$env$covOptions$testCov, list(n = 7L))
    expect_error(setCov(.fit, "testCov", control = list(n = 7L)), "no need to switch")

    # the same options are served from the cache without calling the method
    suppressMessages(setCov(.fit, .method0))
    expect_equal(.fit$cov, .cov0)
    suppressMessages(setCov(.fit, "testCov", control = list(n = 7L)))
    expect_equal(.calls, 1L)
    expect_identical(.fit$covMethod, "testCov")

    # different options recompute, even the installed method
    suppressMessages(setCov(.fit, "testCov", control = list(n = 8L)))
    expect_equal(.calls, 2L)
    expect_equal(.seenControl, list(n = 8L))
    expect_equal(.fit$env$covOptions$testCov, list(n = 8L))
    suppressMessages(setCov(.fit, .method0))
    suppressMessages(setCov(.fit, "testCov"))
    expect_equal(.calls, 3L)
  })

  test_that("a rejected matrix leaves the covariance unchanged", {
    .fit <- .fitOnce()
    .cov0 <- .fit$cov
    .register("testCovBad", function(fit, method, control = NULL, ...) {
      -fit$cov
    })
    .register("testCovNoNames", function(fit, method, control = NULL, ...) {
      unname(fit$cov)
    })
    expect_error(setCov(.fit, "testCovBad"), "left unchanged")
    expect_error(setCov(.fit, "testCovNoNames"), "matching dimnames")
    expect_equal(.fit$cov, .cov0)
  })

  test_that("unknown methods and wrong controls error before computing", {
    .fit <- .fitOnce()
    expect_error(setCov(.fit, "nonesuch"), "not supported; can be one of: .*r,s; ")
    expect_error(setCov(.fit, "r,s (full)", control = saControl()), "rsControl")
    expect_error(setCov(.fit, "sa", control = rsControl()), "saControl")
    expect_error(setCov(.fit, "imp", control = saControl()), "impCovControl")
    expect_error(setCov(.fit, "sa (full)"), "not supported")
    .register("testCovNothing", function(fit, method, ...) NULL)
    expect_error(setCov(.fit, "testCovNothing"), "without installing")
    expect_null(.fit$env$covOptions$testCovNothing)
  })

  test_that("covariance controls hold only their own options", {
    expect_equal(unclass(rsControl()), setNames(list(), character(0)))
    expect_equal(unclass(rsControl(hessEps = 1e-4, gillKcov = 3)), list(hessEps = 1e-4, gillKcov = 3L))
    expect_equal(unclass(saControl(nSaCov = 50)), list(nBurn = 100L, nEm = 100L, nSaCov = 50L, seed = 99L))
    expect_equal(unclass(impCovControl(isample = 10)), list(nIter = 1L, isample = 10L, impSeed = 42L))
    expect_error(rsControl(hessEps = -1))
  })

  test_that("rsControl() options reach the finite-difference refit", {
    .fit <- .fitOnce()
    .seen <- NULL
    local_mocked_bindings(.setCov = function(obj, ...) {
      .seen <<- list(...)
      invisible(NULL)
    })
    setCov(.fit, "s (full)", control = rsControl(hessEps = 1e-4, covSmall = 1e-6))
    expect_equal(.seen, list(covMethod = "s", covFull = TRUE, hessEps = 1e-4, covSmall = 1e-6))
    expect_identical(.fit$covMethod, "s (full)")
  })

  test_that("estimation-time FD covariances match the fit's own options", {
    .fit <- .fitOnce()
    expect_identical(.fit$covMethod, "r,s (full)")
    .n <- 0L
    local_mocked_bindings(.setCov = function(obj, ...) {
      .n <<- .n + 1L
      invisible(NULL)
    })
    # the fit's own hessEps is the estimation-time option, so the cache is used
    setCov(.fit, "r,s", control = rsControl(hessEps = .fit$foceiControl$hessEps))
    expect_equal(.n, 0L)
    expect_identical(.fit$covMethod, "r,s")
    expect_error(setCov(.fit, "r,s"), "no need to switch")
    # a different option recomputes the cached shape
    setCov(.fit, "r,s (full)", control = rsControl(hessEps = 1e-4))
    expect_equal(.n, 1L)
    expect_equal(.fit$env$covOptions[["r,s (full)"]]$hessEps, 1e-4)
  })

  test_that("the print separates the other calculated covariances with ';'", {
    .fit <- .fitOnce()
    .out <- capture.output(print(.fit))
    .line <- grep("other calculated covs", .out, value = TRUE)
    expect_length(.line, 1L)
    expect_match(crayon::strip_style(.line), "r,s; ", fixed = TRUE)
  })

  test_that("estimation records its covariance options before settings change", {
    .fit <- .fitOnce()
    expect_true(all(c("r,s (full)", "r,s") %in% names(.fit$env$covOptions)))
    expect_equal(.fit$env$covOptions[["r,s"]]$hessEps, .fit$foceiControl$hessEps)
    .n <- 0L
    local_mocked_bindings(.setCov = function(obj, ...) {
      .n <<- .n + 1L
      invisible(NULL)
    })
    # the cached "r,s" was computed with the original hessEps, not the new one
    .env <- .fit$env
    .fc <- get("foceiControl0", envir = .env)
    .fc$hessEps <- 1e-3
    assign("foceiControl0", .fc, envir = .env)
    expect_equal(.fit$foceiControl$hessEps, 1e-3)
    setCov(.fit, "r,s")
    expect_equal(.n, 1L)
    expect_equal(.fit$env$covOptions[["r,s"]]$hessEps, 1e-3)
  })

  test_that("saControl()/impCovControl() options reach the engine and key the cache", {
    .fit <- .fitOnce()
    .ctl <- list()
    local_mocked_bindings(.covRecomputeNative = function(fit, est, control, useEtaMat = TRUE) {
      .ctl[[length(.ctl) + 1L]] <<- control
      .cov <- get("cov", envir = fit$env)
      list(cov = .cov, covMethod = if (est == "saem") "sa" else "imp", mixRotated = TRUE)
    })
    # a positional control is recorded just as a named one
    suppressMessages(setCov(.fit, "sa", saControl(nSaCov = 50)))
    expect_length(.ctl, 1L)
    expect_equal(.ctl[[1]]$nSaCov, 50L)
    expect_identical(.ctl[[1]]$covMethod, "sa")
    expect_equal(.fit$env$covOptions$sa, unclass(saControl(nSaCov = 50)))
    suppressMessages(setCov(.fit, "imp", control = impCovControl(isample = 20)))
    expect_length(.ctl, 2L)
    expect_equal(.ctl[[2]]$isample, 20L)
    # same options reuse the cache; the default options differ and recompute
    suppressMessages(setCov(.fit, "sa", control = saControl(nSaCov = 50)))
    expect_length(.ctl, 2L)
    expect_identical(.fit$covMethod, "sa")
    suppressMessages(setCov(.fit, "sa"))
    expect_length(.ctl, 3L)
    expect_equal(.ctl[[3]]$nSaCov, 500L)
    suppressMessages(setCov(.fit, "imp", control = impCovControl(isample = 20)))
    expect_length(.ctl, 3L)
    expect_identical(.fit$covMethod, "imp")
    expect_false("imp" %in% names(.fit$env$covList))
  })

  test_that("a method without a control ignores one and keeps its cache", {
    .fit <- .fitOnce()
    .n <- 0L
    local_mocked_bindings(.setCovAnalytic = function(fit, env, method) {
      .n <<- .n + 1L
      .covInstallResult(env, list(cov = get("cov", envir = env), covMethod = method, mixRotated = TRUE))
      invisible(TRUE)
    })
    suppressMessages(setCov(.fit, "analytic", control = rsControl(hessEps = 1e-3)))
    expect_equal(.n, 1L)
    expect_equal(.fit$env$covOptions$analytic, list())
    expect_error(setCov(.fit, "analytic", control = rsControl()), "no need to switch")
    suppressMessages(setCov(.fit, "r,s (full)"))
    suppressMessages(setCov(.fit, "analytic"))
    expect_equal(.n, 1L)
    expect_identical(.fit$covMethod, "analytic")
  })

  test_that("a self-installing method does not leave a stale cached copy", {
    .fit <- .fitOnce()
    .register("testCovSelf", function(fit, method, control = list(k = 1L), ...) {
      .env <- fit$env
      .covInstallResult(
        .env,
        list(cov = get("cov", envir = .env) * control$k, covMethod = "testCovSelf", mixRotated = TRUE)
      )
      NULL
    })
    suppressMessages(setCov(.fit, "testCovSelf"))
    suppressMessages(setCov(.fit, "r,s (full)"))
    expect_true("testCovSelf" %in% names(.fit$env$covList))
    # the method's own install stashes nothing under its name, so the old copy
    # must be dropped by setCov()
    local_mocked_bindings(.covInstallResult = function(env, r) {
      assign("cov", r$cov, envir = env)
      assign("covMethod", r$covMethod, envir = env)
      TRUE
    })
    suppressMessages(setCov(.fit, "testCovSelf", control = list(k = 4L)))
    expect_false("testCovSelf" %in% names(.fit$env$covList))
    expect_equal(.fit$env$covOptions$testCovSelf, list(k = 4L))
  })

  test_that("a method's mixRotated attribute is honored and not installed", {
    .fit <- .fitOnce()
    .seen <- list()
    local_mocked_bindings(.covInstallResult = function(env, r) {
      .seen[[length(.seen) + 1L]] <<- r
      FALSE
    })
    .register("testCovRot", function(fit, method, ...) {
      structure(get("cov", envir = fit$env), mixRotated = TRUE)
    })
    .register("testCovRaw", function(fit, method, ...) {
      get("cov", envir = fit$env)
    })
    expect_error(setCov(.fit, "testCovRot"), "left unchanged")
    expect_error(setCov(.fit, "testCovRaw"), "left unchanged")
    expect_true(.seen[[1]]$mixRotated)
    expect_null(attr(.seen[[1]]$cov, "mixRotated"))
    expect_false(.seen[[2]]$mixRotated)
    expect_null(.fit$env$covOptions$testCovRot)
  })

  test_that("a scoped name dispatches to its base method with a control", {
    .fit <- .fitOnce()
    suppressWarnings(suppressMessages(
      setCov(.fit, "s (full)", control = rsControl(hessEps = 1e-4))
    ))
    expect_identical(.fit$covMethod, "s (full)")
    .se <- .fit$parFixedDf$SE
    names(.se) <- rownames(.fit$parFixedDf)
    .d <- sqrt(diag(.fit$cov))
    .n <- intersect(names(.se), names(.d))
    expect_true(length(.n) > 0L)
    expect_equal(unname(.se[.n]), unname(.d[.n]))
  })

  test_that("setCovOptions() puts fit state in a method's cache key", {
    .fit <- .fitOnce()
    .method0 <- .fit$covMethod
    .calls <- 0L
    .register("testCovSeeded", function(fit, method, control = structure(list(k = 2L), class = "testSeedCtl"), ...) {
      .calls <<- .calls + 1L
      fit$cov * control$k
    })
    .ns <- asNamespace("nlmixr2est")
    registerS3method(
      "setCovOptions",
      "testSeedCtl",
      function(control, fit, ...) {
        .cm <- fit$covMethod
        # the seed is the installed covariance, or the one recorded when it is ours
        if (identical(.cm, "testCovSeeded")) {
          .cm <- fit$env$covOptions$testCovSeeded$seed
        }
        c(unclass(control), list(seed = .cm))
      },
      envir = .ns
    )
    withr::defer(rm(list = "setCovOptions.testSeedCtl", envir = .ns[[".__S3MethodsTable__."]]))
    suppressMessages(setCov(.fit, "testCovSeeded"))
    expect_equal(.calls, 1L)
    expect_equal(.fit$env$covOptions$testCovSeeded, list(k = 2L, seed = .method0))
    expect_error(setCov(.fit, "testCovSeeded"), "no need to switch")
    # the same seed reuses the cache
    suppressMessages(setCov(.fit, .method0))
    suppressMessages(setCov(.fit, "testCovSeeded"))
    expect_equal(.calls, 1L)
    # a different seed recomputes
    suppressMessages(setCov(.fit, "r,s"))
    suppressMessages(setCov(.fit, "testCovSeeded"))
    expect_equal(.calls, 2L)
    expect_equal(.fit$env$covOptions$testCovSeeded$seed, "r,s")
  })

  test_that("setCovOptions() defaults to the control as a plain list", {
    expect_equal(setCovOptions(saControl(), NULL), unclass(saControl()))
    expect_equal(setCovOptions(list(a = 1), NULL), list(a = 1))
    expect_equal(.covOptionsResolve(NULL, NULL), list())
  })

  test_that("setCov(fit) <- matrix installs, records options and is swappable", {
    .fit <- .fitOnce()
    .cov0 <- .fit$cov
    .method0 <- .fit$covMethod
    .se0 <- .fit$parFixedDf$SE
    setCov(.fit) <- .cov0 * 4
    expect_identical(.fit$covMethod, "user")
    expect_equal(.fit$cov, .cov0 * 4)
    expect_equal(.fit$parFixedDf$SE, .se0 * 2)
    expect_equal(.fit$env$covOptions$user, list())
    expect_true(.method0 %in% names(.fit$env$covList))
    setCov(.fit, "doubled") <- .cov0 * 4
    expect_identical(.fit$covMethod, "doubled")
    expect_true("user" %in% names(.fit$env$covList))
    suppressMessages(setCov(.fit, .method0))
    expect_equal(.fit$cov, .cov0)
    expect_error(setCov(.fit) <- -.cov0, "left unchanged")
    expect_equal(.fit$cov, .cov0)
    expect_error(setCov(.fit) <- "a", "cannot install")
  })

  test_that("setCovValue() methods carry their own name and options", {
    .fit <- .fitOnce()
    .ns <- asNamespace("nlmixr2est")
    registerS3method(
      "setCovValue",
      "testCovResult",
      function(value, fit, method = NULL, ...) {
        list(
          cov = value$cov,
          method = if (is.null(method)) "testRes" else method,
          options = list(n = value$n),
          extra = list(testResObj = value)
        )
      },
      envir = .ns
    )
    withr::defer(rm(list = "setCovValue.testCovResult", envir = .ns[[".__S3MethodsTable__."]]))
    .register("testRes", function(fit, method, control = list(n = 3L), ...) {
      stop("should come from the installed result")
    })
    setCov(.fit) <- structure(list(cov = .fit$cov * 9, n = 3L), class = "testCovResult")
    expect_identical(.fit$covMethod, "testRes")
    expect_equal(.fit$env$covOptions$testRes, list(n = 3L))
    expect_s3_class(.fit$env$testResObj, "testCovResult")
    # a result that fails to install stores nothing
    rm("testResObj", envir = .fit$env)
    expect_error(setCov(.fit) <- structure(list(cov = -.fit$cov, n = 3L), class = "testCovResult"), "left unchanged")
    expect_false(exists("testResObj", envir = .fit$env, inherits = FALSE))
    # the installed result counts as the method's default computation
    expect_error(setCov(.fit, "testRes"), "no need to switch")
  })
})
