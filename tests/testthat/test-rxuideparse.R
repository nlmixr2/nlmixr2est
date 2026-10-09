nmTest({
  test_that("foceiControl() deparse", {
    expect_equal(
      rxUiDeparse.foceiControl(
        foceiControl(
          innerOpt = "lbfgsb3c",
          scaleType = "norm",
          normType = "std",
          derivMethod = "central",
          covDerivMethod = "forward",
          covMethod = "s",
          diagXform = "identity",
          addProp = "combined1",
          eventType = "forward",
          optimHessType = "forward"
        ),
        "ctl"
      ),
      quote(ctl <- foceiControl(derivMethod = "central", covDerivMethod = "forward",
                                           covMethod = "s", diagXform = "identity", optimHessType = "forward",
                                           innerOpt = "lbfgsb3c", scaleType = "norm", normType = "std",
                                           eventType = "forward", addProp = "combined1"))
    )

    expect_equal(
      rxUiDeparse.foceiControl(foceiControl(eventType = "forward"), "ctl"),
      quote(ctl <- foceiControl(eventType = "forward"))
    )

    expect_equal(
      rxUiDeparse.foceiControl(foceiControl(warm = "save"), "ctl"),
      quote(ctl <- foceiControl(warm = "save"))
    )

    expect_equal(rxUiDeparse.foceiControl(foceiControl(), "ctl"), quote(ctl <- foceiControl()))

    expect_warning(rxUiDeparse.foceiControl(foceiControl(outerOpt = optim), "ctl"), "reset")
  })

  test_that("foceiControl() deparse keeps a deferred covariance and a skipped one", {
    for (.m in c("sa", "imp")) {
      .call <- rxUiDeparse.foceiControl(foceiControl(covMethod = .m), "ctl")
      expect_equal(.call, str2lang(paste0("ctl <- foceiControl(covMethod = \"", .m, "\")")))
      expect_identical(eval(.call[[3]])$covMethodDeferred, .m)
    }
    # an analytic control whose covariance step was turned off (slot 0)
    .ctl <- foceiControl(covMethod = "analytic")
    .ctl$covMethod <- 0L
    expect_equal(rxUiDeparse.foceiControl(.ctl, "ctl"), quote(ctl <- foceiControl(covMethod = "")))
  })

  test_that("saemControl() deparse", {
    expect_equal(
      rxUiDeparse.saemControl(saemControl(nBurn = 2, nEm = 2, nmc = 7, nu = c(3, 3, 3)), "ctl"),
      quote(ctl <- saemControl(nBurn = 2, nEm = 2, nmc = 7,
                                          nu = c(3, 3, 3)))
    )
    expect_equal(rxUiDeparse.saemControl(saemControl(), "ctl"), quote(ctl <- saemControl()))
  })

  test_that("bobyqaControl()", {
    expect_equal(rxUiDeparse.bobyqaControl(bobyqaControl(), "var"), quote(var <- bobyqaControl()))

    expect_equal(
      rxUiDeparse.bobyqaControl(bobyqaControl(scaleType = "multAdd"), "var"),
      quote(var <- bobyqaControl(scaleType = "multAdd"))
    )
  })

  test_that("lbfgsb3cControl()", {
    expect_equal(rxUiDeparse.lbfgsb3cControl(lbfgsb3cControl(), "var"), quote(var <- lbfgsb3cControl()))

    expect_equal(
      rxUiDeparse.lbfgsb3cControl(lbfgsb3cControl(normType = "len"), "var"),
      quote(var <- lbfgsb3cControl(normType = "len"))
    )
  })

  test_that("n1qn1Control()", {
    expect_equal(rxUiDeparse.n1qn1Control(n1qn1Control(), "var"), quote(var <- n1qn1Control()))

    expect_equal(
      rxUiDeparse.n1qn1Control(n1qn1Control(covMethod = "n1qn1"), "var"),
      quote(var <- n1qn1Control(covMethod = "n1qn1"))
    )
  })

  test_that("newuoaControl()", {
    expect_equal(rxUiDeparse.newuoaControl(newuoaControl(), "var"), quote(var <- newuoaControl()))
    expect_equal(
      rxUiDeparse.newuoaControl(newuoaControl(addProp = "combined1"), "var"),
      quote(var <- newuoaControl(addProp = "combined1"))
    )
  })

  test_that("nlmeControl()", {
    expect_equal(rxUiDeparse.nlmeControl(nlmeControl(), "var"), quote(var <- nlmeControl()))

    expect_equal(rxUiDeparse.nlmeControl(nlmeControl(opt = "nlm"), "var"), quote(var <- nlmeControl(opt = "nlm")))
  })

  test_that("nlminbControl()", {
    expect_equal(rxUiDeparse.nlminbControl(nlminbControl(), "var"), quote(var <- nlminbControl()))

    expect_equal(
      rxUiDeparse.nlminbControl(nlminbControl(solveType = "grad"), "var"),
      quote(var <- nlminbControl(solveType = "grad"))
    )
    # covMethod's default follows solveType
    expect_equal(
      rxUiDeparse.nlminbControl(nlminbControl(solveType = "grad", covMethod = "nlminb"), "var"),
      quote(var <- nlminbControl(solveType = "grad", covMethod = "nlminb"))
    )
    expect_equal(
      rxUiDeparse.nlminbControl(nlminbControl(covMethod = "r"), "var"),
      quote(var <- nlminbControl(covMethod = "r"))
    )
  })

  test_that("nlmControl()", {
    expect_equal(rxUiDeparse.nlmControl(nlmControl(), "var"), quote(var <- nlmControl()))
    expect_equal(rxUiDeparse.nlmControl(nlmControl(covMethod = "r"), "var"), quote(var <- nlmControl(covMethod = "r")))
    # covMethod's default follows solveType: "fun" defaults to "r"
    expect_equal(rxUiDeparse.nlmControl(nlmControl(solveType = "fun"), "var"), quote(var <- nlmControl(solveType = "fun")))
  })

  test_that("nlsControl()", {
    expect_equal(rxUiDeparse.nlsControl(nlsControl(), "var"), quote(var <- nlsControl()))
    expect_equal(
      rxUiDeparse.nlsControl(nlsControl(algorithm = "port"), "var"),
      quote(var <- nlsControl(algorithm = "port"))
    )
  })

  test_that("optimControl()", {
    expect_equal(rxUiDeparse.optimControl(optimControl(), "var"), quote(var <- optimControl()))
    expect_equal(
      rxUiDeparse.optimControl(optimControl(method = "L-BFGS-B", covMethod = "optim"), "var"),
      quote(var <- optimControl(method = "L-BFGS-B", covMethod="optim"))
    )

    expect_equal(
      rxUiDeparse.optimControl(optimControl(eventType = "forward"), "var"),
      quote(var <- optimControl(eventType = "forward"))
    )
  })

  test_that("uobyqaControl()", {
    expect_equal(rxUiDeparse.uobyqaControl(uobyqaControl(), "var"), quote(var <- uobyqaControl()))
    expect_equal(rxUiDeparse.uobyqaControl(uobyqaControl(scaleTo = 4), "var"), quote(var <- uobyqaControl(scaleTo = 4)))
  })

  test_that("tableControl()", {
    expect_equal(rxUiDeparse.tableControl(tableControl(), "var"), quote(var <- tableControl()))
    expect_equal(
      rxUiDeparse.tableControl(tableControl(censMethod = "epred"), "var"),
      quote(var <- tableControl(censMethod = "epred"))
    )
  })

  test_that("a deparsed control evaluates back to the identical control", {
    .ctors <- c(
      "agqControl",
      "bobyqaControl",
      "emviControl",
      "fbviControl",
      "foceiControl",
      "iagqControl",
      "ilaplaceControl",
      "laplaceControl",
      "lbfgsb3cControl",
      "magqControl",
      "mlaplaceControl",
      "n1qn1Control",
      "newuoaControl",
      "nlmControl",
      "nlmeControl",
      "nlminbControl",
      "nlsControl",
      "optimControl",
      "saemControl",
      "trustControl",
      "uobyqaControl",
      "vaeControl",
      "impmapControl",
      "impControl",
      "npagControl",
      "npbControl"
    )
    for (.c in .ctors) {
      .f <- get(.c)
      for (.a in list(list(), list(sigdig = 4), list(print = 0L), list(sigdig = 5, print = 0L))) {
        .x <- suppressWarnings(do.call(.f, .a))
        .e <- rxode2::rxUiDeparse(.x, "ctl")
        expect_identical(eval(.e[[3]]), .x, info = paste(.c, deparse1(.e)))
      }
    }
  })

  test_that("sigdig and print deparse as themselves", {
    expect_equal(
      rxUiDeparse.foceiControl(foceiControl(sigdig = 4, print = 0L), "ctl"),
      quote(ctl <- foceiControl(print = 0L, sigdig = 4))
    )
    expect_equal(
      rxode2::rxUiDeparse(saemControl(sigdig = 4, nBurn = 5L, nEm = 5L), "ctl"),
      quote(ctl <- saemControl(sigdig = 4, nBurn = 5L, nEm = 5L))
    )
    for (.x in list(
      foceiControl(lbfgsFactr = 1 / 3),
      saemControl(nu = c(2, 2, 2)),
      saemControl(nmc = 3L),
      saemControl(sigdig = 6, sigdigTable = 3),
      saemControl(sigdig = 4, sigdigTable = 3, tol = 1e-5),
      saemControl(trace = 1),
      saemControl(nBurn = c(burn = 200)),
      foceiControl(outerOpt = "bobyqa"),
      foceiControl(fast = TRUE),
      foceiControl(resetEtaP = 0),
      foceiControl(fdChartrand = FALSE),
      foceiControl(maxInnerIterations = 100, n1qn1nsim = 10001),
      optimControl(method = "BFGS", covMethod = "r"),
      impmapControl(proposal = "mixture", propMixScale = c(1, 4, 9), propMixWeight = c(1, 6, 15)),
      foceiControl(rxControl = foceiControl()$rxControl),
      saemControl(rxControl = saemControl()$rxControl),
      impmapControl(ctol = 0.01),
      impControl(isample = 500L, sirSample = 30L),
      impmapControl(gammaRule = "floor", nConvWindow = 20L)
    )) {
      .e <- rxode2::rxUiDeparse(.x, "ctl")
      expect_identical(eval(.e[[3]]), .x, info = deparse1(.e))
    }
  })

  test_that("the imp and np controls deparse through their own constructor", {
    expect_equal(rxode2::rxUiDeparse(impmapControl(), "ctl"), quote(ctl <- impmapControl()))
    expect_equal(rxode2::rxUiDeparse(impControl(isample = 500L), "ctl"), quote(ctl <- impControl(isample = 500L)))
    expect_equal(
      rxode2::rxUiDeparse(npagControl(cycles = 3L, cores = 2L), "ctl"),
      quote(ctl <- npagControl(cycles = 3L, cores = 2L))
    )
    expect_equal(rxode2::rxUiDeparse(npbControl(points = 20L), "ctl"), quote(ctl <- npbControl(points = 20L)))
    ## a nu the fit raised itself (nuAuto) is redone on a refit, not written
    .x <- saemControl(print = 100)
    .x$mcmc$nu <- c(4, 4, 4)
    expect_equal(rxode2::rxUiDeparse(.x, "ctl"), quote(ctl <- saemControl(print = 100L)))
    ## a fit resolves gammaMethod = "auto" and keeps the request
    .x <- impControl()
    .x$gammaMethod <- "global"
    .x$gammaMethodUser <- "auto"
    expect_equal(rxode2::rxUiDeparse(.x, "ctl"), quote(ctl <- impControl()))
    expect_equal(
      rxode2::rxUiDeparse(impControl(covMethod = "r,s"), "ctl"),
      quote(ctl <- impControl(covMethod = "r,s"))
    )
  })

  test_that("every scalar constructor argument survives the deparse", {
    ## each argument set away from its default, one at a time; emviControl's
    ## resume holds a previous fit and is not deparsed
    .alt <- function(d) {
      if (is.null(d)) {
        return(list(2L, 0.5, 1e-5, TRUE))
      }
      if (is.call(d) && identical(d[[1]], as.name("c"))) {
        .v <- eval(d)
        return(if (is.character(.v)) as.list(.v[-1]) else list())
      }
      if (is.logical(d) && length(d) == 1L && !is.na(d)) {
        return(list(!d))
      }
      if (is.integer(d) && length(d) == 1L) {
        return(list(d + 1L, 0L))
      }
      if (is.numeric(d) && length(d) == 1L && is.finite(d)) {
        return(list(d * 2 + 0.1, 0, 1 / 3))
      }
      list()
    }
    for (.c in c(
      "agqControl",
      "bobyqaControl",
      "emviControl",
      "foceiControl",
      "laplaceControl",
      "lbfgsb3cControl",
      "n1qn1Control",
      "newuoaControl",
      "nlmControl",
      "nlmeControl",
      "nlminbControl",
      "nlsControl",
      "optimControl",
      "saemControl",
      "trustControl",
      "uobyqaControl",
      "vaeControl",
      "impmapControl",
      "npagControl",
      "npbControl",
      "tableControl"
    )) {
      .f <- get(.c)
      .fm <- formals(.f)
      for (.a in setdiff(names(.fm), c("...", "rxControl", "gamma", "df", "print", "sigdig", "resume"))) {
        for (.v in .alt(.fm[[.a]])) {
          .x <- tryCatch(
            suppressWarnings(suppressMessages(do.call(.f, stats::setNames(list(.v), .a)))),
            error = function(e) NULL
          )
          if (is.null(.x)) {
            next
          }
          .e <- rxode2::rxUiDeparse(.x, "ctl")
          expect_identical(suppressWarnings(eval(.e[[3]])), .x, info = paste(.c, .a, deparse1(.v)))
        }
      }
    }
  })
})

# Every exported `*Control()` that builds with no arguments and returns a
# classed object, keyed by the object's class.  A control class added to the
# package appears here without being listed by hand.
.deparseControlConstructors <- function() {
  .ns <- asNamespace("nlmixr2est")
  .fn <- sort(grep("Control$", getNamespaceExports(.ns), value = TRUE))
  .ret <- list()
  for (.f in .fn) {
    .ctl <- tryCatch(get(.f, envir = .ns)(), error = function(e) NULL)
    .cls <- class(.ctl)[1]
    if (is.null(.ctl) || !grepl("Control$", .cls)) {
      next
    }
    .ret[[.f]] <- list(cls = .cls, ctl = .ctl)
  }
  .ret
}

# Settings that differ from the defaults, per constructor, so the round trip is
# not only the empty call.  A constructor not named here is checked at its defaults.
.deparseControlSettings <- list(
  foceControl = list(maxOuterIterations = 7L, covMethod = "s"),
  focepControl = list(maxOuterIterations = 7L, covMethod = "s"),
  foControl = list(maxOuterIterations = 7L, posthoc = FALSE),
  foiControl = list(maxOuterIterations = 7L, posthoc = FALSE),
  posthocControl = list(covMethod = "r"),
  mfoceiControl = list(maxOuterIterations = 7L, covMethod = "s"),
  ifoceiControl = list(maxOuterIterations = 7L, covMethod = "s"),
  mfoceControl = list(maxOuterIterations = 7L, covMethod = "s"),
  ifoceControl = list(maxOuterIterations = 7L, covMethod = "s"),
  mfocepControl = list(maxOuterIterations = 7L, covMethod = "s"),
  ifocepControl = list(maxOuterIterations = 7L, covMethod = "s"),
  impmapControl = list(isample = 123L, proposal = "t", df = 5, covMethod = "r", gammaRule = "floor"),
  impControl = list(isample = c(100L, 200L), covMethod = "analytic"),
  qrpemControl = list(isample = 128L),
  npagControl = list(cycles = 7L, cores = 2L, gridBounds = "ini", covMethod = "s"),
  npbControl = list(points = 20L, nsamp = 11L, cores = 2L, covMethod = "s"),
  rsControl = list(hessEps = 1e-4, gillKcov = 3L, covGillF = TRUE, smatNorm = FALSE),
  saControl = list(nBurn = 11L, seed = 5L),
  impCovControl = list(nIter = 3L, isample = 77L),
  iterPrintControl = list(every = 5L, simple = TRUE)
)

test_that("every control class has an rxUiDeparse method that rebuilds it", {
  .all <- .deparseControlConstructors()
  # a class listed here has gone missing, or the scan stopped finding constructors
  expect_gte(length(.all), 40L)
  for (.f in names(.all)) {
    .cls <- .all[[.f]]$cls
    expect_false(is.null(utils::getS3method("rxUiDeparse", .cls, optional = TRUE)), info = paste(.f, .cls))
    .call <- rxUiDeparse(.all[[.f]]$ctl, "ctl")
    expect_identical(.call[[2]], quote(ctl), info = .f)
    expect_identical(eval(.call[[3]]), .all[[.f]]$ctl, info = .f)
    if (!is.null(.deparseControlSettings[[.f]])) {
      .ctl <- do.call(get(.f, envir = asNamespace("nlmixr2est")), .deparseControlSettings[[.f]])
      .call <- rxUiDeparse(.ctl, "ctl")
      expect_identical(eval(.call[[3]]), .ctl, info = paste(.f, "with settings"))
    }
  }
  # every constructor named in the settings table exists
  expect_identical(setdiff(names(.deparseControlSettings), names(.all)), character(0))
})

test_that("every focei-family control deparses as a call that names its settings", {
  for (.cls in c("foceControl", "focepControl", "mfoceiControl", "ifocepControl")) {
    .ctl <- do.call(.cls, list(maxOuterIterations = 7L, covMethod = "s"))
    expect_equal(
      rxUiDeparse(.ctl, "ctl"),
      str2lang(paste0("ctl <- ", .cls, "(maxOuterIterations = 7L, covMethod = \"s\")")),
      info = .cls
    )
  }
  expect_equal(
    rxUiDeparse(saControl(nBurn = 11L), "ctl"),
    quote(ctl <- saControl(nBurn = 11L)),
    ignore_attr = TRUE
  )
  expect_equal(rxUiDeparse(rsControl(), "ctl"), quote(ctl <- rsControl()))
  expect_equal(rxUiDeparse(impmapControl(), "ctl"), quote(ctl <- impmapControl()))
})
