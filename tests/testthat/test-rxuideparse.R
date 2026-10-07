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
  })

  test_that("nlmControl()", {
    expect_equal(rxUiDeparse.nlmControl(nlmControl(), "var"), quote(var <- nlmControl()))
    expect_equal(rxUiDeparse.nlmControl(nlmControl(covMethod = "r"), "var"), quote(var <- nlmControl(covMethod = "r")))
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
