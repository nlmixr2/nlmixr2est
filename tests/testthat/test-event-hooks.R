nmTest({
  .hasBus <- exists("rxEventEmit", envir = asNamespace("rxode2"), inherits = FALSE)
  .rec <- new.env()
  .start <- function(env = parent.frame()) {
    .rec$ev <- list()
    rxode2::rxEventListen("nlmixr2est-test", function(event, ...) {
      .p <- list(...)
      .rec$ev[[length(.rec$ev) + 1L]] <- list(event = event, p = .p)
    })
    withr::defer(rxode2::rxEventUnlisten("nlmixr2est-test"), envir = env)
  }
  .events <- function() vapply(.rec$ev, function(e) e$event, character(1))
  .reset <- function() .rec$ev <- list()

  one.cmt <- function() {
    ini({
      tka <- log(1.57)
      tcl <- log(2.72)
      tv <- log(31.5)
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

  test_that("nlmixrUpdateObject only rebinds a single bound name", {
    fit <- .nlmixr2est_fit_posthoc <- suppressMessages(nlmixr2(one.cmt, nlmixr2data::theo_sd, est = "posthoc"))
    e <- new.env()
    e$f <- fit
    clone <- nlmixrClone(fit)
    expect_true(nlmixrUpdateObject(clone, "f", e, fit$env))
    expect_identical(e$f$env, clone$env)
    expect_false(nlmixrUpdateObject(clone, c("[[", "fits", "1"), e))
    expect_false(nlmixrUpdateObject(clone, "notBound", e))
    expect_false(nlmixrUpdateObject(clone, NULL, e))
  })

  test_that("addCwres/addTable on non-name expressions do not error", {
    fit <- suppressMessages(nlmixr2(one.cmt, nlmixr2data::theo_sd, est = "posthoc"))
    fits <- list(fit)
    env <- new.env()
    env$fit <- fit
    expect_error(suppressMessages(addCwres(fits[[1]])), NA)
    expect_error(suppressMessages(addCwres(env$fit)), NA)
    expect_error(suppressMessages(addTable(fits[[1]], updateObject = TRUE)), NA)
    expect_error(suppressMessages(addTable(env$fit, updateObject = TRUE)), NA)
  })

  test_that("a fit emits exactly one fitComplete and no solveComplete", {
    skip_if_not(.hasBus, "rxode2 has no event bus")
    .start()
    fit <- suppressMessages(nlmixr2(one.cmt, nlmixr2data::theo_sd, est = "focei",
                                    control = foceiControl(print = 0)))
    expect_identical(.events(), "fitComplete")
    expect_identical(.rec$ev[[1]]$p$objName, "one.cmt")
    expect_true(inherits(.rec$ev[[1]]$p$fit, "nlmixr2FitCore"))
    expect_true(is.finite(sum(unlist(.rec$ev[[1]]$p$fit$time))))
    .reset()
    nlmixr2(one.cmt)
    foceiControl()
    saemControl()
    expect_length(.rec$ev, 0L)
    .reset()
    fit2 <- suppressMessages(nlmixr2(fit, est = "posthoc"))
    expect_identical(.events(), "fitComplete")
    expect_true(inherits(.rec$ev[[1]]$p$object, "nlmixr2FitCore"))
    .reset()
    expect_error(suppressMessages(nlmixr2(one.cmt, nlmixr2data::theo_sd, est = "notAnEst")))
    expect_length(.rec$ev, 0L)
    expect_equal(rxode2::rxEventDepth(), 0L)
  })

  test_that("simulate/predict/rxSolve on a fit emit one solveComplete with the fit", {
    skip_if_not(.hasBus, "rxode2 has no event bus")
    fit <- suppressMessages(nlmixr2(one.cmt, nlmixr2data::theo_sd, est = "posthoc"))
    .start()
    suppressMessages(simulate(fit))
    suppressMessages(predict(fit, nlmixr2data::theo_sd))
    suppressMessages(predict(fit, nlmixr2data::theo_sd, level = "individual"))
    expect_identical(.events(), rep("solveComplete", 3))
    for (.e in .rec$ev) expect_true(inherits(.e$p$object, "nlmixr2FitCore"))
    for (.e in .rec$ev) expect_lt(nchar(paste(deparse(.e$p$call), collapse = "")), 200)
  })

  test_that("updates emit one fitUpdate with the right inPlace; no-ops emit nothing", {
    skip_if_not(.hasBus, "rxode2 has no event bus")
    fit <- suppressMessages(nlmixr2(one.cmt, nlmixr2data::theo_sd, est = "saem",
                                    control = saemControl(print = 0)))
    expect_false("CWRES" %in% names(fit))
    .start()
    fit <- suppressMessages(addCwres(fit))
    expect_identical(.events(), "fitUpdate")
    expect_true(.rec$ev[[1]]$p$inPlace)
    expect_identical(.rec$ev[[1]]$p$what, "cwres")
    .reset()
    suppressMessages(addCwres(fit))
    expect_length(.rec$ev, 0L)
    fits <- list(fit)
    suppressMessages(suppressWarnings(addNpde(fits[[1]])))
    expect_identical(.events(), "fitUpdate")
    expect_false(.rec$ev[[1]]$p$inPlace)
    expect_null(.rec$ev[[1]]$p$name)
    .reset()
    setOfv(fit, "focei")
    expect_identical(.events(), "fitUpdate")
    expect_true(.rec$ev[[1]]$p$inPlace)
  })

  test_that("vpcSim and augPred emit their own solveComplete, nothing nested", {
    skip_if_not(.hasBus, "rxode2 has no event bus")
    fit <- suppressMessages(nlmixr2(one.cmt, nlmixr2data::theo_sd, est = "posthoc"))
    .start()
    suppressMessages(vpcSim(fit, n = 3))
    suppressMessages(augPred(fit))
    expect_identical(.events(), c("solveComplete", "solveComplete"))
    expect_identical(vapply(.rec$ev, function(e) e$p$kind, ""), c("vpcSim", "augPred"))
  })

  test_that("reading a deferred saem objective emits one fitUpdate", {
    skip_if_not(.hasBus, "rxode2 has no event bus")
    fit <- suppressMessages(nlmixr2(one.cmt, nlmixr2data::theo_sd, est = "saem",
                                    control = saemControl(print = 0)))
    skip_if_not(is.na(get("objective", fit$env)), "objective not deferred")
    .start()
    suppressMessages(fit$objf)
    expect_identical(.events(), "fitUpdate")
    expect_identical(.rec$ev[[1]]$p$what, "ofv")
    .reset()
    fit$objf
    expect_length(.rec$ev, 0L)
  })
})
