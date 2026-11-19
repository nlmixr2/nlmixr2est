test_that(".simModelCacheKey() only keys off a fit environment", {
  expect_null(.simModelCacheKey(NULL))
  expect_null(.simModelCacheKey(list(ui = "model")))
  .env <- new.env(parent = emptyenv())
  # an environment that holds no model is not a fit environment
  expect_null(.simModelCacheKey(.env))
  assign("ui", "model", envir = .env)
  expect_identical(.simModelCacheKey(.env), "model")
  # an uncompressed (environment) model can be edited in place, so it is not
  # a key the cache can trust
  assign("ui", new.env(parent = emptyenv()), envir = .env)
  expect_null(.simModelCacheKey(.env))
  assign("ui", "model", envir = .env)
  # a rxUi is an environment too, and it must not pick up a `ui` from a parent
  expect_null(.simModelCacheKey(new.env(parent = .env)))
  withr::with_options(list(nlmixr2.simModelCache = FALSE), {
    expect_null(.simModelCacheKey(.env))
  })
})

test_that("the simulation model cache round trips, replaces and evicts", {
  .simModelCacheReset()
  on.exit(.simModelCacheReset())

  expect_null(.simModelCacheGet("a"))
  .simModelCacheSet("a", "modelA")
  expect_identical(.simModelCacheGet("a"), "modelA")
  # keys are compared by value, so an equal key hits the same entry
  expect_identical(.simModelCacheGet(paste0("", "a")), "modelA")

  # re-setting a key replaces its entry rather than adding a second one
  .simModelCacheSet("a", "modelA2")
  expect_identical(.simModelCacheGet("a"), "modelA2")
  expect_equal(length(.simModelCache$entries), 1L)

  # filling the cache drops the least recently used entry ("a")
  for (.i in seq_len(.simModelCacheMax)) {
    .simModelCacheSet(paste0("k", .i), paste0("m", .i))
  }
  expect_equal(length(.simModelCache$entries), .simModelCacheMax)
  expect_null(.simModelCacheGet("a"))
  # reading "k1" makes it the most recently used, so the next insert keeps it
  expect_identical(.simModelCacheGet("k1"), "m1")
  .simModelCacheSet("k99", "m99")
  expect_equal(length(.simModelCache$entries), .simModelCacheMax)
  expect_identical(.simModelCacheGet("k1"), "m1")
  expect_null(.simModelCacheGet("k2"))

  .simModelCacheReset()
  expect_equal(length(.simModelCache$entries), 0L)
  expect_null(.simModelCacheGet("k1"))
})

nmTest({
  test_that("a fit lowers to one shared simulation model (nlmixr2/rxode2#1289)", {
    skip_if(is.null(one.compartment.fit.saem))
    skip_if(is.null(one.compartment.fit.focei))
    .fit <- one.compartment.fit.saem
    .simModelCacheReset()
    on.exit(.simModelCacheReset())

    .mod <- .fit$simulationModel
    expect_true(is.environment(.mod))
    # the same object, not merely an equal one: nothing is lowered twice
    expect_identical(.fit$simulationModel, .mod)
    expect_equal(length(.simModelCache$entries), 1L)

    # solving the fit reuses that same model instead of lowering per call
    .ev <- rxode2::et(amt = 320, cmt = "depot")
    .ev <- rxode2::et(.ev, seq(0, 24, by = 4))
    suppressMessages(rxode2::rxSolve(.fit, .ev))
    expect_equal(length(.simModelCache$entries), 1L)
    expect_identical(.fit$simulationModel, .mod)

    # a different fitted model is a different key, so it gets its own entry
    expect_false(identical(one.compartment.fit.focei$simulationModel, .mod))
    expect_equal(length(.simModelCache$entries), 2L)

    # and the cache can be turned off, which restores the per-call lowering
    withr::with_options(list(nlmixr2.simModelCache = FALSE), {
      expect_false(identical(.fit$simulationModel, .fit$simulationModel))
    })
  })

  test_that("the cached model solves the same as a freshly lowered one", {
    skip_if(is.null(one.compartment.fit.saem))
    .fit <- one.compartment.fit.saem
    .ev <- rxode2::et(amt = 320, cmt = "depot")
    .ev <- rxode2::et(.ev, seq(0, 24, by = 4))
    .simModelCacheReset()
    on.exit(.simModelCacheReset())

    ## rxSetSeed(), not set.seed(), is what pins the between-subject draws:
    ## rxode2 seeds its own engine from R's RNG at solve entry, so with
    ## set.seed() alone the etas depend on how much R RNG was consumed in
    ## between -- and the uncached arm lowers a fresh model while the cached one
    ## skips lowering entirely, so the two do not consume the same amount.  That
    ## made this comparison pass or fail on an accident of the lowering path
    ## (it failed on CI and in long local sessions, passed in isolation).
    .solve <- function() {
      rxode2::rxSetSeed(42)
      set.seed(42)
      suppressMessages(rxode2::rxSolve(.fit, .ev))
    }
    withr::with_options(list(nlmixr2.simModelCache = FALSE), {
      .uncached <- .solve()
    })
    .cached <- .solve()
    .cachedAgain <- .solve()

    expect_equal(as.data.frame(.cached), as.data.frame(.uncached))
    # a shared model must not accumulate state between solves either
    expect_equal(as.data.frame(.cachedAgain), as.data.frame(.uncached))
  })

  test_that("a cached model still solves after its dll is unloaded", {
    skip_if(is.null(one.compartment.fit.saem))
    .fit <- one.compartment.fit.saem
    .ev <- rxode2::et(amt = 320, cmt = "depot")
    .ev <- rxode2::et(.ev, seq(0, 24, by = 4))
    .simModelCacheReset()
    on.exit(.simModelCacheReset())

    suppressMessages(rxode2::rxSolve(.fit, .ev))
    rxode2::rxUnloadAll()
    # the cached model outlives the unload, so it has to reload itself rather
    # than solve against a dll that is gone
    expect_no_error(suppressMessages(rxode2::rxSolve(.fit, .ev)))
  })
})
