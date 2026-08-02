# npag/npb accept an impmapControl, so the importance-sampling control surface
# used to be accepted and then never read -- silently.  These tests pin the
# rejection, and they matter precisely because the bug they replace was SILENCE:
# nothing errored, nothing warned, and the fit simply ignored the request.
#
# Control-level only, no fits, so this file stays in the push/PR subset.
nmTest({

  test_that("inapplicable importance-sampling controls are rejected", {
    for (.e in c("npag", "npb")) {
      .ctl <- if (.e == "npag") npagControl else npbControl
      for (.n in c("isample", "df", "auto", "iaccept", "qr", "sir",
                   "iscaleMin", "autoDfPatience", "mapIter")) {
        expect_error(do.call(.ctl, setNames(list(1), .n)), regexp = .n,
                     info = paste(.e, .n))
      }
    }
  })

  test_that("controls with an np counterpart name it in the message", {
    expect_error(npagControl(nIter = 200), "cycles")
    expect_error(npagControl(ctol = 1e-6), "rhoend")
    expect_error(npbControl(nIter = 200), "cycles")
    # npag has no seed of its own (Sobol-deterministic), npb does
    expect_error(npbControl(impSeed = 7), "seed")
    expect_error(npagControl(impSeed = 7), "importance-sampling")
  })

  test_that("partial matching cannot smuggle an inert control in", {
    # `gamma` is a prefix of the real formal `gammaOptimize`, so R bound it there
    # and npagControl(gamma = 2) set gammaOptimize = isTRUE(2) = FALSE -- turning
    # the assay-error optimisation OFF rather than being ignored.  match.call()
    # normalises partial matching away, so the check has to read sys.call().
    expect_error(npagControl(gamma = 2), "gamma")
    # ... while the real formal still works, and still means what it says
    expect_false(npagControl(gammaOptimize = FALSE)$gammaOptimize)
    expect_true(npagControl()$gammaOptimize)
  })

  test_that("npb's inert FORMALS are rejected too", {
    # cycles and gammaOptimize are formals of npbControl documented unused for
    # npb, so they never reach ... and a names(list(...)) check cannot see them
    expect_error(npbControl(cycles = 500), "cycles")
    expect_error(npbControl(gammaOptimize = TRUE), "gammaOptimize")
    # they are real npag controls, so npag must still take them
    expect_equal(npagControl(cycles = 50)$cycles, 50L)
    expect_true(npagControl(gammaOptimize = TRUE)$gammaOptimize)
  })

  test_that("np controls and shared scaffolding are still accepted", {
    expect_silent(npagControl(points = 40, cycles = 50, gridWidth = 3))
    expect_silent(npbControl(points = 40, alpha = 2, seed = 7, nchains = 2))
    # impCov is shared FOCEI-family scaffolding, not a proposal control
    expect_true(is.list(npagControl(covMethod = "")))
  })

  test_that("a rebuild is exempt -- the validators round-trip", {
    # The supported round-trip is getValidNlmixrCtl, NOT do.call on the
    # constructor: the control stamps npCores while the formal is cores, so
    # do.call(npagControl, npagControl()) has never worked.  Pin what IS
    # supported so a future rejection cannot break it.
    expect_silent(getValidNlmixrCtl.npag(list(npagControl())))
    expect_silent(getValidNlmixrCtl.npb(list(npbControl())))
    expect_equal(getValidNlmixrCtl.npag(list(npagControl(cycles = 7)))$cycles, 7L)
  })

  test_that("a raw list bypassing the constructor is still checked", {
    # getValidNlmixrCtl.npag routes a bare list through .npValidCtl ->
    # getValidNlmixrCtl.impmap, never entering npagControl(), so a check living
    # only in the constructor would miss this
    expect_error(getValidNlmixrCtl.npag(list(list(isample = 500))), "isample")
    expect_error(getValidNlmixrCtl.npb(list(list(df = 8))), "df")
  })

})
