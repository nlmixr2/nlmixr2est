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

  test_that("the five review-found bypasses stay closed", {
    # Each of these passed an earlier implementation.  They are kept as named
    # regressions because every one of them was SILENT: the control was accepted
    # and then ignored, which is the whole failure mode being fixed.
    .foo <- function(...) npagControl(...)
    .bar <- function(...) npbControl(...)
    # 1. an internal field alongside an inert one no longer exempts the call
    expect_error(npagControl(isample = 500, impCov = TRUE), "isample")
    # 2. partial matching through a forwarding wrapper, where sys.call() cannot
    #    see the literal name -- caught by the explicit inert formals
    expect_error(.foo(gamma = 2), "gamma")
    expect_error(.foo(df = 8), "df")
    expect_error(.bar(gamma = 2), "gamma")
    # 3. a pre-built impmapControl handed to an np validator
    expect_error(getValidNlmixrCtl.npag(list(impmapControl(isample = 500))), "isample")
    # 4. npb's inert formals arriving as a raw list, bypassing npbControl()
    expect_error(getValidNlmixrCtl.npb(list(list(cycles = 500))), "cycles")
    # 5. the mu/irls sugar engines are normalised, so impSeed still remaps
    expect_error(getValidNlmixrCtl.mnpb(list(list(impSeed = 7))), "seed")
    # and a name typed with its DEFAULT value is still a request
    expect_error(npagControl(mapIter = 1), "mapIter")
  })

  test_that("a fit-stamped control re-validates", {
    # .impmapFamilyFit RESOLVES gammaMethod ("auto" -> "global"/"individual") and
    # STAMPS autoNonNormal onto the runtime control, so a real fit's control
    # legitimately differs from what the constructor returns.  Two separate bugs
    # met here: the validator read those as user requests, and impmapControl()
    # never stripped autoNonNormal (so do.call on a fit's own control had always
    # died with "unused argument"), which is why this path had never worked.
    for (.e in c("npag", "npb")) {
      .ctl <- if (.e == "npag") npagControl() else npbControl()
      .ctl$gammaMethod <- "global"
      .ctl$autoNonNormal <- TRUE
      .v <- if (.e == "npag") getValidNlmixrCtl.npag else getValidNlmixrCtl.npb
      .out <- .v(list(.ctl))
      expect_true(.out$autoNonNormal, info = .e)
    }
    # impmapControl itself round-trips a stamped control now
    .ic <- impmapControl(); .ic$autoNonNormal <- TRUE
    expect_true(do.call(impmapControl, .ic)$autoNonNormal)
    # but an explicitly TYPED gammaMethod is still a request, and still rejected
    expect_error(npagControl(gammaMethod = "global"), "gammaMethod")
  })

  test_that("an inert FORMAL is rejected on presence, not on value", {
    # gamma/df exist as formals only to be rejected, so their PRESENCE is the
    # request.  Comparing values would let a wrapper through merely because the
    # value happened to equal impmap's default for that control.
    .foo <- function(...) npagControl(...)
    .bar <- function(...) npbControl(...)
    expect_error(.foo(gamma = 1), "gamma")     # 1 IS impmap's default gamma
    expect_error(.foo(df = 0), "df")           # 0 IS impmap's default df
    expect_error(.bar(gamma = 1), "gamma")
    # npag's real controls are untouched -- test-npag-error-models.R relies on
    # passing cycles/gammaOptimize to npagControl
    expect_equal(npagControl(cycles = 30L, gammaOptimize = TRUE)$cycles, 30L)
    expect_true(npagControl(gammaOptimize = TRUE)$gammaOptimize)
    # npbControl(cycles = 100L) IS rejected although 100L is its own default:
    # cycles is documented unused for npb, so typing it is a request for
    # something that does nothing.  Deliberate, and pinned so it stays that way.
    expect_error(npbControl(cycles = 100L), "cycles")
  })

  test_that("an inert control that TOOK EFFECT is caught however it got there", {
    # A forwarding wrapper hides the literal names from sys.call(), so the
    # name-based check cannot see them; checking the BUILT control closes that.
    #
    # Abbreviated names need no help: impmapControl() declares isample/iaccept
    # AFTER `...` in its signature, and R requires an exact match for those, so
    # `isampl` is already an "unused argument" error.  Pinned anyway, because
    # moving a formal ahead of `...` would silently re-open partial matching.
    expect_error(npagControl(isampl = 500))
    expect_error(npagControl(iaccep = 0.2))
    expect_error(npbControl(isampl = 500))
    .foo <- function(...) npagControl(...)
    expect_error(.foo(isample = 500), "isample")
    # ... and npb gained the `df` formal npag already had, so a wrapper cannot
    # smuggle df through partial matching there either
    .bar <- function(...) npbControl(...)
    expect_error(.bar(df = 8), "df")
    expect_error(.bar(df = 0), "df")
    # A value equal to the default changed nothing, so through a wrapper -- where
    # the literal name is invisible -- there is nothing to report.  A DIRECT call
    # still errors, because the name was typed.
    expect_true(is.list(.foo(mapIter = 1)))
    expect_error(npagControl(mapIter = 1), "mapIter")
  })

})
