nmTest({
  # A FOCEi-family control is built by foceiControl() and then reclassed to its
  # own class -- the mu-referenced and method variants (ifocei, mfocei, foce,
  # ...) do NOT keep "foceiControl" in their class vector.  The estimation
  # restart path re-validates the environment with .nlmixrCheckFoceiEnvironment,
  # which used to require inherits(control, "foceiControl") and so aborted any
  # restart of a mu-referenced fit with
  #   "focei$control must be a focei control object"
  # even though the fit itself was set up from a valid control.  The check now
  # accepts the whole family via .nlmixrIsFoceiFamilyControl().

  test_that(".nlmixrIsFoceiFamilyControl accepts every FOCEi-family control", {
    .family <- c(
      "focei",
      "foce",
      "focep",
      "fo",
      "foi",
      "mfocei",
      "ifocei",
      "mfoce",
      "ifoce",
      "mfocep",
      "ifocep",
      "agq",
      "magq",
      "iagq",
      "laplace",
      "mlaplace",
      "ilaplace"
    )
    for (.m in .family) {
      .ctl <- do.call(paste0(.m, "Control"), list())
      expect_true(
        .nlmixrIsFoceiFamilyControl(.ctl),
        info = paste0(.m, "Control (class ", paste(class(.ctl), collapse = ","), ")")
      )
    }
  })

  test_that(".nlmixrIsFoceiFamilyControl rejects non-FOCEi controls", {
    for (.m in c("saem", "nlme", "nlm")) {
      .ctl <- do.call(paste0(.m, "Control"), list())
      expect_false(.nlmixrIsFoceiFamilyControl(.ctl), info = paste0(.m, "Control"))
    }
  })

  test_that(".nlmixrCheckFoceiEnvironment does not reject a mu-referenced control", {
    # a minimal environment with the fields the check inspects; the control is
    # an ifoceiControl, which pre-fix tripped the class assertion
    .env <- new.env()
    .env$dataSav <- data.frame(ID = 1L, TIME = 0, DV = 1)
    .env$thetaIni <- c(1, 2)
    .env$skipCov <- NULL
    .env$rxInv <- structure(list(), class = "rxSymInvCholEnv")
    .env$lower <- c(-Inf, -Inf)
    .env$upper <- c(Inf, Inf)
    .env$etaMat <- NA
    .env$control <- ifoceiControl()
    expect_error(.nlmixrCheckFoceiEnvironment(.env), NA)
    # and it still rejects an unrelated control
    .env$control <- saemControl()
    expect_error(.nlmixrCheckFoceiEnvironment(.env), "focei control object")
  })
})

test_that("foce, focep, laplace and agq convert a family control by what its caller set", {
  .targets <- c(foce = "foceControl", focep = "focepControl", laplace = "laplaceControl", agq = "agqControl")
  .identity <- c("fo", "interaction", "nAGQ", "foce")
  for (.est in names(.targets)) {
    .default <- do.call(.targets[[.est]], list())
    for (.src in c("foControl", "foiControl")) {
      for (.posthoc in c(TRUE, FALSE)) {
        .in <- do.call(.src, list(maxOuterIterations = 7L, posthoc = .posthoc))
        expect_message(
          .ctl <- getValidNlmixrControl(.in, .est),
          paste0("converting ", .src, " to ", .targets[[.est]]),
          fixed = TRUE
        )
        expect_s3_class(.ctl, .targets[[.est]], exact = TRUE)
        expect_identical(.ctl$maxOuterIterations, 7L)
        # posthoc is a field of foControl()/foiControl() only
        expect_false("posthoc" %in% names(.ctl))
        # the target method's own settings, not FO's
        expect_identical(.ctl[.identity], .default[.identity])
      }
    }
  }
  # the settings that make a method what it is come from the target, so agq
  # given a foceControl() runs AGQ, not FOCE, and laplace given an
  # agqControl() runs the Laplace approximation
  .agq <- suppressMessages(getValidNlmixrControl(foceControl(maxOuterIterations = 7L), "agq"))
  expect_identical(.agq[.identity], agqControl()[.identity])
  expect_identical(.agq$maxOuterIterations, 7L)
  .lap <- suppressMessages(getValidNlmixrControl(foceControl(), "laplace"))
  expect_identical(.lap[.identity], laplaceControl()[.identity])
  .lap <- suppressMessages(getValidNlmixrControl(agqControl(), "laplace"))
  expect_identical(.lap[.identity], laplaceControl()[.identity])
  # an explicit setting is kept, and a foceiControl() is passed through as is
  .lap <- suppressMessages(getValidNlmixrControl(agqControl(nAGQ = 5), "laplace"))
  expect_equal(.lap$nAGQ, 5)
  .agq <- suppressMessages(getValidNlmixrControl(foceiControl(nAGQ = 3), "agq"))
  expect_equal(.agq$nAGQ, 3)
  # fo and foi still take their own posthoc
  .fo <- suppressMessages(getValidNlmixrControl(foiControl(posthoc = FALSE), "fo"))
  expect_false(.fo$posthoc)
})
