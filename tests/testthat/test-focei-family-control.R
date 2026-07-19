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
    .family <- c("focei", "foce", "focep", "fo", "foi",
                 "mfocei", "ifocei", "mfoce", "ifoce",
                 "mfocep", "ifocep",
                 "agq", "magq", "iagq",
                 "laplace", "mlaplace", "ilaplace")
    for (.m in .family) {
      .ctl <- do.call(paste0(.m, "Control"), list())
      expect_true(nlmixr2est:::.nlmixrIsFoceiFamilyControl(.ctl),
                  info = paste0(.m, "Control (class ",
                                paste(class(.ctl), collapse = ","), ")"))
    }
  })

  test_that(".nlmixrIsFoceiFamilyControl rejects non-FOCEi controls", {
    for (.m in c("saem", "nlme", "nlm")) {
      .ctl <- do.call(paste0(.m, "Control"), list())
      expect_false(nlmixr2est:::.nlmixrIsFoceiFamilyControl(.ctl),
                   info = paste0(.m, "Control"))
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
    expect_error(nlmixr2est:::.nlmixrCheckFoceiEnvironment(.env), NA)
    # and it still rejects an unrelated control
    .env$control <- saemControl()
    expect_error(nlmixr2est:::.nlmixrCheckFoceiEnvironment(.env),
                 "focei control object")
  })
})
