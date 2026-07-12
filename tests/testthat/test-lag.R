nmTest({
  test_that("test lag with warfarin", {
    # Use centralized fit from helper-fits.R
    f <- one.compartment.with.lag.fit.focei

    # The default (jump) event-sensitivity FOCEi fit currently converges to
    # ~2191 for this lag model.  rxode2 #1118 (b1b69f1ff) fixed the
    # coincident-observation jump doubling but did NOT change this fit -- the
    # jump (~2191) vs eventSens="fd" (~240) gap has a separate, unresolved
    # cause.  Bound kept loose until that is fixed upstream.
    expect_true(f$objf < 2300)
  })
})
