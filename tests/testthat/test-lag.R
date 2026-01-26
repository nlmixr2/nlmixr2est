nmTest({
  test_that("test lag with warfarin", {
    # Use centralized model from helper-models.R (one.compartment.with.lag is aliased as KA1Lode)
    KA1Lode <- one.compartment.with.lag

    d <- nlmixr2data::warfarin |>
      dplyr::filter(dvid=="cp")

    f <- .nlmixr(KA1Lode, d, "focei")

    expect_true(f$objf < 500)
  })
})
