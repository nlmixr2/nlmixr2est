nmTest({
  test_that("model piping", {

    # Use centralized model from helper-models.R (one.compartment.with.lag is aliased as KA1Lode)
    KA1Lode <- one.compartment.with.lag

    d <- nlmixr2data::warfarin |>
      dplyr::filter(dvid=="cp")

    f <- .nlmixr(KA1Lode, data = d, est = "saem", control = saemControlFast)

    # General piping model updates work
    suppressMessages(expect_error(
      fUpV <-
        f |>
        model(v <- exp(lv)),
      NA
    ))
    expect_equal(
      methods::functionBody(as.function(fUpV))[[3]][[2]][[5]],
      str2lang("v <- exp(lv)")
    )

    # piping model updates work with append
    suppressMessages(expect_error(
      fUpFoo <-
        f |>
        model(foo <- exp(lv), append = TRUE),
      NA
    ))
    expect_equal(
      methods::functionBody(as.function(fUpFoo))[[3]][[2]][[11]],
      str2lang("foo <- exp(lv)")
    )
  })
})
