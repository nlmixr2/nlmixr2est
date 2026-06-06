nmTest({
  test_that("model piping", {

    # Use centralized fit from helper-fits.R
    f <- one.compartment.with.lag.fit.saem

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
