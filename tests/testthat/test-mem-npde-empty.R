test_that("npdeCalc with zero-row input errors gracefully (no crash)", {
  # empty input previously caused a SIGSEGV reading id[0]/simId[0]
  empty_npdesim <- list(
    simId = integer(0), id = integer(0), ipred = numeric(0),
    pred = numeric(0), lambda = numeric(0), yj = numeric(0),
    low = numeric(0), hi = numeric(0)
  )
  expect_error(
    .Call("_nlmixr2est_npdeCalc", empty_npdesim, numeric(0),
          NULL, NULL, NULL, list(), PACKAGE = "nlmixr2est"),
    "zero rows"
  )
})
