test_that("npdeCalc with zero-row input errors gracefully (no crash)", {
  # Before fix: getSimIdLoc reads id[0]/simId[0] from empty arma vectors,
  # causing a SIGSEGV (UBSan reports load of misaligned address 0x000000000001
  # at npde.cpp:11).
  # After fix: returns informative error instead of crashing.
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
