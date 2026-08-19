test_that("indLin() is available in the post-fit ODE-method fallback list (#858)", {
  # "dop853"/"liblsoda"/"lsoda" primaries hard-override the candidate list to a
  # single named fallback, so "indLin" never appears there regardless.
  expect_equal(.residOdeFallbackMethods("dop853"), list("dop853", "liblsoda"))
  expect_equal(.residOdeFallbackMethods("liblsoda"), list("liblsoda", "dop853"))
  expect_equal(.residOdeFallbackMethods("lsoda"), list("lsoda", "dop853"))

  # A primary of "indLin" itself always tries "indLin" first; it can never
  # also appear later in its own fallback list (setdiff), so this is
  # unaffected by the #858 change -- confirm the invariant directly.
  .indLinCandidates <- .residOdeFallbackMethods("indLin")
  expect_equal(.indLinCandidates[[1]], "indLin")
  expect_equal(sum(unlist(.indLinCandidates) == "indLin"), 1L)

  # Any other explicitly-requested ODE method (e.g. an obscure integrator set
  # via rxControl(method=)) falls through to the "use all the methods"
  # branch -- this is the only branch #858's setdiff removal affects, and
  # "indLin" must now be a candidate there.
  .rk4Candidates <- .residOdeFallbackMethods("rk4")
  expect_equal(.rk4Candidates[[1]], "rk4")
  expect_true("indLin" %in% unlist(.rk4Candidates))
})
