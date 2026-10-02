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

test_that("the fit's ODE method code is read back by name for any rxode2 method", {
  .codes <- rxode2::odeMethodToInt(NULL)
  # the four methods the old hard-coded factor knew
  expect_equal(.residOdeMethodName(.codes[["dop853"]]), "dop853")
  expect_equal(.residOdeMethodName(.codes[["lsoda"]]), "lsoda")
  expect_equal(.residOdeMethodName(.codes[["liblsoda"]]), "liblsoda")
  expect_equal(.residOdeMethodName(.codes[["indLin"]]), "indLin")
  # newer ones were a "malformed factor" error
  for (.m in intersect(c("cvode", "lsode", "bdf", "dop5"), names(.codes))) {
    .n <- .residOdeMethodName(.codes[[.m]])
    expect_equal(.codes[[.n]], .codes[[.m]])
  }
  # the rxControl() form (named integer) and a name pass through
  expect_equal(.residOdeMethodName(rxode2::rxControl(method = "liblsoda")$method), "liblsoda")
  expect_equal(.residOdeMethodName("cvode"), "cvode")
  expect_error(.residOdeMethodName(99999L), "unknown rxode2 ODE method code")
})

test_that("a fit solved with a newer rxode2 ODE method still gets its table", {
  skip_on_cran()
  .codes <- rxode2::odeMethodToInt(NULL)
  skip_if_not("cvode" %in% names(.codes))
  one.cmt <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }
  .fit <- suppressMessages(suppressWarnings(
    nlmixr2(
      one.cmt,
      nlmixr2data::theo_sd,
      "focei",
      foceiControl(maxOuterIterations = 0, print = 0, rxControl = rxode2::rxControl(method = "cvode"))
    )
  ))
  expect_equal(.residOdeMethodName(.fit$methodOde), "cvode")
  # addTable() recalculates the table with the fit's own ODE method
  .tab <- suppressMessages(addTable(.fit))
  expect_true(inherits(.tab, "nlmixr2FitData"))
  expect_true(all(c("IPRED", "PRED") %in% names(.tab)))
})
