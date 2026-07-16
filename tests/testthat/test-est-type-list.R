test_that("tagged estimation-method list (issue #750)", {

  # built-in methods carry type/description attributes
  expect_equal(attr(nlmixr2Est.focei, "type"), "Linearized")
  expect_equal(attr(nlmixr2Est.focei, "description"), "FOCE with Interaction")
  expect_equal(attr(nlmixr2Est.saem, "type"), "Stochastic EM")
  expect_equal(attr(nlmixr2Est.agq, "type"), "Integral approximation")

  # nlmixr2AllEstType() returns a data frame grouped by category
  .df <- nlmixr2AllEstType()
  expect_s3_class(.df, "data.frame")
  expect_true(all(c("est", "type", "description") %in% names(.df)))
  expect_true(all(c("focei", "saem", "nlm") %in% .df$est))
  # canonical category order is preserved (Linearized first)
  expect_equal(.df$type[1], "Linearized")
  # every listed method has a non-empty description
  expect_false(any(.df$description == ""))

  # display lines are grouped and highlight the requested (bad) method
  .lines <- .nlmixr2EstTypeLines(current = "focei")
  expect_true(length(.lines) > length(unique(.df$type)))
  expect_true(any(grepl("Linearized", .lines)))
})

test_that("unsupported est= prints the grouped list", {
  one.cmt <- function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45; add.sd <- 0.7
    })
    model({
      ka <- exp(tka); cl <- exp(tcl); v <- exp(tv)
      linCmt() ~ add(add.sd)
    })
  }
  expect_error(
    nlmixr2(one.cmt, nlmixr2data::theo_sd, est = "notARealMethod"),
    "not supported"
  )
})
