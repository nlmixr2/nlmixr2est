# A linCmt() model whose carry-eligible etas need more carry columns than the
# loaded rxode2 has is solved as ODEs (#1318).  The fits are in
# test-focei-lincmt-carry-iov.R.

test_that("a model past the carry cap is translated to ODEs before a FOCEi fit", {
  skip_if_not(.rxFoceiLinCmtCarryCapable())
  local_mocked_bindings(.rxFoceiLinCmtCarryMaxPairs = .carryCap8)
  local_mocked_bindings(.rxFoceiLinCmtCarryColumnsNeeded = function(ui, cap) 10L)
  .ui <- suppressMessages(nlmixr2est::nlmixr2(.carryModCov))
  expect_warning(
    .ret <- .uiApplyLinCmtCarryCap(.ui, "focei", NULL, nlmixr2est::foceiControl()),
    "linCmt() solved as ODEs: etas need 10 of 8 carry columns",
    fixed = TRUE
  )
  expect_true(rxode2::.rxLinNcmt(.ret$ui)["numLin"] == 0L)
  expect_true(length(.ret$ui$state) > 0L)
  # every member of the FOCEi family, but not the others
  expect_false(is.null(suppressWarnings(.uiApplyLinCmtCarryCap(.ui, "laplace", NULL, NULL))))
  expect_null(.uiApplyLinCmtCarryCap(.ui, "saem", NULL, NULL))
  expect_null(.uiApplyLinCmtCarryCap(.ui, "nlm", NULL, NULL))
  # nothing is carried with the carry off or linear interpolation
  expect_null(.uiApplyLinCmtCarryCap(
    .ui,
    "focei",
    NULL,
    nlmixr2est::foceiControl(linCmtSensCarry = "none")
  ))
  expect_null(.uiApplyLinCmtCarryCap(
    .ui,
    "focei",
    NULL,
    nlmixr2est::foceiControl(rxControl = rxode2::rxControl(covsInterpolation = "linear"))
  ))
})

test_that("a model within the carry cap keeps linCmt()", {
  skip_if_not(.rxFoceiLinCmtCarryCapable())
  local_mocked_bindings(.rxFoceiLinCmtCarryMaxPairs = .carryCap8)
  local_mocked_bindings(.rxFoceiLinCmtCarryColumnsNeeded = function(ui, cap) 8L)
  .ui <- suppressMessages(nlmixr2est::nlmixr2(.carryModCov))
  expect_null(.uiApplyLinCmtCarryCap(.ui, "focei", NULL, nlmixr2est::foceiControl()))
})

test_that("the carry columns are only counted for a model that can fill them", {
  # no covariate, or too few etas for 8 columns: no symbolic check is needed
  local_mocked_bindings(.foceiLinCmtCarryPairs = function(...) stop("not needed"))
  expect_identical(
    .rxFoceiLinCmtCarryColumnsNeeded(suppressMessages(nlmixr2est::nlmixr2(.carryModNoCov)), 8L),
    0L
  )
  expect_identical(
    .rxFoceiLinCmtCarryColumnsNeeded(suppressMessages(nlmixr2est::nlmixr2(.carryModCov)), 8L),
    0L
  )
})

test_that("pairs dropped past the cap warn that their gradient is approximate", {
  # eta.cl and eight occasion etas in cl's slot, with the slot's tracker
  .pairs <- data.frame(
    slot = 2L,
    etaName = c("eta.cl", paste0("rx.iov.cl.", 1:8)),
    shape = c("mult", rep("multCov", 8L)),
    fD = NA_character_,
    lagD = NA_character_
  )
  local_mocked_bindings(.rxFoceiLinCmtCarryMaxPairs = .carryCap8)
  nlmixr2global$carryCapWarnings <- NULL
  on.exit(nlmixr2global$carryCapWarnings <- NULL, add = TRUE)
  expect_warning(
    .kept <- .rxFoceiLinCmtCarryCapPairs(.pairs),
    "rx.iov.cl.7 and 1 more have an approximate gradient (8 carry columns)",
    fixed = TRUE
  )
  expect_identical(.kept$etaName, c("eta.cl", paste0("rx.iov.cl.", 1:6)))
})
