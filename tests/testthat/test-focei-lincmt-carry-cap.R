# A linCmt() model whose carry-eligible etas need more carry columns than the
# loaded rxode2 has is solved as ODEs (#1318).  The fits are in
# test-focei-lincmt-carry-iov.R.

test_that("a model past the carry cap is translated to ODEs before a FOCEi fit", {
  skip_if_not(.rxFoceiLinCmtCarryCapable())
  local_mocked_bindings(.rxFoceiLinCmtCarryMaxPairs = .carryCap8)
  local_mocked_bindings(.rxFoceiLinCmtCarryColumnsNeeded = function(ui, cap, data = NULL) 10L)
  .ui <- suppressMessages(nlmixr2est::nlmixr2(.carryModCov))
  expect_warning(
    .ret <- .uiApplyOdeLinCmtCarryCap(.ui, "focei", NULL, nlmixr2est::foceiControl()),
    "linCmt() solved as ODEs: etas need 10 of 8 carry columns",
    fixed = TRUE
  )
  expect_true(rxode2::.rxLinNcmt(.ret$ui)["numLin"] == 0L)
  expect_true(length(.ret$ui$state) > 0L)
  # every member of the FOCEi family, but not the others
  for (.est in c("laplace", "agqf", "flaplace", "vae")) {
    expect_false(is.null(suppressWarnings(.uiApplyOdeLinCmtCarryCap(.ui, .est, NULL, NULL))), label = .est)
  }
  expect_null(.uiApplyOdeLinCmtCarryCap(.ui, "saem", NULL, NULL))
  expect_null(.uiApplyOdeLinCmtCarryCap(.ui, "nlm", NULL, NULL))
  # nothing is carried with the carry off or linear interpolation
  expect_null(.uiApplyOdeLinCmtCarryCap(
    .ui,
    "focei",
    NULL,
    nlmixr2est::foceiControl(linCmtSensCarry = "none")
  ))
  expect_null(.uiApplyOdeLinCmtCarryCap(
    .ui,
    "focei",
    NULL,
    nlmixr2est::foceiControl(rxControl = rxode2::rxControl(covsInterpolation = "linear"))
  ))
})

test_that("a model within the carry cap keeps linCmt()", {
  skip_if_not(.rxFoceiLinCmtCarryCapable())
  local_mocked_bindings(.rxFoceiLinCmtCarryMaxPairs = .carryCap8)
  local_mocked_bindings(.rxFoceiLinCmtCarryColumnsNeeded = function(ui, cap, data = NULL) 8L)
  .ui <- suppressMessages(nlmixr2est::nlmixr2(.carryModCov))
  expect_null(.uiApplyOdeLinCmtCarryCap(.ui, "focei", NULL, nlmixr2est::foceiControl()))
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

# four etas and a covariate: enough to need the symbolic count
.carryCapMod <- function() {
  ini({
    tka <- 0.45
    tcl <- 1
    tv <- 3.45
    tq <- 1
    tvp <- 4
    add.sd <- 0.7
    eta.ka ~ 0.6
    eta.cl ~ 0.3
    eta.v ~ 0.1
    eta.q ~ 0.1
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl) * (wt / 70)
    v <- exp(tv + eta.v)
    q <- exp(tq + eta.q)
    vp <- exp(tvp)
    linCmt() ~ add(add.sd)
  })
}

# eta.cl and eight occasion etas in cl's slot (wt on cl), with the slot's
# tracker: 10 columns, of which the last two occasion etas are dropped at 8
.carryCapPairs <- function(modelOk = TRUE) {
  .p <- data.frame(
    slot = 2L,
    etaName = c("eta.cl", paste0("rx.iov.cl.", 1:8)),
    covs = c("wt", rep("wt,occ", 8L)),
    shape = c("mult", rep("multCov", 8L)),
    fD = NA_character_,
    lagD = NA_character_
  )
  attr(.p, "modelOk") <- modelOk
  .p
}

test_that("a model is translated only when a dropped pair varies within a subject", {
  local_mocked_bindings(.foceiLinCmtCarryPairs = function(...) .carryCapPairs())
  .ui <- suppressMessages(nlmixr2est::nlmixr2(.carryCapMod))
  .d <- data.frame(id = rep(1:2, each = 3), wt = 70, occ = 1)
  # constant within every subject: the dropped pairs' gradient is exact
  expect_identical(.rxFoceiLinCmtCarryColumnsNeeded(.ui, 8L, .d), 0L)
  .d$occ <- rep(1:3, 2)
  expect_identical(.rxFoceiLinCmtCarryColumnsNeeded(.ui, 8L, .d), 10L)
  # without data, the covariates may vary
  expect_identical(.rxFoceiLinCmtCarryColumnsNeeded(.ui, 8L, NULL), 10L)
  # within the cap
  expect_identical(.rxFoceiLinCmtCarryColumnsNeeded(.ui, 10L, NULL), 0L)
  # a parameterization the carry cannot handle is never carried
  local_mocked_bindings(.foceiLinCmtCarryPairs = function(...) .carryCapPairs(FALSE))
  expect_identical(.rxFoceiLinCmtCarryColumnsNeeded(.ui, 8L, NULL), 0L)
})

test_that("a model the build cannot fit within the cap is translated whatever the data", {
  local_mocked_bindings(.foceiLinCmtCarryPairs = function(...) {
    .p <- .carryCapPairs()[1:5, ]
    .p$slot <- 1:5
    .p$shape <- "mult"
    .p
  })
  .ui <- suppressMessages(nlmixr2est::nlmixr2(.carryCapMod))
  .d <- data.frame(id = rep(1:2, each = 3), wt = 70, occ = 1)
  expect_identical(.rxFoceiLinCmtCarryColumnsNeeded(.ui, 8L, .d), 10L)
})

test_that("a numeric cmt keeps naming the same compartment after linToOde()", {
  .f <- function() {
    ini({
      tka <- 0
      tcl <- 1
      tv <- 3
      tq <- 1
      tvp <- 3
      tef <- 1
      add.sd <- 1
      ef.sd <- 1
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl)
      v <- exp(tv)
      q <- exp(tq)
      vp <- exp(tvp)
      cp <- linCmt()
      cp ~ add(add.sd)
      ef <- tef
      ef ~ add(ef.sd)
    })
  }
  .ui <- suppressMessages(rxode2::rxode2(.f))
  .ode <- suppressMessages(rxode2::linToOde(.ui))
  .d <- data.frame(id = 1, time = 0:4, CMT = c(1, 2, 3, 4, -1))
  .r <- .linCmtCarryRemapCmt(.d, .ui, .ode)
  expect_identical(
    .linCmtCarryCmtNames(.ode)[abs(.r$CMT)],
    .linCmtCarryCmtNames(.ui)[abs(.d$CMT)]
  )
  expect_identical(.r$CMT[5], -1)
  # a named cmt is left alone
  .d$CMT <- c("depot", "central", "cp", "ef", "depot")
  expect_identical(.linCmtCarryRemapCmt(.d, .ui, .ode), .d)
})

test_that("pairs dropped past the cap warn that their gradient is approximate", {
  local_mocked_bindings(.rxFoceiLinCmtCarryMaxPairs = .carryCap8)
  nlmixr2global$carryCapWarnings <- NULL
  on.exit(nlmixr2global$carryCapWarnings <- NULL, add = TRUE)
  expect_warning(
    .kept <- .rxFoceiLinCmtCarryCapPairs(.carryCapPairs()),
    "rx.iov.cl.7 and 1 more have an approximate gradient (8 carry columns)",
    fixed = TRUE
  )
  expect_identical(.kept$etaName, c("eta.cl", paste0("rx.iov.cl.", 1:6)))
})
