# Work-arounds for rxode2 features newer than the oldest rxode2 nlmixr2est
# supports (DESCRIPTION).  Off CRAN, each block warns once the installed rxode2
# is a release that has its feature, and names what to remove then; before
# that it is quiet.  A new work-around adds a block here.

# The warning a work-around's block gives for the installed rxode2 version, or
# NULL while rxode2 is `dev` (the development version that added the feature)
# or older
.depWorkaroundWarning <- function(installed, dev, remove) {
  if (installed <= dev) {
    return(NULL)
  }
  paste0(
    "rxode2 ",
    installed,
    " is a release after ",
    dev,
    ", so this work-around can go: ",
    remove,
    "; and raise DESCRIPTION to rxode2 (>= ",
    installed,
    ")"
  )
}

.depIntervalCovsRemove <- paste(
  "the linCmt() carry's end-record fallback for an rxode2 without interval covariates",
  "(rxode2 PR 1517): remove R/foceiLinCmtCarryInterval.R and its calls in",
  ".rxFoceiLinCmtCarryMidCovs() and .rxFoceiLinCmtCarryCovInterval()",
  "(R/foceiLinCmtCarryCompose.R), .carryTwinInterps() (helper-lincmt-carry.R), the",
  "probe's skips in test-focei-lincmt-carry-fit.R and its check in",
  "test-focei-lincmt-carry-emit.R, and this block"
)

test_that("the work-around warning is given only for a release after the development version", {
  expect_null(.depWorkaroundWarning(numeric_version("5.1.8"), "5.1.8.9000", "x"))
  expect_null(.depWorkaroundWarning(numeric_version("5.1.8.9000"), "5.1.8.9000", "x"))
  expect_match(
    .depWorkaroundWarning(numeric_version("5.1.9"), "5.1.8.9000", "x"),
    "rxode2 5.1.9 is a release after 5.1.8.9000.*x; and raise DESCRIPTION to rxode2 \\(>= 5.1.9\\)"
  )
})

test_that("interval-covariate work-around: the linCmt() carry's end-record fallback", {
  skip_on_cran()
  .w <- .depWorkaroundWarning(.rxode2Version(), "5.1.8.9000", .depIntervalCovsRemove)
  if (!is.null(.w)) {
    warning(.w, call. = FALSE)
  }
  succeed()
})

.depCarryCapRemove <- paste(
  "the linCmt() carry's 8-column fallback for an rxode2 without",
  "linCmtCarryMaxPairs() (its RX_LINCMT_CARRY_MAXPAIRS of 32): have",
  ".rxFoceiLinCmtCarryMaxPairs() (R/foceiLinCmtCarryGates.R) call",
  "rxode2::linCmtCarryMaxPairs(), drop the skip_if() on",
  ".rxFoceiLinCmtCarryMaxPairs() < 10L in test-focei-lincmt-carry-iov.R's",
  "nine-pair test, and remove this block and the probe test after it"
)

test_that("carry-cap work-around: 8 carry columns without linCmtCarryMaxPairs()", {
  .ns <- asNamespace("rxode2")
  .has <- exists("linCmtCarryMaxPairs", envir = .ns, inherits = FALSE)
  if (.has) {
    expect_identical(.rxFoceiLinCmtCarryMaxPairs(), get("linCmtCarryMaxPairs", envir = .ns)())
  } else {
    expect_identical(.rxFoceiLinCmtCarryMaxPairs(), 8L)
  }
  skip_on_cran()
  # a release after the development version that still lacks the function
  # still needs the fallback
  .w <- if (.has) .depWorkaroundWarning(.rxode2Version(), "5.1.8.9000", .depCarryCapRemove)
  if (!is.null(.w)) {
    warning(.w, call. = FALSE)
  }
  succeed()
})

test_that("the carry-cap probe gives 8 when rxode2's count fails or is no count", {
  skip_if_not(exists("linCmtCarryMaxPairs", envir = asNamespace("rxode2"), inherits = FALSE))
  local_mocked_bindings(linCmtCarryMaxPairs = .carryCapError, .package = "rxode2")
  expect_identical(.rxFoceiLinCmtCarryMaxPairs(), 8L)
  local_mocked_bindings(linCmtCarryMaxPairs = .carryCap0, .package = "rxode2")
  expect_identical(.rxFoceiLinCmtCarryMaxPairs(), 8L)
  local_mocked_bindings(linCmtCarryMaxPairs = .carryCap32, .package = "rxode2")
  expect_identical(.rxFoceiLinCmtCarryMaxPairs(), 32L)
})
