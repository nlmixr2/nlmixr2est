# Work-arounds for rxode2 features newer than the oldest rxode2 nlmixr2est
# supports (DESCRIPTION).  Off CRAN, each block warns once the installed rxode2
# is a release that has its feature, and names what to remove then; before
# that it is quiet.  A new work-around adds a block here.

# The warning a work-around's block gives for the installed version of `pkg`,
# or NULL while it is `dev` (the development version that added the feature,
# or the last one without it) or older; `raise` adds the DESCRIPTION minimum
# to raise
.depWorkaroundWarning <- function(installed, dev, remove, pkg = "rxode2", raise = TRUE) {
  if (installed <= dev) {
    return(NULL)
  }
  paste0(
    pkg,
    " ",
    installed,
    " is a release after ",
    dev,
    ", so this work-around can go: ",
    remove,
    if (raise) paste0("; and raise DESCRIPTION to ", pkg, " (>= ", installed, ")")
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

.depResetRemove <- paste(
  "the linCmt() carry's reset twins wait for an rxode2lincmt that restarts the",
  "carried sensitivities at a reset (rxode2lincmt PR 4): if this one does, remove",
  ".carrySkipReset() (helper-lincmt-carry.R), its calls in",
  "test-focei-lincmt-carry-fit.R, and this block"
)

test_that("the work-around warning names the package, and the minimum only when asked", {
  expect_match(
    .depWorkaroundWarning(numeric_version("0.1.2"), "0.1.1", "x", pkg = "rxode2lincmt", raise = FALSE),
    "^rxode2lincmt 0.1.2 is a release after 0.1.1, so this work-around can go: x$"
  )
})

test_that("rxode2lincmt reset work-around: the carry's reset twins", {
  skip_on_cran()
  .w <- .depWorkaroundWarning(.rxode2lincmtVersion(), "0.1.1", .depResetRemove, pkg = "rxode2lincmt", raise = FALSE)
  if (!is.null(.w)) {
    warning(.w, call. = FALSE)
  }
  succeed()
})

test_that("interval-covariate work-around: the linCmt() carry's end-record fallback", {
  skip_on_cran()
  .w <- .depWorkaroundWarning(.rxode2Version(), "5.1.8.9000", .depIntervalCovsRemove)
  if (!is.null(.w)) {
    warning(.w, call. = FALSE)
  }
  succeed()
})
