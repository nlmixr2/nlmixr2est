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

.depSlot104Remove <- paste(
  "rxode2's setIndNewind() (slot 104): remove nmSetIndNewind from src/rxode2Optional.h",
  "and its row of rxode2OptionalSlots in src/inner.cpp, the NULL check in lhsWalkNewind()",
  "(src/inner.h, which then calls setIndNewind()), the skip in test-lhs-walk-newind.R,",
  "and this block (and the rxode2Optional machinery if no entry is left)"
)

test_that("the work-around warning is given only for a release after the development version", {
  expect_null(.depWorkaroundWarning(numeric_version("5.1.8"), "5.1.8.9000", "x"))
  expect_null(.depWorkaroundWarning(numeric_version("5.1.8.9000"), "5.1.8.9000", "x"))
  expect_match(
    .depWorkaroundWarning(numeric_version("5.1.9"), "5.1.8.9000", "x"),
    "rxode2 5.1.9 is a release after 5.1.8.9000.*x; and raise DESCRIPTION to rxode2 \\(>= 5.1.9\\)"
  )
})

test_that("slot-104 work-around: rxode2's setIndNewind()", {
  skip_on_cran()
  .w <- .depWorkaroundWarning(.rxode2Version(), "5.1.8.9000", .depSlot104Remove)
  if (!is.null(.w)) {
    warning(.w, call. = FALSE)
  }
  succeed()
})
