# Work-arounds for rxode2 features newer than the oldest rxode2 nlmixr2est
# supports (DESCRIPTION).  Off CRAN, each block warns once the installed rxode2
# is a release that has its feature, and names what to remove then; before
# that it is quiet.  A new work-around adds a block here, and a row to
# .depBlocks so the files it touches are checked against what it names.

# The warning a work-around's block gives, or NULL unless the block's probe
# found the feature (`present`) and the installed version is a release after
# `dev`, the development version that added it.  A development build
# (x.y.z.9000 and later) is never a release, and a release without the
# feature still needs the work-around.
.depWorkaroundWarning <- function(installed, dev, remove, present) {
  if (!isTRUE(present) || !.depIsRelease(installed) || installed <= dev) {
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

# A release rather than a development version (a fourth component of 9000 or
# more)
.depIsRelease <- function(version) {
  .p <- unlist(numeric_version(version))
  length(.p) < 4L || .p[[4L]] < 9000L
}

# Does the file use any of `symbols`?
.depFileUses <- function(file, symbols) {
  .l <- readLines(file, warn = FALSE)
  any(vapply(symbols, .depLinesUse, logical(1), lines = .l))
}

.depLinesUse <- function(symbol, lines) {
  any(grepl(symbol, lines, fixed = TRUE))
}

# The files under R/ and src/ (when the source tree is there) and
# tests/testthat that use any of `symbols`, other than this file
.depFilesUsing <- function(symbols) {
  .root <- testthat::test_path("..", "..")
  .dirs <- c(file.path(.root, c("R", "src")), testthat::test_path())
  .f <- unlist(lapply(
    .dirs[dir.exists(.dirs)],
    list.files,
    pattern = "\\.(R|c|cpp|h|hpp)$",
    full.names = TRUE
  ))
  .f <- .f[basename(.f) != "test-dependency-workarounds.R"]
  basename(.f[vapply(.f, .depFileUses, logical(1), symbols = symbols)])
}

# stand-ins for the version and probe mocks
.depVersion519 <- function() numeric_version("5.1.9")
.depVersionDev <- function() numeric_version("5.1.9.9000")
.depTrue <- function() TRUE
.depFalse <- function() FALSE

.depIntervalCovsRemove <- paste(
  "the linCmt() carry's end-record fallback for an rxode2 without interval covariates",
  "(rxode2 PR 1517): remove R/foceiLinCmtCarryInterval.R and its calls in",
  ".rxFoceiLinCmtCarryAliasCovs() and .rxFoceiLinCmtCarryCovInterval()",
  "(R/foceiLinCmtCarryCompose.R); the linear refusal, .rxFoceiLinCmtCarryLinearRefused()",
  "and .rxFoceiLinCmtCarryLinearRefuse() (R/foceiLinCmtCarryEligible.R), with their calls",
  "in R/focei.R and .rxFoceiLinCmtCarryBuildEnabled() (R/foceiLinCmtCarryGates.R);",
  ".carryTwinInterps() and .carryInterps() (helper-lincmt-carry.R), and the loop over",
  ".carryInterps() in test-focei-lincmt-carry-theta.R; the probe's skips in",
  "test-focei-lincmt-carry-fit.R, test-focei-lincmt-carry-fit-fallback.R and",
  "test-focei-lincmt-carry-emit-fd.R and its check in test-focei-lincmt-carry-emit.R;",
  "and this block"
)

# The interval-covariate block's warning for the installed rxode2
.depIntervalCovsWarning <- function() {
  .depWorkaroundWarning(
    .rxode2Version(),
    "5.1.8.9000",
    .depIntervalCovsRemove,
    present = .rxFoceiLinCmtIntervalCovs()
  )
}

# Each block: what its warning names to remove, the symbols of its probe,
# skips and fallbacks, and a file the scan must find (anchor), so a scan run
# from the wrong root fails instead of passing on what little it saw
.depBlocks <- list(
  intervalCovs = list(
    remove = .depIntervalCovsRemove,
    anchor = "foceiLinCmtCarryInterval.R",
    symbols = c(
      ".rxFoceiLinCmtIntervalCovs",
      ".rxFoceiLinCmtIntervalProbe",
      ".carryIntervalEnv",
      ".carryTwinInterps",
      ".carryInterps",
      ".rxFoceiLinCmtCarryLinearRefused",
      ".rxFoceiLinCmtCarryLinearRefuse"
    )
  )
)

test_that("the work-around warning is given only for a release with the feature", {
  expect_null(.depWorkaroundWarning(numeric_version("5.1.8"), "5.1.8.9000", "x", present = TRUE))
  expect_null(.depWorkaroundWarning(numeric_version("5.1.8.9000"), "5.1.8.9000", "x", present = TRUE))
  # a development build after the development version is not a release
  expect_null(.depWorkaroundWarning(numeric_version("5.1.8.9001"), "5.1.8.9000", "x", present = TRUE))
  expect_null(.depWorkaroundWarning(numeric_version("5.1.9.9000"), "5.1.8.9000", "x", present = TRUE))
  # a release without the feature still needs the work-around
  expect_null(.depWorkaroundWarning(numeric_version("5.1.9"), "5.1.8.9000", "x", present = FALSE))
  expect_null(.depWorkaroundWarning(numeric_version("5.1.9"), "5.1.8.9000", "x", present = NA))
  expect_match(
    .depWorkaroundWarning(numeric_version("5.1.9"), "5.1.8.9000", "x", present = TRUE),
    "rxode2 5.1.9 is a release after 5.1.8.9000.*x; and raise DESCRIPTION to rxode2 \\(>= 5.1.9\\)"
  )
  expect_true(.depIsRelease("5.1.9"))
  expect_true(.depIsRelease("6.0.0.1"))
  expect_false(.depIsRelease("5.1.8.9000"))
})

test_that("each work-around names every file that uses its probe, skips and fallbacks", {
  # tests/ alone (an installed-package check) cannot show what R/ uses
  skip_if_not(
    dir.exists(file.path(testthat::test_path("..", ".."), "R")),
    "no R/ beside tests/: the file-name check needs the source tree"
  )
  for (.b in names(.depBlocks)) {
    .files <- .depFilesUsing(.depBlocks[[.b]]$symbols)
    expect_true(.depBlocks[[.b]]$anchor %in% .files, label = paste(.b, "scan finds its anchor"))
    for (.f in .files) {
      expect_true(grepl(.f, .depBlocks[[.b]]$remove, fixed = TRUE), label = paste(.b, "names", .f))
    }
  }
})

test_that("interval-covariate work-around: warns only for a release with interval covariates", {
  local_mocked_bindings(.rxode2Version = .depVersion519, .rxFoceiLinCmtIntervalCovs = .depFalse)
  expect_null(.depIntervalCovsWarning())
  local_mocked_bindings(.rxFoceiLinCmtIntervalCovs = .depTrue)
  expect_match(.depIntervalCovsWarning(), "rxode2 5.1.9 is a release after 5.1.8.9000", fixed = TRUE)
  local_mocked_bindings(.rxode2Version = .depVersionDev)
  expect_null(.depIntervalCovsWarning())
})

test_that("interval-covariate work-around: the linCmt() carry's end-record fallback", {
  skip_on_cran()
  .w <- .depIntervalCovsWarning()
  if (!is.null(.w)) {
    warning(.w, call. = FALSE)
  }
  succeed()
})
