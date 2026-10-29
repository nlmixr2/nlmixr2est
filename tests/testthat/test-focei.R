# x with a fastmatch hash attribute, as fastmatch::fmatch() leaves on names
.stripHashed <- function(x) {
  attr(x, ".match.hash") <- "hash"
  x
}

test_that(".stripFastmatchItem() keeps a list's own attributes", {
  .m <- matrix(1, 1, 1, dimnames = list(.stripHashed("eta.v"), .stripHashed("eta.v")))
  .clean <- matrix(1, 1, 1, dimnames = list("eta.v", "eta.v"))
  # a nested omega is a "lotri" list, which rxode2 reads by its class
  .om <- structure(list(id = .m, occ = .m), class = c("lotri", "lotriFix", "list"))
  .ret <- .stripFastmatchItem(.om)
  expect_identical(class(.ret), c("lotri", "lotriFix", "list"))
  expect_identical(names(.ret), c("id", "occ"))
  expect_identical(.ret$occ, .clean)
  # a plain list of matrices is stripped as before
  expect_identical(.stripFastmatchItem(list(a = .m)), list(a = .clean))
  expect_identical(.stripFastmatchItem(.m), .clean)
})

# a .foceiLagTxt() that stops at the first d(pred)/d(eta) row, and one that
# stops at the first d(R)/d(eta) row, as a row of the build that fails does
.foceiOrigLagTxt <- .foceiLagTxt
.foceiStopPredRow <- function(lag, d, dfe) {
  stop("the d(pred)/d(eta) row stopped", call. = FALSE)
}
.foceiStopRRow <- function(lag, d, dfe) {
  if (grepl("rx_r__BY_", dfe, fixed = TRUE)) {
    stop("the d(R)/d(eta) row stopped", call. = FALSE)
  }
  .foceiOrigLagTxt(lag, d, dfe)
}

.foceiProgressModel <- function() {
  ini({
    tka <- 0.45
    tcl <- 1
    tv <- 3.45
    eta.ka ~ 0.6
    eta.cl ~ 0.3
    prop.sd <- 0.1
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv)
    d / dt(depot) <- -ka * depot
    d / dt(center) <- ka * depot - cl / v * center
    cp <- center / v
    cp ~ prop(prop.sd)
  })
}

# both builds show a progress bar; the error that stops a row reaches the
# caller, not the bar's generic abort message
test_that("an error in the d(pred)/d(eta) build is the one raised", {
  .ui <- rxode2::rxode2(.foceiProgressModel)
  local_mocked_bindings(.foceiLagTxt = .foceiStopPredRow)
  expect_error(
    suppressMessages(rxUiGet.foceiHdEta(list(.ui))),
    "the d(pred)/d(eta) row stopped",
    fixed = TRUE
  )
})

test_that("an error in the d(R)/d(eta) build is the one raised", {
  .ui <- rxode2::rxode2(.foceiProgressModel)
  local_mocked_bindings(.foceiLagTxt = .foceiStopRRow)
  expect_error(
    suppressMessages(rxUiGet.foceiEnv(list(.ui))),
    "the d(R)/d(eta) row stopped",
    fixed = TRUE
  )
})
