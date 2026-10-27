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
