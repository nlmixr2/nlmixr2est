# Core dispatch scaffold for the nonparametric engines (npag/npb) and their
# mu-referenced sugar (mnpag/inpag, mnpb/inpb).  Always-run (kept out of
# .slowBatches): no model fit, just control construction and S3 registration.

test_that("nonparametric control constructors return valid controls", {
  .npag <- npagControl()
  expect_true(inherits(.npag, "impmapControl"))
  expect_equal(.npag$est, "npag")
  .npb <- npbControl()
  expect_true(inherits(.npb, "impmapControl"))
  expect_equal(.npb$est, "npb")
})

test_that("npag auto grid size scales with the number of dimensions", {
  # default points is auto (NA sentinel), resolved from the eta count in setup
  expect_true(is.na(npagControl()$points))
  expect_equal(npagControl(points = 300L)$points, 300L)   # explicit override survives
  # the rule: floor at the Pmetrics default (2028), then 512 per eta beyond
  .f <- nlmixr2est:::.npAutoPoints
  expect_equal(.f(1L), 2028L)      # low dim -> floor
  expect_equal(.f(3L), 2028L)      # theo (3 etas) -> matches Pmetrics
  expect_equal(.f(4L), 2048L)      # just past the floor
  expect_equal(.f(8L), 4096L)      # warfarin (8 etas) -> denser high-dim grid
  expect_true(.f(20L) > .f(8L))    # monotone increasing
  # the NA sentinel survives the control validator round-trip
  expect_true(is.na(nlmixr2est:::.npValidCtl(list(npagControl()), "npag")$points))
})

test_that("npag/npb and their mu sugar are registered est methods", {
  .methods <- as.character(utils::methods("nlmixr2Est"))
  for (.est in c("npag", "npb", "mnpag", "inpag", "mnpb", "inpb")) {
    expect_true(paste0("nlmixr2Est.", .est) %in% .methods,
                info = paste0("nlmixr2Est.", .est, " registered"))
    # control validator exists (getValidNlmixrCtl.default errors otherwise)
    expect_silent(getValidNlmixrControl(NULL, .est))
  }
})

test_that("mu sugar forces the mu-referencing variant on its control", {
  expect_equal(getValidNlmixrControl(NULL, "mnpag")$muModel, "lin")
  expect_equal(getValidNlmixrControl(NULL, "inpag")$muModel, "irls")
  expect_equal(getValidNlmixrControl(NULL, "mnpb")$muModel, "lin")
  expect_equal(getValidNlmixrControl(NULL, "inpb")$muModel, "irls")
})

test_that("mu sugar methods are flagged as mu-referenced", {
  for (.est in c("mnpag", "inpag", "mnpb", "inpb")) {
    .mu <- attr(utils::getS3method("nlmixr2Est", .est), "mu")
    expect_true(is.function(.mu) && isTRUE(.mu(NULL)),
                info = paste0(.est, " is mu-referenced"))
  }
})

test_that("npb rejects generalized (non-normal) likelihoods; npag supports them", {
  pois <- function() {
    ini({ tlam <- log(2); eta.lam ~ 0.1 })
    model({ lam <- exp(tlam + eta.lam); y ~ pois(lam) })
  }
  norm <- function() {
    ini({ tv <- log(32); add.sd <- 0.7; eta.v ~ 0.1 })
    model({ v <- exp(tv + eta.v); cp <- v; cp ~ add(add.sd) })
  }
  # npb (Gibbs) still rejects a non-normal endpoint ...
  expect_error(.npAssertNormal(rxode2::assertRxUi(pois), "npb"), "generalized")
  expect_error(.npAssertNormal(rxode2::assertRxUi(norm), "npb"), NA)   # normal ok
  # ... but npag now supports generalized / ll() likelihoods (no error)
  expect_error(.npAssertNormal(rxode2::assertRxUi(pois), "npag"), NA)
})
