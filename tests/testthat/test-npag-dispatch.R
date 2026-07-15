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
