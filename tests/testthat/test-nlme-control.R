test_that("nlmeControl() accepts the common 'print' alias (issue: list(print=) with est='nlme')", {
  # nlme prints through 'verbose'; 'print' is the shared nlmixr control alias, so
  # a generic list(print=0) must not error but map to quiet output.
  expect_false(nlmeControl(print=0)$verbose)
  expect_true(nlmeControl(print=1)$verbose)
  expect_true(nlmeControl(print=5)$verbose)

  # default and an explicit 'verbose' (without 'print') are unchanged
  expect_true(nlmeControl()$verbose)
  expect_false(nlmeControl(verbose=FALSE)$verbose)

  # 'print' is validated and genuinely unknown args are still rejected
  expect_error(nlmeControl(print=-1))
  expect_error(nlmeControl(bogus=1), "unused argument")

  # the shared control validator (what nlmixr2(..., 'nlme', list(...)) uses) now
  # builds a valid nlmeControl from a bare list(print=)
  .ctl <- getValidNlmixrControl(list(print=0), "nlme")
  expect_s3_class(.ctl, "nlmeControl")
  expect_false(.ctl$verbose)
})
