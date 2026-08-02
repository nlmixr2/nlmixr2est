# The "*f" convenience methods: each is its base method with fast=TRUE default.
#
# Registration / control checks only, so this stays in the push/PR subset.  The
# fit-based checks live in test-focei-fast-methods-fit.R (weekly-batched).

nmTest({
  test_that("all nine *f methods are registered and default to fast=TRUE", {
    .fm <- c("focef", "focepf", "foceif", "mfocef", "mfocepf", "mfoceif",
             "ifocef", "ifocepf", "ifoceif")
    expect_true(all(.fm %in% nlmixr2AllEst()))
    # each *f control forces fast=TRUE, even from an empty control
    for (m in .fm) {
      .ctl <- getValidNlmixrCtl(structure(list(NULL), class = c(m, "getValidNlmixrControl")))
      expect_true(isTRUE(.ctl$fast), info = m)
    }
    # a derivative-free outerOpt still downgrades fast (runs in the base constructor)
    expect_warning(
      .ctl <- getValidNlmixrCtl(structure(list(foceiControl(outerOpt = "bobyqa")),
                                          class = c("foceif", "getValidNlmixrControl"))),
      "derivative-free")
    expect_false(isTRUE(.ctl$fast))
  })
})
