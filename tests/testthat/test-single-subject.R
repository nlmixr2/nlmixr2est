nmTest({

  # A model whose only random effect is fixed to zero collapses to a
  # fixed-effect / single-subject ("N of 1") model once the zero eta is dropped
  # (see issue #493).
  one.cmt <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      eta.ka ~ fix(0)
      add.err <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl)
      v <- exp(tv)
      linCmt() ~ add(add.err)
    })
  }

  test_that("methods that require random effects give an actionable error", {
    for (.est in c("fo", "foi", "saem", "fsaem")) {
      .err <- tryCatch(.nlmixr(one.cmt, nlmixr2data::theo_sd, .est, list(print = 0)),
                       error = function(e) conditionMessage(e))
      # keeps the original model name (not the internal '.mod') ...
      expect_true(grepl("'one.cmt'", .err, fixed = TRUE))
      # ... and points to methods that can fit a no-random-effect model
      expect_true(grepl("focei", .err, fixed = TRUE))
      expect_true(grepl("nlminb", .err, fixed = TRUE))
    }
  })

  test_that("focei / foce fit a single-subject model", {
    for (.est in c("focei", "foce")) {
      .fit <- .nlmixr(one.cmt, nlmixr2data::theo_sd, .est,
                      list(print = 0, calcTables = FALSE, maxOuterIterations = 0))
      expect_s3_class(.fit, "nlmixr2FitCore")
    }
  })

  test_that("nlm-family fits a single-subject model", {
    for (.est in c("nlminb", "bobyqa")) {
      .fit <- .nlmixr(one.cmt, nlmixr2data::theo_sd, .est, list(print = 0))
      expect_s3_class(.fit, "nlmixr2FitCore")
    }
  })

})
