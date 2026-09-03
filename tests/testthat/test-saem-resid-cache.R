# The residual M-step caches the transformed predictions/observations it scores
# (ensureSaemFixedTransformCache, src/saem.cpp).  That cache used to be keyed on
# the ADDRESS of the ysb/fsb buffers, which are per-iteration locals: freed and
# reallocated every M-step, so the allocator routinely hands back the same
# address with different CONTENTS.  The guard then reported "unchanged" and the
# optimizer scored later iterations against the FIRST iteration's predictions --
# which still carry all the between-subject variability as error, so the
# residual stayed near its warm start and the variance components collapsed to
# compensate.  Whether the address is reused is not deterministic, so the same
# fit gave different answers from one run to the next.
#
# Only residual models that need the optimizer are affected; pure add() and
# pure prop() have closed-form M-steps.

test_that("saem repeats a combined-error fit exactly (residual transform cache)", {
  skip_on_cran()

  .d <- nlmixr2data::theo_md

  # nolint start: object_usage_linter. rxode2 ini()/model() are NSE blocks:
  # every assignment here is model specification, consumed by rxode2, not a
  # local variable lintr can see used.
  .mod <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      add.sd <- 0.7
      prop.sd <- 0.1
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      cp <- linCmt()
      cp ~ add(add.sd) + prop(prop.sd) + combined1()
    })
  }
  # nolint end

  .ctl <- saemControl(nBurn = 60, nEm = 60, seed = 42L, print = 0L,
                      covMethod = "", calcTables = FALSE)
  .f1 <- suppressWarnings(nlmixr2(.mod(), .d, est = "saem", control = .ctl))
  .f2 <- suppressWarnings(nlmixr2(.mod(), .d, est = "saem", control = .ctl))

  # same model, same data, same seed, same session: bit-identical
  expect_identical(.f1$fixef[["add.sd"]], .f2$fixef[["add.sd"]])
  expect_identical(.f1$fixef[["prop.sd"]], .f2$fixef[["prop.sd"]])
  expect_identical(.f1$fixef[["tka"]], .f2$fixef[["tka"]])
  expect_identical(.f1$fixef[["tcl"]], .f2$fixef[["tcl"]])
  expect_identical(.f1$fixef[["tv"]], .f2$fixef[["tv"]])
  .o1 <- if (is.list(.f1$omega)) .f1$omega$id else .f1$omega
  .o2 <- if (is.list(.f2$omega)) .f2$omega$id else .f2$omega
  expect_identical(diag(.o1), diag(.o2))
})
