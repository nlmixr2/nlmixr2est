# saemControl(iovMethod=) decides whether the shared IOV pre-processing rewrite
# runs at all.  The switch is the "iov" attribute on nlmixr2Est.saem, which
# .isIovMethod() already supports as a function(control).
test_that("iovMethod selects the IOV handling", {
  expect_equal(saemControl()$iovMethod, "theta")
  expect_equal(saemControl(iovMethod = "twoLevel")$iovMethod, "twoLevel")
  expect_error(saemControl(iovMethod = "nope"))

  # the shared rewrite applies for "theta" and stands down for "twoLevel"
  expect_true(.isIovMethod("saem", saemControl()))
  expect_true(.isIovMethod("saem", saemControl(iovMethod = "theta")))
  expect_false(.isIovMethod("saem", saemControl(iovMethod = "twoLevel")))

  # every other method keeps the rewrite unconditionally
  expect_true(.isIovMethod("focei", foceiControl()))
})

test_that("the shared rewrite runs, or does not, according to iovMethod", {
  .theoIov <- nlmixr2data::theo_md
  .theoIov$occ <- 1
  .theoIov$occ[.theoIov$TIME >= 144] <- 2

  .mod <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      add.sd <- 0.7
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      iov.cl ~ 0.1 | occ
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl + iov.cl)
      v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }
  .ui <- rxode2::rxUiDecompress(.mod())

  # "theta": iov.cl becomes a magnitude theta plus one unit-variance eta per
  # observed occasion
  .rw <- .uiApplyIov(.ui, "saem", .theoIov, saemControl(iovMethod = "theta"))
  expect_true(is.list(.rw))
  .ini <- .rw$ui$iniDf
  expect_true("iov.cl" %in% .ini$name[is.na(.ini$neta1)])
  expect_true(all(c("rx.iov.cl.1", "rx.iov.cl.2") %in% .ini$name))

  # "twoLevel": nothing is rewritten, so the `| occ` variance component survives
  # into the saem build
  expect_null(.uiApplyIov(.ui, "saem", .theoIov, saemControl(iovMethod = "twoLevel")))
  expect_equal(.ui$iniDf$condition[.ui$iniDf$name == "iov.cl"], "occ")
  expect_true(is.list(.ui$omega))
  expect_true(all(c("id", "occ") %in% names(.ui$omega)))
})
