# saemControl(iovMethod=) decides whether the shared IOV pre-processing rewrite
# runs at all.  The switch is the "iov" attribute on nlmixr2Est.saem, which
# .isIovMethod() already supports as a function(control).
test_that("iovMethod selects the IOV handling", {
  expect_equal(saemControl()$iovMethod, "twoLevel")
  expect_equal(saemControl(iovMethod = "theta")$iovMethod, "theta")
  expect_error(saemControl(iovMethod = "nope"))

  # the shared rewrite applies for "theta" and stands down for "twoLevel"
  expect_false(.isIovMethod("saem", saemControl()))
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

# A FOCEi family method that honours a repeated (`same()`) occasion block must
# also declare the "iov" attribute -- otherwise .uiApplyIov() stands down and
# nothing expands the occasion parameters, which is how est="ifocei" and
# est="mfocei" came to error on every IOV model (#1083).
test_that("every .iovSameMethods method declares the 'iov' attribute", {
  .missing <- .iovSameMethods[
    !vapply(.iovSameMethods, function(.e) .isIovMethod(.e, foceiControl()), logical(1), USE.NAMES = FALSE)
  ]
  expect_equal(.missing, character(0))
  # the two that were missing it, named so a regression is unambiguous
  expect_true(.isIovMethod("ifocei", foceiControl()))
  expect_true(.isIovMethod("mfocei", foceiControl()))
})

test_that("the shared IOV rewrite runs for ifocei/mfocei (#1083)", {
  .d <- nlmixr2data::theo_sd
  .d$occ <- 1 + (.d$TIME >= 5)
  .mod <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      add.sd <- 0.7
      eta.ka ~ 0.6
      iov.cl ~ 0.1 | occ
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + iov.cl)
      v <- exp(tv)
      linCmt() ~ add(add.sd)
    })
  }
  for (.est in c("ifocei", "mfocei", "focei")) {
    .ui <- rxode2::rxUiDecompress(.mod())
    .rw <- .uiApplyIov(.ui, .est, .d, foceiControl())
    expect_true(is.list(.rw), info = .est)
    .ini <- .rw$ui$iniDf
    expect_true(all(c("rx.iov.cl.1", "rx.iov.cl.2") %in% .ini$name), info = .est)
  }
})

# The "*f" (R/foceiFast.R) and "f*" (R/foceiFull.R) methods are thin delegates:
# nlmixr2Est.<delegate> calls nlmixr2Est.<base> with a different control default.
# Their IOV capabilities must therefore match their base's, in BOTH registries --
# the "iov" attribute (.isIovMethod) and .iovSameMethods (.isIovSameMethod).  The
# two drifted apart twice: ifocei/mfocei lacked "iov" (#1083), and the six "full"
# delegates were missing from .iovSameMethods, so est="flaplace" refused a
# correlated occasion block that est="laplace" fits.
test_that("thin delegates match their base method's IOV capabilities", {
  .delegates <- c(
    foceif = "focei",
    focef = "foce",
    focepf = "focep",
    mfoceif = "mfocei",
    mfocef = "mfoce",
    mfocepf = "mfocep",
    ifoceif = "ifocei",
    ifocef = "ifoce",
    ifocepf = "ifocep",
    agqf = "agq",
    magqf = "magq",
    iagqf = "iagq",
    flaplace = "laplace",
    mflaplace = "mlaplace",
    iflaplace = "ilaplace",
    fagq = "agq",
    mfagq = "magq",
    ifagq = "iagq"
  )
  for (.d in names(.delegates)) {
    .b <- .delegates[[.d]]
    expect_equal(.isIovMethod(.d, foceiControl()), .isIovMethod(.b, foceiControl()), info = .d)
    expect_equal(.isIovSameMethod(.d), .isIovSameMethod(.b), info = .d)
  }
})

# A correlated occasion block needs iovMethod="omega" (one estimated block
# repeated per occasion).  est="flaplace" used to refuse it outright while the
# est="laplace" it delegates to fitted the same model.
test_that("a correlated occasion block survives a full-Laplace delegate", {
  .d <- nlmixr2data::theo_sd
  .d$occ <- 1 + (.d$TIME >= 5)
  .mod <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      add.sd <- 0.7
      eta.ka ~ 0.6
      iov.cl + iov.v ~ c(0.1, 0.03, 0.2) | occ
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + iov.cl)
      v <- exp(tv + iov.v)
      linCmt() ~ add(add.sd)
    })
  }
  for (.est in c("laplace", "flaplace", "fagq", "iflaplace")) {
    .ui <- rxode2::rxUiDecompress(.mod())
    .rw <- .uiApplyIov(.ui, .est, .d, foceiControl())
    expect_true(is.list(.rw), info = .est)
    # the repeated block is how the correlation is carried
    expect_true(any(grepl(":same:", .rw$ui$iniDf$condition, fixed = TRUE)), info = .est)
  }
})
