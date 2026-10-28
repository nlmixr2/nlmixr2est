# The IOV rewrite reads an occasion variable as a covariate (`occ == 1`).  It
# declares locf() for it, so an occasion steps at its first record whatever
# covsInterpolation is, as rxode2 steps the occasion variables of a nested
# omega.

# an IOV model, with `first` as the first line of its model block
.iovLocfMod <- function(first = "") {
  eval(parse(
    text = sprintf(
      "function() {
        ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; eta.ka ~ 0.6; iov.cl ~ 0.1 | occ; add.sd <- 0.7 })
        model({
          %s
          ka <- exp(tka + eta.ka)
          cl <- exp(tcl + iov.cl)
          v <- exp(tv)
          d/dt(depot) <- -ka * depot
          d/dt(center) <- ka * depot - cl / v * center
          cp <- center / v
          cp ~ add(add.sd)
        })
      }",
      first
    )
  ))
}

.iovLocfData <- function() {
  .d <- nlmixr2data::theo_sd
  .d$occ <- 1 + (.d$TIME >= 5)
  .d
}

.iovLocfApply <- function(mod, ...) {
  .uiApplyIov(rxode2::rxode2(mod), "focei", .iovLocfData(), foceiControl(...))$ui
}

# the number of locf(occ) lines of a ui's model
.iovLocfCount <- function(ui) {
  sum(vapply(ui$lstExpr, identical, logical(1), quote(locf(occ))))
}

# the objective of a focei fit at the initial estimates under `ci`
.iovLocfObjf <- function(ci) {
  .fit <- suppressWarnings(suppressMessages(nlmixr2(
    .iovLocfMod(),
    .iovLocfData(),
    "focei",
    control = foceiControl(
      print = 0,
      maxOuterIterations = 0L,
      covMethod = "",
      rxControl = rxode2::rxControl(covsInterpolation = ci, atol = 1e-10, rtol = 1e-10)
    )
  )))
  .fit$objf
}

# an IOV model from its ini and model statements
.iovLocfModOf <- function(ini, body) {
  eval(parse(text = sprintf("function() {\n ini({ %s })\n model({ %s })\n}", ini, body)))
}

# a focei fit at the initial estimates under `ci`
.iovLocfFit <- function(mod, data, ci) {
  suppressWarnings(suppressMessages(nlmixr2(
    mod,
    data,
    "focei",
    control = foceiControl(
      print = 0,
      maxOuterIterations = 0L,
      covMethod = "",
      rxControl = rxode2::rxControl(covsInterpolation = ci, atol = 1e-10, rtol = 1e-10)
    )
  )))
}

# a short saem fit with `iovMethod` under `ci`
.iovLocfSaem <- function(mod, iovMethod, ci) {
  suppressWarnings(suppressMessages(nlmixr2(
    mod,
    .iovLocfData(),
    "saem",
    control = saemControl(
      print = 0,
      nBurn = 2,
      nEm = 2,
      nmc = 1,
      seed = 42,
      iovMethod = iovMethod,
      rxControl = rxode2::rxControl(covsInterpolation = ci)
    )
  )))
}

.iovLocfOde <- "d/dt(depot) <- -ka * depot\n d/dt(center) <- ka * depot - cl / v * center\n cp <- center / v\n cp ~ add(add.sd)"

nmTest({
  test_that("the IOV rewrite declares locf() for each occasion variable", {
    for (.m in c("theta", "omega")) {
      .ui <- .iovLocfApply(.iovLocfMod(), iovMethod = .m)
      expect_identical(.ui$lstExpr[[1]], quote(locf(occ)), label = .m)
      expect_identical(.iovLocfCount(.ui), 1L, label = .m)
      expect_identical(as.character(.ui$mv0$interp[["occ"]]), "locf", label = .m)
    }
    # a model that declares it already keeps its one line
    .own <- .iovLocfMod("locf(occ)")
    expect_identical(.iovLocfCount(rxode2::rxode2(.own)), 1L)
    expect_identical(.iovLocfCount(.iovLocfApply(.own)), 1L)
  })

  test_that("an occasion variable declared with another interpolation stops", {
    # rxode2 from 5.1.8.9000 refuses the model itself, with the same message
    .e <- expect_error(.iovLocfApply(.iovLocfMod("nocb(occ)")))
    expect_identical(conditionMessage(.e), "occasion variable(s) 'occ' (nocb) must use 'locf' interpolation")
  })

  test_that("the fit's model keeps a user's own locf(occ) and not the rewrite's", {
    expect_identical(.iovLocfCount(.iovLocfFit(.iovLocfMod(), .iovLocfData(), "nocb")$ui), 0L)
    expect_identical(.iovLocfCount(.iovLocfFit(.iovLocfMod("locf(occ)"), .iovLocfData(), "nocb")$ui), 1L)
  })

  test_that("two occasion variables each get a locf() line and step", {
    .m <- .iovLocfModOf(
      "tka <- 0.45; tcl <- 1; tv <- 3.45; eta.ka ~ 0.6; iov.cl ~ 0.1 | occ; iov.v ~ 0.1 | occ2; add.sd <- 0.7",
      paste("ka <- exp(tka + eta.ka)\n cl <- exp(tcl + iov.cl)\n v <- exp(tv + iov.v)\n", .iovLocfOde)
    )
    .d <- .iovLocfData()
    .d$occ2 <- 1 + (.d$TIME >= 8)
    .ui <- .uiApplyIov(rxode2::rxode2(.m), "focei", .d, foceiControl())$ui
    expect_identical(.ui$lstExpr[1:2], list(quote(locf(occ)), quote(locf(occ2))))
    .locf <- .iovLocfFit(.m, .d, "locf")$objf
    expect_equal(.iovLocfFit(.m, .d, "nocb")$objf, .locf, tolerance = 1e-8)
  })

  test_that("an occasion the model also reads as a covariate steps there too", {
    .m <- .iovLocfModOf(
      "tka <- 0.45; tcl <- 1; tv <- 3.45; eta.ka ~ 0.6; iov.cl ~ 0.1 | occ; add.sd <- 0.7",
      paste("ka <- exp(tka + eta.ka)\n cl <- exp(tcl + iov.cl + 0.1 * occ)\n v <- exp(tv)\n", .iovLocfOde)
    )
    expect_identical(.iovLocfCount(.uiApplyIov(rxode2::rxode2(.m), "focei", .iovLocfData(), foceiControl())$ui), 1L)
    .locf <- .iovLocfFit(.m, .iovLocfData(), "locf")$objf
    expect_equal(.iovLocfFit(.m, .iovLocfData(), "nocb")$objf, .locf, tolerance = 1e-8)
  })

  test_that("a user's own locf(occ) beside nocb(WT) keeps both", {
    .m <- .iovLocfModOf(
      "tka <- 0.45; tcl <- 1; tv <- 3.45; eta.ka ~ 0.6; iov.cl ~ 0.1 | occ; add.sd <- 0.7",
      paste(
        "locf(occ)\n nocb(WT)\n ka <- exp(tka + eta.ka)\n cl <- exp(tcl + iov.cl + 0.01 * (WT - 70))\n v <- exp(tv)\n",
        .iovLocfOde
      )
    )
    .d <- .iovLocfData()
    .d$WT <- .d$WT + 5 * (.d$TIME >= 10)
    .f <- .iovLocfFit(.m, .d, "locf")
    expect_identical(.iovLocfCount(.f$ui), 1L)
    expect_true(any(vapply(.f$ui$lstExpr, identical, logical(1), quote(nocb(WT)))))
    # WT steps at the next record and occ at its own, under either global setting
    expect_equal(.iovLocfFit(.m, .d, "nocb")$objf, .f$objf, tolerance = 1e-8)
  })

  test_that("saem's twoLevel and collapsed expansions step the occasion too", {
    .m <- .iovLocfModOf(
      "tka <- 0.45; tcl <- 1; tv <- 3.45; eta.ka ~ 0.6; eta.cl ~ 0.3; iov.cl ~ 0.1 | occ; add.sd <- 0.7",
      paste("ka <- exp(tka + eta.ka)\n cl <- exp(tcl + eta.cl + iov.cl)\n v <- exp(tv)\n", .iovLocfOde)
    )
    for (.meth in c("twoLevel", "collapsed")) {
      .f <- .iovLocfSaem(.m, .meth, "locf")
      # the expansion's locf(occ) is not in the fit's model
      expect_identical(.iovLocfCount(.f$ui), 0L, label = .meth)
      for (.ci in c("nocb", "midpoint")) {
        .g <- .iovLocfSaem(.m, .meth, .ci)
        expect_equal(.g$IPRED, .f$IPRED, tolerance = 1e-10, label = paste(.meth, .ci))
        expect_equal(.g$objf, .f$objf, tolerance = 1e-10, label = paste(.meth, .ci))
      }
    }
  })

  test_that("an IOV fit under nocb or midpoint steps its occasions as under locf", {
    .locf <- .iovLocfObjf("locf")
    expect_equal(.iovLocfObjf("nocb"), .locf, tolerance = 1e-8)
    expect_equal(.iovLocfObjf("midpoint"), .locf, tolerance = 1e-8)
  })
})
