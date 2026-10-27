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

  test_that("an IOV fit under nocb or midpoint steps its occasions as under locf", {
    .locf <- .iovLocfObjf("locf")
    expect_equal(.iovLocfObjf("nocb"), .locf, tolerance = 1e-8)
    expect_equal(.iovLocfObjf("midpoint"), .locf, tolerance = 1e-8)
  })
})
