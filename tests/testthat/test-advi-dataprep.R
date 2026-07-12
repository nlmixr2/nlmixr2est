# .adviDataPrep + .adviClassifyPars: eta<->theta map and the mu-ref vs non-mu
# structural vs sigma classification that drives the ADVI gradient path.
# Always-run core unit test (no fit).

test_that("all-mu-referenced model has an empty theta-sensitivity set", {
  ## every structural theta carries an eta and is mu-referenced -> no theta-sens
  allMu <- function() {
    ini({ tka <- 0.45; tcl <- 1; tv <- 3.45
      eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
      add.sd <- 0.7 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
      d/dt(depot) <- -ka*depot; d/dt(center) <- ka*depot - cl/v*center
      cp <- center/v; cp ~ add(add.sd) })
  }
  ui <- rxode2::rxode2(allMu)
  cls <- .adviClassifyPars(ui)
  expect_length(cls$etaNames, 3L)
  expect_true(all(!is.na(cls$muRefThetaIdx)))   # every eta mu-references a theta
  expect_length(cls$struct, 0L)                 # no non-mu structural thetas
  expect_length(cls$sigma, 1L)                  # add.sd is a sigma theta
  expect_setequal(cls$thetaSensIdx, cls$sigma)

  prep <- .adviDataPrep(ui, nlmixr2data::theo_sd)
  expect_length(prep$muRefThetaIdx, 3L)
  expect_false(any(prep$isFree))
  expect_length(prep$structIdx, 0L)
})

test_that("fixed-effect-only thetas become non-mu structural sensitivities", {
  ## tcl, tv have no eta -> non-mu structural thetas needing theta sensitivities
  mixedMu <- function() {
    ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; eta.ka ~ 0.6; add.sd <- 0.7 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv)
      d/dt(depot) <- -ka*depot; d/dt(center) <- ka*depot - cl/v*center
      cp <- center/v; cp ~ add(add.sd) })
  }
  ui <- rxode2::rxode2(mixedMu)
  cls <- .adviClassifyPars(ui)
  expect_length(cls$etaNames, 1L)
  ## tcl (ntheta 2) and tv (ntheta 3) are non-mu structural
  expect_setequal(cls$struct, c(2L, 3L))
  expect_true(all(cls$struct %in% cls$thetaSensIdx))
  expect_true(all(cls$sigma %in% cls$thetaSensIdx))

  prep <- .adviDataPrep(ui, nlmixr2data::theo_sd)
  expect_gt(length(prep$structIdx), 0L)
  expect_equal(prep$N, length(unique(nlmixr2data::theo_sd$ID)))
  expect_gt(prep$Nobs, 0L)
})
