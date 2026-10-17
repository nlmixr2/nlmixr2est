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
  expect_true(all(!is.na(cls$muRefThetaIdx))) # every eta mu-references a theta
  expect_length(cls$struct, 0L) # no non-mu structural thetas
  expect_length(cls$sigma, 1L) # add.sd is a sigma theta
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
})

# a model with an eta on every structural parameter, for the data tests below
.vdMod <- function() {
  ini({
    tka <- 0.45
    tcl <- 1
    tv <- 3.45
    eta.ka ~ 0.6
    eta.cl ~ 0.3
    eta.v ~ 0.1
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    d/dt(depot) <- -ka * depot
    d/dt(center) <- ka * depot - cl / v * center
    cp <- center / v
    cp ~ add(add.sd)
  })
}

test_that(".adviDataPrep drops a subject with no observations, as the FOCEi inner problem does", {
  ui <- rxode2::rxode2(.vdMod)
  .ids <- unique(nlmixr2data::theo_sd$ID)
  # neither a missing DV nor MDV = 1 makes a record an observation
  .noDv <- nlmixr2data::theo_sd
  .noDv$DV[.noDv$ID == 12 & .noDv$EVID == 0] <- NA
  .mdv <- nlmixr2data::theo_sd
  .mdv$MDV <- as.integer(.mdv$EVID != 0 | .mdv$ID == 12)
  for (.d in list(.noDv, .mdv)) {
    expect_warning(.p <- .adviDataPrep(ui, .d), "IDs without observations dropped: 12", fixed = TRUE)
    expect_identical(.p$N, 11L)
    expect_identical(.p$ids, .ids[.ids != 12])
  }
  .d <- nlmixr2data::theo_sd
  .d$ID[5] <- NA
  expect_error(.adviDataPrep(ui, .d), "^'id' cannot have NA values$")
})

test_that(".adviDataPrep takes an MDV column as the EVID when there is none, as FOCEi does", {
  ui <- rxode2::rxode2(.vdMod)
  .d <- nlmixr2data::theo_sd
  .d$MDV <- as.integer(.d$EVID != 0)
  .d$EVID <- NULL
  # MDV 0 keeps a record an observation even with an AMT; from AMT alone every
  # record of subject 12 would be a dose and the subject dropped
  .d$AMT[.d$ID == 12 & .d$MDV == 0] <- 50
  .p <- .adviDataPrep(ui, .d)
  expect_identical(.p$N, 12L)
  expect_identical(.p$ids, unique(.d$ID))
})

test_that(".adviDataPrep takes dose-free data with no AMT or EVID column", {
  ui <- rxode2::rxode2(.vdMod)
  .d <- nlmixr2data::theo_sd
  .d <- .d[.d$EVID == 0, c("ID", "TIME", "DV")]
  .p <- .adviDataPrep(ui, .d)
  expect_identical(.p$N, 12L)
  expect_identical(.p$ids, unique(.d$ID))
})
