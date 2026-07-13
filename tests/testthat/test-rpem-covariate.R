# mu2 covariates (D22): RPEM estimates a covariate coefficient on a mu-referenced
# (eta) param via the regression M-step, and must match FOCEI.

test_that("RPEM estimates a mu2 covariate coefficient (matches FOCEI)", {
  skip_on_cran()

  struct <- rxode2::rxode2({
    ka <- exp(tka + kawt * WT + eta); cl <- exp(tcl); v <- exp(tv); cp <- linCmt()
  })
  set.seed(5)
  nsub <- 40L; trueTka <- 0.4; trueKawt <- 0.5; addSd <- 0.15
  obsT <- seq(0.25, 24, by = 1.5)
  WT <- round(rnorm(nsub, 0, 1), 3)                 # centered covariate
  etasTrue <- rnorm(nsub, 0, sqrt(0.25))
  dat <- do.call(rbind, lapply(seq_len(nsub), function(i) {
    ev <- rxode2::et(amt = 100, cmt = "depot"); ev <- rxode2::et(ev, obsT)
    s <- rxode2::rxSolve(struct,
                         params = c(tka = trueTka, kawt = trueKawt, tcl = 1, tv = 3.45,
                                    WT = WT[i], eta = etasTrue[i]),
                         events = ev, returnType = "data.frame", addDosing = FALSE,
                         atol = 1e-8, rtol = 1e-8)
    d <- as.data.frame(ev); d$id <- i; d$WT <- WT[i]; d$DV <- 0
    o <- d$evid == 0; d$DV[o] <- s$cp + rnorm(nrow(s), 0, addSd); d
  }))

  # FOCEI reference (hold tcl, tv fixed as RPEM does)
  mod <- function() {
    ini({ tka <- 0.2; kawt <- 0.3; tcl <- fix(1.0); tv <- fix(3.45); add.sd <- 0.3; eta.ka ~ 1.0 })
    model({ ka <- exp(tka + kawt * WT + eta.ka); cl <- exp(tcl); v <- exp(tv)
            cp <- linCmt(); cp ~ add(add.sd) })
  }
  fit <- suppressMessages(nlmixr2(mod, dat, est = "focei",
                                  control = foceiControl(print = 0, calcTables = FALSE)))
  fTka <- fit$parFixedDf["tka", "Estimate"]; fKawt <- fit$parFixedDf["kawt", "Estimate"]

  # RPEM (no fix(); tcl/tv held at ini)
  rmod <- function() {
    ini({ tka <- 0.2; kawt <- 0.3; tcl <- fix(1.0); tv <- fix(3.45); add.sd <- 0.3; eta.ka ~ 0.5 })
    model({ ka <- exp(tka + kawt * WT + eta.ka); cl <- exp(tcl); v <- exp(tv)
            cp <- linCmt(); cp ~ add(add.sd) })
  }
  ui <- rxode2::rxode2(rmod)
  rf <- .rpemFit(ui, dat, rpemControl(nGauss = 300L, nMH = 80000L, mhBurn = 8000L,
                                      niter = 30L, collect = 12L, seed = 2000L))

  expect_true("kawt" %in% names(rf$covCoef))
  expect_equal(unname(rf$mu["tka"]), fTka, tolerance = 0.05)
  expect_equal(unname(rf$covCoef["kawt"]), fKawt, tolerance = 0.05)
})

# Time-varying vs non-time-varying covariate split (same detection mechanism as
# the mu-referenced SAEM family): a time-varying covariate cannot be a per-subject
# mu regressor, so it is routed to the structural beta re-solve; a non-time-varying
# covariate stays in the mu2 regression.

test_that("RPEM classifier routes time-varying covariates to the structural set", {
  set.seed(3); nsub <- 8L; obsT <- c(0.5, 2, 6, 12)
  dat <- do.call(rbind, lapply(seq_len(nsub), function(i) {
    ntv <- round(rnorm(1), 2)
    rbind(data.frame(id = i, time = 0, evid = 1, amt = 100, cmt = 1, DV = 0, NTV = ntv, TV = 0),
          data.frame(id = i, time = obsT, evid = 0, amt = 0, cmt = 1, DV = 1, NTV = ntv,
                     TV = round(rnorm(length(obsT)), 2)))
  }))
  mod <- function() {
    ini({ tka <- 0.3; tcl <- fix(1.0); lv <- fix(3.45); b_tv <- 0.1; b_ntv <- 0.1
          add.sd <- 0.2; eta.ka ~ 0.5 })
    model({ ka <- exp(tka + b_tv * TV + b_ntv * NTV + eta.ka); cl <- exp(tcl); v <- exp(lv)
            cp <- linCmt(); cp ~ add(add.sd) })
  }
  ui <- rxode2::rxode2(mod)
  tv <- .nlmixrTimeVaryingCovariates(dat, ui, rxode2::rxControl())
  expect_true("TV" %in% tv)         # TV varies within subject -> detected
  expect_false("NTV" %in% tv)       # NTV is constant per subject -> not detected

  cl <- .rpemClassify(ui, tv)
  # non-time-varying stays in the mu regression
  expect_true("b_ntv" %in% cl$covCoefNames)
  expect_false("b_tv" %in% cl$covCoefNames)
  # time-varying coefficient falls into the structural (beta re-solve) set
  expect_true("b_tv" %in% cl$thetaNames[cl$structIdx + 1L])

  # without the split both would be mu regressors (regression M-step could not
  # represent a within-subject covariate)
  cl0 <- .rpemClassify(ui, character(0))
  expect_true(all(c("b_tv", "b_ntv") %in% cl0$covCoefNames))
})

test_that("RPEM recovers time-varying (structural) and non-time-varying (mu) covariates", {
  skip_on_cran()
  simMod <- rxode2::rxode2({
    ka <- exp(0.45 + 0.3 * NTV + 0.25 * TV + eka); cl <- exp(1); v <- exp(3.45); cp <- linCmt()
  })
  set.seed(5); nsub <- 70L; obsT <- seq(0.5, 24, by = 1.5)
  dat <- do.call(rbind, lapply(seq_len(nsub), function(i) {
    ntv <- rnorm(1, 0, 1); eka <- rnorm(1, 0, sqrt(0.15)); n <- length(obsT)
    ev <- data.frame(id = i, time = c(0, obsT), evid = c(1, rep(0, n)), amt = c(100, rep(0, n)),
                     cmt = 1, NTV = ntv, TV = c(0, rnorm(n, 0, 1)), eka = eka)
    s <- rxode2::rxSolve(simMod, ev, returnType = "data.frame", addDosing = FALSE)
    ev$DV <- 0; o <- ev$evid == 0; ev$DV[o] <- s$cp + rnorm(sum(o), 0, 0.1)
    ev$eka <- NULL; ev
  }))
  rmod <- function() {
    ini({ tka <- 0.3; tcl <- fix(1.0); lv <- fix(3.45); b_ntv <- 0.1; b_tv <- 0.1
          add.sd <- 0.2; eta.ka ~ 0.3 })
    model({ ka <- exp(tka + b_ntv * NTV + b_tv * TV + eta.ka); cl <- exp(tcl); v <- exp(lv)
            cp <- linCmt(); cp ~ add(add.sd) })
  }
  rf <- .rpemFit(rxode2::rxode2(rmod),
                 dat, rpemControl(nGauss = 300L, nMH = 60000L, mhBurn = 6000L, niter = 25L,
                                  collect = 10L, seed = 1L))
  expect_equal(unname(rf$covCoef["b_ntv"]), 0.3, tolerance = 0.2)   # mu regression
  expect_equal(unname(rf$struct["b_tv"]), 0.25, tolerance = 0.2)    # structural beta
})
