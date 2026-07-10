# Fisher-score standard errors (design/rpem/08): at the converged estimates RPEM
# forms each subject's marginal-likelihood score from the stored samples (Fisher
# identity, importance-weighted complete-data scores), the empirical Fisher
# information I = sum_i s_i s_i^T, and reports SE = sqrt(diag(I^-1)).  The fit's
# covMethod is "fisher" and the SEs must be finite, close to FOCEI's, and shrink
# ~1/sqrt(n) as n grows.

test_that("RPEM reports Fisher-score SEs (close to FOCEI, shrinking ~1/sqrt(n))", {
  skip_on_cran()
  skip_on_ci()  # heavy: several FOCEI + RPEM fits

  struct <- rxode2::rxode2({ ka <- exp(tka + eta); cl <- exp(tcl); v <- exp(tv); cp <- linCmt() })
  simDat <- function(nsub, seed) {
    set.seed(seed); etasTrue <- rnorm(nsub, 0, sqrt(0.3)); obsT <- seq(0.5, 24, by = 1.5)
    do.call(rbind, lapply(seq_len(nsub), function(i) {
      ev <- rxode2::et(amt = 100, cmt = "depot"); ev <- rxode2::et(ev, obsT)
      s <- rxode2::rxSolve(struct, params = c(tka = 0.45, tcl = 1, tv = 3.45, eta = etasTrue[i]),
                           events = ev, returnType = "data.frame", addDosing = FALSE)
      d <- as.data.frame(ev); d$id <- i; o <- d$evid == 0; d$DV <- 0
      d$DV[o] <- s$cp + rnorm(nrow(s), 0, 0.1)
      d
    }))
  }
  mod <- function() {
    ini({ tka <- 0.3; tcl <- fix(1.0); tv <- fix(3.45); add.sd <- 0.2; eta.ka ~ 0.6 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv); cp <- linCmt(); cp ~ add(add.sd) })
  }
  ctl <- rpemControl(nGauss = 300L, nMH = 60000L, mhBurn = 6000L, niter = 30L,
                     collect = 12L, seed = 100L, cores = 4L)

  dat <- simDat(40L, 8)
  rf <- suppressMessages(nlmixr2(mod, dat, est = "rpem", control = ctl))
  ff <- suppressMessages(nlmixr2(mod, dat, est = "focei",
                                 control = foceiControl(print = 0, calcTables = FALSE)))

  # the fit records the Fisher-score covariance method
  expect_match(as.character(rf$covMethod), "fisher")
  # SEs are finite and positive for the estimated parameters
  .se <- rf$parFixedDf[c("tka", "add.sd"), "SE"]
  expect_true(all(is.finite(.se) & .se > 0))
  # and close to FOCEI's covariance SEs
  .fse <- ff$parFixedDf[c("tka", "add.sd"), "SE"]
  expect_equal(unname(.se[1]), unname(.fse[1]), tolerance = 0.25)   # tka SE
  expect_equal(unname(.se[2]), unname(.fse[2]), tolerance = 0.25)   # add.sd SE
  # the omega-variance SE is reported in $cov (named om.<eta>), finite and positive
  expect_true("om.eta.ka" %in% rownames(rf$cov))
  expect_true(is.finite(rf$cov["om.eta.ka", "om.eta.ka"]) &&
                rf$cov["om.eta.ka", "om.eta.ka"] > 0)

  # SE shrinks ~1/sqrt(n): 4x subjects -> ~2x smaller tka SE
  rf4 <- suppressMessages(nlmixr2(mod, simDat(160L, 9), est = "rpem", control = ctl))
  .ratio <- rf$parFixedDf["tka", "SE"] / rf4$parFixedDf["tka", "SE"]
  expect_gt(.ratio, 1.4)
  expect_lt(.ratio, 2.6)
})

test_that("RPEM reports multi-eta Fisher-score SEs (diagonal Omega)", {
  skip_on_cran()
  skip_on_ci()  # heavy: FOCEI + RPEM fits

  sim <- rxode2::rxode2({ ka <- exp(tka + eka); cl <- exp(tcl + ecl); v <- exp(tv + ev)
                          cp <- linCmt() })
  set.seed(3); nsub <- 100L; obsT <- seq(0.25, 24, by = 1.5); omT <- c(0.3, 0.1, 0.06)
  dat <- do.call(rbind, lapply(seq_len(nsub), function(i) {
    e <- c(rnorm(1, 0, sqrt(omT[1])), rnorm(1, 0, sqrt(omT[2])), rnorm(1, 0, sqrt(omT[3])))
    ev <- rxode2::et(amt = 100, cmt = "depot"); ev <- rxode2::et(ev, obsT)
    s <- rxode2::rxSolve(sim, params = c(tka = 0.45, tcl = 1, tv = 3.45,
                                         eka = e[1], ecl = e[2], ev = e[3]),
                         events = ev, returnType = "data.frame", addDosing = FALSE)
    d <- as.data.frame(ev); d$id <- i; o <- d$evid == 0; d$DV <- 0
    d$DV[o] <- s$cp + rnorm(nrow(s), 0, 0.1)
    d
  }))
  mod <- function() {
    ini({ tka <- 0.3; tcl <- 0.9; tv <- 3.3; add.sd <- 0.2
          eta.ka ~ 0.3; eta.cl ~ 0.1; eta.v ~ 0.05 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
            cp <- linCmt(); cp ~ add(add.sd) })
  }
  ff <- suppressMessages(nlmixr2(mod, dat, "focei", control = foceiControl(print = 0, calcTables = FALSE)))
  rf <- suppressMessages(nlmixr2(mod, dat, "rpem",
                                 control = rpemControl(nGauss = 300L, nMH = 60000L, mhBurn = 6000L,
                                                       niter = 35L, collect = 12L, seed = 100L, cores = 4L)))

  expect_match(as.character(rf$covMethod), "fisher")
  # a typical-value SE per eta + one omega-variance SE per eta, all finite/positive
  .th <- c("tka", "tcl", "tv")
  .se <- rf$parFixedDf[.th, "SE"]
  expect_true(all(is.finite(.se) & .se > 0))
  .om <- c("om.eta.ka", "om.eta.cl", "om.eta.v")
  expect_true(all(.om %in% rownames(rf$cov)))
  expect_true(all(is.finite(diag(rf$cov)[.om]) & diag(rf$cov)[.om] > 0))
  # the well-identified structural SEs agree with FOCEI at this sample size
  expect_equal(unname(rf$parFixedDf["tcl", "SE"]), unname(ff$parFixedDf["tcl", "SE"]), tolerance = 0.3)
  expect_equal(unname(rf$parFixedDf["tv",  "SE"]), unname(ff$parFixedDf["tv",  "SE"]), tolerance = 0.3)
})
