# mu2 covariates (D22): RPEM estimates a covariate coefficient on a mu-referenced
# (eta) param via the regression M-step, and must match FOCEI.

test_that("RPEM estimates a mu2 covariate coefficient (matches FOCEI)", {
  skip_on_cran()
  skip_on_ci()  # heavy: FOCEI fit + multi-iteration RPEM loop

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
