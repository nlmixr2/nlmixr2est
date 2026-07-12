# Censoring / BLQ (M3): below-limit observations contribute the CENS probability
# instead of the density.  The E-step uses the censored normal log-likelihood
# (doCensNormal1) and the residual M-step maximizes it (design/rpem/).  Must match
# FOCEI's M3, and the fit object must report the censoring type that was performed.

test_that("RPEM supports M3 BLQ censoring (matches FOCEI) and reports the censoring type", {
  skip_on_cran()
  skip_on_ci()  # heavy: FOCEI fit + multi-iteration RPEM loop

  struct <- rxode2::rxode2({ ka <- exp(tka + eta); cl <- exp(tcl); v <- exp(tv); cp <- linCmt() })
  set.seed(7); nsub <- 40L; obsT <- seq(0.5, 24, by = 1.5); addT <- 0.15; loq <- 0.5
  etasTrue <- rnorm(nsub, 0, sqrt(0.3))
  dat <- do.call(rbind, lapply(seq_len(nsub), function(i) {
    ev <- rxode2::et(amt = 100, cmt = "depot"); ev <- rxode2::et(ev, obsT)
    s <- rxode2::rxSolve(struct, params = c(tka = 0.45, tcl = 1, tv = 3.45, eta = etasTrue[i]),
                         events = ev, returnType = "data.frame", addDosing = FALSE)
    d <- as.data.frame(ev); d$id <- i; o <- d$evid == 0; d$DV <- 0
    d$DV[o] <- s$cp + rnorm(nrow(s), 0, addT)
    d$cens <- 0L; bl <- o & d$DV < loq; d$cens[bl] <- 1L; d$DV[bl] <- loq   # BLQ -> DV = LOQ
    d
  }))
  expect_gt(sum(dat$cens), 30)   # a meaningful fraction censored

  mod <- function() {
    ini({ tka <- 0.3; tcl <- fix(1.0); tv <- fix(3.45); add.sd <- 0.2; eta.ka ~ 1.0 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv); cp <- linCmt(); cp ~ add(add.sd) })
  }
  fit <- suppressMessages(nlmixr2(mod, dat, est = "focei",
                                  control = foceiControl(print = 0, calcTables = FALSE)))
  fTka <- fit$parFixedDf["tka", "Estimate"]; fAdd <- fit$parFixedDf["add.sd", "Estimate"]

  rmod <- function() {
    ini({ tka <- 0.3; tcl <- fix(1.0); tv <- fix(3.45); add.sd <- 0.2; eta.ka ~ 0.6 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv); cp <- linCmt(); cp ~ add(add.sd) })
  }
  rf <- suppressMessages(nlmixr2(rmod, dat, est = "rpem",
                                 control = rpemControl(nGauss = 300L, nMH = 80000L, mhBurn = 8000L,
                                                       niter = 30L, collect = 12L, seed = 100L, cores = 4L)))
  # estimates match FOCEI's M3
  expect_equal(unname(rf$parFixedDf["tka", "Estimate"]), fTka, tolerance = 0.05)
  expect_equal(unname(rf$parFixedDf["add.sd", "Estimate"]), fAdd, tolerance = 0.03)
  # the fit object records the censoring method that was performed
  expect_match(as.character(rf$censInformation), "M3")
})

test_that("RPEM M3 censoring runs in the C++ cLoop and matches the R loop", {
  skip_on_cran()
  skip_on_ci()  # heavy: two RPEM fits

  struct <- rxode2::rxode2({ ka <- exp(tka + eta); cl <- exp(tcl); v <- exp(tv); cp <- linCmt() })
  set.seed(7); nsub <- 40L; obsT <- seq(0.5, 24, by = 1.5); addT <- 0.15; loq <- 0.5
  etasTrue <- rnorm(nsub, 0, sqrt(0.3))
  dat <- do.call(rbind, lapply(seq_len(nsub), function(i) {
    ev <- rxode2::et(amt = 100, cmt = "depot"); ev <- rxode2::et(ev, obsT)
    s <- rxode2::rxSolve(struct, params = c(tka = 0.45, tcl = 1, tv = 3.45, eta = etasTrue[i]),
                         events = ev, returnType = "data.frame", addDosing = FALSE)
    d <- as.data.frame(ev); d$id <- i; o <- d$evid == 0; d$DV <- 0
    d$DV[o] <- s$cp + rnorm(nrow(s), 0, addT)
    d$cens <- 0L; bl <- o & d$DV < loq; d$cens[bl] <- 1L; d$DV[bl] <- loq; d
  }))
  rmod <- function() {
    ini({ tka <- 0.3; tcl <- fix(1.0); tv <- fix(3.45); add.sd <- 0.2; eta.ka ~ 0.6 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv); cp <- linCmt(); cp ~ add(add.sd) })
  }
  ui <- rxode2::rxUiDecompress(rxode2::rxode2(rmod))
  ctl <- function(cl) rpemControl(nGauss = 300L, nMH = 80000L, mhBurn = 8000L, niter = 30L,
                                  collect = 12L, seed = 100L, cores = 4L, cLoop = cl)
  rfR <- .rpemFit(ui, dat, ctl(FALSE))
  rfC <- .rpemFit(ui, dat, ctl(TRUE))
  # the C++ loop's censored residual M-step matches the R loop and de-biases add.sd
  expect_equal(unname(rfC$mu["tka"]), unname(rfR$mu["tka"]), tolerance = 0.03)
  expect_equal(rfC$addSd, rfR$addSd, tolerance = 0.02)
  expect_equal(rfC$addSd, 0.15, tolerance = 0.03)   # censoring de-biases the residual sd
})
