# est="rpem" dispatch: nlmixr2(model, data, est="rpem") must resolve the control,
# run the E-M loop, and return a fit object with sensible estimates (K=1 core).

test_that("est='rpem' dispatch runs end-to-end (K=1)", {
  skip_on_cran()

  struct <- rxode2::rxode2({ ka <- exp(tka+eta); cl <- exp(tcl); v <- exp(tv); cp <- linCmt() })
  set.seed(3)
  nsub <- 25L; addSd <- 0.15; obsT <- seq(0.25, 24, by=1.5)
  etasTrue <- rnorm(nsub, 0, sqrt(0.3))
  dat <- do.call(rbind, lapply(seq_len(nsub), function(i) {
    ev <- rxode2::et(amt=100, cmt="depot"); ev <- rxode2::et(ev, obsT)
    s <- rxode2::rxSolve(struct, params=c(tka=0.45, tcl=1, tv=3.45, eta=etasTrue[i]),
                         events=ev, returnType="data.frame", addDosing=FALSE)
    d <- as.data.frame(ev); d$id <- i; d$DV <- 0; o <- d$evid == 0
    d$DV[o] <- s$cp + rnorm(nrow(s), 0, addSd); d
  }))

  rmod <- function() {
    ini({ tka <- 0.3; tcl <- fix(1.0); tv <- fix(3.45); add.sd <- 0.3; eta.ka ~ 0.6 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv); cp <- linCmt(); cp ~ add(add.sd) })
  }

  fit <- nlmixr2(rmod, dat, est="rpem",
                 control=rpemControl(nGauss=200L, nMH=50000L, mhBurn=5000L,
                                     niter=20L, collect=8L, seed=100L))

  # est="rpem" returns a full nlmixr2FitData (population fit assembled via the
  # eval-only FOCEI path at the RPEM estimates).
  expect_s3_class(fit, "nlmixr2FitData")
  expect_true(is.data.frame(fit))
  expect_true(all(c("IPRED", "PRED", "CWRES", "eta.ka") %in% names(fit)))
  expect_equal(nrow(fit), sum(dat$evid == 0))
  # parameter table has the RPEM estimates (finite, plausible ranges)
  .pf <- fit$parFixedDf
  expect_true(is.finite(.pf["tka", "Estimate"]))
  expect_gt(.pf["tka", "Estimate"], 0.1);  expect_lt(.pf["tka", "Estimate"], 1.2)
  expect_gt(fit$omega[1, 1], 0.05); expect_lt(fit$omega[1, 1], 0.8)
  expect_true(is.finite(fit$objDf[1, "OBJF"]))
  # standard errors (FOCEI covariance at the RPEM estimates)
  expect_true(is.finite(.pf["tka", "SE"]) && .pf["tka", "SE"] > 0)
  expect_true(is.finite(.pf["add.sd", "SE"]) && .pf["add.sd", "SE"] > 0)
})

test_that("est='rpem' runs with a bounded mu-referenced typical value (centered-eta demotion)", {
  skip_on_cran()

  # A bounded typical value with an eta (tcl <- log(c(0, 2.7, 100))) is rewritten by the
  # bounded transform into a structural theta (rxBoundedTr.tcl) + model-computed tcl,
  # demoting eta.cl to a centered (non-mu-referenced) eta.  The C++ loop takes a per-eta
  # muIdx with -1 for centered etas; this used to fail with
  # 'muIdx / mu0 / omDiag0 must have nEta entries'.
  one.cmt <- function() {
    ini({
      tka <- 0.45
      tcl <- log(c(0, 2.7, 100))
      tv <- 3.45
      eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }
  fit <- suppressMessages(suppressWarnings(nlmixr2(one.cmt, nlmixr2data::theo_sd, est = "rpem",
    control = rpemControl(nGauss = 100L, nMH = 20000L, mhBurn = 2000L,
                          niter = 10L, collect = 4L, seed = 42L, print = 0L))))
  expect_s3_class(fit, "nlmixr2FitData")
  .pf <- fit$parFixedDf
  # tcl is reported back-transformed on the natural scale, inside its declared box
  expect_gt(.pf["tcl", "Estimate"], log(1e-8)); expect_lt(.pf["tcl", "Estimate"], log(100))
  # all three omegas (including the centered eta.cl) are estimated, finite and positive
  .om <- diag(fit$omega)
  expect_equal(length(.om), 3L)
  expect_true(all(is.finite(.om) & .om > 0))
})
