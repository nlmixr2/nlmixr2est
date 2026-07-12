# Opt-in mode-centered importance sampling (impInflate).  Default (impInflate=0)
# uses the paper's prior sampling; impInflate>=1 draws eta ~ N(EBE, c*Omega) and
# importance-weights, improving coverage of the largest Omega in multi-eta models
# (see design/rpem/04).  This is a partial mitigation, not a full fix.

test_that("RPEM impInflate opt-in raises the largest omega toward FOCEI/SAEM", {
  skip_on_cran()
  skip_on_ci()  # heavy: two multi-iteration RPEM loops

  struct <- rxode2::rxode2({ ka <- exp(tka + eka); cl <- exp(tcl + ecl); v <- exp(tv + ev); cp <- linCmt() })
  set.seed(11); nsub <- 80L; obsT <- seq(0.5, 24, by = 2)
  trueAdd <- 0.1; omKa <- 0.3; omCl <- 0.1; omV <- 0.06
  eKa <- rnorm(nsub, 0, sqrt(omKa)); eCl <- rnorm(nsub, 0, sqrt(omCl)); eV <- rnorm(nsub, 0, sqrt(omV))
  dat <- do.call(rbind, lapply(seq_len(nsub), function(i) {
    ev <- rxode2::et(amt = 100, cmt = "depot"); ev <- rxode2::et(ev, obsT)
    s <- rxode2::rxSolve(struct, params = c(tka = 0.5, tcl = 1, tv = 3.45,
                                            eka = eKa[i], ecl = eCl[i], ev = eV[i]),
                         events = ev, returnType = "data.frame", addDosing = FALSE)
    dose <- data.frame(id = i, time = 0, DV = 0, evid = 1L, amt = 100, cmt = 1L)
    obs <- data.frame(id = i, time = obsT, DV = s$cp + rnorm(length(s$cp), 0, trueAdd),
                      evid = 0L, amt = 0, cmt = 1L)
    rbind(dose, obs)
  }))
  dat <- dat[order(dat$id, dat$time), ]

  rmod <- function() {
    ini({ tka <- 0.5; tcl <- 1.0; tv <- 3.45; add.sd <- 0.2
          eka ~ 0.2; ecl ~ 0.08; ev ~ 0.05 })
    model({ ka <- exp(tka + eka); cl <- exp(tcl + ecl); v <- exp(tv + ev); cp <- linCmt()
            cp ~ add(add.sd) })
  }
  ui <- rxode2::rxode2(rmod)

  ctl <- function(ii) rpemControl(nGauss = 400L, nMH = 60000L, mhBurn = 6000L,
                                  niter = 30L, collect = 12L, seed = 77L, impInflate = ii)
  # reset the global threefry stream so the R-loop impInflate M-step MH (which draws from
  # it) is reproducible regardless of test order
  rxode2::rxSetSeed(77)
  rf0 <- .rpemFit(ui, dat, ctl(0))
  rf4 <- .rpemFit(ui, dat, ctl(4))

  omKa0 <- unname(rf0$omega["eka"])
  omKa4 <- unname(rf4$omega["eka"])

  # Partial mitigation (not a full fix): prior sampling under-estimates the
  # largest omega (~0.13 here vs truth 0.30); opt-in pulls it substantially up
  # but not all the way.
  expect_true(omKa4 > omKa0)
  expect_true(omKa4 > 0.18)
  # smaller omegas recovered by both (wide band: a stochastic MCPEM omega estimate whose
  # in-suite value drifts a little with the global-RNG state left by earlier tests);
  # residual stays sane (inflation bumps it a bit)
  expect_equal(unname(rf4$omega["ecl"]), omCl, tolerance = 0.25)
  expect_equal(unname(rf4$omega["ev"]), omV, tolerance = 0.25)
  expect_true(abs(rf4$addSd - trueAdd) < 0.07)
})
