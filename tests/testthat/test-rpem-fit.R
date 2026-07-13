# .rpemFit drives the whole K=1 E-M loop from a mu-referenced UI (parameter
# classification + engine).  It must recover the FOCEI estimates on the same
# well-identified data (design/rpem/11 Bar 2).

test_that(".rpemFit recovers FOCEI mu/Omega/add.sd via the packaged loop (K=1)", {
  skip_on_cran()

  struct <- rxode2::rxode2({ ka <- exp(tka+eta); cl <- exp(tcl); v <- exp(tv); cp <- linCmt() })
  set.seed(7)
  nsub <- 40L; trueTka <- 0.45; addSd <- 0.15; obsT <- seq(0.25, 24, by=1.5)
  etasTrue <- rnorm(nsub, 0, sqrt(0.3))
  dat <- do.call(rbind, lapply(seq_len(nsub), function(i) {
    ev <- rxode2::et(amt=100, cmt="depot"); ev <- rxode2::et(ev, obsT)
    s <- rxode2::rxSolve(struct, params=c(tka=trueTka, tcl=1, tv=3.45, eta=etasTrue[i]),
                         events=ev, returnType="data.frame", addDosing=FALSE)
    d <- as.data.frame(ev); d$id <- i; d$DV <- 0; o <- d$evid == 0
    d$DV[o] <- s$cp + rnorm(nrow(s), 0, addSd); d
  }))

  # FOCEI reference (hold tcl, tv fixed at truth -- RPEM holds them fixed too).
  mod <- function() {
    ini({ tka <- 0.1; tcl <- fix(1.0); tv <- fix(3.45); add.sd <- 0.3; eta.ka ~ 1.0 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv); cp <- linCmt(); cp ~ add(add.sd) })
  }
  fit <- suppressMessages(nlmixr2(mod, dat, est="focei", control=foceiControl(print=0)))
  fTka <- fit$parFixedDf["tka", "Estimate"]; fOm <- fit$omega[1, 1]
  fSd  <- fit$parFixedDf["add.sd", "Estimate"]

  # RPEM through the packaged loop (no fix(); tcl/tv held fixed at their ini).
  rmod <- function() {
    ini({ tka <- 0.1; tcl <- fix(1.0); tv <- fix(3.45); add.sd <- 0.3; eta.ka ~ 0.6 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv); cp <- linCmt(); cp ~ add(add.sd) })
  }
  ui <- rxode2::rxode2(rmod)
  rfit <- .rpemFit(ui, dat, rpemControl(nGauss=300L, nMH=80000L, mhBurn=8000L,
                                        niter=30L, collect=12L, seed=2000L))

  expect_equal(unname(rfit$mu["tka"]), fTka, tolerance = 0.04)
  expect_equal(unname(rfit$omega["eta.ka"]), fOm, tolerance = 0.10)
  expect_equal(rfit$addSd, fSd, tolerance = 0.10)

  # Per-subject EBEs (posterior-mean etas) should track the true simulated etas
  # (correlation is shift-invariant, so it holds despite the mu offset).
  expect_equal(nrow(rfit$ebe), nsub)
  expect_gt(stats::cor(rfit$ebe[, 1], etasTrue), 0.8)
})
