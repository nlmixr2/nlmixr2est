# Full C++ E-M loop (rpemControl(cLoop = TRUE), design/rpem/12 M5): the additive /
# proportional, diagonal-omega, mu-referenced core runs the whole E-M loop in C++,
# drawing etas from rxode2's per-thread threefry engine with a deterministic
# per-(iteration, subject) seed.  It must recover the same estimates as the R loop,
# be reproducible run-to-run, and produce a full nlmixr2FitData.

test_that("est=rpem cLoop recovers, matches the R loop, and is reproducible", {
  skip_on_cran()
  skip_on_ci()  # heavy: two RPEM fits

  struct <- rxode2::rxode2({ ka <- exp(tka + eta); cl <- exp(tcl); v <- exp(tv); cp <- linCmt() })
  set.seed(8); nsub <- 40L; etas <- rnorm(nsub, 0, sqrt(0.3)); obsT <- seq(0.5, 24, by = 1.5)
  dat <- do.call(rbind, lapply(seq_len(nsub), function(i) {
    ev <- rxode2::et(amt = 100, cmt = "depot"); ev <- rxode2::et(ev, obsT)
    s <- rxode2::rxSolve(struct, params = c(tka = 0.45, tcl = 1, tv = 3.45, eta = etas[i]),
                         events = ev, returnType = "data.frame", addDosing = FALSE)
    d <- as.data.frame(ev); d$id <- i; o <- d$evid == 0; d$DV <- 0
    d$DV[o] <- s$cp + rnorm(nrow(s), 0, 0.1)
    d
  }))
  mod <- function() {
    ini({ tka <- 0.3; tcl <- fix(1.0); tv <- fix(3.45); add.sd <- 0.2; eta.ka ~ 0.6 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv); cp <- linCmt(); cp ~ add(add.sd) })
  }
  ui <- rxode2::rxode2(mod)
  ctlR <- rpemControl(nGauss = 300L, nMH = 60000L, mhBurn = 6000L, niter = 25L,
                      collect = 10L, seed = 42L, cores = 4L, cLoop = FALSE)
  ctlC <- rpemControl(nGauss = 300L, nMH = 60000L, mhBurn = 6000L, niter = 25L,
                      collect = 10L, seed = 42L, cores = 4L, cLoop = TRUE)

  fitC <- suppressMessages(nlmixr2(mod, dat, est = "rpem", control = ctlC))
  expect_s3_class(fitC, "nlmixr2FitData")
  # recovers the truth (tka 0.45, add.sd 0.1)
  expect_equal(unname(fitC$parFixedDf["tka", "Estimate"]), 0.45, tolerance = 0.15)
  expect_equal(unname(fitC$parFixedDf["add.sd", "Estimate"]), 0.1, tolerance = 0.1)

  # matches the R loop closely (both target the same posterior; small Monte-Carlo gap)
  rfR <- .rpemFit(ui, dat, ctlR)
  rfC <- .rpemFit(ui, dat, ctlC)
  expect_equal(unname(rfC$mu["tka"]), unname(rfR$mu["tka"]), tolerance = 0.06)
  expect_equal(rfC$addSd, rfR$addSd, tolerance = 0.06)

  # reproducible run-to-run for a fixed core count (deterministic threefry seeds)
  rfC2 <- .rpemFit(ui, dat, ctlC)
  expect_equal(rfC$mu, rfC2$mu)
  expect_equal(rfC$omega, rfC2$omega)
  expect_equal(rfC$addSd, rfC2$addSd)
})
