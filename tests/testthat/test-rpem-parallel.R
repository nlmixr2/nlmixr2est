# Parallel E-step: cores > 1 solves each Monte Carlo sample's subjects via
# rxode2's par_solve (per-thread ODE workspace).  The eta draw is done in R, so it
# is independent of the core count; a multi-core fit must agree with a single-core
# fit to Monte-Carlo precision (design/rpem/06, D21).

test_that("RPEM parallel E-step (cores>1) agrees with serial", {
  skip_on_cran()
  skip_on_ci()  # heavy: two multi-iteration RPEM fits

  struct <- rxode2::rxode2({ ka <- exp(tka + eta); cl <- exp(tcl); v <- exp(tv); cp <- linCmt() })
  set.seed(3); nsub <- 30L; obsT <- seq(0.5, 24, by = 1.5)
  etasTrue <- rnorm(nsub, 0, sqrt(0.3))
  dat <- do.call(rbind, lapply(seq_len(nsub), function(i) {
    ev <- rxode2::et(amt = 100, cmt = "depot"); ev <- rxode2::et(ev, obsT)
    s <- rxode2::rxSolve(struct, params = c(tka = 0.45, tcl = 1, tv = 3.45, eta = etasTrue[i]),
                         events = ev, returnType = "data.frame", addDosing = FALSE)
    d <- as.data.frame(ev); d$id <- i; d$DV <- 0; o <- d$evid == 0
    d$DV[o] <- s$cp + rnorm(nrow(s), 0, 0.1)
    d
  }))
  rmod <- function() {
    ini({ tka <- 0.4; tcl <- fix(1.0); tv <- fix(3.45); add.sd <- 0.15; eta.ka ~ 0.4 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv); cp <- linCmt(); cp ~ add(add.sd) })
  }
  ui <- rxode2::rxode2(rmod)
  ctl <- function(cores) rpemControl(nGauss = 300L, nMH = 60000L, mhBurn = 6000L,
                                     niter = 20L, collect = 8L, seed = 7L, cores = cores)
  rf1 <- .rpemFit(ui, dat, ctl(1L))
  rf2 <- .rpemFit(ui, dat, ctl(2L))

  # same core count is exactly reproducible; across cores, agreement is to
  # Monte-Carlo precision: the MH is a Markov chain, so a single accept/reject flip
  # from the ~1e-8 ind_solve-vs-par_solve difference makes the chain diverge to a
  # statistically-equivalent trajectory (well within FOCEI-comparison tolerances).
  expect_equal(unname(rf2$mu["tka"]), unname(rf1$mu["tka"]), tolerance = 0.05)
  expect_equal(rf2$addSd, rf1$addSd, tolerance = 0.05)
})
