## Regression test for parallel FOCEI determinism (issue 641 follow-up).
##
## With sortIds(rx, 0) in inner.cpp's post-parallel reorder, the next outer
## iteration would re-dispatch subjects in wall-time-decreasing order, which
## made the optimizer trajectory depend on system load and produce
## different fitted parameters across runs even with the same seed.
##
## Switching to sortIds(rx, 2) (deterministic iota order) plus rxode2's
## per-thread tolerance arrays closes that source of non-determinism.
## This test pins cores=2 and runs the same fit twice — they must agree
## bit-for-bit.

test_that("issue 641: parallel FOCEI is deterministic across runs", {
  skip_on_cran()
  skip_if(rxode2::getRxThreads() < 2L,
          "test requires at least 2 threads available")

  emax_fun <- function(x, e0, emax, ex50, hill = 1) {
    e0 + emax * x^hill / (ex50^hill + x^hill)
  }

  set.seed(641)
  d_mod <- data.frame(
    WEEK = rep(0:12, 5),
    ID = rep(1:5, each = 13),
    COVV = 1
  )
  d_mod$emax_val <- rnorm(nrow(d_mod), mean = 20, sd = 5) + d_mod$ID
  d_mod$DV <- emax_fun(d_mod$WEEK, e0 = 0, emax = d_mod$emax_val, ex50 = 4)
  d_mod$TIME <- d_mod$WEEK

  mod <- function() {
    ini({
      e0 <- fixed(0)
      tvemax <- -40
      covValueEmax <- fixed(0)
      let50 <- log(3)
      covValueET50 <- fixed(0)
      iivemax ~ 10
      iivet50 ~ 0
      addSd <- 5
    })
    model({
      emax <- tvemax + iivemax + covValueEmax * COVV
      et50 <- exp(let50 + iivet50 + covValueET50 * COVV)
      eff <- e0 + emax * WEEK / (et50 + WEEK)
      eff ~ add(addSd)
    })
  }

  prevThreads <- rxode2::getRxThreads()
  on.exit(rxode2::setRxThreads(prevThreads), add = TRUE)
  rxode2::setRxThreads(2L)

  doFit <- function() {
    nlmixr(mod, data = d_mod, est = "focei",
           control = foceiControl(print = 0, covMethod = "",
                                  calcTables = FALSE))
  }

  fitA <- doFit()
  fitB <- doFit()

  expect_identical(fitA$objf, fitB$objf)
  expect_identical(as.numeric(fitA$parFixedDf$Estimate),
                   as.numeric(fitB$parFixedDf$Estimate))
})
