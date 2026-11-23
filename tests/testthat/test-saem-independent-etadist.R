test_that("SAEM estimates the parameters of an independent gamma ETA", {
  skip_on_cran()
  skip_if_not(exists("rxUiEtaDists", asNamespace("rxode2"), inherits = FALSE))
  skip_if_not(exists("rxEtaDistExpand", asNamespace("rxode2"), inherits = FALSE))
  mod <- function() {
    ini({
      lmean <- log(4)
      lrv <- log(0.3)
      eta.cl ~ 1
      add.sd <- 0.1
      dist(eta.cl) ~ dgamma(shape = exp(-lrv), rate = exp(-lrv - lmean))
    })
    model({
      cl <- eta.cl
      v <- 8
      linCmt() ~ add(add.sd)
    })
  }

  threads <- rxode2::getRxThreads()
  on.exit(rxode2::setRxThreads(threads = threads))
  rxode2::setRxThreads(threads = 1L)
  data <- withr::with_seed(59, {
    cl <- rgamma(60, shape = 4, rate = 4/4.8)
    time <- c(0.15, 0.4, 0.8, 1.5, 3, 5)
    do.call(rbind, lapply(seq_along(cl), function(i) {
      pred <- 50/8 * exp(-cl[i]/8 * time)
      data.frame(ID = i, TIME = c(0, time), EVID = c(1, rep(0, length(time))),
                 CMT = 1, AMT = c(50, rep(0, length(time))),
                 DV = c(NA_real_, pred + rnorm(length(time), sd = 0.08)))
    }))
  })

  fit <- suppressWarnings(nlmixr2(mod, data, "saem",
    saemControl(nBurn = 45, nEm = 30, nmc = 3L, seed = 99L,
                print = 0L, covMethod = "", calcTables = FALSE)))
  .est <- fixef(fit)[c("lmean", "lrv")]
  expect_true(all(is.finite(.est)))
  expect_lt(abs(unname(.est["lmean"]) - log(4.8)), 0.3)
  expect_gt(abs(unname(.est["lmean"]) - log(4)), 0.02)
})
