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

test_that("SAEM expands a declared correlated block with a Gaussian copula", {
  skip_if_not(exists("rxUiEtaDists", asNamespace("rxode2"), inherits = FALSE))
  skip_if_not(exists("rxEtaDistExpand", asNamespace("rxode2"), inherits = FALSE))
  mod <- function() {
    ini({
      lcl <- log(4)
      lv <- log(8)
      eta.cl + eta.v ~ c(1, 0.35, 1)
      dist(eta.cl) ~ dgamma(shape = 4, rate = 4/exp(lcl))
      dist(eta.v) ~ dgamma(shape = 3, rate = 3/exp(lv))
      add.sd <- 0.1
    })
    model({
      cl <- eta.cl
      v <- eta.v
      linCmt() ~ add(add.sd)
    })
  }

  expanded <- .preProcessEtaDist(rxode2::rxode2(mod), "saem", NULL, saemControl())
  corRow <- expanded$ui$iniDf[grepl("^rxCor[.]", expanded$ui$iniDf$name), ]
  expect_equal(nrow(corRow), 1L)
  expect_equal(tanh(corRow$est), 0.35)
  expect_equal(unname(as.numeric(expanded$ui$omega)), as.numeric(diag(2L)))
  expect_identical(expanded$control$nonMuThetaOpt, "bhhh")
  for (method in c("focei", "imp", "impmap", "vae")) {
    alt <- .preProcessEtaDist(rxode2::rxode2(mod), method, NULL, NULL)
    expect_null(alt$control)
    expect_equal(nrow(rxode2::rxUiEtaDists(alt$ui)), 0L)
  }
})
