nmTest({

  ## #915: a logitNorm()/probitNorm() endpoint with an unfixed residual sd
  ## errored during fit finalization with "missing value where TRUE/FALSE
  ## needed" -- .getSaemTheta() never copied the saem-estimated residual sd
  ## for these two error types into ui$iniDf$est, so the FOCEI scaleC setup
  ## later hit an NA when reentering to build the fit table.

  dose <- 320; V <- 70; times <- c(0.5, 1, 2, 4, 7, 12, 24)
  LO <- 0; HI <- 10
  lgt <- function(y) log((y - LO) / (HI - y))
  set.seed(11); n <- 12; e <- rnorm(n, 0, sqrt(0.09))
  d <- do.call(rbind, lapply(seq_len(n), function(i) {
    f <- dose / V * exp(-exp(log(4) + e[i]) / V * times)
    rbind(data.frame(ID = i, TIME = 0, DV = NA_real_, AMT = dose, EVID = 1),
          data.frame(ID = i, TIME = times,
                     DV = LO + (HI - LO) * plogis(lgt(f) + rnorm(length(times), 0, 0.25)),
                     AMT = 0, EVID = 0))
  }))

  .ctl <- saemControl(nBurn = 20, nEm = 20, print = 0L, nmc = 2, calcTables = FALSE)

  test_that("SAEM logitNorm() with explicit bounds runs and estimates the residual sd", {
    m <- function() {
      ini({ tcl <- log(4); eta.cl ~ 0.09; logit.sd <- 0.5 })
      model({ cl <- exp(tcl + eta.cl); v <- 70; linCmt() ~ logitNorm(logit.sd, 0, 10) })
    }
    fit <- suppressWarnings(.nlmixr(m, d, est = "saem", control = .ctl))
    expect_true(all(is.finite(fit$parFixedDf$Estimate)))
    expect_true(!is.na(fit$theta["logit.sd"]))
    expect_true(fit$theta["logit.sd"] > 0.1 && fit$theta["logit.sd"] < 0.5)
  })

  test_that("SAEM probitNorm() with explicit bounds runs and estimates the residual sd", {
    m <- function() {
      ini({ tcl <- log(4); eta.cl ~ 0.09; probit.sd <- 0.5 })
      model({ cl <- exp(tcl + eta.cl); v <- 70; linCmt() ~ probitNorm(probit.sd, 0, 10) })
    }
    fit <- suppressWarnings(.nlmixr(m, d, est = "saem", control = .ctl))
    expect_true(all(is.finite(fit$parFixedDf$Estimate)))
    expect_true(!is.na(fit$theta["probit.sd"]))
  })

})
