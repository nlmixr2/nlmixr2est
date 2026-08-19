nmTest({

  ## #915: a logitNorm()/probitNorm() endpoint with an unfixed residual sd
  ## errored during fit finalization with "missing value where TRUE/FALSE
  ## needed" -- .getSaemTheta() never copied the saem-estimated residual sd
  ## for these two error types (and, the same gap, propF()/powF()/powF2())
  ## into ui$iniDf$est, so the FOCEI scaleC setup later hit an NA when
  ## reentering to build the fit table.  While investigating, an antigravity
  ## review also caught a second, independent gap: rxUiGet.saemAres/Bres/Cres
  ## (R/saemRxUiGet.R) likewise never matched these error types, so the SAEM
  ## kernel silently *started* from the hardcoded fallback (10/1) instead of
  ## the user's ini() estimate.

  dose <- 320; V <- 70; times <- c(0.5, 1, 2, 4, 7, 12, 24)
  LO <- 0; HI <- 10

  .ctl <- saemControl(nBurn = 20, nEm = 20, print = 0L, nmc = 2, calcTables = FALSE)

  test_that("saemAres/saemBres honor the user's ini() estimate for the new error types", {
    m1 <- function() {
      ini({ tcl <- log(4); eta.cl ~ 0.09; logit.sd <- 0.37 })
      model({ cl <- exp(tcl + eta.cl); v <- 70; linCmt() ~ logitNorm(logit.sd, 0, 10) })
    }
    expect_equal(rxode2::rxode2(m1)$saemAres, 0.37)

    m2 <- function() {
      ini({ tcl <- log(4); eta.cl ~ 0.09; probit.sd <- 0.41 })
      model({ cl <- exp(tcl + eta.cl); v <- 70; linCmt() ~ probitNorm(probit.sd, 0, 10) })
    }
    expect_equal(rxode2::rxode2(m2)$saemAres, 0.41)

    m3 <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45
            eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1; prop.sd <- 0.42 })
      model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
              f2 <- 1
              linCmt() ~ propF(prop.sd, f2) })
    }
    expect_equal(rxode2::rxode2(m3)$saemBres, 0.42)

    m4 <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45
            eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1; pow.sd <- 0.53; pw <- 1 })
      model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
              f2 <- 1
              linCmt() ~ powF(pow.sd, pw, f2) })
    }
    expect_equal(rxode2::rxode2(m4)$saemBres, 0.53)
  })

  test_that("SAEM logitNorm() with explicit bounds runs and estimates the residual sd", {
    lgt <- function(y) log((y - LO) / (HI - y))
    set.seed(11); n <- 12; e <- rnorm(n, 0, sqrt(0.09))
    d <- do.call(rbind, lapply(seq_len(n), function(i) {
      f <- dose / V * exp(-exp(log(4) + e[i]) / V * times)
      rbind(data.frame(ID = i, TIME = 0, DV = NA_real_, AMT = dose, EVID = 1),
            data.frame(ID = i, TIME = times,
                       DV = LO + (HI - LO) * plogis(lgt(f) + rnorm(length(times), 0, 0.25)),
                       AMT = 0, EVID = 0))
    }))
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
    pbt <- function(y) qnorm((y - LO) / (HI - LO))
    set.seed(13); n <- 12; e <- rnorm(n, 0, sqrt(0.09))
    d <- do.call(rbind, lapply(seq_len(n), function(i) {
      f <- dose / V * exp(-exp(log(4) + e[i]) / V * times)
      rbind(data.frame(ID = i, TIME = 0, DV = NA_real_, AMT = dose, EVID = 1),
            data.frame(ID = i, TIME = times,
                       DV = LO + (HI - LO) * pnorm(pbt(f) + rnorm(length(times), 0, 0.25)),
                       AMT = 0, EVID = 0))
    }))
    m <- function() {
      ini({ tcl <- log(4); eta.cl ~ 0.09; probit.sd <- 0.5 })
      model({ cl <- exp(tcl + eta.cl); v <- 70; linCmt() ~ probitNorm(probit.sd, 0, 10) })
    }
    fit <- suppressWarnings(.nlmixr(m, d, est = "saem", control = .ctl))
    expect_true(all(is.finite(fit$parFixedDf$Estimate)))
    expect_true(!is.na(fit$theta["probit.sd"]))
    expect_true(fit$theta["probit.sd"] > 0.1 && fit$theta["probit.sd"] < 0.5)
  })

  test_that("SAEM propF() runs and estimates the residual sd", {
    f <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45
            eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1; prop.sd <- 0.3 })
      model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
              f2 <- 1
              linCmt() ~ propF(prop.sd, f2) })
    }
    fit <- suppressWarnings(.nlmixr(f, nlmixr2data::theo_sd, est = "saem", control = .ctl))
    expect_true(all(is.finite(fit$parFixedDf$Estimate)))
    expect_true(!is.na(fit$theta["prop.sd"]))
    expect_true(fit$theta["prop.sd"] > 0.01 && fit$theta["prop.sd"] < 1)
  })

  test_that("SAEM powF() runs and estimates the residual sd", {
    ## note: a plain pow(pow.sd, pw) on theo_sd is itself not tightly
    ## identified at this nBurn/nEm (a pre-existing, unrelated instability,
    ## reproducible with no "F" variant involved) -- match the "finite"
    ## convention test-saem-resid-lambda.R uses for the pow family rather
    ## than pinning a numeric range; the exact defect this PR fixes (the
    ## wrong *starting* value) is pinned deterministically above via saemBres.
    f <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45
            eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1; pow.sd <- 0.3; pw <- 0.8 })
      model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
              f2 <- 1
              linCmt() ~ powF(pow.sd, pw, f2) })
    }
    fit <- suppressWarnings(.nlmixr(f, nlmixr2data::theo_sd, est = "saem", control = .ctl))
    expect_true(all(is.finite(fit$parFixedDf$Estimate)))
    expect_true(!is.na(fit$theta["pow.sd"]))
  })

})
