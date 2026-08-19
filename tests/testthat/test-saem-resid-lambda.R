nmTest({

  ## B2 regression: SAEM lambda (Box-Cox) combined error models under-sized xmin, overflowing the Nelder-Mead buffer.

  .ctl <- saemControl(nBurn = 20, nEm = 20, print = 0L, nmc = 2)

  test_that("SAEM pow + boxCox (rmPowLam) runs with finite estimates", {
    f <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45
            eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1; prop.sd <- 0.3; pw <- 0.8; lm <- 0.8 })
      model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
              linCmt() ~ pow(prop.sd, pw) + boxCox(lm) })
    }
    fit <- .nlmixr(f, theo_sd, est = "saem", control = .ctl)
    expect_true(all(is.finite(fit$parFixedDf$Estimate)))
  })

  test_that("SAEM add + prop + boxCox (rmAddPropLam) runs with finite estimates", {
    f <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45
            eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1; add.sd <- 0.3; prop.sd <- 0.2; lm <- 0.8 })
      model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
              linCmt() ~ add(add.sd) + prop(prop.sd) + boxCox(lm) })
    }
    fit <- .nlmixr(f, theo_sd, est = "saem", control = .ctl)
    expect_true(all(is.finite(fit$parFixedDf$Estimate)))
  })

  test_that("SAEM add + pow + boxCox (rmAddPowLam, n=4) runs with finite estimates", {
    f <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45
            eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1; add.sd <- 0.3; prop.sd <- 0.2; pw <- 0.8; lm <- 0.8 })
      model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
              linCmt() ~ add(add.sd) + pow(prop.sd, pw) + boxCox(lm) })
    }
    fit <- .nlmixr(f, theo_sd, est = "saem", control = .ctl)
    expect_true(all(is.finite(fit$parFixedDf$Estimate)))
  })

  # #914: a fixed boxCox lambda never reached the transform (saem.cpp's `lambda`
  # member was never synced from `lres`, so every kernel _powerD() call ran at
  # lambda=1) and the pure add+lambda residual model (rmAddLam) never zeroed
  # bres, leaving g = ares + 1*|ft| instead of g = ares.  Both defects are silent
  # -- the fit still converges, just to badly biased estimates -- so recovery
  # against simulated truth is asserted here rather than just finiteness.
  test_that("SAEM add + boxCox (rmAddLam) with a fixed lambda recovers truth", {
    set.seed(914)
    .dose <- 320; .v <- 70; .times <- c(0.5, 1, 2, 4, 7, 12, 24)
    .bcF <- function(y, l) (y^l - 1) / l
    .bcI <- function(x, l) (l * x + 1)^(1 / l)
    .n <- 20; .lam <- 0.5
    .eta <- rnorm(.n, 0, sqrt(0.09))
    .d <- do.call(rbind, lapply(seq_len(.n), function(i) {
      .f <- .dose / .v * exp(-exp(log(4) + .eta[i]) / .v * .times)
      rbind(data.frame(ID = i, TIME = 0, DV = NA_real_, AMT = .dose, EVID = 1),
            data.frame(ID = i, TIME = .times,
                       DV = .bcI(.bcF(.f, .lam) + rnorm(length(.times), 0, 0.15), .lam),
                       AMT = 0, EVID = 0))
    }))
    f <- function() {
      ini({ tcl <- log(4); eta.cl ~ 0.09; add.sd <- 0.15; lambda <- fixed(0.5) })
      model({ cl <- exp(tcl + eta.cl); v <- 70; linCmt() ~ add(add.sd) + boxCox(lambda) })
    }
    fit <- .nlmixr(f, .d, est = "saem",
                   control = saemControl(nBurn = 100, nEm = 100, print = 0L, nmc = 2))
    # the fixed lambda must reach transMat, not stay at the identity default
    expect_equal(unname(fit$saem$transMat[1, 1]), 0.5)
    expect_lt(abs(fit$theta[["tcl"]] - log(4)), 0.3)
    expect_lt(fit$parFixedDf["add.sd", "Estimate"], 0.3)
  })

})
