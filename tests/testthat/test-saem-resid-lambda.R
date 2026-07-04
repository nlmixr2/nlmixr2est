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

})
