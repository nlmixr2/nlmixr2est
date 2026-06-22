nmTest({

  ## D2: handleCensNpdeCdf computed j = trunc(pd*K) before pd was clamped, so a
  ## (right-)censored point whose pd rounds to 1.0 read curRow[K]/curRow[K+1]
  ## out of bounds (and with unsigned j, j-1 wrapped for K==1).  This exercises
  ## the censored cdf npde path; the OOB read is caught by valgrind/ASAN CI.
  test_that("npde cdf censoring on BLQ data yields finite NPDE (D2 OOB regression)", {
    d <- theo_sd
    lloq <- 2.0
    d$CENS <- ifelse(d$EVID == 0 & d$DV < lloq, 1L, 0L)
    d$DV[d$CENS == 1] <- lloq
    skip_if_not(sum(d$CENS) > 0)
    one.cmt <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45
            eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1; add.sd <- 0.7 })
      model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
              linCmt() ~ add(add.sd) })
    }
    fit <- suppressMessages(suppressWarnings(
      .nlmixr(one.cmt, d, est = "focei",
              control = foceiControl(print = 0L),
              table = tableControl(censMethod = "cdf", npde = TRUE))))
    expect_true("NPDE" %in% names(fit))
    expect_true(all(is.finite(fit$NPDE)))
  })

  ## B1: do_mcmc passed the raw and transformed predictions to handleF in the
  ## wrong order (vs saem.cpp:942), so for propT/transformed-error endpoints the
  ## residual SD used in the MCMC acceptance was computed from the un-transformed
  ## prediction.  Exercise the transformed-prediction do_mcmc path.
  test_that("SAEM prop+boxCox (transformed-prediction do_mcmc path) runs finite (B1)", {
    f <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45
            eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1; prop.sd <- 0.3; lm <- 0.8 })
      model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
              linCmt() ~ prop(prop.sd) + boxCox(lm) })
    }
    fit <- suppressMessages(suppressWarnings(
      .nlmixr(f, theo_sd, est = "saem",
              control = saemControl(nBurn = 30, nEm = 30, print = 0L, nmc = 3))))
    expect_true(all(is.finite(fit$parFixedDf$Estimate)))
  })

})
