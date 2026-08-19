nmTest({

  ## #972 regression: saem's S-step (do_mcmc/arDYFhyp family) ignored the
  ## pow() exponent (cres), evaluating every MCMC-acceptance proposal as if
  ## the residual model were a plain prop() (exponent fixed at 1).  With the
  ## exponent estimated, this created a destructive feedback loop: pw
  ## collapsed toward 0 while prop.sd inflated to compensate.  Reproducer and
  ## truth values are from the issue.

  test_that("SAEM recovers an estimated pow() exponent (rmPow, #972)", {
    .testSeed(42)
    dose <- 320; V <- 70; times <- c(0.5, 1, 2, 4, 7, 12, 24)
    n <- 60; eta <- rnorm(n, 0, sqrt(0.09))
    pwTrue <- 0.5; sdTrue <- 0.3
    d <- do.call(rbind, lapply(seq_len(n), function(i) {
      f <- dose / V * exp(-exp(log(4) + eta[i]) / V * times)
      g <- sdTrue * f^pwTrue
      rbind(data.frame(ID = i, TIME = 0, DV = NA_real_, AMT = dose, EVID = 1),
            data.frame(ID = i, TIME = times, DV = f + rnorm(length(times), 0, g),
                       AMT = 0, EVID = 0))
    }))

    mod <- function() {
      ini({ tcl <- log(4); eta.cl ~ 0.09; prop.sd <- 0.3; pw <- 0.5 })
      model({ cl <- exp(tcl + eta.cl); v <- 70; linCmt() ~ pow(prop.sd, pw) })
    }
    fit <- .nlmixr(mod, d, est = "saem")
    .pw <- fit$parFixedDf["pw", "Estimate"]
    .propSd <- fit$parFixedDf["prop.sd", "Estimate"]
    ## Before the fix pw collapsed to ~0.01-0.04 (truth 0.5) and prop.sd
    ## inflated to ~0.6-3.8 (truth 0.3) to compensate; these loose bounds
    ## pass with the fix and fail without it.
    expect_true(.pw > 0.3 && .pw < 0.7)
    expect_true(.propSd > 0.15 && .propSd < 0.45)
  })

})
