## External likelihood-contribution API across the methods that do NOT read the
## objective straight off fInd->llik:
##
##   - npag / npb build the conditional density by summing llikObs, so likInner0
##     folds the contributed LL into llikObs[kk] as well.  Without that the
##     contribution is silently dropped for the nonparametric methods.
##   - the analytic VAE decoder never enters likInner0 at all; it cycles the same
##     registry inside vaeDecoderPxzCore.
##
## Same probe as test-lik-contrib.R: the test contributor adds a constant c to
## every observation's log-likelihood, which cannot move any optimum, so the
## objective must shift by exactly -2*c*nObs (-2LL) or -c*nObs (the ELBO loss,
## which is a plain negative log-likelihood plus KL).

.likContribModel <- function() {
  ini({
    tka <- fix(0.45); tcl <- fix(1); tv <- fix(3.45)
    eta.ka ~ fix(0.6); eta.cl ~ fix(0.3); eta.v ~ fix(0.1)
    add.sd <- fix(0.7)
  })
  model({
    ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
    d/dt(depot) <- -ka * depot
    d/dt(center) <- ka * depot - cl / v * center
    cp <- center / v
    cp ~ add(add.sd)
  })
}

test_that("likelihood contributions reach the nonparametric objective (npag/npb)", {
  skip_on_cran()

  .old <- rxode2::getRxThreads(); on.exit(rxode2::setRxThreads(.old), add = TRUE)
  rxode2::setRxThreads(1L)   # test contributor uses global accumulators

  .nObs <- sum(theo_sd$EVID == 0)
  cc <- 0.01
  ## a fit that errors must not leave the hook registered for the rest of the
  ## suite -- removal is idempotent, so an unconditional on.exit is safe
  on.exit(.Call("_nlmixr2est_removeTestContrib", PACKAGE = "nlmixr2est"), add = TRUE)

  for (.est in c("npag", "npb")) {
    .ctl <- do.call(paste0(.est, "Control"), list(print = 0L))
    f0 <- .nlmixr(.likContribModel, theo_sd, est = .est, control = .ctl)

    .Call("_nlmixr2est_registerTestContrib", PACKAGE = "nlmixr2est")
    .Call("_nlmixr2est_setTestContribAddLL", cc, PACKAGE = "nlmixr2est")
    f1 <- .nlmixr(.likContribModel, theo_sd, est = .est, control = .ctl)
    res <- .Call("_nlmixr2est_getTestContrib", PACKAGE = "nlmixr2est")
    .Call("_nlmixr2est_removeTestContrib", PACKAGE = "nlmixr2est")

    ## the llikObs fold-in, not just fInd->llik: a constant c per observation
    ## shifts -2LL by exactly -2*c*nObs
    expect_equal(f1$objf - f0$objf, -2 * cc * .nObs, tolerance = 1e-6,
                 info = .est)
    ## and the hook actually ran, in balanced subject brackets
    expect_gt(res[[1]], 0)                 # nObs
    expect_equal(res[[5]], res[[6]])       # nBegin == nEnd
  }
})

nmTest({
  test_that("likelihood contributions reach the analytic VAE decoder ELBO", {
    skip_on_cran()

    theo <- function() {
      ini({
        lka <- log(1.8); lke <- log(0.086); lV <- log(32)
        eta.ka ~ 0.3; eta.ke ~ 0.03; eta.V ~ 0.03
        add.err <- 0.7
      })
      model({
        ka <- exp(lka + eta.ka); ke <- exp(lke + eta.ke); V <- exp(lV + eta.V)
        d/dt(depot) = -ka * depot
        d/dt(central) = ka * depot - ke * central
        cp <- central / V
        cp ~ add(add.err)
      })
    }

    .old <- rxode2::getRxThreads(); on.exit(rxode2::setRxThreads(.old), add = TRUE)
    rxode2::setRxThreads(1L)

    ui <- rxode2::assertRxUi(theo)
    prep <- .vaeDataPrep(ui, nlmixr2data::theo_sd)
    am <- .vaeDecoderModel(ui)
    zDim <- prep$zDim; hDim <- 6L; nCov <- ncol(prep$covIn)
    .testSeed(7)
    params <- .vaeEncoderInitParams(zDim, hDim, nCov, prep$zPop, rep(0.1, zDim))
    .testSeed(123); eps <- matrix(rnorm(prep$N * zDim), prep$N, zDim)
    nObsTot <- sum(vapply(prep$subj, function(s) length(s$y), integer(1)))

    st0 <- .vaeElboStep(params, prep, am, prep$zPop, prep$omega, prep$a, 0.7, eps,
                        withGrad = FALSE)

    cc <- 0.01
    .Call("_nlmixr2est_registerTestContrib", PACKAGE = "nlmixr2est")
    on.exit(.Call("_nlmixr2est_removeTestContrib", PACKAGE = "nlmixr2est"), add = TRUE)
    .Call("_nlmixr2est_setTestContribAddLL", cc, PACKAGE = "nlmixr2est")
    st1 <- .vaeElboStep(params, prep, am, prep$zPop, prep$omega, prep$a, 0.7, eps,
                        withGrad = FALSE)
    res <- .Call("_nlmixr2est_getTestContrib", PACKAGE = "nlmixr2est")
    .Call("_nlmixr2est_removeTestContrib", PACKAGE = "nlmixr2est")

    ## pxz is the per-subject negative log-likelihood, so +c per obs lowers the
    ## ELBO loss by exactly c*nObs
    expect_equal(st1$loss - st0$loss, -cc * nObsTot, tolerance = 1e-8)
    ## the decoder threads id/eta through, so begin/end bracket every subject once
    expect_equal(res[[1]], nObsTot)
    expect_equal(res[[5]], prep$N)
    expect_equal(res[[6]], prep$N)
  })
})
