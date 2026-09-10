## VAE C++ inner-likelihood driver: set up the FOCEi inner problem once
## (.vaeInnerSetup -> vaeInnerSetup_) then evaluate the per-subject objective and
## eta-gradient through likInner0/lpInner in parallel (.vaeInnerEval ->
## vaeInnerLik). Must match the R-exported likInner/foceiInnerLp, and for a
## mixture model must give distinct per-component objectives (hard assignment).

nmTest({
  test_that("vae inner driver matches likInner/foceiInnerLp (theophylline)", {
    theo <- function() {
      ini({ lka <- log(1.8); lke <- log(0.086); lV <- log(32)
        eta.ka ~ 0.3; eta.ke ~ 0.03; eta.V ~ 0.03; add.err <- 0.7 })
      model({ ka <- exp(lka + eta.ka); ke <- exp(lke + eta.ke); V <- exp(lV + eta.V)
        d/dt(depot) = -ka * depot; d/dt(central) = ka * depot - ke * central
        cp <- central / V; cp ~ add(add.err) })
    }
    ui <- rxode2::assertRxUi(theo)
    ctl <- vaeControl()
    N <- length(unique(nlmixr2data::theo_sd$ID))
    .testSeed(1); etaMat <- matrix(rnorm(N * 3, 0, 0.1), N, 3)

    .vaeInnerSetup(ui, nlmixr2data::theo_sd, etaMat, ctl)
    on.exit(.vaeInnerFree(), add = TRUE)
    r <- .vaeInnerEval(etaMat, ctl, grad = TRUE)
    expect_equal(length(r$obj), N)
    ## match the single-subject R-exported inner functions
    expect_lt(abs(r$obj[1] - likInner(etaMat[1, ], 1L)), 1e-6)
    expect_lt(max(abs(r$lp[1, ] - foceiInnerLp(etaMat[1, ], 1L))), 1e-6)
    expect_lt(abs(r$obj[5] - likInner(etaMat[5, ], 5L)), 1e-6)
  })

  test_that("vae likelihood choices map to the right FOCEi inner control", {
    ## focei -> interaction; foce -> NONMEM FOCE; focep -> FOCE+ (R at live eta)
    expect_identical(formals(vaeControl)$likelihood,
                     quote(c("focei", "foce", "focep", "laplace")))
    fi <- .vaeInnerFoceiControl(vaeControl(likelihood = "focei"))
    fe <- .vaeInnerFoceiControl(vaeControl(likelihood = "foce"))
    fp <- .vaeInnerFoceiControl(vaeControl(likelihood = "focep"))
    expect_equal(fi$interaction, 1L)
    expect_equal(fe$interaction, 0L)
    expect_equal(fp$interaction, 0L)
    expect_equal(fe$foceType, 0L)   # NONMEM FOCE
    expect_equal(fp$foceType, 1L)   # FOCE+
  })

  test_that("vaeInnerUpdatePar_ fast path matches the full re-setup path", {
    ## The per-gradient-step fast path (vaeInnerUpdatePar_: updateTheta on the
    ## cached reduced par vector) replaces the full re-setup (.vaeInnerUpdate ->
    ## rxSymInvCholCreate + foceiSetup_).  At the SETUP parameter values the two
    ## agree to machine precision; away from them a tiny (<1e-2) discrepancy
    ## remains because foceiSetup_ recomputes the internal scaling-normalization
    ## reference (c1/c2/scaleTo, a function of initPar) each call whereas the fast
    ## path holds the setup's -- a change of internal reparameterization only,
    ## which does not affect the inner likelihood beyond that tolerance.
    theo <- function() {
      ini({ lka <- log(1.8); lke <- log(0.086); lV <- log(32)
        eta.ka ~ 0.3; eta.ke ~ 0.03; eta.V ~ 0.03; add.err <- 0.7 })
      model({ ka <- exp(lka + eta.ka); ke <- exp(lke + eta.ke); V <- exp(lV + eta.V)
        d/dt(depot) = -ka * depot; d/dt(central) = ka * depot - ke * central
        cp <- central / V; cp ~ add(add.err) })
    }
    ui <- rxode2::assertRxUi(theo)
    ctl <- vaeControl()
    N <- length(unique(nlmixr2data::theo_sd$ID))
    .testSeed(3); etaMat <- matrix(rnorm(N * 3, 0, 0.1), N, 3)
    prep <- .vaeDataPrep(ui, nlmixr2data::theo_sd)
    env <- .vaeInnerSetup(ui, nlmixr2data::theo_sd, etaMat, ctl)
    on.exit(.vaeInnerFree(), add = TRUE)
    ## exact at the setup parameter values
    vaeInnerUpdatePar_(as.numeric(prep$th), as.numeric(prep$omega))
    rFast0 <- .vaeInnerEval(etaMat, ctl, grad = TRUE)
    .vaeInnerUpdate(env, prep$th, prep$omega, etaMat)
    rRef0 <- .vaeInnerEval(etaMat, ctl, grad = TRUE)
    expect_equal(rFast0$obj, rRef0$obj, tolerance = 1e-10)
    expect_equal(rFast0$lp, rRef0$lp, tolerance = 1e-10)
    ## near the setup values, agreement to a small tolerance
    for (i in 1:3) {
      .testSeed(i)
      th <- prep$th * (1 + rnorm(length(prep$th), 0, 0.1))
      om <- prep$omega * exp(rnorm(3, 0, 0.3))
      vaeInnerUpdatePar_(as.numeric(th), as.numeric(om))
      rFast <- .vaeInnerEval(etaMat, ctl, grad = TRUE)
      .vaeInnerUpdate(env, th, om, etaMat)
      rRef <- .vaeInnerEval(etaMat, ctl, grad = TRUE)
      expect_lt(max(abs(rFast$obj - rRef$obj)), 1e-2)
      expect_lt(max(abs(rFast$lp - rRef$lp)), 1e-2)
    }
  })

  test_that("vaeInnerUpdatePar_ fast path matches re-setup for a CORRELATED omega", {
    ## The correlated branch packs chol(Omega^-1) onto the omega block by the
    ## position list from the model structure; a row/column swap or a dropped
    ## diagonal sqrt() would leave the fast path disagreeing with the full
    ## rxSymInvCholCreate + foceiSetup_ re-setup.  A 3-eta block exercises a
    ## packing order the diagonal case cannot.
    theoCor <- function() {
      ini({ lka <- log(1.8); lke <- log(0.086); lV <- log(32)
        eta.ka + eta.ke + eta.V ~ c(0.3,
                                    0.01, 0.03,
                                    0.02, 0.005, 0.03)
        add.err <- 0.7 })
      model({ ka <- exp(lka + eta.ka); ke <- exp(lke + eta.ke); V <- exp(lV + eta.V)
        d/dt(depot) = -ka * depot; d/dt(central) = ka * depot - ke * central
        cp <- central / V; cp ~ add(add.err) })
    }
    ui <- rxode2::assertRxUi(theoCor)
    ctl <- vaeControl()
    N <- length(unique(nlmixr2data::theo_sd$ID))
    .testSeed(3); etaMat <- matrix(rnorm(N * 3, 0, 0.1), N, 3)
    prep <- .vaeDataPrep(ui, nlmixr2data::theo_sd)
    expect_true(.omegaHasOffDiag(prep$omegaMat))
    env <- .vaeInnerSetup(ui, nlmixr2data::theo_sd, etaMat, ctl)
    on.exit(.vaeInnerFree(), add = TRUE)
    ## exact at the setup parameter values
    vaeInnerUpdatePar_(as.numeric(prep$th), prep$omegaMat)
    rFast0 <- .vaeInnerEval(etaMat, ctl, grad = TRUE)
    .vaeInnerUpdate(env, prep$th, prep$omegaMat, etaMat)
    rRef0 <- .vaeInnerEval(etaMat, ctl, grad = TRUE)
    expect_equal(rFast0$obj, rRef0$obj, tolerance = 1e-10)
    expect_equal(rFast0$lp, rRef0$lp, tolerance = 1e-10)
    ## and away from them, with the correlation itself perturbed
    for (i in 1:3) {
      .testSeed(i)
      th <- prep$th * (1 + rnorm(length(prep$th), 0, 0.1))
      om <- prep$omegaMat
      diag(om) <- diag(om) * exp(rnorm(3, 0, 0.2))
      om[1L, 2L] <- om[2L, 1L] <- om[1L, 2L] * (1 + rnorm(1, 0, 0.2))
      vaeInnerUpdatePar_(as.numeric(th), om)
      rFast <- .vaeInnerEval(etaMat, ctl, grad = TRUE)
      .vaeInnerUpdate(env, th, om, etaMat)
      rRef <- .vaeInnerEval(etaMat, ctl, grad = TRUE)
      expect_lt(max(abs(rFast$obj - rRef$obj)), 1e-2)
      expect_lt(max(abs(rFast$lp - rRef$lp)), 1e-2)
    }
    ## a wrongly sized omega is rejected, not read out of bounds
    expect_error(vaeInnerUpdatePar_(as.numeric(prep$th), matrix(1, 1, 1)),
                 "expected")
  })

  test_that("vae inner driver selects mixture components per id", {
    mixmod <- function() {
      ini({ lka <- log(1.8); lke1 <- log(0.15); lke2 <- log(0.04); lV <- log(32); p1 <- 0.6
        eta.ka ~ 0.3; eta.ke ~ 0.03; eta.V ~ 0.03; add.err <- 0.7 })
      model({ ka <- exp(lka + eta.ka)
        ke <- mix(exp(lke1 + eta.ke), p1, exp(lke2 + eta.ke)); V <- exp(lV + eta.V)
        d/dt(depot) = -ka * depot; d/dt(central) = ka * depot - ke * central
        cp <- central / V; cp ~ add(add.err) })
    }
    ui <- rxode2::assertRxUi(mixmod)
    expect_equal(ui$saemNMix, 2L)
    ctl <- vaeControl()
    N <- length(unique(nlmixr2data::theo_sd$ID)); nMix <- ui$saemNMix
    etaSetup <- matrix(0, N, 3)

    .vaeInnerSetup(ui, nlmixr2data::theo_sd, etaSetup, ctl)
    on.exit(.vaeInnerFree(), add = TRUE)
    ## nSub*nMix ids, component-major (encoder etas repeated per component)
    etaEval <- do.call(rbind, rep(list(etaSetup), nMix))
    r <- .vaeInnerEval(etaEval, ctl)
    expect_equal(length(r$obj), N * nMix)
    o1 <- r$obj[1:N]; o2 <- r$obj[(N + 1):(2 * N)]
    ## the two components give distinct per-subject objectives (mixture selection);
    ## the very slow component (ke=0.04) at eta=0 can occasionally fail to solve, so
    ## assert on the subjects that solved rather than requiring every one
    ok <- is.finite(o1) & is.finite(o2)
    expect_gt(sum(ok), N / 2)
    expect_true(all(abs(o1[ok] - o2[ok]) > 1e-6))
    ## hard-assignment mixnum is a valid component per subject
    mixnum <- ifelse(o1[ok] <= o2[ok], 1L, 2L)
    expect_true(all(mixnum %in% c(1L, 2L)))
  })

  test_that("vae inner-driver training recovers mixture assignment", {
    skip_on_cran()
    ## two well-separated subpopulations (fast vs slow ke); training via the inner
    ## driver must recover the population parameters AND assign every subject to its
    ## true component (mixnum) -- the FOCEi feature-parity milestone.
    sim <- function() {
      ini({ lka <- log(1.5); lV <- log(32) })
      model({ ka <- exp(lka + eta.ka); ke <- KE * exp(eta.ke); V <- exp(lV + eta.V)
        d/dt(depot) = -ka * depot; d/dt(central) = ka * depot - ke * central
        cp <- central / V })
    }
    .testSeed(42)
    nPer <- 20L
    ev <- rxode2::et(amt = 320, cmt = "depot") %>% rxode2::et(seq(0.5, 24, length.out = 8))
    mkGroup <- function(ke, ids) {
      d <- rxode2::rxSolve(sim, rxode2::et(ev, id = ids),
                           params = c(lka = log(1.5), lV = log(32), KE = ke),
                           omega = lotri::lotri(eta.ka ~ 0.04, eta.ke ~ 0.02, eta.V ~ 0.02))
      d <- as.data.frame(d)[, c("id", "time", "cp")]
      names(d) <- c("ID", "TIME", "DV"); d$ID <- d$ID + (ids[1] - 1)
      d$DV <- d$DV + stats::rnorm(nrow(d), 0, 0.25); d
    }
    dat <- rbind(mkGroup(0.15, 1:nPer), mkGroup(0.04, (nPer + 1):(2 * nPer)))
    dose <- data.frame(ID = unique(dat$ID), TIME = 0, DV = 0, EVID = 1, AMT = 320, CMT = 1)
    dat$EVID <- 0; dat$AMT <- 0; dat$CMT <- 2
    dat <- rbind(dose, dat); dat <- dat[order(dat$ID, dat$TIME, -dat$EVID), ]
    trueGrp <- ifelse(unique(dat$ID) <= nPer, 1L, 2L)

    mixmod <- function() {
      ini({ lka <- log(1.5); lke1 <- log(0.15); lke2 <- log(0.04); lV <- log(32); p1 <- 0.5
        eta.ka ~ 0.04; eta.ke ~ 0.02; eta.V ~ 0.02; add.err <- 0.25 })
      model({ ka <- exp(lka + eta.ka)
        ke <- mix(exp(lke1 + eta.ke), p1, exp(lke2 + eta.ke)); V <- exp(lV + eta.V)
        d/dt(depot) = -ka * depot; d/dt(central) = ka * depot - ke * central
        cp <- central / V; cp ~ add(add.err) })
    }
    ui <- rxode2::assertRxUi(mixmod)
    ctl <- vaeControl(itersBurnIn = 40L, iters = 100L, klWarmup = 30L, gammaIter = 60L,
                      nGradStep = 4L, covariateSelection = FALSE, seed = 1L)
    prep <- .vaeDataPrep(ui, dat)
    nMix <- as.integer(ui$saemNMix)
    ## prep$th holds the mixture slots on the MLOGIT scale (the scale the inner
    ## problem reads them on); mexpit back, exactly as .vaeFitModel does
    mixProb <- .getMixFromLog(prep$th, ui$thetaMixIndex)
    expect_equal(mixProb, c(0.5, 0.5))
    innerEnv <- .vaeInnerSetup(ui, dat, matrix(0, prep$N, prep$zDim), ctl)
    on.exit(.vaeInnerFree(), add = TRUE)
    fit <- .vaeTrain(prep, innerEnv, ctl, nMix, mixProb)

    expect_true(all(is.finite(fit$zPop)) && is.finite(fit$a) && fit$a > 0)
    expect_equal(fit$zPop[2], 0)                        # mixture eta stays centered at 0
    ## component labels are arbitrary; agreement is max(match, 1-match)
    agree <- mean(fit$mixnum == trueGrp)
    expect_gt(max(agree, 1 - agree), 0.9)
  })

  test_that("vaeElboStepCpp_ ELBO step: shapes, and encoder gradient matches FD", {
    ## The C++ ELBO core (vaeElboStepCpp_, the same one vaeTrainCpp_ drives) exposed
    ## to R via .vaeElboStepInner.  Validate its structure and that the analytic
    ## encoder-parameter gradient it returns matches a finite-difference gradient of
    ## the ELBO loss.
    theo <- function() {
      ini({ lka <- log(1.8); lke <- log(0.086); lV <- log(32)
        eta.ka ~ 0.3; eta.ke ~ 0.03; eta.V ~ 0.03; add.err <- 0.7 })
      model({ ka <- exp(lka + eta.ka); ke <- exp(lke + eta.ke); V <- exp(lV + eta.V)
        d/dt(depot) = -ka * depot; d/dt(central) = ka * depot - ke * central
        cp <- central / V; cp ~ add(add.err) })
    }
    ui <- rxode2::assertRxUi(theo)
    # sigdig = 6 is REQUIRED for the finite-difference check below, not cosmetic.
    # The loss contains an ODE solve, so it carries a relative noise floor set by
    # the solver tolerances (sigdig drives rtol/atol).  A central difference's
    # roundoff error is noise/h, so at the default sigdig the h = 1e-5 step below
    # is far inside the noise and the FD "reference" is meaningless -- measured
    # relative error against the analytic gradient, same seed and step:
    #
    #   sigdig default   h=1e-3 1.0e-4 | h=1e-4 1.0e-3 | h=1e-5 1.0e-2   (~1e-7/h)
    #   sigdig = 6       h=1e-3 1.1e-5 | h=1e-4 9.7e-8 | h=1e-5 5.8e-8
    #
    # i.e. the error scales as 1/h at the default (roundoff-dominated) and becomes
    # a clean U-shape with its minimum at this step once the solve is tight enough.
    # The analytic gradient is right either way -- it agrees to 6e-8 here.
    ctl <- vaeControl(sigdig = 6)
    prep <- .vaeDataPrep(ui, nlmixr2data::theo_sd)
    N <- prep$N; zDim <- prep$zDim; hDim <- 12L
    innerEnv <- .vaeInnerSetup(ui, nlmixr2data::theo_sd, matrix(0, N, zDim), ctl)
    on.exit(.vaeInnerFree(), add = TRUE)
    .testSeed(1)
    ## nCov from the prep: the encoder head is [hDim + nCov] wide because the
    ## encoder is conditioned on the covariates
    params <- .vaeEncoderInitParams(zDim, hDim, ncol(prep$covIn), prep$zPop, rep(0.1, zDim))
    eps <- matrix(rnorm(N * zDim), N, zDim)
    st <- .vaeElboStepInner(params, prep, innerEnv, prep$zPop, prep$omega, prep$a, 1, eps, ctl)

    expect_true(is.finite(st$loss) && is.finite(st$pxz) && is.finite(st$DKL))
    expect_equal(st$loss, st$pxz + st$DKL)              # alphaKL = 1
    expect_equal(dim(st$mu), c(N, zDim))
    expect_equal(dim(st$z), c(N, zDim))
    expect_equal(dim(st$L), c(zDim, zDim, N))
    expect_setequal(names(st$grads), c("Wih", "Whh", "bih", "bhh", "fcW", "fcB"))
    expect_length(st$preds, N)
    expect_true(all(st$mixnum == 1L))                  # single component

    ## finite-difference check of dLoss/d(fcB) (a small, well-conditioned block)
    Lf <- function(p) .vaeElboStepInner(p, prep, innerEnv, prep$zPop, prep$omega,
                                        prep$a, 1, eps, ctl, withGrad = FALSE)$loss
    h <- 1e-5
    fd <- vapply(seq_along(params$fcB), function(j) {
      pp <- params; pp$fcB[j] <- pp$fcB[j] + h
      pm <- params; pm$fcB[j] <- pm$fcB[j] - h
      (Lf(pp) - Lf(pm)) / (2 * h)
    }, numeric(1))
    expect_lt(max(abs(fd - st$grads$fcB)) / max(abs(st$grads$fcB)), 1e-3)
  })

  test_that("the mixture ELBO is the marginal -2LL, at the right temperature", {
    ## obj is -log p(y_i, eta_i | m) at 1x scale, so the marginal is
    ## -sum_i log sum_m pi_m exp(-obj_im).  The code used exp(-0.5*obj) and
    ## negated twice, i.e. it marginalized a SQUARE ROOT likelihood and carried
    ## a spurious -log(pi_best).  Both errors cancel exactly when nMix == 1, and
    ## also when the components are identical AND pi is uniform -- which is why
    ## p1 here is 0.7 and not 0.5.  Nothing else in the suite could see this.
    mixmod <- function() {
      ini({ lka <- log(1.5); lke1 <- log(0.15); lke2 <- log(0.04); lV <- log(32)
            p1 <- 0.7
            eta.ka ~ 0.04; eta.ke ~ 0.02; eta.V ~ 0.02; add.err <- 0.25 })
      model({ ka <- exp(lka + eta.ka)
        ke <- mix(exp(lke1 + eta.ke), p1, exp(lke2 + eta.ke)); V <- exp(lV + eta.V)
        d/dt(depot) = -ka * depot; d/dt(central) = ka * depot - ke * central
        cp <- central / V; cp ~ add(add.err) })
    }
    ui <- rxode2::assertRxUi(mixmod)
    ctl <- vaeControl(itersBurnIn = 2L, iters = 2L, covariateSelection = FALSE, seed = 1L)
    prep <- .vaeDataPrep(ui, nlmixr2data::theo_sd)
    N <- prep$N; zDim <- prep$zDim
    nMix <- as.integer(ui$saemNMix)
    mixProb <- .getMixFromLog(prep$th, ui$thetaMixIndex)
    expect_equal(mixProb, c(0.7, 0.3))

    innerEnv <- .vaeInnerSetup(ui, nlmixr2data::theo_sd, matrix(0, N, zDim), ctl)
    on.exit(.vaeInnerFree(), add = TRUE)
    .testSeed(7)
    ## the encoder is conditioned on the component, so the head carries nMix
    ## extra one-hot inputs alongside the covariates
    params <- .vaeEncoderInitParams(zDim, 12L, ncol(prep$covIn) + nMix, prep$zPop,
                                    rep(0.1, zDim))
    eps <- matrix(rnorm(N * zDim), N, zDim)
    st <- .vaeElboStepInner(params, prep, innerEnv, prep$zPop, prep$omega, prep$a,
                            1, eps, ctl, nMix, mixProb, withGrad = FALSE)

    ## pxz = jointTot - sum(pzI); rebuild pzI from the returned z to recover the
    ## mixture term the step actually computed
    eta <- sweep(st$z, 2, prep$zPop, "-")
    .om <- if (is.matrix(prep$omega)) prep$omega else
      diag(as.numeric(prep$omega), nrow = length(prep$omega))
    pzI <- 0.5 * (rowSums((eta %*% solve(.om)) * eta) +
                    as.numeric(determinant(.om, logarithm = TRUE)$modulus) +
                    zDim * log(2 * pi))
    .jointTot <- st$pxz + sum(pzI)

    ## the same quantity, recomputed in R from the raw per-component objectives
    .obj <- matrix(.vaeInnerEval(do.call(rbind, rep(list(eta), nMix)), ctl)$obj, nrow = N)
    .ll <- sweep(-.obj, 2, log(mixProb), "+")
    .mm <- apply(.ll, 1, max)
    .want <- -sum(.mm + log(rowSums(exp(.ll - .mm))))
    expect_equal(.jointTot, .want, tolerance = 1e-8)

    ## Every (subject, component) pair gets its OWN posterior from the encoder.
    ## The components used to be scored at a single shared eta tiled across them
    ## -- an eta fitted to none of them.  This is the assertion that would
    ## regress if the tiling came back.
    expect_equal(dim(st$muAll), c(N * nMix, zDim))
    expect_equal(dim(st$mu), c(N, zDim))
    .m1 <- st$muAll[seq_len(N), , drop = FALSE]
    .m2 <- st$muAll[N + seq_len(N), , drop = FALSE]
    expect_false(isTRUE(all.equal(.m1, .m2)))
    ## and what is reported per subject is the SELECTED component's row
    .sel <- (st$mixnum - 1L) * N + seq_len(N)
    expect_equal(st$mu, st$muAll[.sel, , drop = FALSE])
    ## responsibilities are a distribution over components
    expect_length(st$mixW, nMix)
    expect_equal(sum(st$mixW), 1)

    ## The mixture proportion is estimated on the MLOGIT scale through its own
    ## analytic gradient -- the same chain focei uses (mixGrad): the per-subject
    ## responsibility difference against the last (non-free) component, times
    ## the dmexpit Jacobian.  Check it against finite differences of the term it
    ## is the gradient of.
    .stepAt <- function(prp, alpha = 0) {
      .vaeElboStepInner(params, prp, innerEnv, prp$zPop, prp$omega, prp$a,
                        alpha, eps, ctl, nMix, mixProb, withGrad = TRUE)
    }
    .g <- .stepAt(prep)
    expect_length(.g$gMixTheta, nMix - 1L)
    .h <- 1e-5
    .pp <- prep; .pp$th[ui$thetaMixIndex] <- prep$th[ui$thetaMixIndex] + .h
    .pm <- prep; .pm$th[ui$thetaMixIndex] <- prep$th[ui$thetaMixIndex] - .h
    .fd <- (.stepAt(.pp)$pxz - .stepAt(.pm)$pxz) / (2 * .h)
    expect_equal(as.numeric(.g$gMixTheta), .fd, tolerance = 1e-4)
    ## and the proportions the step used are the ones the inner problem holds,
    ## on the simplex
    expect_equal(sum(.g$mixProb), 1)
    expect_equal(as.numeric(.g$mixProb), mixProb, tolerance = 1e-5)

    ## The encoder gradient under a mixture must be the gradient of the loss the
    ## step actually reports -- the marginal data term, the SELECTED component's
    ## prior and KL.  Finite-difference it, which settles at once whether the
    ## prior correction and the KL are applied at the right rows and scaled the
    ## right way.
    .lossAt <- function(pp) {
      .vaeElboStepInner(pp, prep, innerEnv, prep$zPop, prep$omega, prep$a,
                        1, eps, ctl, nMix, mixProb, withGrad = FALSE)$loss
    }
    .withG <- .vaeElboStepInner(params, prep, innerEnv, prep$zPop, prep$omega,
                                prep$a, 1, eps, ctl, nMix, mixProb, withGrad = TRUE)
    .hh <- 1e-6
    .anaB <- as.numeric(.withG$grads$fcB)
    .fdB <- vapply(seq_along(.anaB), function(j) {
      .pp <- params; .pp$fcB[j] <- params$fcB[j] + .hh
      .pm <- params; .pm$fcB[j] <- params$fcB[j] - .hh
      (.lossAt(.pp) - .lossAt(.pm)) / (2 * .hh)
    }, numeric(1))
    expect_equal(.anaB, .fdB, tolerance = 1e-3)

    ## and it is NOT the square-root marginalization the code used to compute
    .ll2 <- sweep(-0.5 * .obj, 2, log(mixProb), "+")
    .mm2 <- apply(.ll2, 1, max)
    .old <- -2 * sum(.mm2 + log(rowSums(exp(.ll2 - .mm2))))
    expect_false(isTRUE(all.equal(.want, .old)))
    expect_false(isTRUE(all.equal(.jointTot, .old)))
  })

  test_that("vae estimates the mixture proportion", {
    skip_on_cran()
    ## Two subpopulations in a 3:1 split, started from 0.5.  The proportion was
    ## previously either held at its ini() value or moved by the bobyqa regress
    ## step against (-Inf, Inf) on the wrong scale; it is now estimated on the
    ## mlogit scale by its own analytic gradient, consumed by the same Adam
    ## loop that trains the encoder.
    sim <- function() {
      ini({ lka <- log(1.5); lV <- log(32) })
      model({ ka <- exp(lka + eta.ka); ke <- KE * exp(eta.ke); V <- exp(lV + eta.V)
        d/dt(depot) = -ka * depot; d/dt(central) = ka * depot - ke * central
        cp <- central / V })
    }
    .testSeed(42)
    nFast <- 30L; nSlow <- 10L
    ev <- rxode2::et(amt = 320, cmt = "depot") %>% rxode2::et(seq(0.5, 24, length.out = 8))
    mkGroup <- function(ke, ids) {
      d <- rxode2::rxSolve(sim, rxode2::et(ev, id = ids),
                           params = c(lka = log(1.5), lV = log(32), KE = ke),
                           omega = lotri::lotri(eta.ka ~ 0.04, eta.ke ~ 0.02, eta.V ~ 0.02))
      d <- as.data.frame(d)[, c("id", "time", "cp")]
      names(d) <- c("ID", "TIME", "DV"); d$ID <- d$ID + (ids[1] - 1)
      d$DV <- d$DV + stats::rnorm(nrow(d), 0, 0.25); d
    }
    dat <- rbind(mkGroup(0.15, 1:nFast), mkGroup(0.04, (nFast + 1):(nFast + nSlow)))
    dose <- data.frame(ID = unique(dat$ID), TIME = 0, DV = 0, EVID = 1, AMT = 320, CMT = 1)
    dat$EVID <- 0; dat$AMT <- 0; dat$CMT <- 2
    dat <- rbind(dose, dat); dat <- dat[order(dat$ID, dat$TIME, -dat$EVID), ]
    trueGrp <- ifelse(unique(dat$ID) <= nFast, 1L, 2L)

    mixmod <- function() {
      ini({ lka <- log(1.5); lke1 <- log(0.15); lke2 <- log(0.04); lV <- log(32)
            p1 <- 0.5                       # deliberately NOT the truth (0.75)
            eta.ka ~ 0.04; eta.ke ~ 0.02; eta.V ~ 0.02; add.err <- 0.25 })
      model({ ka <- exp(lka + eta.ka)
        ke <- mix(exp(lke1 + eta.ke), p1, exp(lke2 + eta.ke)); V <- exp(lV + eta.V)
        d/dt(depot) = -ka * depot; d/dt(central) = ka * depot - ke * central
        cp <- central / V; cp ~ add(add.err) })
    }
    ui <- rxode2::assertRxUi(mixmod)
    ctl <- vaeControl(itersBurnIn = 40L, iters = 100L, klWarmup = 30L, gammaIter = 60L,
                      nGradStep = 4L, covariateSelection = FALSE, seed = 1L)
    prep <- .vaeDataPrep(ui, dat)
    nMix <- as.integer(ui$saemNMix)
    mixProb <- .getMixFromLog(prep$th, ui$thetaMixIndex)
    expect_equal(mixProb, c(0.5, 0.5))
    innerEnv <- .vaeInnerSetup(ui, dat, matrix(0, prep$N, prep$zDim), ctl)
    on.exit(.vaeInnerFree(), add = TRUE)
    fit <- .vaeTrain(prep, innerEnv, ctl, nMix, mixProb)

    ## the gradient actually ran -- a held proportion would report 0 steps and
    ## still look plausible
    expect_true(fit$nMixThetaStep > 0L)
    expect_length(fit$mixProb, nMix)
    expect_equal(sum(fit$mixProb), 1, tolerance = 1e-8)
    expect_true(all(fit$mixProb > 0 & fit$mixProb < 1))
    ## it MOVED off its start, and toward the truth
    expect_false(isTRUE(all.equal(as.numeric(fit$mixProb), c(0.5, 0.5))))
    .agree <- mean(fit$mixnum == trueGrp)
    .p1 <- if (.agree >= 0.5) fit$mixProb[1] else fit$mixProb[2]   # allow label swap
    expect_equal(.p1, nFast / (nFast + nSlow), tolerance = 0.1)
    expect_gt(max(.agree, 1 - .agree), 0.9)

    ## a fix()ed proportion is NOT estimated -- it stays exactly where ini() put
    ## it, even though the gradient machinery still runs
    mixFixed <- function() {
      ini({ lka <- log(1.5); lke1 <- log(0.15); lke2 <- log(0.04); lV <- log(32)
            p1 <- fix(0.5)
            eta.ka ~ 0.04; eta.ke ~ 0.02; eta.V ~ 0.02; add.err <- 0.25 })
      model({ ka <- exp(lka + eta.ka)
        ke <- mix(exp(lke1 + eta.ke), p1, exp(lke2 + eta.ke)); V <- exp(lV + eta.V)
        d/dt(depot) = -ka * depot; d/dt(central) = ka * depot - ke * central
        cp <- central / V; cp ~ add(add.err) })
    }
    .uiF <- rxode2::assertRxUi(mixFixed)
    .ctlF <- vaeControl(itersBurnIn = 5L, iters = 10L, klWarmup = 3L, gammaIter = 6L,
                        nGradStep = 2L, covariateSelection = FALSE, seed = 1L, print = 0L)
    .prepF <- .vaeDataPrep(.uiF, dat)
    .mpF <- .getMixFromLog(.prepF$th, .uiF$thetaMixIndex)
    .envF <- .vaeInnerSetup(.uiF, dat, matrix(0, .prepF$N, .prepF$zDim), .ctlF)
    .fitF <- .vaeTrain(.prepF, .envF, .ctlF, nMix, .mpF)
    expect_equal(as.numeric(.fitF$mixProb), c(0.5, 0.5))
  })

  test_that("the mixture proportion gradient is exact beyond two components", {
    ## mexpit's Jacobian is DENSE -- d(pi_m)/d(theta_l) = pi_m (delta_ml - pi_l)
    ## -- so carrying d/d(pi) and applying one scalar per free parameter (which
    ## is what mixGrad does) drops the off-diagonal terms.  That is exact at
    ## nMix == 2, where there is a single free parameter, and wrong by 25% at
    ## nMix == 3.  Three components is the smallest case that can see it.
    mixmod3 <- function() {
      ini({ lka <- log(1.5); lke1 <- log(0.15); lke2 <- log(0.08); lke3 <- log(0.04)
            lV <- log(32); p1 <- 0.5; p2 <- 0.3
            eta.ka ~ 0.04; eta.ke ~ 0.02; eta.V ~ 0.02; add.err <- 0.25 })
      model({ ka <- exp(lka + eta.ka)
        ke <- mix(exp(lke1 + eta.ke), p1, exp(lke2 + eta.ke), p2, exp(lke3 + eta.ke))
        V <- exp(lV + eta.V)
        d/dt(depot) = -ka * depot; d/dt(central) = ka * depot - ke * central
        cp <- central / V; cp ~ add(add.err) })
    }
    ui <- rxode2::assertRxUi(mixmod3)
    ctl <- vaeControl(itersBurnIn = 2L, iters = 2L, covariateSelection = FALSE, seed = 1L)
    prep <- .vaeDataPrep(ui, nlmixr2data::theo_sd)
    N <- prep$N; zDim <- prep$zDim
    nMix <- as.integer(ui$saemNMix)
    expect_equal(nMix, 3L)
    mixProb <- .getMixFromLog(prep$th, ui$thetaMixIndex)
    expect_equal(sum(mixProb), 1)

    innerEnv <- .vaeInnerSetup(ui, nlmixr2data::theo_sd, matrix(0, N, zDim), ctl)
    on.exit(.vaeInnerFree(), add = TRUE)
    .testSeed(11)
    params <- .vaeEncoderInitParams(zDim, 12L, ncol(prep$covIn) + nMix, prep$zPop,
                                    rep(0.1, zDim))
    eps <- matrix(rnorm(N * zDim), N, zDim)
    .stepAt <- function(prp) {
      .vaeElboStepInner(params, prp, innerEnv, prp$zPop, prp$omega, prp$a,
                        0, eps, ctl, nMix, mixProb, withGrad = TRUE)
    }
    .g <- .stepAt(prep)
    expect_length(.g$gMixTheta, nMix - 1L)

    .idx <- ui$thetaMixIndex
    .h <- 1e-5
    .fd <- vapply(seq_along(.idx), function(j) {
      .pp <- prep; .pp$th[.idx[j]] <- prep$th[.idx[j]] + .h
      .pm <- prep; .pm$th[.idx[j]] <- prep$th[.idx[j]] - .h
      (.stepAt(.pp)$pxz - .stepAt(.pm)$pxz) / (2 * .h)
    }, numeric(1))
    expect_equal(as.numeric(.g$gMixTheta), .fd, tolerance = 1e-4)

    ## the gradient is N*(pi_l - mean responsibility), with no Jacobian involved
    expect_equal(as.numeric(.g$gMixTheta),
                 N * (as.numeric(.g$mixProb)[seq_along(.idx)] -
                        as.numeric(.g$mixW)[seq_along(.idx)]),
                 tolerance = 1e-8)
  })
})
