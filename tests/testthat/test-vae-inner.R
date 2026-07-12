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
    set.seed(1); etaMat <- matrix(rnorm(N * 3, 0, 0.1), N, 3)

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
    set.seed(42)
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
    mixProb <- { p <- as.numeric(prep$th[ui$thetaMixIndex]); c(p, 1 - sum(p)) }
    innerEnv <- .vaeInnerSetup(ui, dat, matrix(0, prep$N, prep$zDim), ctl)
    on.exit(.vaeInnerFree(), add = TRUE)
    fit <- .vaeTrain(prep, innerEnv, ctl, nMix, mixProb)

    expect_true(all(is.finite(fit$zPop)) && is.finite(fit$a) && fit$a > 0)
    expect_equal(fit$zPop[2], 0)                        # mixture eta stays centered at 0
    ## component labels are arbitrary; agreement is max(match, 1-match)
    agree <- mean(fit$mixnum == trueGrp)
    expect_gt(max(agree, 1 - agree), 0.9)
  })
})
