## Phase 7: error-model parity + general (non-normal) likelihood support. Because
## training drives the FOCEi inner likelihood (likInner0 / foceiInnerLp), every
## error structure and endpoint distribution the inner problem supports works
## through the same interface; only the closed-form error-param M-step is
## error-type specific (.vaeUpdateErr: additive / proportional / combined).

nmTest({
  .vaeErrData <- function(seed = 7, N = 30L) {
    .testSeed(seed)
    simMod <- function() {
      ini({ lka <- log(1.5); lke <- log(0.1); lV <- log(32) })
      model({ ka <- exp(lka + eta.ka); ke <- exp(lke + eta.ke); V <- exp(lV + eta.V)
        d/dt(depot) = -ka * depot; d/dt(central) = ka * depot - ke * central
        cp <- central / V })
    }
    ev <- rxode2::et(amt = 320, cmt = "depot") %>%
      rxode2::et(seq(0.5, 24, length.out = 8)) %>% rxode2::et(id = 1:N)
    base <- rxode2::rxSolve(simMod, ev, params = c(lka = log(1.5), lke = log(0.1), lV = log(32)),
                            omega = lotri::lotri(eta.ka ~ 0.05, eta.ke ~ 0.02, eta.V ~ 0.02))
    base <- as.data.frame(base)[, c("id", "time", "cp")]; names(base) <- c("ID", "TIME", "DV")
    mk <- function(dv) {
      d <- base; d$DV <- dv
      dose <- data.frame(ID = unique(d$ID), TIME = 0, DV = 0, EVID = 1, AMT = 320, CMT = 1)
      d$EVID <- 0; d$AMT <- 0; d$CMT <- 2
      d <- rbind(dose, d); d[order(d$ID, d$TIME, -d$EVID), ]
    }
    f <- base$DV
    list(f = f, mk = mk)
  }

  .vaeRunErr <- function(ui, dat, nMix = 1L, mixProb = 1) {
    ctl <- vaeControl(itersBurnIn = 40L, iters = 90L, klWarmup = 30L, gammaIter = 55L,
                      nGradStep = 4L, covariateSelection = FALSE, seed = 1L)
    prep <- .vaeDataPrep(ui, dat)
    innerEnv <- .vaeInnerSetup(ui, dat, matrix(0, prep$N, prep$zDim), ctl)
    on.exit(.vaeInnerFree(), add = TRUE)
    .vaeTrain(prep, innerEnv, ctl, nMix, mixProb)
  }

  test_that("vae recovers a proportional error model", {
    skip_on_cran()
    d <- .vaeErrData()
    dat <- d$mk(d$f + stats::rnorm(length(d$f), 0, 0.15 * abs(d$f)))
    propMod <- function() {
      ini({ lka <- log(1.5); lke <- log(0.1); lV <- log(32)
        eta.ka ~ 0.05; eta.ke ~ 0.02; eta.V ~ 0.02; prop.err <- 0.1 })
      model({ ka <- exp(lka + eta.ka); ke <- exp(lke + eta.ke); V <- exp(lV + eta.V)
        d/dt(depot) = -ka * depot; d/dt(central) = ka * depot - ke * central
        cp <- central / V; cp ~ prop(prop.err) })
    }
    fit <- .vaeRunErr(rxode2::assertRxUi(propMod), dat)
    expect_equal(names(fit$a), "prop.err")
    expect_lt(abs(fit$a[["prop.err"]] - 0.15), 0.04)
    expect_lt(abs(fit$zPop[3] - log(32)), 0.1)
  })

  test_that("vae recovers a combined error model", {
    skip_on_cran()
    d <- .vaeErrData()
    dat <- d$mk(d$f + stats::rnorm(length(d$f), 0, sqrt(0.2^2 + (0.1 * d$f)^2)))
    combMod <- function() {
      ini({ lka <- log(1.5); lke <- log(0.1); lV <- log(32)
        eta.ka ~ 0.05; eta.ke ~ 0.02; eta.V ~ 0.02; add.err <- 0.2; prop.err <- 0.1 })
      model({ ka <- exp(lka + eta.ka); ke <- exp(lke + eta.ke); V <- exp(lV + eta.V)
        d/dt(depot) = -ka * depot; d/dt(central) = ka * depot - ke * central
        cp <- central / V; cp ~ add(add.err) + prop(prop.err) })
    }
    fit <- .vaeRunErr(rxode2::assertRxUi(combMod), dat)
    expect_setequal(names(fit$a), c("add.err", "prop.err"))
    expect_lt(abs(fit$a[["add.err"]] - 0.2), 0.1)
    expect_lt(abs(fit$a[["prop.err"]] - 0.1), 0.05)
  })

  test_that("vae supports a general (Poisson) likelihood via the inner problem", {
    skip_on_cran()
    .testSeed(11)
    N <- 40L; times <- seq(1, 10, by = 1); lb <- log(8)
    rows <- do.call(rbind, lapply(seq_len(N), function(i) {
      lam <- exp(lb + stats::rnorm(1, 0, 0.5)) * exp(-0.05 * times)
      data.frame(ID = i, TIME = times, DV = stats::rpois(length(times), lam), EVID = 0)
    }))
    poisMod <- function() {
      ini({ lb <- log(5); eta.b ~ 0.1 })
      model({ lambda <- exp(lb + eta.b) * exp(-0.05 * time); cp <- lambda; cp ~ pois(lambda) })
    }
    ui <- rxode2::assertRxUi(poisMod)
    expect_true(all(ui$predDfFocei$distribution != "norm"))  # general likelihood
    ctl <- vaeControl(itersBurnIn = 30L, iters = 70L, klWarmup = 25L, gammaIter = 45L,
                      nGradStep = 4L, covariateSelection = FALSE, seed = 1L)
    prep <- .vaeDataPrep(ui, rows)
    expect_length(prep$a, 0L)                                # no residual error param
    innerEnv <- .vaeInnerSetup(ui, rows, matrix(0, prep$N, prep$zDim), ctl)
    on.exit(.vaeInnerFree(), add = TRUE)
    ## the inner problem supplies a finite objective + eta-gradient for the count LL
    r <- .vaeInnerEval(matrix(stats::rnorm(N, 0, 0.1), N, 1), ctl, grad = TRUE)
    expect_true(all(is.finite(r$obj)) && all(is.finite(r$lp)))
    fit <- .vaeTrain(prep, innerEnv, ctl, 1L, 1)
    expect_true(is.finite(fit$zPop[1]) && is.finite(fit$omega[1]) && fit$omega[1] > 0)
    expect_lt(abs(fit$zPop[1] - lb), 0.25)                   # recovers the rate
  })
})
