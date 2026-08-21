nmTest({
  # Phase 2 (see the SAEM general-likelihood theta plan): rxUiGet.saemPhi1Inner
  # reuses FOCEi's own inner (eta-sensitivity-only) model for SAEM's general-
  # likelihood phi1 theta step -- NOT the augmented outer sensitivity model,
  # which differentiates w.r.t. theta directions too (needed only for FOCEi's
  # own analytic outer GRADIENT; SAEM's phi1 step optimizes theta with bobyqa,
  # derivative-free, so no theta-gradient is ever needed). What IS needed, for
  # each theta bobyqa proposes: push theta and a candidate eta into this inner
  # model to read the log-density (rx_pred_) and its 1st/2nd-order eta
  # sensitivities, used to find the conditional mode and its Laplace
  # correction.  This file verifies those sensitivities against finite
  # differences -- matching output alone is not sufficient evidence the
  # analytic path is correct, so a "ran" check (non-NULL, right shape) is
  # paired with a numeric comparison.
  mLl <- function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45
      lsd <- fixed(log(0.7))
      eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
    })
    model({
      ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      sd <- exp(lsd)
      ll(err) ~ -lsd - 0.5 * log(2 * pi) - 0.5 * ((DV - cp) / sd)^2
    })
  }

  test_that("saemPhi1Inner is NULL for a normal (non-general-lik) endpoint", {
    mNorm <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45
        add.sd <- 0.7
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        d / dt(depot) <- -ka * depot
        d / dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }
    .ui <- rxode2::rxUiDecompress(rxode2::rxode2(mNorm))
    expect_null(.ui$saemPhi1Inner)
  })

  test_that("saemPhi1Inner does not build innerHess2 by default (phi1Hessian=FALSE)", {
    # saemControl(phi1Hessian=FALSE) is the default -- an ablation check
    # found the Laplace log|H| correction was not what fixed a diverging
    # Gaussian twin, and it can dominate/diverge for a heavy-tailed
    # t()/cauchy() endpoint (nlmixr2/nlmixr2est#999), so building the
    # (expensive, foceiControl(fast=TRUE)-only) innerHess2 model is wasted
    # work unless a caller opts in.
    .ui <- rxode2::rxUiDecompress(rxode2::rxode2(mLl))
    assign("control", saemControl(nBurn = 5, nEm = 5), envir = .ui)
    .res <- .ui$saemPhi1Inner
    expect_false(is.null(.res))
    expect_s3_class(.res$inner, "rxode2")
    expect_null(.res$innerHess2)
    expect_s3_class(.res$predNoLhs, "rxode2")
  })

  test_that("saemPhi1Inner builds innerHess2 and preserves the ui's saemControl when phi1Hessian=TRUE (#Phase2)", {
    .ui <- rxode2::rxUiDecompress(rxode2::rxode2(mLl))
    assign("control", saemControl(nBurn = 5, nEm = 5, phi1Hessian = TRUE), envir = .ui)
    .res <- .ui$saemPhi1Inner
    expect_false(is.null(.res))
    expect_s3_class(.res$inner, "rxode2")
    expect_s3_class(.res$innerHess2, "rxode2")
    expect_s3_class(.res$predNoLhs, "rxode2")
    # the accessor must leave the ui's ORIGINAL control installed, not the
    # temporary foceiControl(fast=TRUE) it needs internally to build innerHess2
    # (saemControl() packs nBurn/nEm into $mcmc$niter, not a top-level field)
    expect_s3_class(.ui$control, "saemControl")
    expect_equal(.ui$control$mcmc$niter, c(5, 5))
  })

  test_that("saemPhi1Inner's eta gradient/Hessian match finite differences (#Phase2)", {
    .ui <- rxode2::rxUiDecompress(rxode2::rxode2(mLl))
    assign("control", saemControl(phi1Hessian = TRUE), envir = .ui)
    .res <- .ui$saemPhi1Inner
    .inner <- .res$inner
    .innerHess2 <- .res$innerHess2

    .ev <- rxode2::et(amt = 320, time = 0) %>% rxode2::et(seq(0.5, 24, by = 4))
    .ev$DV <- 5
    .pars <- c("THETA[1]" = 0.45, "THETA[2]" = 1, "THETA[3]" = 3.45,
              "THETA[4]" = log(0.7), "ETA[1]" = 0.1, "ETA[2]" = -0.05, "ETA[3]" = 0.02)
    .h <- 1e-4

    .solveAt <- function(mod, eta1) {
      .p <- .pars; .p["ETA[1]"] <- eta1
      .s <- suppressWarnings(rxode2::rxSolve(mod, .p, .ev, returnType = "data.frame"))
      sum(.s$rx_pred_)
    }
    .fdGrad <- (.solveAt(.inner, .pars[["ETA[1]"]] + .h) -
                  .solveAt(.inner, .pars[["ETA[1]"]] - .h)) / (2 * .h)
    .s0 <- suppressWarnings(rxode2::rxSolve(.inner, .pars, .ev, returnType = "data.frame"))
    .analyticGrad <- sum(.s0$rx__sens_rx_pred__BY_ETA_1___)
    expect_equal(.analyticGrad, .fdGrad, tolerance = 1e-3)

    .solveAtH2 <- function(eta1) {
      .p <- .pars; .p["ETA[1]"] <- eta1
      .s <- suppressWarnings(rxode2::rxSolve(.innerHess2, .p, .ev, returnType = "data.frame"))
      sum(.s$rx__sens_rx_pred__BY_ETA_1___)
    }
    .fdHess <- (.solveAtH2(.pars[["ETA[1]"]] + .h) - .solveAtH2(.pars[["ETA[1]"]] - .h)) / (2 * .h)
    .sH0 <- suppressWarnings(rxode2::rxSolve(.innerHess2, .pars, .ev, returnType = "data.frame"))
    .analyticHess <- sum(.sH0$rx__d2pred_1_1__)
    expect_equal(.analyticHess, .fdHess, tolerance = 1e-3)
  })

  test_that("saemPhi1Inner selects innerLlik (true log-density) for t()/cauchy() sugar (#Phase2)", {
    # add()/dt()/cauchy() sugar's FIRST (non-Llik) inner model is FOCEi's own
    # closed-form Gaussian approximation (what FOCEi's Newton search actually
    # fits); only the "Llik" round is the true log-density -- this is what
    # SAEM's general-likelihood phi1 step must use.
    mT <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45
        add.sd <- fixed(0.7)
        nu <- fixed(30)
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        d / dt(depot) <- -ka * depot
        d / dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd) + dt(nu)
      })
    }
    .ui <- rxode2::rxUiDecompress(rxode2::rxode2(mT))
    .res <- .ui$saemPhi1Inner
    expect_false(is.null(.res))
    .mv <- paste(rxode2::rxModelVars(.res$inner)$model, collapse = "\n")
    expect_true(grepl("llikT", .mv))
  })

  test_that("saemPhi1TargetMap resolves THETA[k]/ETA[k] -> phi column (#Phase4)", {
    .ui <- rxode2::rxUiDecompress(rxode2::rxode2(mLl))
    .res <- .ui$saemPhi1Inner
    expect_true(isTRUE(.res$ok))
    # tka, tcl, tv are mu-referenced (phi1, in declaration order); lsd is
    # fixed() but still a phi0 column (SAEM's own fixedIx0 pins its value,
    # matching refinePhi0Lik's own fixed-theta handling -- it is not dropped
    # from the theta vector, so it must not be dropped from this map either)
    expect_equal(.res$thetaKind, c(1L, 1L, 1L, 0L))
    expect_equal(.res$thetaCol, c(0L, 1L, 2L, 0L))
    expect_equal(.res$etaCol, c(0L, 1L, 2L))
  })

  test_that("saemPhi1TargetMap declines a covariate on a mu-ref theta (#Phase4)", {
    mCov <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45
        tka.wt <- 0.1
        lsd <- fixed(log(0.7))
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
      })
      model({
        ka <- exp(tka + tka.wt * WT + eta.ka)
        cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        d / dt(depot) <- -ka * depot
        d / dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        sd <- exp(lsd)
        ll(err) ~ -lsd - 0.5 * log(2 * pi) - 0.5 * ((DV - cp) / sd)^2
      })
    }
    .ui <- rxode2::rxUiDecompress(rxode2::rxode2(mCov))
    .res <- .ui$saemPhi1Inner
    expect_false(isTRUE(.res$ok))
  })
})
