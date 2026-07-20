## ADVI inner/outer setup wiring: .adviInnerSetup sets up the reused FOCEi inner
## problem (per-subject log-joint + eta-gradient via likInner0/lpInner) and, for
## a model with non-mu structural/sigma thetas, builds the impmap
## theta-sensitivity model (reused for the outer population gradient).

nmTest({
  test_that("advi inner driver returns finite obj/lp and matches likInner/foceiInnerLp", {
    theo <- function() {
      ini({ lka <- log(1.8); lke <- log(0.086); lV <- log(32)
        eta.ka ~ 0.3; eta.ke ~ 0.03; eta.V ~ 0.03; add.err <- 0.7 })
      model({ ka <- exp(lka + eta.ka); ke <- exp(lke + eta.ke); V <- exp(lV + eta.V)
        d/dt(depot) = -ka * depot; d/dt(central) = ka * depot - ke * central
        cp <- central / V; cp ~ add(add.err) })
    }
    ui <- rxode2::assertRxUi(theo)
    ctl <- adviControl()
    N <- length(unique(nlmixr2data::theo_sd$ID))
    .testSeed(1); etaMat <- matrix(rnorm(N * 3, 0, 0.1), N, 3)

    .adviInnerSetup(ui, nlmixr2data::theo_sd, etaMat, ctl)
    on.exit(.adviInnerFree(), add = TRUE)
    r <- .adviInnerEval(etaMat, ctl, grad = TRUE)
    expect_equal(length(r$obj), N)
    expect_true(all(is.finite(r$obj)))
    expect_true(all(is.finite(r$lp)))
    ## match the single-subject R-exported inner functions (reused verbatim)
    expect_lt(abs(r$obj[1] - likInner(etaMat[1, ], 1L)), 1e-6)
    expect_lt(max(abs(r$lp[1, ] - foceiInnerLp(etaMat[1, ], 1L))), 1e-6)
  })

  test_that("advi builds the theta-sensitivity model for non-mu structural thetas", {
    ## tcl, tv are fixed-effect-only (non-mu) -> theta-sensitivity model required
    mixedMu <- function() {
      ini({ lka <- log(1.8); lcl <- log(0.086); lV <- log(32)
        eta.ka ~ 0.3; add.err <- 0.7 })
      model({ ka <- exp(lka + eta.ka); cl <- exp(lcl); V <- exp(lV)
        d/dt(depot) = -ka * depot; d/dt(central) = ka * depot - cl / V * central
        cp <- central / V; cp ~ add(add.err) })
    }
    ui <- rxode2::assertRxUi(mixedMu)
    ctl <- adviControl()
    N <- length(unique(nlmixr2data::theo_sd$ID))
    env <- .adviInnerSetup(ui, nlmixr2data::theo_sd, matrix(0, N, 1L), ctl)
    on.exit(.adviInnerFree(), add = TRUE)
    ## the sensitivity model compiled and is attached (reused outer-gradient source)
    expect_false(is.null(env$model$thetaSens))
    ## it still evaluates a finite inner objective/gradient
    r <- .adviInnerEval(matrix(0, N, 1L), ctl, grad = TRUE)
    expect_true(all(is.finite(r$obj)))
    expect_true(all(is.finite(r$lp)))
  })
})
