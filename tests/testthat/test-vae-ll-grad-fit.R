## est="vae" with a log-likelihood endpoint: the two M-step gaps left open when
## the analytic fast=TRUE ll() gradient landed.
##
## residOptimize="twoStage" stage 2 was gated on "has a slot in `a`", and an
## ll() model has no error-parameter rows, so stage 2 never ran and "twoStage"
## degraded to the joint "optimize" solve.  These fits assert stage 2 runs and
## actually moves the log-density-only theta.
##
## (The nonMuTheta="grad" ll() fits that belong here are NOT present: that path is
## gated off in .vaeGradInScope because the "grad" pool arrangement corrupts the
## heap for an ll() model.  Add them when that is fixed.)
##
## Fit-based and multi-iteration: weekly batch, not the push/PR subset.

nmTest({
  ## Gaussian likelihood written out as ll(), which makes lsd a PLAIN theta
  ## (no `err` row) that appears ONLY in the log-density -- the stage-2 case.
  ## lka/lcl feed d/dt, so they stay in stage 1.
  .llMod <- function() {
    ini({ lka <- 0.4; lcl <- -0.9; lb <- log(3); lsd <- log(1.2); eta.b ~ 0.1 })
    model({
      ka <- exp(lka); cl <- exp(lcl)
      d / dt(depot) <- -ka * depot
      d / dt(central) <- ka * depot - cl * central
      mu <- exp(lb + eta.b) * central
      sd <- exp(lsd)
      ll(cp) ~ -0.5 * log(2 * pi) - log(sd) - 0.5 * ((DV - mu) / sd)^2
    })
  }

  .mkDat <- function(seed = 7L, N = 24L) {
    .testSeed(seed)
    .t <- c(0.5, 1, 2, 4, 6, 8, 12, 24)
    do.call(rbind, lapply(seq_len(N), function(i) {
      .e <- data.frame(ID = i, TIME = 0, DV = NA_real_, EVID = 1, AMT = 100)
      .b <- exp(log(3) + stats::rnorm(1, 0, 0.3))
      .ka <- exp(0.4); .cl <- exp(-0.9)
      .cen <- 100 * .ka / (.ka - .cl) * (exp(-.cl * .t) - exp(-.ka * .t))
      .o <- data.frame(ID = i, TIME = .t,
                       DV = .b * .cen + stats::rnorm(length(.t), 0, 1.2),
                       EVID = 0, AMT = 0)
      rbind(.e, .o)
    }))
  }

  .ctl <- function(...) {
    vaeControl(print = 0L, calcTables = FALSE, returnVae = TRUE,
               covariateSelection = FALSE, seed = 3L,
               itersBurnIn = 12L, iters = 34L, klWarmup = 6L, gammaIter = 24L, ...)
  }

  test_that("twoStage stage 2 runs and moves the log-density-only theta", {
    skip_on_cran()
    .dat <- .mkDat()
    .r <- suppressWarnings(suppressMessages(
      nlmixr2(.llMod(), .dat, est = "vae",
              control = .ctl(nonMuTheta = "regress", residOptimize = "twoStage"))))
    ## the mechanism assertion: stage 2 was entered.  Before the eligibility
    ## generalization an ll() model had NO stage-2 parameter, so this was 0 and
    ## "twoStage" was indistinguishable from the joint "optimize" solve.
    expect_true(.r$nStage2 > 0L)
    ## and it actually moved lsd off its ini() value, toward the simulated 1.2
    expect_false(isTRUE(all.equal(unname(.r$regressTheta[["lsd"]]), log(1.2),
                                  tolerance = 1e-8)))
    expect_lt(abs(.r$regressTheta[["lsd"]] - log(1.2)), 0.6)
    ## the structural thetas stayed in stage 1 and are still sane
    expect_true(all(is.finite(.r$regressTheta[c("lka", "lcl")])))
  })

  test_that("a Gaussian twoStage fit is unchanged by the eligibility rule", {
    skip_on_cran()
    ## regression guard for the reclassification: an err parameter is stage 2 and
    ## a structural theta that reaches d/dt is not -- exactly as before.
    .gMod <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- c(2, 3.45, 5); add.sd <- c(0, 0.7, 5)
        eta.ka ~ 0.6; eta.cl ~ 0.3 })
      model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv)
        d / dt(depot) <- -ka * depot
        d / dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd) })
    }
    .ui <- rxode2::assertRxUi(.gMod())
    ## tv is structural (feeds d/dt through v) and add.sd is the error param
    expect_equal(.vaeRegressStage2(.ui, c("tv", "add.sd"), c(-1L, 0L)),
                 c(0L, 1L))
    .r <- suppressWarnings(suppressMessages(
      nlmixr2(.gMod(), nlmixr2data::theo_sd, est = "vae",
              control = .ctl(nonMuTheta = "regress", residOptimize = "twoStage"))))
    expect_true(.r$nStage2 > 0L)
    expect_true(is.finite(.r$a[["add.sd"]]))
  })
})
