## Phase 5 (Milestone B): BICc-ELBO covariate selection on theophylline. The VAE
## must auto-discover WT and select it on ka and V (not ke), matching the paper,
## and the WT effect on ka must pull omega_ka down from its no-covariate value.
## (Slow: runs a moderate training schedule; skipped on CRAN.)

nmTest({
  test_that("vae selects WT on ka and V (not ke) for theophylline", {
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
    ui <- rxode2::assertRxUi(theo)
    ## Pin the 7.0.1 search explicitly (one shape family, mean centering) so the
    ## paper-matching numbers below keep their original meaning; the multi-shape
    ## default is covered separately.
    ctl <- vaeControl(itersBurnIn = 80L, klWarmup = 40L, gammaIter = 120L,
                      iters = 160L, hiddenDim = 25L, seed = 1L, covariateSelection = TRUE,
                      print = 0L, shapes = "power", covCenterType = "mean")
    prep <- .vaeDataPrep(ui, nlmixr2data::theo_sd, ctl)
    expect_equal(prep$covNames, "WT_power")
    expect_equal(prep$covRaw, "WT")
    expect_equal(prep$covType, "continuous")
    ## the design is the historic log(WT/mean(WT)) column
    .wt <- vapply(unique(nlmixr2data::theo_sd$ID),
                  function(i) nlmixr2data::theo_sd$WT[nlmixr2data::theo_sd$ID == i][1],
                  numeric(1))
    expect_equal(prep$covPop, mean(.wt))
    expect_equal(unname(prep$covMat[, 1]), unname(log(.wt / mean(.wt))))

    innerEnv <- .vaeInnerSetup(ui, nlmixr2data::theo_sd,
                               matrix(0, prep$N, prep$zDim), ctl)
    on.exit(.vaeInnerFree(), add = TRUE)
    fit <- rxode2::rxWithSeed(1L, .vaeTrain(prep, innerEnv, ctl))

    ## rows = params (ka, ke, V), single column WT
    expect_true(fit$selected[1, 1])   # WT -> ka
    expect_false(fit$selected[2, 1])  # WT -> ke NOT selected
    expect_true(fit$selected[3, 1])   # WT -> V
    ## WT->ka effect is large and positive (paper ~2.55)
    expect_gt(fit$beta[1, 1], 1.5)
    ## omega_ka pulled down toward the with-covariate value (< no-covariate ~0.61)
    expect_lt(sqrt(fit$omega[1]), 0.58)
    ## fixed effects sane
    expect_lt(abs(exp(fit$zPop[2]) - 0.0867) / 0.0867, 0.05)   # ke
    expect_lt(abs(exp(fit$zPop[3]) - 31.97) / 31.97, 0.10)     # V
  })

  ## End-to-end regression for the two covariate-output bugs, exercised together
  ## via a fixed residual parameter (literalFix=TRUE default), which is what
  ## triggered the parFixedDf drop:
  ##   Bug 1 -- the selected covariate coefficients (beta.*) must appear in the
  ##            population-parameter table, not just in $theta/$cov.
  ##   Bug 2 -- covariate-bearing mu-parameters must back-transform (exp) rather
  ##            than print the raw log-scale estimate.
  test_that("est=vae keeps covariate betas and back-transforms with a fixed parameter", {
    skip_on_cran()
    theoFix <- function() {
      ini({
        lka <- log(1.8); lke <- log(0.086); lV <- log(32)
        eta.ka ~ 0.3; eta.ke ~ 0.03; eta.V ~ 0.03
        add.err <- fix(0.7)
      })
      model({
        ka <- exp(lka + eta.ka); ke <- exp(lke + eta.ke); V <- exp(lV + eta.V)
        d/dt(depot) = -ka * depot
        d/dt(central) = ka * depot - ke * central
        cp <- central / V
        cp ~ add(add.err)
      })
    }
    ## the 7.0.1 search, so this stays a test of the OUTPUT bugs rather than of
    ## which shape the expanded search happens to pick
    ctl <- vaeControl(itersBurnIn = 80L, klWarmup = 40L, gammaIter = 120L,
                      iters = 160L, hiddenDim = 25L, seed = 1L,
                      covariateSelection = TRUE, print = 0L,
                      shapes = "power", covCenterType = "mean")
    fit <- suppressMessages(suppressWarnings(
      nlmixr2(theoFix, nlmixr2data::theo_sd, est = "vae", control = ctl)))

    pf <- fit$parFixedDf
    ## Bug 1: covariate coefficients present in the population-parameter table
    ## (WT selected on ka and V for theophylline).  Coefficients carry the shape
    ## they were written in.
    expect_true(all(c("beta.lka.WT.power", "beta.lV.WT.power") %in% rownames(pf)))
    expect_true(all(rownames(pf) %in% rownames(fit$cov) |
                      rownames(pf) == "add.err"))

    ## Bug 2: the covariate-bearing mu-parameters back-transform (exp), so the
    ## Back-transformed value differs from the raw log-scale Estimate
    expect_equal(pf["lka", "Back-transformed"], exp(pf["lka", "Estimate"]),
                 tolerance = 1e-6)
    expect_equal(pf["lV", "Back-transformed"], exp(pf["lV", "Estimate"]),
                 tolerance = 1e-6)
    ## the covariate-free lke also back-transforms
    expect_equal(pf["lke", "Back-transformed"], exp(pf["lke", "Estimate"]),
                 tolerance = 1e-6)
    ## the fixed residual parameter is on the natural scale (no exp)
    expect_equal(pf["add.err", "Estimate"], 0.7, tolerance = 1e-6)
    ## covariate coefficients are reported raw (not back-transformed)
    expect_equal(pf["beta.lka.WT.power", "Back-transformed"],
                 pf["beta.lka.WT.power", "Estimate"], tolerance = 1e-6)
  })

  ## The multi-shape default: BICc now arbitrates between the log and linear
  ## families rather than being forced into log(cov/center).  Assert on the
  ## COVARIATE and on exclusivity, not on which family happens to win.
  test_that("the prep list carries hockey blocks through to the encoder input", {
    skip_on_cran()
    theoHk <- function() {
      ini({ lka <- log(1.8); lke <- log(0.086); lV <- log(32)
        eta.ka ~ 0.3; eta.ke ~ 0.03; eta.V ~ 0.03; add.err <- 0.7 })
      model({ ka <- exp(lka + eta.ka); ke <- exp(lke + eta.ke); V <- exp(lV + eta.V)
        d/dt(depot) = -ka * depot; d/dt(central) = ka * depot - ke * central
        cp <- central / V; cp ~ add(add.err) })
    }
    ui <- rxode2::assertRxUi(theoHk)
    ctl <- vaeControl(shapes = c("lin", "hockey"), print = 0L, covMethod = "")
    prep <- .vaeDataPrep(ui, nlmixr2data::theo_sd, ctl)
    ## one block id per column, and the arms share theirs
    expect_length(prep$covBlock, ncol(prep$covMat))
    expect_equal(prep$covShape, c("lin", "hockeyLow", "hockeyHi"))
    expect_equal(prep$covGroup, c(1L, 1L, 1L))
    expect_equal(prep$covBlock, c(1L, 2L, 2L))
    ## the encoder input keeps one BLOCK per group, so a selected hockey reaches
    ## the encoder whole rather than truncated at the knot
    expect_equal(ncol(prep$covIn), 1L)
    hkOnly <- .vaeDataPrep(ui, nlmixr2data::theo_sd,
                          vaeControl(shapes = "hockey", print = 0L,
                                     covMethod = ""))
    expect_equal(hkOnly$covShape, c("hockeyLow", "hockeyHi"))
    expect_equal(hkOnly$covBlock, c(1L, 1L))
    expect_equal(ncol(hkOnly$covIn), 2L)
    ## and with no hockey requested every column is its own block, as before
    plain <- .vaeDataPrep(ui, nlmixr2data::theo_sd,
                          vaeControl(shapes = c("power", "lin"), print = 0L,
                                     covMethod = ""))
    expect_false(any(grepl("hockey", plain$covShape)))
    expect_equal(plain$covBlock, seq_along(plain$covShape))
  })

  test_that("the default search picks one shape per covariate end-to-end", {
    skip_on_cran()
    theo <- function() {
      ini({ lka <- log(1.8); lke <- log(0.086); lV <- log(32)
        eta.ka ~ 0.3; eta.ke ~ 0.03; eta.V ~ 0.03; add.err <- 0.7 })
      model({ ka <- exp(lka + eta.ka); ke <- exp(lke + eta.ke); V <- exp(lV + eta.V)
        d/dt(depot) = -ka * depot; d/dt(central) = ka * depot - ke * central
        cp <- central / V; cp ~ add(add.err) })
    }
    ui <- rxode2::assertRxUi(theo)
    ctl <- vaeControl(itersBurnIn = 80L, klWarmup = 40L, gammaIter = 120L,
                      iters = 160L, hiddenDim = 25L, seed = 1L,
                      covariateSelection = TRUE, print = 0L, covMethod = "")
    prep <- .vaeDataPrep(ui, nlmixr2data::theo_sd, ctl)
    ## WT contributes every shape family, all sharing one exclusion group; the
    ## hockey arms are two columns of one block within it
    expect_equal(unique(prep$covRaw), "WT")
    expect_equal(prep$covShape, c("power", "lin", "hockeyLow", "hockeyHi"))
    expect_equal(length(unique(prep$covGroup)), 1L)
    expect_equal(length(unique(prep$covBlock)), 3L)

    fit <- suppressWarnings(rxode2::rxWithSeed(
      1L, nlmixr2(ui, nlmixr2data::theo_sd, est = "vae", control = ctl)))
    sel <- fit$vae$selected
    ## WT still lands on ka and V but not ke, whichever shape won
    expect_true(any(sel[1, ]));  expect_false(any(sel[2, ])); expect_true(any(sel[3, ]))
    ## exclusivity: never two shapes of one covariate on one parameter
    for (k in seq_len(nrow(sel))) expect_lte(sum(sel[k, ]), 1L)
    ## exactly one coefficient per selected parameter, named for the shape used
    bn <- grep("^beta\\.", fit$ui$iniDf$name, value = TRUE)
    expect_equal(length(bn), sum(sel))
    expect_true(all(grepl("\\.(power|lin|log|identity|center)$", bn)))
    ## and the emitted model still re-parses with the mu-ref exp() intact
    expect_equal(fit$ui$muRefCurEval$curEval[fit$ui$muRefCurEval$parameter == "lka"],
                 "exp")
  })

  ## The L0-penalty warmup ramp (covSelectAlpha) is a distinct step in the
  ## iteration table -- it must be labeled "CovSel ramp" on the ramp iterations,
  ## which are exactly the ones that evaluate and print the ELBO objective.
  test_that("est=vae shows the CovSel ramp step in the iteration print", {
    skip_on_cran()
    theo <- function() {
      ini({ lka <- log(1.8); lke <- log(0.086); lV <- log(32)
        eta.ka ~ 0.3; eta.ke ~ 0.03; eta.V ~ 0.03; add.err <- 0.7 })
      model({ ka <- exp(lka + eta.ka); ke <- exp(lke + eta.ke); V <- exp(lV + eta.V)
        d/dt(depot) = -ka * depot; d/dt(central) = ka * depot - ke * central
        cp <- central / V; cp ~ add(add.err) })
    }
    ctl <- vaeControl(itersBurnIn = 2L, iters = 6L, klWarmup = 4L, gammaIter = 5L,
                      nGradStep = 2L, covariateSelection = TRUE, covSelectAlpha = 2,
                      print = 1L)
    tc <- textConnection("msgs", "w", local = TRUE)
    sink(tc, type = "message")
    ## guarantee the message sink and connection are released even if the fit errors
    withr::defer({ sink(type = "message"); close(tc) })
    out <- capture.output(
      suppressWarnings(nlmixr2(theo, nlmixr2data::theo_sd, est = "vae", control = ctl)),
      type = "output")
    allout <- c(out, msgs)
    ## the ramp iterations (it < klWarmup) are labeled, and the legend documents it
    expect_true(any(grepl("CovSel ramp", allout)))
    expect_true(any(grepl("covariate-selection L0-penalty warmup", allout)))
  })

  test_that("covariateSelection=FALSE estimates a model-declared coefficient (end-to-end)", {
    skip_on_cran()
    ## nonMuTheta='none' would formerly hold cl.wt fixed; it must now be regressed,
    ## move off its init, and appear as a parameter-history column
    cov <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; cl.wt <- 0.1; add.err <- 0.7; eta.cl ~ 0.1 })
      model({ ka <- exp(tka); cl <- exp(tcl + cl.wt * log(WT / 70) + eta.cl); v <- exp(tv)
        d/dt(depot) <- -ka * depot; d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v; cp ~ add(add.err) })
    }
    ctl <- vaeControl(covariateSelection = FALSE, nonMuTheta = "none", itersBurnIn = 3L,
                      klWarmup = 5L, iters = 15L, seed = 1L, print = 0L, calcTables = FALSE)
    fit <- suppressWarnings(nlmixr2(cov, nlmixr2data::theo_sd, est = "vae", control = ctl))
    est <- fit$parFixedDf["cl.wt", "Estimate"]
    expect_true(is.finite(est))
    expect_false(isTRUE(all.equal(est, 0.1)))        # moved off the init
    expect_true("cl.wt" %in% colnames(fit$parHist)) # the regress M-step ran on it
  })

  ## pinCovariates=TRUE restricts the BICc search to the model-declared pairs: a
  ## strong declared covariate is kept on ITS parameter only; a noise declared
  ## covariate is dropped and written back as 0; no non-declared cell is selected.
  test_that("pinCovariates restricts selection to declared pairs and zeros dropped ones", {
    skip_on_cran()
    d <- nlmixr2data::theo_sd
    ids <- unique(d$ID)
    ## subject-constant, positive, pure-noise covariate (no PK signal)
    set.seed(42)
    nz <- stats::setNames(stats::runif(length(ids), 0.5, 1.5), as.character(ids))
    d$NZ <- nz[as.character(d$ID)]
    theo <- function() {
      ini({
        lka <- log(1.8); lke <- log(0.086); lV <- log(32)
        beta_lka_WT <- 0.1
        beta_lke_NZ <- 0.1
        eta.ka ~ 0.3; eta.ke ~ 0.03; eta.V ~ 0.03
        add.err <- 0.7
      })
      model({
        ka <- exp(lka + beta_lka_WT * log(WT / 70) + eta.ka)
        ke <- exp(lke + beta_lke_NZ * log(NZ / 1) + eta.ke)
        V <- exp(lV + eta.V)
        d/dt(depot) = -ka * depot
        d/dt(central) = ka * depot - ke * central
        cp <- central / V
        cp ~ add(add.err)
      })
    }
    ui <- rxode2::assertRxUi(theo)
    ctl <- vaeControl(itersBurnIn = 80L, klWarmup = 40L, gammaIter = 120L, iters = 160L,
                      seed = 1L, print = 0L, covMethod = "")
    fit <- suppressWarnings(rxode2::rxWithSeed(1L, nlmixr2(ui, d, est = "vae", control = ctl)))

    ## The centered covariates are mu2 references, so the search runs on the
    ## derived nlmixrMuDerCov# columns; assert on the user-facing OUTPUT model
    ## (robust to the internal column naming) rather than the selected matrix.
    thetaForEta <- .foceiEtaThetaMap(ui)$thetaForEta
    kV <- match("lV", thetaForEta)
    ## restriction: V declares no covariate, so nothing is ever selected on its row
    expect_true(all(!fit$vae$selected[kV, ]))
    idf <- fit$ui$iniDf
    ## the noise covariate (NZ on ke) is dropped to exactly 0; the strong one
    ## (WT on ka) is kept with a non-trivial coefficient
    expect_equal(idf$est[idf$name == "beta_lke_NZ"], 0)
    expect_gt(abs(idf$est[idf$name == "beta_lka_WT"]), 0.5)
    ## pinning never adds an undeclared covariate: only the two declared beta_
    ## coefficients exist in the model
    expect_setequal(grep("^beta_", idf$name, value = TRUE), c("beta_lka_WT", "beta_lke_NZ"))
    ## mu2 rewriting is undone: original centered forms restored, derived cols gone
    mfun <- deparse(fit$ui$fun)
    expect_true(any(grepl("log\\(WT/70\\)|log\\(WT / 70\\)", mfun)))
    expect_false(any(grepl("nlmixrMuDerCov", mfun)))
    expect_false(any(grepl("nlmixrMuDerCov", names(fit))))
    ## the pinned note reached $runInfo
    expect_true(any(grepl("pinned to model-specified", fit$runInfo)))
  })

  ## issue #801: a covariate reaching the coefficient line only through an
  ## intermediate model variable (wt70 <- WT/70) was mis-classified as a plain
  ## non-mu structural theta -- frozen under nonMuTheta="none" and, under "eta",
  ## an eta was injected into the mu-referenced expression, ERRORING the fit with
  ## "2+ single population parameters in a single mu-referenced expression".  It
  ## must now estimate the declared coefficient in every nonMuTheta mode.
  test_that("est=vae estimates a covariate behind an intermediate var (nonMuTheta=eta)", {
    skip_on_cran()
    indirect <- function() {
      ini({
        lka <- log(1.8); lke <- log(0.086); lV <- log(32)
        beta.ka <- 0.1
        eta.ka ~ 0.3; eta.ke ~ 0.03; eta.V ~ 0.03
        add.err <- 0.7
      })
      model({
        wt70 <- WT / 70
        ka <- exp(lka + beta.ka * log(wt70) + eta.ka)
        ke <- exp(lke + eta.ke)
        V <- exp(lV + eta.V)
        d/dt(depot) = -ka * depot
        d/dt(central) = ka * depot - ke * central
        cp <- central / V
        cp ~ add(add.err)
      })
    }
    ctl <- vaeControl(itersBurnIn = 80L, klWarmup = 40L, gammaIter = 120L, iters = 160L,
                      seed = 1L, print = 0L, covMethod = "", nonMuTheta = "eta")
    fit <- suppressWarnings(rxode2::rxWithSeed(1L,
      nlmixr2(indirect, nlmixr2data::theo_sd, est = "vae", control = ctl)))
    bk <- fit$ui$iniDf$est[fit$ui$iniDf$name == "beta.ka"]
    ## not frozen at the 0.1 ini value; pulled to the theophylline WT-on-ka effect
    expect_gt(abs(bk - 0.1), 0.5)
    expect_gt(bk, 1)
    ## the coefficient stays in the reported model (no injected eta collapsed it)
    expect_true("beta.ka" %in% fit$ui$iniDf$name)
  })
})
