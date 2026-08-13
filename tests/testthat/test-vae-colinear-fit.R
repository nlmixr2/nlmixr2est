## End-to-end colinearity handling in the VAE covariate M-step: the hysteresis
## and near-tie passes only exist inside a real training run, so these need
## actual fits.  Weekly, not essential (see .slowBatches in tests/testthat.R);
## the pure helpers are covered by test-vae-colinear.R.

nmTest({

  ## CL genuinely depends on WT; LBM is a ~0.99 correlated decoy that carries no
  ## independent signal.  This is the situation the feature exists for: the two
  ## score within noise of one another, so the winner is near-arbitrary.
  .colinData <- function(nid = 40L, seed = 9L) {
    .testSeed(seed)
    wt <- round(stats::runif(nid, 50, 100), 1)
    lbm <- round(0.75 * wt + stats::rnorm(nid, 0, 1.2), 1)
    ctr <- stats::median(wt)
    cl <- 2.7 * exp(0.75 * log(wt / ctr) + stats::rnorm(nid, 0, 0.15))
    v <- 31 * exp(stats::rnorm(nid, 0, 0.1))
    ka <- 1.5 * exp(stats::rnorm(nid, 0, 0.3))
    tms <- c(0.25, 0.5, 1, 2, 4, 6, 8, 12, 24)
    do.call(rbind, lapply(seq_len(nid), function(i) {
      ke <- cl[i] / v[i]
      f <- 320 / v[i] * ka[i] / (ka[i] - ke) * (exp(-ke * tms) - exp(-ka[i] * tms))
      rbind(data.frame(ID = i, TIME = 0, AMT = 320, EVID = 1, DV = 0,
                       WT = wt[i], LBM = lbm[i]),
            data.frame(ID = i, TIME = tms, AMT = 0, EVID = 0,
                       DV = f + stats::rnorm(length(tms), 0, 0.25),
                       WT = wt[i], LBM = lbm[i]))
    }))
  }

  .colinModel <- function() {
    ini({
      tka <- log(1.5); tcl <- log(2.7); tv <- log(31)
      eta.ka ~ 0.3; eta.cl ~ 0.2; eta.v ~ 0.1
      add.sd <- 0.3
    })
    model({
      ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }

  .ctl <- function(...) {
    vaeControl(iters = 60L, itersBurnIn = 15L, calcTables = FALSE, ...)
  }

  .fitQuiet <- function(d, ...) {
    suppressMessages(suppressWarnings(
      nlmixr2(.colinModel, d, est = "vae", control = .ctl(...))))
  }

  test_that("the colinear pass runs and reports near-tied alternates", {
    skip_on_cran()
    d <- .colinData()
    ## the design really is clustered, so a counter reaching zero would mean the
    ## mechanism never engaged rather than that there was nothing to do
    cv <- vaeCovariates(d, warn = FALSE)
    expect_true(nlmixr2est:::.vaeClusterBinds(cv$cluster, cv$group))

    ## a wide window so the report is exercised deterministically; the default is
    ## much tighter and correctly stays quiet on this fit
    f <- .fitQuiet(d, covSelectAltTol = 20)
    cd <- f$vae$colinear

    ## mechanism ran
    expect_gt(cd$nColinTest, 0L)
    ## and produced the diagnostic
    alt <- cd$alternates
    expect_s3_class(alt, "data.frame")
    expect_gt(nrow(alt), 0L)
    expect_identical(names(alt), c("param", "covariate", "alternate", "delta"))
    ## every reported alternate is a DIFFERENT covariate.  Alternate shapes of
    ## one covariate live in the same cluster and are near-tied on almost every
    ## model; admitting them would make this fire universally.
    expect_false(any(sub("_.*", "", alt$covariate) == sub("_.*", "", alt$alternate)))
    ## the true covariate is what got selected, and the decoy is the alternative
    expect_true(all(grepl("^WT", alt$covariate)))
    expect_true(all(grepl("^LBM", alt$alternate)))
    ## and the user is told, rather than having to go looking
    expect_match(f$runInfo, "near-tied colinear", all = FALSE)
  })

  test_that("the default window stays quiet on a well-determined fit", {
    skip_on_cran()
    f <- .fitQuiet(.colinData())
    ## a 0.99 decoy still fits measurably worse when the effect is strong, so at
    ## the default tolerance there is nothing to report and no note is raised
    expect_identical(nrow(f$vae$colinear$alternates), 0L)
    expect_false(any(grepl("near-tied colinear", f$runInfo)))
    ## the pass still ran -- silence here means "nothing was close", not "the
    ## mechanism was never reached"
    expect_gt(f$vae$colinear$nColinTest, 0L)
  })

  test_that("covSelectColinear=FALSE makes the mechanism unreachable", {
    skip_on_cran()
    d <- .colinData()
    fOff <- .fitQuiet(d, covSelectColinear = FALSE, covSelectAltTol = 20)
    ## not merely "reported nothing" -- never scored a single swap
    expect_identical(fOff$vae$colinear$nColinTest, 0L)
    expect_identical(fOff$vae$colinear$nColinHold, 0L)
    expect_identical(nrow(fOff$vae$colinear$alternates), 0L)
    expect_false(any(grepl("near-tied colinear", fOff$runInfo)))
  })

  test_that("hysteresis does not change a fit that has nothing to hold", {
    skip_on_cran()
    ## theo_sd offers one covariate group, so no cluster can bind and the whole
    ## colinear path is unreachable -- the selection must be bit-identical with
    ## the mechanism on and off
    skip_if_not_installed("nlmixr2data")
    d <- nlmixr2data::theo_sd
    mod <- .colinModel
    on <- suppressMessages(suppressWarnings(
      nlmixr2(mod, d, est = "vae", control = .ctl())))
    off <- suppressMessages(suppressWarnings(
      nlmixr2(mod, d, est = "vae", control = .ctl(covSelectColinear = FALSE))))
    expect_identical(on$vae$selected, off$vae$selected)
    expect_equal(on$vae$zPop, off$vae$zPop)
    expect_identical(on$vae$colinear$nColinTest, 0L)
  })
})
