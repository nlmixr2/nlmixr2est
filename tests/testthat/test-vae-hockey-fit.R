## End-to-end hockey-stick covariate selection.  The unit-level pieces (shape
## expressions, block feasibility, the branch-and-bound) live in
## test-vae-cov-shapes.R / test-vae-covariate-accessor.R / test-vae-cov-groups.R;
## this
## file checks the whole path on data with a KNOWN answer.  (Slow: trains a
## moderate schedule, so it runs in a weekly batch rather than on every push.)

nmTest({
  ## One-compartment data whose ka carries a piecewise-linear WT effect knotted
  ## at the median.  `slopeLow`/`slopeHi` equal makes it a plain straight line,
  ## which is how the no-false-positive case is built from the same generator.
  .hockeyData <- function(slopeLow, slopeHi, nid = 60L, seed = 11L) {
    set.seed(seed)
    wt <- round(stats::runif(nid, 40, 140), 1)
    ctr <- stats::median(wt)
    lka <- log(1.5) + slopeLow * (wt < ctr) * (wt - ctr) +
      slopeHi * (wt >= ctr) * (wt - ctr)
    ka <- exp(lka + stats::rnorm(nid, 0, 0.15))
    ke <- 0.09; v <- 32
    tms <- c(0.25, 0.5, 1, 2, 4, 6, 8, 12, 24)
    d <- do.call(rbind, lapply(seq_len(nid), function(i) {
      cp <- 320 / v * ka[i] / (ka[i] - ke) * (exp(-ke * tms) - exp(-ka[i] * tms))
      rbind(data.frame(ID = i, TIME = 0, DV = NA_real_, AMT = 320, EVID = 1,
                       WT = wt[i]),
            data.frame(ID = i, TIME = tms,
                       DV = cp + stats::rnorm(length(tms), 0, 0.25),
                       AMT = 0, EVID = 0, WT = wt[i]))
    }))
    list(data = d, knot = ctr)
  }

  .hockeyModel <- function() {
    ini({ lka <- log(1.5); lke <- log(0.09); lV <- log(32)
      eta.ka ~ 0.1; add.err <- 0.3 })
    model({ ka <- exp(lka + eta.ka); ke <- exp(lke); V <- exp(lV)
      d/dt(depot) = -ka * depot
      d/dt(central) = ka * depot - ke * central
      cp <- central / V; cp ~ add(add.err) })
  }

  .hockeyCtl <- function() {
    vaeControl(itersBurnIn = 60L, klWarmup = 30L, gammaIter = 100L, iters = 130L,
               seed = 1L, print = 0L, covMethod = "")
  }

  ## train once and hand back everything the assertions need
  .runHockey <- function(sim) {
    ui <- rxode2::assertRxUi(.hockeyModel)
    ctl <- .hockeyCtl()
    prep <- .vaeDataPrep(ui, sim$data, ctl)
    inner <- .vaeInnerSetup(ui, sim$data, matrix(0, prep$N, prep$zDim), ctl)
    on.exit(.vaeInnerFree(), add = TRUE)
    fit <- rxode2::rxWithSeed(1L, .vaeTrain(prep, inner, ctl))
    list(ui = ui, prep = prep, fit = fit)
  }

  test_that("a true kink is found, written and reproduced exactly", {
    skip_on_cran()
    sim <- .hockeyData(slopeLow = -0.010, slopeHi = 0.020)
    r <- suppressWarnings(.runHockey(sim))
    prep <- r$prep; fit <- r$fit

    ## the search picks the hockey BLOCK -- both arms, never one, and never the
    ## tie-equivalent `lin + one arm`
    sel <- prep$covNames[which(fit$selected[1, ])]
    expect_equal(sel, c("WT_hockeyLow", "WT_hockeyHi"))

    ui2 <- suppressMessages(suppressWarnings(.vaeUpdateModel(r$ui, fit)))
    ini2 <- ui2$iniDf
    ## the coefficients are named for the arms and recover the simulated slopes
    lo <- ini2$est[ini2$name == "beta.lka.WT.hockey.low"]
    hi <- ini2$est[ini2$name == "beta.lka.WT.hockey.hi"]
    expect_length(lo, 1L); expect_length(hi, 1L)
    expect_equal(lo, -0.010, tolerance = 0.3)
    expect_equal(hi, 0.020, tolerance = 0.3)

    ## the written text is the agreed form, knotted at the median, and the
    ## mu-referenced exp(theta + ... + eta) shape survives so the theta still
    ## back-transforms
    txt <- paste(deparse(ui2$lstExpr[[1]]), collapse = " ")
    expect_match(txt, "beta.lka.WT.hockey.low", fixed = TRUE)
    expect_match(txt, "beta.lka.WT.hockey.hi", fixed = TRUE)
    expect_match(txt, paste0("WT < ", signif(sim$knot, 12)), fixed = TRUE)
    expect_match(txt, paste0("WT >= ", signif(sim$knot, 12)), fixed = TRUE)
    expect_match(txt, "eta.ka", fixed = TRUE)
    expect_equal(ui2$muRefCurEval$curEval[ui2$muRefCurEval$parameter == "lka"],
                 "exp")

    ## ROUND TRIP -- the strongest check there is.  Evaluating the WRITTEN model
    ## on the data must reproduce the M-step's own population prediction for
    ## every subject; any centering, sign or intercept error shows up here.
    tv <- ini2$est[ini2$name == "lka"]
    pred <- rep(tv, prep$N)
    for (j in which(fit$selected[1, ])) {
      .nm <- paste0("beta.lka.", prep$covRaw[j], ".",
                    .vaeShapeCoefTag(prep$covShape[j]))
      pred <- pred + ini2$est[ini2$name == .nm] * prep$covMat[, j]
    }
    expect_equal(pred, unname(fit$zPopMat[, 1]), tolerance = 1e-10)
  })

  test_that("a straight-line covariate does not buy the second coefficient", {
    skip_on_cran()
    ## same generator, equal slopes: hockey spans this exactly, but it costs two
    ## coefficients against lin's one, so BICc must reject it
    sim <- .hockeyData(slopeLow = 0.012, slopeHi = 0.012)
    r <- suppressWarnings(.runHockey(sim))
    sel <- r$prep$covNames[which(r$fit$selected[1, ])]
    expect_false(any(grepl("hockey", sel)))
    expect_true(any(grepl("WT_", sel)))
  })
})
