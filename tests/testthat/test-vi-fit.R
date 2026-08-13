## est="emvi" end-to-end (mean-field, point-estimate): the raw C++ loop (ELBO
## trend, reproducibility) and the finalized nlmixr2FitData (objective, tables,
## covariance) assembled via nlmixr2CreateOutputFromUi at the ADVI estimates.

nmTest({
  one.cmt <- function() {
    ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; eta.ka ~ 0.6; add.sd <- 0.7 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv)
      d/dt(depot) <- -ka*depot; d/dt(center) <- ka*depot - cl/v*center
      cp <- center/v; cp ~ add(add.sd) })
  }

  test_that("est='emvi' raw loop: ELBO increases and is reproducible", {
    ctl <- emviControl(iters = 150L, seed = 7L, print = 0L, returnVi = TRUE)
    res <- suppressMessages(suppressWarnings(
      nlmixr2(one.cmt, nlmixr2data::theo_sd, est = "emvi", control = ctl)))
    expect_s3_class(res, "nlmixr2vi")
    e <- res$elbo; n <- length(e); d <- max(1L, n %/% 10L)
    expect_gt(mean(e[(n - d + 1L):n]), mean(e[1:d]))   # ELBO trend up
    expect_true(all(is.finite(res$theta)))
    expect_true(all(res$popOmega > 0))

    ## same seed -> identical
    res2 <- suppressMessages(suppressWarnings(
      nlmixr2(one.cmt, nlmixr2data::theo_sd, est = "emvi", control = ctl)))
    expect_identical(res$theta, res2$theta)
    expect_identical(res$elbo, res2$elbo)
    expect_identical(res$mu, res2$mu)
  })

  test_that("est='advi' assembles a full nlmixr2FitData (objf + tables + cov)", {
    fit <- suppressMessages(suppressWarnings(
      nlmixr2(one.cmt, nlmixr2data::theo_sd, est = "emvi",
              control = emviControl(iters = 120L, print = 0L))))
    expect_s3_class(fit, "nlmixr2FitData")
    expect_true(is.finite(fit$objf))
    expect_true(all(c("IPRED", "CWRES") %in% names(fit)))
    ## population parameter table present with finite estimates
    expect_true(is.data.frame(fit$parFixedDf))
    expect_true(all(is.finite(fit$parFixedDf$Estimate)))
    ## ADVI artifacts carried on the fit env
    expect_false(is.null(fit$env$viElbo))
    expect_false(is.null(fit$env$viState))
    ## the optimization walk is standard parHistData (captured even with print=0)
    .ph <- fit$parHist
    expect_true(is.data.frame(.ph))
    expect_true(all(c("iter", "tka", "add.sd", "o(eta.ka)") %in% names(.ph)))
    expect_gte(max(.ph$iter), 120L)
    ## the finalize reuses the loop's compiled models (no symengine rebuild)
    expect_false(is.null(fit$env$foceiModel))
  })

  test_that("emviControl(tol=) stops early on ELBO convergence", {
    ## The mechanism must be OBSERVABLE, not just "the fit still works": with a
    ## loose tolerance the loop must stop before `iters`, and with tol=0 it must
    ## run every iteration.  Before this was wired up, `tol` was documented but
    ## never reached the C++ loop, so both runs were identical -- a test that
    ## only checked the fit succeeded would not have caught that.
    .base <- function(tol) emviControl(iters = 400L, seed = 7L, print = 0L,
                                       returnVi = TRUE, tol = tol,
                                       evalElbo = 25L)
    rOff <- suppressMessages(suppressWarnings(
      nlmixr2(one.cmt, nlmixr2data::theo_sd, est = "emvi", control = .base(0))))
    expect_equal(rOff$itRun, 400L)
    expect_false(isTRUE(rOff$tolStopped))
    expect_equal(length(rOff$elbo), 400L)

    rOn <- suppressMessages(suppressWarnings(
      nlmixr2(one.cmt, nlmixr2data::theo_sd, est = "emvi", control = .base(1))))
    expect_true(isTRUE(rOn$tolStopped))      # rel change < 1 always -> first check
    expect_lt(rOn$itRun, 400L)
    ## the reported trace is truncated to what actually ran, not zero-padded
    expect_equal(length(rOn$elbo), rOn$itRun)
    expect_true(all(is.finite(rOn$elbo)))
  })

  test_that("prior tempering does not leak into the reported omega or eta search", {
    ## A run that ENDS INSIDE the warm-up never reaches the iteration that
    ## restores the tempering scale.  Because the reported population omega is
    ## built through the same helper that applies that scale, the fit would
    ## otherwise report an INFLATED omega -- silently, and only for short runs.
    r <- suppressMessages(suppressWarnings(
      nlmixr2(one.cmt, nlmixr2data::theo_sd, est = "emvi",
              control = emviControl(iters = 10L, klWarmup = 50L, temperInit = 10,
                                    seed = 7L, print = 0L, returnVi = TRUE,
                                    tol = 0))))
    expect_equal(unname(diag(r$popOmegaMat)), unname(r$popOmega), tolerance = 1e-12)
    expect_true(all(r$popOmega > 0))

    ## The adaptEta search must score candidates on the UNTEMPERED objective.
    ##
    ## Assert the per-candidate SCORES, not the selected etaScale: the winner is
    ## one of five discrete values, so equality of the winner can hold while the
    ## scoring is tempered (both reviewers showed the earlier version of this
    ## assertion passing against the bug it was meant to catch).  The scores are
    ## continuous and every one of them moves under tempering.
    ##
    ## klWarmup must also EXCEED nAdapt (= min(iters, 75)); candidates are scored
    ## on their last nAdapt/3 iterations, so a smaller klWarmup leaves that
    ## window past the warm-up and untempered even without the fix.
    .fit <- function(kl) suppressMessages(suppressWarnings(
      nlmixr2(one.cmt, nlmixr2data::theo_sd, est = "emvi",
              control = emviControl(iters = 60L, klWarmup = kl, temperInit = 100,
                                    seed = 7L, print = 0L, returnVi = TRUE,
                                    tol = 0))))
    .a <- .fit(0L); .b <- .fit(80L)
    expect_equal(length(.a$etaScores), 5L)          # the search actually ran
    expect_true(all(is.finite(.a$etaScores)))
    expect_equal(.a$etaScores, .b$etaScores)        # untempered either way
    expect_equal(.a$etaScale, .b$etaScale)
  })

  test_that("a step-size search at the etaCandidates edge is reported", {
    ## Read $runInfo, NOT a warning handler: warnings raised during estimation
    ## are COLLECTED into the fit's $runInfo rather than propagated, so
    ## withCallingHandlers(warning=) sees nothing and would pass whatever the
    ## code did.
    ##
    ## The message must appear when the winner is at an EDGE and be absent when
    ## no search runs -- a report that always appears is as useless as one that
    ## never does.
    .info <- function(cand) {
      .f <- suppressMessages(suppressWarnings(
        nlmixr2(one.cmt, nlmixr2data::theo_sd, est = "emvi",
                control = emviControl(iters = 60L, seed = 7L, print = 0L,
                                      tol = 0, etaCandidates = cand))))
      as.character(.f$runInfo)
    }
    ## every entry of this grid is tiny, so the search runs to its top edge
    expect_true(any(grepl("top of etaCandidates", .info(c(1e-4, 5e-4, 1e-3)))))
    ## a single candidate means no search ran, so there is no edge to report
    expect_false(any(grepl("etaCandidates", .info(0.05))))
  })

  test_that("a resumed correlated fit reproduces the single long run", {
    ## emviControl(resume=) promises a resumed run is identical to one fresh run
    ## of the combined length.  perNoCor broke that: nbCorrel was recomputed from
    ## the RESUMED call's iters and compared against a CONTINUING global index, so
    ## the resumed run re-applied a hold the first run had already passed and
    ## restarted the off-diagonal gain in the wrong place.
    ##
    ## Assert the OMEGA BLOCK, not just the thetas: the thetas can agree while the
    ## correlation schedule differs, which is exactly the failure being guarded.
    .cor <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45
        eta.ka + eta.cl ~ c(0.6,
                            0.05, 0.3)
        add.sd <- 0.7 })
      model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v; cp ~ add(add.sd) })
    }
    ## perNoCor is given ABSOLUTELY (90 iterations).  A FRACTION cannot satisfy
    ## this: 0.75 of a 120-iteration run releases at 90, while 0.75 of the first
    ## 60-iteration leg releases at 45, so the two schedules genuinely differ and
    ## no amount of state-passing can reconcile them.  That is the whole reason
    ## the absolute form exists -- a restartable fit needs its schedule points
    ## pinned to absolute iterations rather than to a share of whatever slice is
    ## being run.
    .ctl <- function(n, res = NULL)
      emviControl(iters = n, seed = 7L, print = 0L, tol = 0, resume = res,
                  perNoCor = 90, adaptEta = FALSE, etaCandidates = 0.05)

    one <- suppressMessages(suppressWarnings(
      nlmixr2(.cor, nlmixr2data::theo_sd, est = "emvi", control = .ctl(120L))))
    half <- suppressMessages(suppressWarnings(
      nlmixr2(.cor, nlmixr2data::theo_sd, est = "emvi", control = .ctl(60L))))
    cont <- suppressMessages(suppressWarnings(
      nlmixr2(.cor, nlmixr2data::theo_sd, est = "emvi",
              control = .ctl(60L, res = half$env$viState))))

    expect_equal(unname(cont$omega), unname(one$omega), tolerance = 1e-8)
    expect_equal(unname(cont$theta), unname(one$theta), tolerance = 1e-8)
    ## the release point travelled with the state rather than being recomputed
    expect_equal(half$env$viState$nbCorrel, cont$env$viState$nbCorrel)
    expect_equal(one$env$viState$nbCorrel, 90)   # absolute, not 0.75 * iters
  })

  test_that("est decides the algorithm, whatever the control says", {
    ## a control whose pointEstimate contradicts est must not produce a fit whose
    ## $est misdescribes the algorithm that ran.  This is the direct-dispatch
    ## path (nlmixr2Est.emvi), which does not go through getValidNlmixrCtl.
    .run <- function(est, pe) {
      ctl <- emviControl(iters = 40L, seed = 7L, print = 0L, returnVi = TRUE,
                       pointEstimate = pe)
      suppressMessages(suppressWarnings(
        nlmixr2(one.cmt, nlmixr2data::theo_sd, est = est, control = ctl)))
    }
    expect_true(.run("emvi", FALSE)$pointEstimate)
    expect_false(.run("fbvi", TRUE)$pointEstimate)
  })

  test_that("a fit resumed under the OTHER method is refused", {
    ## emvi state carries sTheta/sLpo, fbvi carries mPop/Lpop; crossing them
    ## reaches adviOptimize_ with the wrong half missing and used to die on an
    ## Rcpp NULL conversion instead of naming the problem
    fE <- suppressMessages(suppressWarnings(
      nlmixr2(one.cmt, nlmixr2data::theo_sd, est = "emvi",
              control = emviControl(iters = 40L, seed = 7L, print = 0L))))
    expect_error(suppressMessages(suppressWarnings(
      nlmixr2(one.cmt, nlmixr2data::theo_sd, est = "fbvi",
              control = emviControl(iters = 40L, seed = 7L, print = 0L, resume = fE)))),
      "cannot resume")
    ## the same-method resume still works
    expect_s3_class(suppressMessages(suppressWarnings(
      nlmixr2(one.cmt, nlmixr2data::theo_sd, est = "emvi",
              control = emviControl(iters = 40L, seed = 7L, print = 0L, resume = fE)))),
      "nlmixr2FitData")
  })

  test_that("a fit does not depend on what ran before it in the same session", {
    ## The advi run state lives in file statics, and forgetting to reset a NEW
    ## one has now shipped three times.  This pins the invariant those bugs each
    ## broke: the same fit must give the same answer whatever preceded it.
    ##
    ## It is a REGRESSION GUARD, not a demonstration of a current bug -- as of
    ## this commit adviOptimize_ happens to rewrite every static from its args,
    ## so the leak is only reachable through the direct-gradient entry points.
    ## The guard is the point: it fails the next time a static is added without
    ## a default, which is exactly how the previous three got in.
    .runB <- function() {
      suppressMessages(suppressWarnings(
        nlmixr2(one.cmt, nlmixr2data::theo_sd, est = "emvi",
                control = emviControl(iters = 60L, seed = 3L, print = 0L,
                                      returnVi = TRUE))))
    }
    b1 <- .runB()
    ## dirty every static a run can touch: tempering on, an early ELBO stop, a
    ## non-default correlation schedule, and the full-Bayes Jacobian tables
    suppressMessages(suppressWarnings(
      nlmixr2(one.cmt, nlmixr2data::theo_sd, est = "fbvi",
              control = fbviControl(iters = 60L, seed = 11L, print = 0L,
                                    returnVi = TRUE, klWarmup = 40L,
                                    temperInit = 10, tol = 1e-1, evalElbo = 10L,
                                    perNoCor = 5))))
    b2 <- .runB()
    expect_identical(b1$theta, b2$theta)
    expect_identical(b1$popOmega, b2$popOmega)
    expect_identical(b1$elbo, b2$elbo)
  })
})
