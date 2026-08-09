nmTest({
  test_that("saemControl fast* options validate and round-trip", {
    .c0 <- saemControl()
    expect_false(.c0$fast)
    expect_equal(.c0$fastKernel, "firstN")
    expect_equal(.c0$fastCov, "auto")
    expect_equal(.c0$fastIter, 20L)
    expect_equal(.c0$fastLik, "focei")

    .c1 <- saemControl(fast=TRUE, fastKernel="throughout", fastCov="hessian",
                       fastIter=5, fastLik="foce")
    expect_true(.c1$fast)
    expect_equal(.c1$fastKernel, "throughout")
    expect_equal(.c1$fastCov, "hessian")
    expect_equal(.c1$fastIter, 5L)
    expect_equal(.c1$fastLik, "foce")

    # feeding a control back through saemControl() (as getValidNlmixrCtl does)
    # preserves the fast settings
    .c1b <- do.call(saemControl, .c1)
    expect_true(.c1b$fast)
    expect_equal(.c1b$fastKernel, "throughout")
    expect_equal(.c1b$fastCov, "hessian")
    expect_equal(.c1b$fastIter, 5L)
    expect_equal(.c1b$fastLik, "foce")

    expect_error(saemControl(fastIter=0), "fastIter")
  })

  test_that("nu may carry a 4th element: the f-SAEM IMH sweep count", {
    expect_equal(saemControl()$mcmc$nu, c(2, 2, 2))
    .c4 <- saemControl(fast = TRUE, nu = c(2, 2, 2, 3))
    expect_equal(.c4$mcmc$nu, c(2, 2, 2, 3))
    # survives the control round-trip getValidNlmixrCtl does
    expect_equal(do.call(saemControl, .c4)$mcmc$nu, c(2, 2, 2, 3))
    # still exactly 3 or 4
    expect_error(saemControl(nu = c(2, 2)))
    expect_error(saemControl(nu = c(2, 2, 2, 3, 4)))
    # nu[4] is the sweep count ONLY -- it does not switch the kernel on
    expect_false(saemControl(nu = c(2, 2, 2, 3))$fast)
  })

  test_that("nu[4] sets the number of IMH sweeps per iteration", {
    skip_if_not_installed("nlmixr2data")
    one.cmt <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    .base <- function(...) saemControl(nBurn = 40, nEm = 10, nmc = 3, seed = 42,
                                       print = 0L, calcTables = FALSE, ...)
    .f5 <- suppressMessages(nlmixr2(one.cmt, nlmixr2data::theo_sd, est = "fsaem",
                                    control = .base()))
    .f2 <- suppressMessages(nlmixr2(one.cmt, nlmixr2data::theo_sd, est = "fsaem",
                                    control = .base(nu = c(2, 2, 2, 2))))
    # same number of MAP+IMH steps, but 2 proposals per subject per chain per
    # step instead of the default 5
    expect_equal(.f5$fsaemDiag$nStep, .f2$fsaemDiag$nStep)
    expect_equal(.f5$fsaemDiag$nProp / .f2$fsaemDiag$nProp, 5 / 2)
    # fewer sweeps still converges to the same place
    expect_lt(max(abs(fixef(.f5) - fixef(.f2))), 0.1)
  })

  test_that("fastFallback validates and is inert when every Gamma is usable", {
    skip_if_not_installed("nlmixr2data")
    expect_equal(saemControl()$fastFallback, "skip")
    expect_equal(saemControl(fastFallback = "prior")$fastFallback, "prior")
    expect_equal(do.call(saemControl, saemControl(fastFallback = "prior"))$fastFallback,
                 "prior")
    expect_error(saemControl(fastFallback = "nope"))

    one.cmt <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    .b <- function(...) saemControl(nBurn = 40, nEm = 10, nmc = 3, seed = 42,
                                    print = 0L, calcTables = FALSE, ...)
    .sk <- suppressMessages(nlmixr2(one.cmt, nlmixr2data::theo_sd, est = "fsaem",
                                    control = .b()))
    .pr <- suppressMessages(nlmixr2(one.cmt, nlmixr2data::theo_sd, est = "fsaem",
                                    control = .b(fastFallback = "prior")))
    # NOTE: no subject on a healthy fit has an unusable Gamma -- the FOCEi inner
    # information carries the prior term (H = J' Sigma^-1 J + Omega^-1), so it is
    # positive definite by construction.  Reaching the rescue needs an inner
    # optimizer failure, which is not reachable from the control surface, so this
    # asserts the option is INERT rather than exercising the rescue itself.
    expect_equal(.sk$fsaemDiag$nBadGamma, 0)
    expect_equal(.pr$fsaemDiag$nBadGamma, 0)
    expect_equal(.pr$fsaemDiag$nPriorFallback, 0)
    # inert means bit-identical, not merely close
    expect_equal(fixef(.pr), fixef(.sk))
  })

  test_that("fastMode='chainMean' centres the proposal on the chain mean", {
    skip_if_not_installed("nlmixr2data")
    expect_equal(saemControl()$fastMode, "map")
    expect_equal(do.call(saemControl, saemControl(fastMode = "chainMean"))$fastMode,
                 "chainMean")
    expect_error(saemControl(fastMode = "nope"))

    one.cmt <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    .b <- function(...) saemControl(nBurn = 100, nEm = 50, nmc = 3, seed = 42,
                                    print = 0L, calcTables = FALSE, ...)
    .map <- suppressMessages(nlmixr2(one.cmt, nlmixr2data::theo_sd, est = "fsaem",
                                     control = .b()))
    .cm <- suppressMessages(nlmixr2(one.cmt, nlmixr2data::theo_sd, est = "fsaem",
                                    control = .b(fastMode = "chainMean")))
    # the centre really moved: same work, different (lower) acceptance, because
    # the chain mean is a worse centre than the mode on a well-behaved posterior
    expect_equal(.cm$fsaemDiag$nStep, .map$fsaemDiag$nStep)
    expect_lt(.cm$fsaemDiag$accRate, .map$fsaemDiag$accRate)
    # still an efficient sampler, and still the same MLE -- a centre only costs
    # acceptance, never correctness, since the ratio uses whatever was proposed from
    expect_gt(.cm$fsaemDiag$accRate, 0.3)
    expect_lt(max(abs(fixef(.cm) - fixef(.map))), 0.05)
  })

  test_that("fastHRefresh reuses the MAP + Gamma between refreshes", {
    skip_if_not_installed("nlmixr2data")
    expect_equal(saemControl()$fastHRefresh, 1L)
    expect_equal(do.call(saemControl, saemControl(fastHRefresh = 5))$fastHRefresh, 5L)
    expect_error(saemControl(fastHRefresh = 0))

    one.cmt <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    .b <- function(...) saemControl(nBurn = 100, nEm = 50, nmc = 3, seed = 42,
                                    print = 0L, calcTables = FALSE, ...)
    .f1 <- suppressMessages(nlmixr2(one.cmt, nlmixr2data::theo_sd, est = "fsaem",
                                    control = .b()))
    .f5 <- suppressMessages(nlmixr2(one.cmt, nlmixr2data::theo_sd, est = "fsaem",
                                    control = .b(fastHRefresh = 5)))
    # default never reuses; k=5 over the 20 fast iterations recomputes at
    # kiter 0/5/10/15 and reuses the other 16
    expect_equal(.f1$fsaemDiag$nMapReuse, 0)
    expect_equal(.f5$fsaemDiag$nStep, .f1$fsaemDiag$nStep)
    expect_equal(.f5$fsaemDiag$nMapReuse, 16)
    # a stale proposal costs acceptance, not correctness
    expect_lt(.f5$fsaemDiag$accRate, .f1$fsaemDiag$accRate)
    expect_gt(.f5$fsaemDiag$accRate, 0.3)
    expect_lt(max(abs(fixef(.f5) - fixef(.f1))), 0.05)
  })

  test_that("fsaemControl forces fast=TRUE and validates as an fsaem control", {
    .fc <- fsaemControl(fast=FALSE, nBurn=7)
    expect_s3_class(.fc, "fsaemControl")
    expect_true(.fc$fast)
    expect_equal(.fc$mcmc$niter[1], 7)

    .gv <- getValidNlmixrControl(NULL, "fsaem")
    expect_s3_class(.gv, "fsaemControl")
    expect_true(.gv$fast)
  })

  test_that("est='fsaem' runs the live IMH kernel and converges to the SAEM MLE", {
    one.cmt <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v  <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    .ctl <- saemControl(nBurn=200, nEm=100, nmc=3, seed=42, print=0L, calcTables=FALSE)
    .fs <- suppressMessages(nlmixr2(one.cmt, nlmixr2data::theo_sd, est="fsaem", control=.ctl))
    .ss <- suppressMessages(nlmixr2(one.cmt, nlmixr2data::theo_sd, est="saem", control=.ctl))
    # The kernel's own counters, snapshotted onto each fit env.  Estimates alone
    # cannot tell a working IMH kernel from a silent degrade to plain SAEM.
    .diag <- .fs$fsaemDiag
    # fast flag flows through to the stored control (mechanism is wired up)
    expect_true(.fs$saemControl$fast)
    expect_false(.ss$saemControl$fast)
    # one MAP+IMH step per fast iteration, and every subject got a proposal
    expect_equal(.diag$nStep, 20)                      # fastIter default
    expect_gt(.diag$nProp, 0)
    expect_equal(.diag$nMapFail, 0)
    expect_equal(.diag$nBadGamma, 0)
    # an independent sampler drawn from the Laplace approximation of the actual
    # posterior accepts most of what it proposes (the paper reports ~0.9); a low
    # rate would mean the proposal was built at the wrong parameterization
    expect_gt(.diag$accRate, 0.5)
    # The IMH move must be followed by a rescore of the predictions / observation
    # loss / per-subject llik, or the random walks below it score their proposals
    # against the PREVIOUS state.  One rescore per fast iteration, and it really
    # is correcting something (measured max ~22 nats on this model).  Do NOT
    # assert on the acceptance rate instead: do_mcmc writes U_y(ind)=Uc_y(ind) on
    # every acceptance, so a stale entry self-heals the first time a subject
    # accepts and the aggregate rate barely moves (measured 0.2369 vs 0.2344
    # without the rescore) -- the symptom is a shifted answer, not over-acceptance.
    expect_equal(.fs$saemDiag$nRescore, .diag$nStep)
    expect_gt(.fs$saemDiag$uYStaleMax, 0)
    expect_equal(.ss$saemDiag$nRescore, 0)             # plain saem never rescores
    # and the composite kernel behaves like plain saem's: adding the IMH must not
    # change how readily the random walks accept
    expect_lt(max(abs(.fs$saemDiag$accRate - .ss$saemDiag$accRate)), 0.1)
    # the fast kernel changes the simulation trajectory (it fired) -- so fsaem is
    # NOT bit-identical to saem, but converges to the same MLE
    expect_false(isTRUE(all.equal(unname(fixef(.fs)), unname(fixef(.ss)))))
    expect_lt(max(abs(fixef(.fs) - fixef(.ss))), 0.05)
  })

  test_that("est='fsaem' reaches the MLE from poor starts at a short burn-in", {
    # Deliberately poor initial estimates and a short burn-in.
    #
    # This asserts CONVERGENCE, not a speedup.  theo_sd is a well-identified
    # 1-cmt model: 20 burn-in iterations already put BOTH methods at the MLE
    # (measured RMSE ~0.01 for each, which is SAEM's own Monte-Carlo noise
    # floor), so a "fsaem beats saem" inequality here compares two draws from
    # the same noise -- it flips with the seed and says nothing about the
    # kernel.  The f-SAEM paper reports exactly this ("no regression") for a
    # well-identified model; the gain is on ill-conditioned / sparse ones, e.g.
    # the collinear 2-cmt N=60 case in issue #845, which is too slow to fit
    # twice in a test.  The kernel's own counters in the test above are what
    # prove it fired.
    poor <- function() {
      ini({
        tka <- 0.1; tcl <- 0.5; tv <- 3.0
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    .d <- nlmixr2data::theo_sd
    .ref <- fixef(suppressMessages(nlmixr2(poor, .d, est = "saem",
      control = saemControl(nBurn = 400, nEm = 150, nmc = 3, seed = 7, print = 0L, calcTables = FALSE))))
    .short <- function(est, nb) {
      .f <- suppressMessages(nlmixr2(poor, .d, est = est,
        control = saemControl(nBurn = nb, nEm = 20, nmc = 3, seed = 7, print = 0L,
                              calcTables = FALSE, fastIter = nb)))
      sqrt(sum((fixef(.f) - .ref)^2))
    }
    .fs <- .short("fsaem", 20L)
    .ss <- .short("saem", 20L)
    # both land at the MLE; the bound is well outside the seed-to-seed spread of
    # either method (measured max ~0.04 over six seeds) and well inside the
    # distance a broken kernel would leave (the poor start is ~0.8 away)
    expect_lt(.fs, 0.1)
    expect_lt(.ss, 0.1)
  })

  test_that("est='fsaem' supports proportional and combined error", {
    .ctl <- saemControl(nBurn = 150, nEm = 100, nmc = 3, seed = 8, print = 0L, calcTables = FALSE)
    prop <- function() {
      ini({ tka<-0.45; tcl<-1; tv<-3.45; eta.ka~0.6; eta.cl~0.3; eta.v~0.1; prop.sd<-0.15 })
      model({ ka<-exp(tka+eta.ka); cl<-exp(tcl+eta.cl); v<-exp(tv+eta.v); linCmt() ~ prop(prop.sd) })
    }
    comb <- function() {
      ini({ tka<-0.45; tcl<-1; tv<-3.45; eta.ka~0.6; eta.cl~0.3; eta.v~0.1; add.sd<-0.3; prop.sd<-0.1 })
      model({ ka<-exp(tka+eta.ka); cl<-exp(tcl+eta.cl); v<-exp(tv+eta.v); linCmt() ~ add(add.sd) + prop(prop.sd) })
    }
    for (.m in list(prop, comb)) {
      .fs <- suppressMessages(nlmixr2(.m, nlmixr2data::theo_sd, est = "fsaem", control = .ctl))
      .ss <- suppressMessages(nlmixr2(.m, nlmixr2data::theo_sd, est = "saem", control = .ctl))
      # fast kernel active (fsaem differs) but converges near the SAEM estimate
      expect_true(.fs$saemControl$fast)
      # Compare only the well-identified STRUCTURAL parameters: the combined
      # add+prop residual split is weakly identified (multi-modal, a large-additive
      # corner optimum), so its two residual params can flip between fsaem and saem
      # -- that is model ill-conditioning, not a fast-kernel disagreement.
      .th <- c("tka", "tcl", "tv")
      expect_lt(max(abs(fixef(.fs)[.th] - fixef(.ss)[.th])), 0.08)
    }
  })

  test_that("est='fsaem' handles a (non-time-varying) covariate model", {
    skip_if_not_installed("nlmixr2data")
    covm <- function() {
      ini({ tka<-0.45; tcl<-1; tv<-3.45; cl.wt<-0.5; eta.ka~0.6; eta.cl~0.3; eta.v~0.1; add.sd<-0.7 })
      model({ ka<-exp(tka+eta.ka); cl<-exp(tcl+eta.cl+cl.wt*log(WT/70)); v<-exp(tv+eta.v); linCmt()~add(add.sd) })
    }
    .ctl <- saemControl(nBurn=200, nEm=100, nmc=3, seed=9, print=0L, calcTables=FALSE)
    .fs <- suppressMessages(nlmixr2(covm, nlmixr2data::theo_sd, est="fsaem", control=.ctl))
    .ss <- suppressMessages(nlmixr2(covm, nlmixr2data::theo_sd, est="saem",  control=.ctl))
    # the covariate-aware fast kernel fires (fsaem is not bit-identical to saem)
    expect_false(isTRUE(all.equal(unname(fixef(.fs)), unname(fixef(.ss)))))
    # and converges to the SAEM estimate, including the covariate coefficient
    expect_true("cl.wt" %in% names(fixef(.fs)))
    expect_lt(max(abs(fixef(.fs) - fixef(.ss))), 0.07)
  })

  test_that("est='fsaem' handles a covariate model with proportional error", {
    skip_if_not_installed("nlmixr2data")
    prop <- function() {
      ini({ tka<-0.45; tcl<-1; tv<-3.45; cl.wt<-0.5; eta.ka~0.6; eta.cl~0.3; eta.v~0.1; prop.sd<-0.15 })
      model({ ka<-exp(tka+eta.ka); cl<-exp(tcl+eta.cl+cl.wt*log(WT/70)); v<-exp(tv+eta.v); linCmt()~prop(prop.sd) })
    }
    .ctl <- saemControl(nBurn=200, nEm=100, nmc=3, seed=9, print=0L, calcTables=FALSE)
    .fs <- suppressMessages(nlmixr2(prop, nlmixr2data::theo_sd, est="fsaem", control=.ctl))
    .ss <- suppressMessages(nlmixr2(prop, nlmixr2data::theo_sd, est="saem",  control=.ctl))
    expect_false(isTRUE(all.equal(unname(fixef(.fs)), unname(fixef(.ss)))))  # kernel fires
    expect_true("cl.wt" %in% names(fixef(.fs)))
    expect_lt(max(abs(fixef(.fs) - fixef(.ss))), 0.07)
  })

  test_that("est='fsaem' handles a TIME-VARYING mu-ref covariate end to end", {
    skip_if_not_installed("nlmixr2data")
    # A time-varying mu-ref covariate is the one covariate shape that is NOT
    # absorbed into the per-subject mprior data -- it stays a beta regressor in
    # the inner and is refreshed from the live Plambda every iteration
    # (.fsaemInnerUpdateCov).  Nothing else fits one end to end, so that refresh
    # could regress silently.
    tvm <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45; cl.tv <- 0.1
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl + cl.tv * TVC)
        v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    .d <- nlmixr2data::theo_sd
    .d$TVC <- as.numeric(scale(.d$TIME))            # varies WITHIN subject
    .ctl <- saemControl(nBurn = 150, nEm = 80, nmc = 3, seed = 5, print = 0L,
                        calcTables = FALSE)
    .fs <- suppressMessages(nlmixr2(tvm, .d, est = "fsaem", control = .ctl))
    .ss <- suppressMessages(nlmixr2(tvm, .d, est = "saem", control = .ctl))
    # the fast kernel really ran on this model rather than degrading
    expect_gt(.fs$fsaemDiag$nStep, 0)
    expect_gt(.fs$fsaemDiag$accRate, 0.3)
    expect_equal(.fs$fsaemDiag$nMapFail, 0)
    # and the time-varying coefficient agrees with plain SAEM's
    expect_true("cl.tv" %in% names(fixef(.fs)))
    expect_true(is.finite(fixef(.fs)[["cl.tv"]]))
    expect_lt(max(abs(fixef(.fs) - fixef(.ss))), 0.1)
  })

  test_that("fastLik selects the inner likelihood family", {
    skip_if_not_installed("nlmixr2data")
    # foceiControl(foce=) is IGNORED when interaction is on, so a fastLik that
    # did not also turn interaction off was a no-op on the hessian path --
    # "focei" and "foce" built an identical inner.  Check each value now reaches
    # a distinct foceiControl.
    .mk <- function(lik, cov = "hessian") {
      .fsaemInnerFoceiControl(
        list(rxControl = rxode2::rxControl(), fastCov = cov, fastLik = lik,
             fastInnerIt = 100L, sumProd = FALSE, optExpression = TRUE,
             literalFix = FALSE, addProp = "combined2", eventSens = "jump",
             indTolRelax = TRUE, maxOdeRecalc = 5L, odeRecalcFactor = 10^0.5))
    }
    expect_equal(.mk("focei")$interaction, 1L)
    expect_equal(.mk("foce")$interaction, 0L)          # was 1L -- the no-op
    expect_equal(.mk("focep")$interaction, 0L)
    expect_equal(.mk("foce")$foce, "nonmem")
    expect_equal(.mk("focep")$foce, "foce+")
    # jacobian is already a no-interaction form, so fastLik only picks the variant
    expect_equal(.mk("focei", "jacobian")$interaction, 0L)

    one.cmt <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    .b <- function(...) saemControl(nBurn = 100, nEm = 50, nmc = 3, seed = 42,
                                    print = 0L, calcTables = FALSE, ...)
    .ss <- suppressMessages(nlmixr2(one.cmt, nlmixr2data::theo_sd, est = "saem",
                                    control = .b()))
    for (.lik in c("foce", "focep")) {
      .f <- suppressMessages(nlmixr2(one.cmt, nlmixr2data::theo_sd, est = "fsaem",
                                     control = .b(fastLik = .lik)))
      expect_equal(.f$saemControl$fastLik, .lik)
      # a different likelihood family still gives a usable proposal and the MLE
      expect_gt(.f$fsaemDiag$accRate, 0.3)
      expect_equal(.f$fsaemDiag$nMapFail, 0)
      expect_lt(max(abs(fixef(.f) - fixef(.ss))), 0.05)
    }
  })

  test_that("est='fsaem' is correct when ini() and eta order disagree", {
    skip_if_not_installed("nlmixr2data")
    # The inner's structural theta is filled from phiPop POSITIONALLY -- phi
    # column j is structural theta j.  That holds because the phi vector follows
    # the ini (ntheta) order, NOT the eta order, and the two can differ:
    # ini(tka, tv, tcl) with etas declared (ka, cl, v) gives
    # ui$muRefDataFrame$theta = tka, tcl, tv.  Routing phiPop through
    # muRefDataFrame "to fix the mismatch" drops acceptance from 0.88 to 0.02,
    # which is what proves the positional read is the correct one.  This pins it.
    .perm <- function() {
      ini({
        tka <- 0.45; tv <- 3.45; tcl <- 1        # deliberately NOT eta order
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    .ui <- rxode2::rxUiDecompress(rxode2::rxode2(.perm))
    .th <- .ui$iniDf[!is.na(.ui$iniDf$ntheta), ]
    .th <- .th[order(.th$ntheta), ]
    # the premise of the test: the two orders really do disagree here
    expect_false(identical(.th$name[is.na(.th$err)], .ui$muRefDataFrame$theta))

    .ctl <- saemControl(nBurn = 100, nEm = 50, nmc = 3, seed = 42, print = 0L,
                        calcTables = FALSE)
    .fs <- suppressMessages(nlmixr2(.perm, nlmixr2data::theo_sd, est = "fsaem",
                                    control = .ctl))
    .ss <- suppressMessages(nlmixr2(.perm, nlmixr2data::theo_sd, est = "saem",
                                    control = .ctl))
    # a transposed inner parameterization shows up as a collapsed acceptance
    # rate long before it shows up in the estimates
    expect_gt(.fs$fsaemDiag$accRate, 0.5)
    expect_equal(.fs$fsaemDiag$nMapFail, 0)
    expect_lt(max(abs(fixef(.fs) - fixef(.ss))), 0.05)
  })

  test_that("est='fsaem' converges with the Hessian proposal (fastCov='hessian')", {
    # the Hessian proposal path (calcEtaHessian, interaction=1) is the general
    # covariance route; on a continuous model it must converge to the same MLE
    # as the Jacobian/linearization path.
    one.cmt <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v  <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    .base <- saemControl(nBurn = 200, nEm = 100, nmc = 3, seed = 42, print = 0L,
                         calcTables = FALSE)
    .ss <- suppressMessages(nlmixr2(one.cmt, nlmixr2data::theo_sd, est = "saem",
                                    control = .base))
    .fh <- suppressMessages(nlmixr2(one.cmt, nlmixr2data::theo_sd, est = "fsaem",
      control = saemControl(nBurn = 200, nEm = 100, nmc = 3, seed = 42, print = 0L,
                            calcTables = FALSE, fastCov = "hessian")))
    # the requested proposal covariance flows through to the stored control
    expect_equal(.fh$saemControl$fastCov, "hessian")
    # the Hessian-proposal fast kernel fires but converges to the SAEM MLE
    expect_false(isTRUE(all.equal(unname(fixef(.fh)), unname(fixef(.ss)))))
    expect_lt(max(abs(fixef(.fh) - fixef(.ss))), 0.05)
  })

  test_that("est='fsaem' fits a modeled-lag model (event-sensitivity path)", {
    skip_if_not_installed("nlmixr2data")
    # A model with modeled dosing (lag/f/rate) is the ONLY kind whose FOCEi inner
    # carries event ("jump") sensitivities -- a plain bolus linCmt() model has
    # eventSensInfo NULL.  Every other fsaem test uses a bolus model, so this is
    # the one that fits the event-sensitivity path end to end.
    #
    # It does NOT pin the odeSwapEsOff() guard, and must not be read as doing so:
    # removing that guard entirely leaves this fit bit-identical (measured), because
    # the fit's first load of the SAEM model rebinds rxode2's event globals anyway.
    # test-fsaem-es.R constructs the re-entry case the guard actually covers.
    lagm <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45; tlag <- -0.7
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        lag(depot) <- exp(tlag)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }
    .ctl <- saemControl(nBurn = 150, nEm = 80, nmc = 3, seed = 11, print = 0L,
                        calcTables = FALSE)
    .fs <- suppressMessages(nlmixr2(lagm, nlmixr2data::theo_sd, est = "fsaem",
                                    control = .ctl))
    .ss <- suppressMessages(nlmixr2(lagm, nlmixr2data::theo_sd, est = "saem",
                                    control = .ctl))
    .diag <- .fs$fsaemDiag
    # the kernel fired on a model that DOES have event sensitivities
    expect_gt(.diag$nStep, 0)
    expect_gt(.diag$accRate, 0.3)
    expect_equal(.diag$nMapFail, 0)
    # the post-IMH rescore runs here too -- elsewhere it is only asserted on a
    # bolus model, so this is the one fit where the event-sensitivity path and
    # the rescore are both live
    expect_equal(.fs$saemDiag$nRescore, .diag$nStep)
    expect_gt(.fs$saemDiag$uYStaleMax, 0)
    # nothing is left installed for the next fit in the session to inherit
    expect_equal(.odeSwapInfo()$esLive, 0L)
    expect_lt(max(abs(fixef(.fs) - fixef(.ss))), 0.1)
  })

  test_that("est='fsaem' handles a fix()ed structural theta", {
    skip_if_not_installed("nlmixr2data")
    # A fix()ed structural theta is NOT in phi1/phi0, so counting it among the
    # structural positions made their number exceed nphi1 + nphi0, tripped the
    # C++ length guard, and filled EVERY structural slot from the fallback -- and
    # the fixed theta itself was left at 0 rather than its own value.
    .fixm <- function() {
      ini({
        tka <- 0.45; tcl <- fix(1); tv <- 3.45
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    .ctl <- saemControl(nBurn = 100, nEm = 50, nmc = 3, seed = 42, print = 0L,
                        calcTables = FALSE)
    .fs <- suppressMessages(nlmixr2(.fixm, nlmixr2data::theo_sd, est = "fsaem",
                                    control = .ctl))
    .ss <- suppressMessages(nlmixr2(.fixm, nlmixr2data::theo_sd, est = "saem",
                                    control = .ctl))
    # the kernel is healthy -- a mis-filled inner theta shows up here first
    expect_gt(.fs$fsaemDiag$accRate, 0.5)
    expect_equal(.fs$fsaemDiag$nMapFail, 0)
    # the fixed value is honored, and the rest agrees with plain SAEM
    expect_equal(unname(fixef(.fs)[["tcl"]]), 1)
    expect_lt(max(abs(fixef(.fs) - fixef(.ss))), 0.05)
  })

  test_that("est='fsaem' rejects a correlated omega from the fast kernel", {
    # The inner rebuilds Omega from the DIAGONAL of Gamma2_phi1, so a declared
    # off-diagonal block would have the IMH scoring against a different prior
    # than the SAEM chain -- the composed kernel would not target the chain's
    # distribution.  .fsaemSupported must degrade instead.
    .corr <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45
        eta.ka + eta.cl ~ c(0.6, 0.01, 0.3)      # correlated block
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    .ui <- rxode2::rxUiDecompress(rxode2::rxode2(.corr))
    .nlmixrSetMuRefTimeVarying(.ui, character(0))
    on.exit(.nlmixrRmMuRefTimeVarying(.ui), add = TRUE)
    expect_false(.fsaemSupported(.ui))

    .diag <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    .u2 <- rxode2::rxUiDecompress(rxode2::rxode2(.diag))
    .nlmixrSetMuRefTimeVarying(.u2, character(0))
    on.exit(.nlmixrRmMuRefTimeVarying(.u2), add = TRUE)
    expect_true(.fsaemSupported(.u2))
  })

  test_that("est='fsaem' rejects unsupported models (mixture) from the fast kernel", {
    # mixture models are outside the fast-kernel envelope: .fsaemSupported must
    # return FALSE so the fit degrades to standard SAEM.
    mixm <- function() {
      ini({
        tka <- 0.45; tcl1 <- 1; tcl2 <- 1.4; tv <- 3.45; p1 <- 0.5
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1; add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- mix(exp(tcl1 + eta.cl), p1, exp(tcl2 + eta.cl))
        v  <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    .ui <- rxode2::rxUiDecompress(rxode2::rxode2(mixm))
    expect_true(length(.ui$mixProbs) > 0L)
    .nlmixrSetMuRefTimeVarying(.ui, character(0))
    on.exit(.nlmixrRmMuRefTimeVarying(.ui), add = TRUE)
    expect_false(.fsaemSupported(.ui))
  })
})
