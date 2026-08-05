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
    # eventSensInfo NULL.  The jump shape is a process global, so setting the
    # inner up leaves it pointing at the inner; the SAEM model has no jump
    # sensitivities and must not be solved under it.  Every other fsaem test uses
    # a bolus model, so this is the one that exercises that path.
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
    # and the SAEM solve was not corrupted by the inner's jump shape: same MLE
    expect_lt(max(abs(fixef(.fs) - fixef(.ss))), 0.1)
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
    nlmixr2est:::.nlmixrSetMuRefTimeVarying(.ui, character(0))
    on.exit(nlmixr2est:::.nlmixrRmMuRefTimeVarying(.ui), add = TRUE)
    expect_false(nlmixr2est:::.fsaemSupported(.ui))
  })
})
