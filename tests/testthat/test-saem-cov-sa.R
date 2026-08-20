nmTest({
  test_that("SAEM covMethod='sa' (stochastic-approximation Louis FIM) full covariance", {
    # The dedicated fixed-theta cov phase must (a) leave the estimate unperturbed,
    # (b) produce a PD full theta + Omega + residual covariance, and (c) agree with
    # the linearized FIM.

    one.cmt <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45; add.sd <- 0.7
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v  <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }

    ctlL <- saemControl(nBurn = 200, nEm = 300, print = 0, seed = 1L, covMethod = "linFim")
    ctlS <- saemControl(nBurn = 200, nEm = 300, print = 0, seed = 1L, covMethod = "sa",
                        nSaCov = 1000)

    fL <- .nlmixr(one.cmt, theo_sd, est = "saem", control = ctlL)
    fS <- .nlmixr(one.cmt, theo_sd, est = "saem", control = ctlS)

    # (a) the cov phase is frozen at theta_hat -- the reported estimate is unchanged
    expect_equal(unname(fS$theta), unname(fL$theta), tolerance = 1e-6)

    # (b) full PD covariance, method label preserved
    expect_equal(fS$covMethod, "sa")
    expect_true(is.matrix(fS$cov) && nrow(fS$cov) >= 4L)
    expect_true(all(is.finite(fS$cov)))
    expect_true(min(eigen(fS$cov, symmetric = TRUE, only.values = TRUE)$values) > 0)

    # Omega diagonal and residual rows are present and named on the eta
    expect_true(all(c("om.eta.ka", "om.eta.cl", "om.eta.v", "add.sd") %in% rownames(fS$cov)))

    # (c) SA and linFim SEs agree (different estimators -> allow Monte-Carlo tolerance)
    .cmn <- intersect(rownames(fS$cov), rownames(fL$cov))
    .seS <- sqrt(diag(fS$cov))[.cmn]
    .seL <- sqrt(diag(fL$cov))[.cmn]
    expect_equal(unname(.seS), unname(.seL), tolerance = 0.25)

    # residual SE surfaced in the parameter table
    expect_true(is.finite(fS$parFixedDf["add.sd", "SE"]))
    expect_gt(fS$parFixedDf["add.sd", "SE"], 0)

    # issue #816: the FORMATTED $parFixed must carry the same residual SE as
    # $parFixedDf / sqrt(diag($cov)) -- it previously printed uninitialized
    # memory (a denormal like 9.4e-323).  A theta with no covariance row must
    # be blank, never garbage.
    for (.f in list(fS, fL)) {
      if ("add.sd" %in% rownames(.f$cov)) {
        expect_equal(unname(.f$parFixedDf["add.sd", "SE"]),
                     unname(sqrt(diag(.f$cov))["add.sd"]))
        .seNum <- suppressWarnings(as.numeric(.f$parFixed["add.sd", "SE"]))
        expect_true(is.finite(.seNum))
        expect_gt(.seNum, 1e-100)
        expect_equal(.seNum, signif(unname(.f$parFixedDf["add.sd", "SE"]), 3),
                     tolerance = 1e-2)
        .rseNum <- suppressWarnings(as.numeric(.f$parFixed["add.sd", "%RSE"]))
        expect_true(is.finite(.rseNum))
        expect_gt(.rseNum, 1e-100)
      } else {
        expect_identical(unname(.f$parFixed["add.sd", "SE"]), "")
      }
    }
  })

  test_that("SAEM covMethod='fim' inverts the (mu-block-corrected) estimation-phase FIM", {
    # Regression: the shared per-iteration integrand omits the deterministic mu-block
    # complete Hessian, so before the correction Ha's theta block was indefinite and
    # covMethod='fim' produced NaN SEs.  With the correction fim is a valid PD full
    # covariance agreeing with the linearized FIM.

    one.cmt <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45; add.sd <- 0.7
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v  <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }

    ctlL <- saemControl(nBurn = 200, nEm = 300, print = 0, seed = 1L, covMethod = "linFim")
    ctlF <- saemControl(nBurn = 200, nEm = 300, print = 0, seed = 1L, covMethod = "fim")

    fL <- .nlmixr(one.cmt, theo_sd, est = "saem", control = ctlL)
    fF <- .nlmixr(one.cmt, theo_sd, est = "saem", control = ctlF)

    expect_equal(fF$covMethod, "fim")
    expect_true(is.matrix(fF$cov) && nrow(fF$cov) >= 4L)
    expect_true(all(is.finite(fF$cov)))
    expect_true(min(eigen(fF$cov, symmetric = TRUE, only.values = TRUE)$values) > 0)
    expect_true(all(c("om.eta.ka", "om.eta.cl", "om.eta.v", "add.sd") %in% rownames(fF$cov)))

    .cmn <- intersect(rownames(fF$cov), rownames(fL$cov))
    expect_equal(unname(sqrt(diag(fF$cov))[.cmn]),
                 unname(sqrt(diag(fL$cov))[.cmn]), tolerance = 0.25)

    expect_true(is.finite(fF$parFixedDf["add.sd", "SE"]))
    expect_gt(fF$parFixedDf["add.sd", "SE"], 0)
  })

  test_that("SAEM iteration history records off-diagonal Omega covariances", {
    # block-Omega model: the off-diagonal covariance trajectory (cov.<eta>.<eta>) is
    # recorded in parHistData, named consistently with the covariance matrix, and its
    # final value matches the fitted Omega off-diagonal.
    blk <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45; add.sd <- 0.3; prop.sd <- 0.1
        eta.ka ~ 0.6
        eta.cl + eta.v ~ c(0.3, 0.05, 0.1)
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd) + prop(prop.sd)
      })
    }
    f <- .nlmixr(blk, theo_sd, est = "saem",
                 control = saemControl(nBurn = 150, nEm = 200, print = 0, seed = 1L,
                                       covMethod = "linFim"))
    ph <- f$parHistData
    expect_true("cov.eta.v.eta.cl" %in% names(ph))
    expect_equal(ph[["cov.eta.v.eta.cl"]][nrow(ph)],
                 unname(f$omega["eta.cl", "eta.v"]), tolerance = 0.05)

    # a diagonal-Omega model must not gain any covariance columns (no regression)
    diagM <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; add.sd <- 0.7
            eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1 })
      model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd) })
    }
    fd <- .nlmixr(diagM, theo_sd, est = "saem",
                  control = saemControl(nBurn = 100, nEm = 120, print = 0, seed = 1L))
    expect_false(any(grepl("^cov\\.", names(fd$parHistData))))
  })

  test_that("fim/sa splice linFim's variance block for off-diagonal Omega / combined residuals", {
    # The analytic FIM cannot reliably do off-diagonal Omega or non-additive residuals;
    # those variance params are spliced from linFim's blocB.  On a block-Omega model the
    # sa covariance must include cov.<eta>.<eta>, and every variance-block SE must equal
    # linFim's varCov computed on the same fit (theta stays simulation-based).
    blk <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45; add.sd <- 0.7
        eta.ka ~ 0.6
        eta.cl + eta.v ~ c(0.3, 0.05, 0.1)
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    fS <- .nlmixr(blk, theo_sd, est = "saem",
                  control = saemControl(nBurn = 200, nEm = 300, print = 0, seed = 1L,
                                        covMethod = "sa", nSaCov = 1000))
    skip_if_not(identical(fS$covMethod, "sa"))   # near-singular fits legitimately fall back
    expect_true("cov.eta.v.eta.cl" %in% rownames(fS$cov))
    .saem <- fS$saem
    attr(.saem, "env") <- fS$env
    .vc <- attr(suppressWarnings(calc.COV(.saem)), "varCov")
    .vn <- colnames(.vc)
    expect_equal(unname(sqrt(diag(fS$cov[.vn, .vn]))),
                 unname(sqrt(diag(.vc))), tolerance = 1e-6)
  })

  test_that("a general log-likelihood endpoint is scored by its own density (#871)", {
    # calc.2LL and calc.COV used to build a Gaussian residual SD for every row.  An
    # ll() endpoint has res.mod == 0, so the residual M-step never runs and ares/bres
    # keep their placeholder 10/1 -- every log-density was scored as a normal mean
    # with SD 10 + |ll|.  Both are checked against an exponential time-to-event model
    # whose marginal likelihood is a one-dimensional quadrature, so the reference is
    # closed form rather than a second fit.
    .mkTte <- function(seed = 1L, n = 150L, meanT = 40) {
      .testSeed(seed)
      do.call(rbind, lapply(seq_len(n), function(i) {
        lami <- meanT * exp(rnorm(1, 0, sqrt(0.15)))
        data.frame(ID = i, TIME = lami * -log(runif(1)), DV = 1, EVID = 0, CMT = 1)
      }))
    }
    expTte <- function() {
      ini({ tlam <- log(25); eta.lam ~ 0.2 })
      model({
        lam <- exp(tlam + eta.lam)
        ll(dv) ~ -log(lam) - time / lam
      })
    }
    .d <- .mkTte(1L)
    # covMethod="sa" is the saemControl default; it must route to linFim here
    .f <- .nlmixr(expTte, .d, est = "saem",
                  control = saemControl(nBurn = 150, nEm = 80, nmc = 3, seed = 1,
                                        print = 0L, calcTables = FALSE,
                                        covMethod = "sa"))

    # the per-observation ll() mask needs res.mod to survive into the trimmed
    # post-fit saem.cfg -- without it neither function can tell the rows apart
    expect_false(is.null(attr(.f$saem, "saem.cfg")$res.mod))
    expect_equal(attr(.f$saem, "saem.cfg")$res.mod, 0L)
    # sa/fim are refused for a general likelihood (the Louis correction has no
    # residual error to anchor it); linFim is what gets reported
    expect_equal(.f$covMethod, "linFim")

    # nb_param (#893): the C++ kernel still accumulates Ha/HaSa every fit
    # regardless of the reported covMethod -- for a general-likelihood model
    # (distribution==4) it must carry NO residual slot at all (theta + eta
    # only), not a spurious always-zero one
    expect_equal(nrow(.f$saem$Ha), 1L + 1L)

    .tlam <- fixef(.f)[["tlam"]]
    .om <- .f$omega[1, 1]
    # exact marginal -2LL and its observed information, by quadrature over the eta
    .nll <- function(p) {
      if (p[2] <= 0) return(1e10)
      .s <- sqrt(p[2])
      -sum(vapply(.d$TIME, function(t)
        log(integrate(function(e) {
          .lam <- exp(p[1] + e)
          exp(-log(.lam) - t / .lam) * dnorm(e, 0, .s)
        }, -8 * .s, 8 * .s, rel.tol = 1e-10)$value), numeric(1)))
    }

    # objective: a fine quadrature must reproduce the exact marginal -2LL.  Before
    # the fix this read 1124 against an exact 1392.
    expect_equal(calc.2LL(.f$saem, nnodes.gq = 13, nsd.gq = 5, .f$phiM),
                 2 * .nll(c(.tlam, .om)), tolerance = 1e-3)

    # standard errors: the linearized FIM must approximate the exact observed
    # information.  Before the fix these were off by roughly the placeholder SD.
    .se <- sqrt(diag(solve(optimHess(c(.tlam, .om), .nll))))
    .got <- sqrt(diag(.f$cov))
    expect_equal(length(.got), 2L)
    expect_equal(unname(.got), unname(.se), tolerance = 0.25)
    # a hard guard: the pre-fix values were an order of magnitude out
    expect_true(all(.got / .se > 0.5 & .got / .se < 2))

    # "fim" is refused for the same reason as "sa"
    .ff <- .nlmixr(expTte, .d, est = "saem",
                   control = saemControl(nBurn = 40, nEm = 20, nmc = 3, seed = 1,
                                         print = 0L, calcTables = FALSE,
                                         covMethod = "fim"))
    expect_equal(.ff$covMethod, "linFim")
  })

  test_that("covMethod='fim' orders the Fisher information by [phi1][phi0], not model order (#906)", {
    # Two models where the SAME two thetas mu-reference (get an eta) and the SAME
    # one does not (phi0) -- only the phi0 theta's *position* in the model changes.
    # Before the fix: (1) the phi0 mu information decays to ~0, so its reported SE
    # collapsed to a nonsense ~1e-6 (defect 1); and (2) the kernel orders its Fisher
    # information [phi1 mu][phi0 mu], not saemParamsToEstimate/iniDf order, so
    # whichever theta actually sits at that phi0 row got its SE silently attributed
    # to a DIFFERENT parameter name whenever phi0 was not already last (defect 2).
    lastPhi0 <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; add.sd <- 0.7
            eta.ka ~ 0.6; eta.cl ~ 0.3 })
      model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv)
              linCmt() ~ add(add.sd) })
    }
    firstPhi0 <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; add.sd <- 0.7
            eta.cl ~ 0.3; eta.v ~ 0.1 })
      model({ ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
              linCmt() ~ add(add.sd) })
    }

    ctlFim <- saemControl(nBurn = 150, nEm = 200, print = 0, seed = 1L, covMethod = "fim")

    for (.mod in list(list(fn = lastPhi0, phi0 = "tv"), list(fn = firstPhi0, phi0 = "tka"))) {
      fFim <- .nlmixr(.mod$fn, theo_sd, est = "saem", control = ctlFim)
      skip_if_not(identical(fFim$covMethod, "fim"))

      # reference: the linearized FIM on the SAME converged fit (exact-match target
      # for the phi0 splice; every structural theta is linearized directly here, so
      # this is unaffected by defects 1/2 and safe to compare row-by-row against)
      .saem <- fFim$saem
      attr(.saem, "env") <- fFim$env
      .cm <- suppressWarnings(calc.COV(.saem))
      .tn <- fFim$ui$saemParamsToEstimate[!fFim$ui$saemFixed]
      dimnames(.cm) <- list(.tn, .tn)
      seLin <- sqrt(diag(.cm))       # theta-only (calc.COV's return has no variance block)
      seFim <- sqrt(diag(fFim$cov))
      .cmn <- c("tka", "tcl", "tv", "add.sd")
      expect_true(all(c("tka", "tcl", "tv") %in% names(seLin)))
      expect_true(all(.cmn %in% names(seFim)))

      # defect 1: no SE has collapsed toward 0 -- in particular the phi0 theta's,
      # which is spliced in from this exact linearized FIM and so must match it
      expect_true(all(seFim[.cmn] > 1e-3))
      expect_equal(unname(seFim[[.mod$phi0]]), unname(seLin[[.mod$phi0]]), tolerance = 1e-8)

      # defect 2: every structural theta's SE is closer to its OWN linearized
      # reference than to any other theta's -- catches a silent row permutation
      # that #906 could not (add.sd is excluded: it is never part of the reordered
      # phi block, so cross-checking it here is unrelated to ordering)
      .thn <- c("tka", "tcl", "tv")
      for (.nm in .thn) {
        .own <- abs(seFim[[.nm]] - seLin[[.nm]])
        .other <- min(abs(seFim[[.nm]] - seLin[setdiff(.thn, .nm)]))
        expect_true(.own < .other,
                    info = sprintf("%s (model with phi0=%s): own diff %.4g not < nearest other %.4g",
                                   .nm, .mod$phi0, .own, .other))
      }
    }
  })

  test_that("covMethod='sa' also splices a real phi0 SE and keeps names in order (#906)", {
    firstPhi0 <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; add.sd <- 0.7
            eta.cl ~ 0.3; eta.v ~ 0.1 })
      model({ ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
              linCmt() ~ add(add.sd) })
    }
    ctlSa <- saemControl(nBurn = 150, nEm = 200, print = 0, seed = 1L, covMethod = "sa",
                         nSaCov = 300)
    fSa <- .nlmixr(firstPhi0, theo_sd, est = "saem", control = ctlSa)
    skip_if_not(identical(fSa$covMethod, "sa"))

    .saem <- fSa$saem
    attr(.saem, "env") <- fSa$env
    .cm <- suppressWarnings(calc.COV(.saem))
    .tn <- fSa$ui$saemParamsToEstimate[!fSa$ui$saemFixed]
    dimnames(.cm) <- list(.tn, .tn)
    seLin <- sqrt(diag(.cm)); seSa <- sqrt(diag(fSa$cov))

    expect_true(all(seSa[c("tka", "tcl", "tv", "add.sd")] > 1e-3))
    expect_equal(unname(seSa[["tka"]]), unname(seLin[["tka"]]), tolerance = 1e-8)
    # tcl must not have picked up tv's (pre-fix, mislabeled) value or vice versa
    expect_true(abs(seSa[["tcl"]] - seLin[["tcl"]]) < abs(seSa[["tcl"]] - seLin[["tv"]]))
    expect_true(abs(seSa[["tv"]] - seLin[["tv"]]) < abs(seSa[["tv"]] - seLin[["tcl"]]))
  })

  test_that("covMethod='fim' drops phi0/fixed rows by raw FIM position, not fixed-filtered position (#906 review)", {
    # Antigravity review of #906 caught a follow-on bug in the fix itself: the
    # kernel's FIM keeps a row for a fix()ed theta (nb_param in src/saem.cpp
    # does not subtract fixed thetas), so computing the phi0 drop position
    # against a fixed-FILTERED name vector silently drops the wrong row
    # whenever a fixed phi1 theta sits ahead of the (correctly identified)
    # phi0 theta in kernel order.  Here tcl (phi1, mu-referenced) is fixed and
    # sits before tka (phi0) in kernel order ([phi1 mu][phi0 mu] = tcl,tv,tka).
    m <- function() {
      ini({ tka <- 0.45; tcl <- fix(1); tv <- 3.45; add.sd <- 0.7
            eta.cl ~ 0.3; eta.v ~ 0.1 })
      model({ ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
              linCmt() ~ add(add.sd) })
    }
    ctlFim <- saemControl(nBurn = 150, nEm = 200, print = 0, seed = 1L, covMethod = "fim")
    fFim <- .nlmixr(m, theo_sd, est = "saem", control = ctlFim)
    skip_if_not(identical(fFim$covMethod, "fim"))

    # tcl is fixed -- no row reported for it at all
    expect_false("tcl" %in% rownames(fFim$cov))
    # tka (phi0) and tv (phi1) must each carry a real, correctly-labeled SE,
    # not a near-zero one (the wrong-row-dropped failure mode kept the
    # degenerate phi0 row in the inverted block instead of removing it)
    .se <- sqrt(diag(fFim$cov))
    expect_true(all(c("tka", "tv") %in% names(.se)))
    expect_true(all(.se[c("tka", "tv")] > 1e-3))

    .saem <- fFim$saem
    attr(.saem, "env") <- fFim$env
    .cm <- suppressWarnings(calc.COV(.saem))
    .tn <- fFim$ui$saemParamsToEstimate[!fFim$ui$saemFixed]
    dimnames(.cm) <- list(.tn, .tn)
    seLin <- sqrt(diag(.cm))
    expect_equal(unname(.se[["tka"]]), unname(seLin[["tka"]]), tolerance = 1e-8)
    expect_true(abs(.se[["tv"]] - seLin[["tv"]]) < abs(.se[["tv"]] - seLin[["tka"]]))
  })

  test_that("covMethod='fim'/'sa' refuses (falls back to linFim) rather than mislabel a phi0 covariate model (#906 review)", {
    # A mu-ref covariate makes saemParamsToEstimate interleave covariate
    # coefficient names with the plain thetas, so it no longer lines up 1:1
    # with the kernel's phi block -- there is no known-good FIM row order for
    # a model with BOTH a covariate and a phi0 theta.  Silently reporting from
    # raw model order would reintroduce #906; the safe behavior is to refuse
    # and fall back, the same way a general likelihood is refused.
    m <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; add.sd <- 0.7
            wt.cl <- 0.1
            eta.cl ~ 0.3; eta.v ~ 0.1 })
      model({
        ka <- exp(tka)   # phi0
        cl <- exp(tcl + eta.cl + wt.cl * WT)
        v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    ctlFim <- saemControl(nBurn = 100, nEm = 100, print = 0, seed = 1L, covMethod = "fim")
    fFim <- .nlmixr(m, theo_sd, est = "saem", control = ctlFim)
    expect_equal(fFim$covMethod, "linFim")
  })

  test_that("multi-endpoint fim/sa: one residual FIM slot per endpoint (#893)", {
    # Before the fix, src/saem.cpp had exactly ONE log-sigma2 slot no matter how
    # many endpoints the model declared, so a multi-endpoint fit's residual score
    # always came from whichever endpoint the per-chain loop processed last,
    # divided by endpoint 0's sigma2 -- and because that slot couples to theta/
    # Omega through the full Fisher information matrix, R's .saemFimToCov() never
    # even surfaced a residual SE for a multi-endpoint model (nrow(.ri) > 1
    # always skipped the residual block entirely).  A pure-additive two-endpoint
    # model must now get a REAL, separate SE for each endpoint's residual.
    pkpd <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45; tslope <- 1
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
        add.sd <- 0.7; add2.sd <- 5
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v  <- exp(tv + eta.v)
        slope <- exp(tslope)
        cp <- linCmt()
        pca <- slope * cp
        cp ~ add(add.sd)
        pca ~ add(add2.sd)
      })
    }
    ctl <- saemControl(nBurn = 150, nEm = 200, print = 0, seed = 1L,
                       covMethod = "sa", nSaCov = 500)
    f <- .nlmixr(pkpd, warfarin, est = "saem", control = ctl)
    skip_if_not(identical(f$covMethod, "sa"))   # near-singular fits legitimately fall back

    # nb_param carries one slot per endpoint: theta(4) + eta(3) + residual(2)
    expect_equal(nrow(f$saem$Ha), 4L + 3L + 2L)

    # both endpoints' additive residual SEs are real (not dropped/NA)
    expect_true(all(c("add.sd", "add2.sd") %in% rownames(f$cov)))
    expect_true(is.finite(f$parFixedDf["add.sd", "SE"]) && f$parFixedDf["add.sd", "SE"] > 0)
    expect_true(is.finite(f$parFixedDf["add2.sd", "SE"]) && f$parFixedDf["add2.sd", "SE"] > 0)

    # a mixed add+prop / add model: the pure-additive endpoint still gets a real
    # analytic slot, the combined endpoint's slot is an exact-zero row dropped
    # before solve() (its SE instead comes from the linFim splice)
    pkpd2 <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45; tslope <- 1
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
        add.sd <- 0.7; prop.sd <- 0.1; add2.sd <- 5
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v  <- exp(tv + eta.v)
        slope <- exp(tslope)
        cp <- linCmt()
        pca <- slope * cp
        cp ~ add(add.sd) + prop(prop.sd)
        pca ~ add(add2.sd)
      })
    }
    f2 <- .nlmixr(pkpd2, warfarin, est = "saem",
                  control = saemControl(nBurn = 150, nEm = 200, print = 0, seed = 1L,
                                        covMethod = "sa", nSaCov = 500))
    .Ha <- f2$saem$HaSa
    expect_equal(nrow(.Ha), 4L + 3L + 2L)
    # the combined endpoint's slot (cp, first endpoint) is exactly zero;
    # the pure-additive endpoint's slot (pca, second) is not
    .zeroRow <- apply(.Ha, 1L, function(r) all(r == 0))
    expect_equal(unname(.zeroRow[8]), TRUE)
    expect_equal(unname(.zeroRow[9]), FALSE)
    if (identical(f2$covMethod, "sa")) {
      expect_true("add2.sd" %in% rownames(f2$cov))
      expect_false("add.sd" %in% rownames(f2$cov))
      expect_false("prop.sd" %in% rownames(f2$cov))
    }
    # the pure-additive endpoint's SE is real regardless (its slot was never
    # dropped); the combined endpoint's SE depends on the linFim splice
    # succeeding, which is a convergence question unrelated to this fix and
    # already covered by the "fim/sa splice" test above
    expect_true(is.finite(f2$parFixedDf["add2.sd", "SE"]) && f2$parFixedDf["add2.sd", "SE"] > 0)
  })

  test_that(".saemLlObsMask refuses to guess rather than mis-score (#871)", {
    .ix <- c(1L, 1L, 2L, 2L)
    # res.mod present: per-observation, res.mod == 0 marks the ll() endpoint
    expect_equal(.saemLlObsMask(list(res.mod = c(0L, 1L), opt = list(distribution = 4)), .ix),
                 c(TRUE, TRUE, FALSE, FALSE))
    expect_equal(.saemLlObsMask(list(res.mod = c(1L, 1L), opt = list(distribution = 1)), .ix),
                 rep(FALSE, 4))
    # res.mod missing (a fit stored before it was kept): a single endpoint can be
    # attributed from the scalar distribution, more than one cannot
    expect_equal(.saemLlObsMask(list(opt = list(distribution = 4)), c(1L, 1L)), c(TRUE, TRUE))
    expect_equal(.saemLlObsMask(list(opt = list(distribution = 1)), .ix), rep(FALSE, 4))
    expect_error(.saemLlObsMask(list(opt = list(distribution = 4)), .ix), "res.mod")
    # a general-likelihood cfg whose res.mod marks no ll() row is inconsistent;
    # scoring those rows as normal is the defect, so fail instead
    expect_error(.saemLlObsMask(list(res.mod = c(1L, 1L), opt = list(distribution = 4)), .ix),
                 "res.mod")
    # no distribution at all (an old cfg) must not error
    expect_equal(.saemLlObsMask(list(opt = list()), .ix), rep(FALSE, 4))
  })

  test_that(".saemResEndpointIdx matches a residual parameter to its ix_endpnt endpoint (#904)", {
    # predDf$cond is in ix_endpnt (endpoint) order; a residual parameter's condition
    # is matched positionally against it, not assumed to line up with .ri's row order
    expect_equal(.saemResEndpointIdx(c("cp", "effect", "cp"), c("cp", "effect")), c(1L, 2L, 1L))
    # a condition naming no endpoint fails rather than being silently attributed
    # to the wrong (or first) one
    expect_equal(.saemResEndpointIdx(c("cp", "bogus"), c("cp", "effect")), c(1L, NA_integer_))
  })
})
