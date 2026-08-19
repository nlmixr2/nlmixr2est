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
