nmTest({
  # Multi-endpoint f-SAEM.  Nothing in the numerics needed changing: likInner0
  # reads each observation's own transform and variance, so the FOCEi inner
  # information H = Gamma_i^-1 is already the sum over that subject's
  # observations whatever endpoints they belong to, and the SAEM chain has always
  # carried per-endpoint ares/bres/yj/lambda indexed by ix_endpnt.  What was
  # single-endpoint was the .fsaemSupported gate and the covariate path's
  # reconstruction of the inner THETA.
  .mkInnerCtl <- function(cov = "jacobian") {
    list(rxControl = rxode2::rxControl(), fastCov = cov, fastLik = "focei",
         fastInnerIt = 100L, sumProd = FALSE, optExpression = TRUE, literalFix = FALSE,
         addProp = "combined2", eventSens = "jump", indTolRelax = TRUE,
         maxOdeRecalc = 5L, odeRecalcFactor = 10^0.5)
  }
  .supports <- function(fn) {
    .ui <- rxode2::rxUiDecompress(rxode2::rxode2(fn))
    .nlmixrSetMuRefTimeVarying(.ui, character(0))
    on.exit(.nlmixrRmMuRefTimeVarying(.ui), add = TRUE)
    .fsaemSupported(.ui)
  }

  test_that(".fsaemSupported accepts several normal endpoints", {
    expect_true(.supports(pkpd.two.endpoint))
    # combined error on one endpoint and additive on the other
    .comb <- function() {
      ini({
        tcl <- -3.2; tv <- -1; tec50 <- 0.5
        eta.cl ~ 0.09; eta.v ~ 0.09
        pk.prop <- 0.1; pk.sd <- 0.4; pd.sd <- 4
      })
      model({
        ka <- exp(0.5); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        e0 <- 100; ec50 <- exp(tec50)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        eff <- e0 * (1 - cp / (ec50 + cp))
        cp ~ prop(pk.prop) + add(pk.sd) | cp
        eff ~ add(pd.sd) | eff
      })
    }
    expect_true(.supports(.comb))
  })

  test_that(".fsaemSupported accepts a normal endpoint beside an ll() one", {
    # SAEM carries one scalar `distribution`, so a mixed model stays on the normal
    # path and cfg$distEp overrides only the ll() observations' loss.
    .mixed <- function() {
      ini({
        tcl <- -3.2; tv <- -1; tlam <- 3
        eta.cl ~ 0.09; eta.v ~ 0.09
        pk.sd <- 0.4
      })
      model({
        ka <- exp(0.5); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        lam <- exp(tlam)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(pk.sd) | cp
        ll(eff) ~ -log(lam) - time / lam
      })
    }
    expect_true(.supports(.mixed))
    .ui <- rxode2::rxUiDecompress(rxode2::rxode2(.mixed))
    # the per-endpoint distribution the kernel dispatches on, and its per-endpoint
    # residual count (an ll() endpoint has no residual theta)
    expect_equal(rxUiGet.saemDistEp(list(.ui)), c(1L, 4L))
    expect_equal(unname(rxUiGet.saemModNumEst(list(.ui))), c(1L, 0L))
    # .saemGeneralLik is all-ll(); .saemAnyGeneralLik is what forces the fast
    # kernel to run throughout and skips the transform-normal assertion
    expect_false(.saemGeneralLik(.ui))
    expect_true(.saemAnyGeneralLik(.ui))
  })

  test_that(".fsaemSupported still rejects what the kernel cannot target", {
    # a transform-normal residual is outside the add/prop envelope on either endpoint
    .tbs <- function() {
      ini({
        tcl <- -3.2; tv <- -1; tec50 <- 0.5
        eta.cl ~ 0.09; eta.v ~ 0.09
        pk.sd <- 0.4; pd.sd <- 4; pd.lam <- 1
      })
      model({
        ka <- exp(0.5); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        e0 <- 100; ec50 <- exp(tec50)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        eff <- e0 * (1 - cp / (ec50 + cp))
        cp ~ add(pk.sd) | cp
        eff ~ add(pd.sd) + boxCox(pd.lam) | eff
      })
    }
    expect_false(.supports(.tbs))
  })

  test_that("the covariate inner gets each endpoint's OWN residual", {
    # The covariate path rebuilds the inner THETA from SAEM's ares/bres and used
    # to take ares[1]/bres[1] for every residual, so endpoint 2 was preconditioned
    # with endpoint 1's error.  A fit cannot pin this: an independent
    # Metropolis-Hastings proposal need not match the target, so a wrong proposal
    # costs acceptance, not correctness (measured 0.95 vs 0.96 on the fit below).
    # The pin is the inner THETA itself.
    .d <- mkPkpdTwoEndpointData(n = 6L)
    set.seed(3)
    .wt <- stats::setNames(stats::runif(6, 50, 100), as.character(1:6))
    .d$WT <- unname(.wt[as.character(.d$id)])
    .covm <- function() {
      ini({
        tcl <- -3.2; tv <- -1; tec50 <- 0.5; cl.wt <- 0.5
        eta.cl ~ 0.09; eta.v ~ 0.09
        pk.sd <- 0.4; pd.sd <- 4
      })
      model({
        ka <- exp(0.5)
        cl <- exp(tcl + eta.cl + cl.wt * log(WT / 70))
        v <- exp(tv + eta.v)
        e0 <- 100; ec50 <- exp(tec50)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        eff <- e0 * (1 - cp / (ec50 + cp))
        cp ~ add(pk.sd) | cp
        eff ~ add(pd.sd) | eff
      })
    }
    .ui <- rxode2::rxUiDecompress(rxode2::rxode2(.covm))
    .nlmixrSetMuRefTimeVarying(.ui, character(0))
    on.exit(.nlmixrRmMuRefTimeVarying(.ui), add = TRUE)
    .N <- length(unique(.d$id))
    .mp0 <- matrix(c(-3.2, -1), .N, 2L, byrow = TRUE)
    .s <- .fsaemInnerSetupCov(.ui, .d, .mp0, .mkInnerCtl())
    on.exit(.fsaemInnerFree(), add = TRUE)
    expect_equal(.s$residEp, c(1L, 2L))

    .s2 <- .fsaemInnerUpdateCov(.s, .mp0, ares = c(0.4, 4), bres = c(0, 0),
                                plambda = numeric(0), omega = c(0.09, 0.09),
                                kiter = 1L)
    .th <- .s2$env$thetaIni
    .it <- .s$innerTheta
    expect_equal(unname(.th[.it$ntheta[match("pk.sd", .it$name)]]), 0.4)
    expect_equal(unname(.th[.it$ntheta[match("pd.sd", .it$name)]]), 4)
  })

  test_that("fsaemInnerMap: the Eq-17 proposal precision spans both endpoints", {
    skip_if_not_installed("numDeriv")
    # The multi-endpoint claim that the counters cannot make: H must equal
    # J' Sigma^-1 J + Omega^-1 with Sigma BLOCK-diagonal -- each observation
    # divided by its OWN endpoint's residual variance.  The two SDs here differ
    # by 10x, so using endpoint 1's residual for both is off by ~100x.
    .d <- mkPkpdTwoEndpointData(n = 8L)
    .ui <- rxode2::rxUiDecompress(rxode2::rxode2(pkpd.two.endpoint))
    .ctl <- .mkInnerCtl("jacobian")
    .N <- length(unique(.d$id))
    .neta <- 2L
    .fsaemInnerSetup(.ui, .d, matrix(0, .N, .neta), .ctl)
    on.exit(.fsaemInnerFree(), add = TRUE)
    .map <- .fsaemInnerMap(.ctl, .neta)
    expect_equal(nrow(.map$eta), .N)
    expect_true(all(.map$ok == 1L))

    .ids <- unique(.d$id)
    .sdOf <- function(i) {
      .obs <- .d[.d$id == .ids[i] & .d$evid == 0, ]
      ifelse(.obs$dvid == "cp", 0.4, 4)
    }
    .OmegaInv <- solve(diag(c(0.09, 0.09)))
    .predOf <- function(i, e) {
      .M <- .map$eta; .M[i, ] <- e
      as.numeric(vaeInnerLik(.M, 1L, FALSE, TRUE)$f[[i]])
    }
    .objOf <- function(i, e) {
      .M <- .map$eta; .M[i, ] <- e
      vaeInnerLik(.M, 1L, FALSE, FALSE)$obj[i]
    }
    .maxRel <- 0
    .maxStep <- 0
    for (i in seq_len(.N)) {
      .J <- numDeriv::jacobian(function(e) .predOf(i, e), .map$eta[i, ])
      .s <- .sdOf(i)
      expect_equal(nrow(.J), length(.s))
      .Hindep <- t(.J) %*% diag(1 / .s^2) %*% .J + .OmegaInv
      .maxRel <- max(.maxRel, max(abs(.map$hess[[i]] - .Hindep) / pmax(abs(.Hindep), 1)))
      # stationarity as a NEWTON step, not a raw gradient: this model's H entries
      # reach ~1e5, so a raw gradient of a couple of nats is already at the
      # optimizer's tolerance
      .g <- numDeriv::grad(function(e) .objOf(i, e), .map$eta[i, ])
      .maxStep <- max(.maxStep, max(abs(solve(.map$hess[[i]], .g))))
    }
    expect_lt(.maxRel, 1e-2)
    expect_lt(.maxStep, 1e-3)
    for (i in seq_len(.N)) {
      expect_true(all(is.finite(.map$gamma[[i]])))
      expect_true(all(eigen(.map$gamma[[i]], only.values = TRUE)$values > 0))
    }
  })
})
