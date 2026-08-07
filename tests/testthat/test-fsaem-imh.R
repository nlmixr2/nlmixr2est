nmTest({
  # f-SAEM independent Metropolis-Hastings kernel (paper Alg 2 / Eq 23), run on
  # the FOCEi inner allocation with likInner0 as the joint-target oracle.
  .mkInnerCtl <- function(cov = "jacobian") {
    list(rxControl = rxode2::rxControl(), fastCov = cov, fastLik = "focei",
         fastInnerIt = 100L, sumProd = FALSE, optExpression = TRUE, literalFix = FALSE,
         addProp = "combined2", eventSens = "jump", indTolRelax = TRUE,
         maxOdeRecalc = 5L, odeRecalcFactor = 10^0.5)
  }

  test_that("linear-Gaussian model: Laplace proposal is exact -> accept w.p. 1", {
    skip_if_not_installed("nlmixr2data")
    # prediction linear in eta + additive constant error -> Gaussian posterior,
    # so the Laplace proposal equals the target and every candidate is accepted
    # (paper Remark 3).
    lin <- function() {
      ini({ tcp <- 5; eta.cp ~ 2; add.sd <- 1 })
      model({ cp <- tcp + eta.cp; cp ~ add(add.sd) })
    }
    .ui <- rxode2::rxUiDecompress(rxode2::rxode2(lin))
    .data <- nlmixr2data::theo_sd
    .N <- length(unique(.data$ID))
    .ctl <- .mkInnerCtl("jacobian")
    .env <- .fsaemInnerSetup(.ui, .data, matrix(0, .N, 1L), .ctl)
    on.exit(.fsaemInnerFree(), add = TRUE)
    .map <- .fsaemInnerMap(.ctl, 1L)

    set.seed(1)
    .nchain <- 5L; .nsweep <- 20L
    .etaCur <- matrix(rnorm(.nchain*.N, 0, 3), .nchain*.N, 1L) # start far from mode
    .imh <- .fsaemImh(.map, .etaCur, .nchain, .nsweep)
    # every proposal accepted, from any starting point
    expect_equal(.imh$nAcc, rep(.nchain*.nsweep, .N))
  })

  test_that("nonlinear model: IMH targets the true conditional posterior", {
    skip_if_not_installed("nlmixr2data")
    nl <- function() {
      ini({ tka<-0.45;tcl<-1;tv<-3.45; eta.ka~0.6;eta.cl~0.3;eta.v~0.1; add.sd<-0.7 })
      model({ ka<-exp(tka+eta.ka); cl<-exp(tcl+eta.cl); v<-exp(tv+eta.v); linCmt()~add(add.sd) })
    }
    .ui <- rxode2::rxUiDecompress(rxode2::rxode2(nl))
    .data <- nlmixr2data::theo_sd
    .N <- length(unique(.data$ID))
    .neta <- 3L
    .ctl <- .mkInnerCtl("jacobian")
    .env <- .fsaemInnerSetup(.ui, .data, matrix(0, .N, .neta), .ctl)
    on.exit(.fsaemInnerFree(), add = TRUE)
    .map <- .fsaemInnerMap(.ctl, .neta)

    set.seed(2)
    .nchain <- 30L; .nsweep <- 40L; .burn <- 15L
    .chol <- t(vapply(seq_len(.N), function(i) as.numeric(t(chol(.map$gamma[[i]]))), numeric(.neta*.neta)))
    .eta <- do.call(rbind, lapply(seq_len(.nchain), function(k) .map$eta)) # start at MAP
    .samp <- vector("list", .N); .acc <- integer(.N)
    for (s in seq_len(.nsweep)) {
      .r <- fsaemImhKernel_(.eta, .map$eta, .chol, .nchain,
                            as.integer(rxode2::getRxThreads()),
                            matrix(0, .N, .neta), numeric(0), numeric(0), integer(0),
                            as.double(s * .N * .nchain), 10L)
      .eta <- .r$eta; .acc <- .acc + .r$nAcc
      if (s > .burn) for (i in seq_len(.N)) {
        .samp[[i]] <- rbind(.samp[[i]], .eta[seq(i, by = .N, length.out = .nchain), , drop = FALSE])
      }
    }
    # efficient independent sampler: high acceptance (paper reports ~0.9)
    expect_gt(mean(.acc/(.nchain*.nsweep)), 0.6)
    # posterior mean ~ MAP and posterior cov ~ Gamma for a couple of subjects
    for (i in c(1L, 6L)) {
      .cm <- colMeans(.samp[[i]])
      expect_lt(max(abs(.cm - .map$eta[i, ])), 0.15)
      .cc <- diag(cov(.samp[[i]]))
      .gg <- diag(.map$gamma[[i]])
      expect_lt(max(abs(.cc - .gg)/pmax(.gg, 1e-3)), 0.6)
    }
  })

  # The two degenerate branches of the kernel.  Neither is reached by a normal
  # fit -- every fsaem fit in the suite reports nBadGamma = 0 and none exercises a
  # bounded proposal at all -- so they are driven directly through the exported
  # kernel and pinned on the diagnostic counters, which is the only evidence the
  # branch actually ran.
  test_that("IMH kernel handles an unusable Gamma and an unsatisfiable bound", {
    skip_if_not_installed("nlmixr2data")
    lin <- function() {
      ini({ tcp <- 5; eta.cp ~ 2; add.sd <- 1 })
      model({ cp <- tcp + eta.cp; cp ~ add(add.sd) })
    }
    .ui <- rxode2::rxUiDecompress(rxode2::rxode2(lin))
    .data <- nlmixr2data::theo_sd
    .N <- length(unique(.data$ID))
    .ctl <- .mkInnerCtl("jacobian")
    .env <- .fsaemInnerSetup(.ui, .data, matrix(0, .N, 1L), .ctl)
    on.exit(.fsaemInnerFree(), add = TRUE)
    .map <- .fsaemInnerMap(.ctl, 1L)
    .nchain <- 2L
    .eta <- matrix(0, .nchain * .N, 1L)
    .cores <- as.integer(rxode2::getRxThreads())

    # (a) a FINITE Cholesky that will not invert (all-zero L) -- the one way to
    # lose a subject that the cholGamma builder cannot see, counted as nBadGamma.
    # A non-finite row is skipped silently instead, with no counter.
    .chol <- matrix(as.numeric(t(chol(.map$gamma[[1L]]))), .N, 1L, byrow = TRUE)
    .chol[1L, ] <- 0                          # finite, singular -> nBadGamma
    .chol[2L, ] <- NA_real_                   # non-finite      -> skipped, no count
    fsaemDiagReset_()
    .r <- fsaemImhKernel_(.eta, .map$eta, .chol, .nchain, .cores,
                          matrix(0, .N, 1L), numeric(0), numeric(0), integer(0),
                          1000, 10L)
    expect_equal(fsaemDiag_()$nBadGamma, 1)   # subject 1 only, not subject 2
    # neither subject has a usable proposal, so both keep their current state
    expect_equal(.r$nAcc[1L], 0L)
    expect_equal(.r$nAcc[2L], 0L)
    expect_equal(.r$eta[1L, 1L], .eta[1L, 1L])

    # (b) a bound no draw can satisfy (phi = mprior + eta must exceed 1e6) -- the
    # retries are exhausted and the state is HELD, which is an exact Metropolis
    # rejection.  Clamping instead would break the acceptance ratio.
    .cholOk <- matrix(as.numeric(t(chol(.map$gamma[[1L]]))), .N, 1L, byrow = TRUE)
    fsaemDiagReset_()
    .r2 <- fsaemImhKernel_(.eta, .map$eta, .cholOk, .nchain, .cores,
                           matrix(0, .N, 1L), rep(1e6, 1L), rep(Inf, 1L), 1L,
                           2000, 3L)
    .d2 <- fsaemDiag_()
    expect_equal(.d2$nBoundFail, .nchain * .N)  # every proposal ran out of retries
    expect_equal(.d2$nAcc, 0)
    expect_true(all(.r2$eta == .eta))           # state held, nothing clamped
  })
})
