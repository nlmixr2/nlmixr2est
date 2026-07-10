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
      .r <- fsaemImhKernel_(.eta, .map$eta, .chol, .nchain, 1L)
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
})
