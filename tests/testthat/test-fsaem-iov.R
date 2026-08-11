nmTest({
  # f-SAEM with non-mu etas (issue #875).  Every IOV eta is a nonMuEta, and a
  # nonMuEta appends a phi column that has no theta -- so the SAEM phi vector
  # stops being a positional copy of the structural thetas and the fast kernel's
  # inner THETA has to be filled by NAME.  fsaemDiag()$lastTheta is the inner
  # THETA the last fast step was built with; the counters alone cannot see a
  # mis-parameterized inner (it costs acceptance, it does not crash).

  .iovData <- function() {
    .d <- nlmixr2data::theo_md
    .d$occ <- 1
    .d$occ[.d$TIME >= 144] <- 2
    .d
  }

  # fastKernel="throughout" so the LAST fast step is the last iteration and
  # lastTheta can be compared against the fit's own estimates; the default
  # "firstN" stops the kernel at iteration 20 and lastTheta is then a snapshot
  # of a much earlier chain.
  .iovCtl <- function() {
    saemControl(nBurn = 60, nEm = 30, nmc = 3, seed = 42, print = 0L,
                calcTables = FALSE, fastKernel = "throughout")
  }

  test_that("est='fsaem' fits an IOV model and keeps the inner's IOV sd live", {
    skip_if_not_installed("nlmixr2data")
    # Every structural theta is mu-referenced here, so the mu-ref slots line up
    # positionally by accident.  The IOV magnitude theta does NOT: it is a phi0,
    # and the positional read handed the inner the prior mean of a nonMuEta phi
    # column -- identically 0 for every iteration, i.e. an inner with no IOV.
    .iov <- function() {
      ini({
        tka <- log(1.2); tcl <- log(0.25); tv <- log(5)
        eta.ka ~ 0.1; eta.cl ~ 0.09; eta.v ~ 0.05
        iov.cl ~ 0.05 | occ
        prop.sd <- 0.15
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl + iov.cl)
        v  <- exp(tv + eta.v)
        d/dt(depot) <- -ka*depot
        d/dt(central) <- ka*depot - cl/v*central
        cp <- central/v
        cp ~ prop(prop.sd)
      })
    }
    .d <- .iovData()
    .fs <- suppressMessages(nlmixr2(.iov, .d, est = "fsaem", control = .iovCtl()))
    .ss <- suppressMessages(nlmixr2(.iov, .d, est = "saem", control = .iovCtl()))

    # the fast kernel really ran, and ran healthily
    expect_gt(.fs$fsaemDiag$nStep, 0)
    expect_gt(.fs$fsaemDiag$accRate, 0.5)
    expect_equal(.fs$fsaemDiag$nMapFail, 0)

    # the inner carries all five thetas (tka, tcl, tv, prop.sd, iov.cl) ...
    .lt <- .fs$fsaemDiag$lastTheta
    expect_equal(length(.lt), 5L)
    # ... the structural + residual ones at the fit's own values ...
    .fe <- fixef(.fs)
    expect_lt(max(abs(.lt[1:4] - .fe[c("tka", "tcl", "tv", "prop.sd")])), 0.05)
    # ... and the IOV magnitude at the chain's live value, not the 0 the
    # nonMuEta phi column would have supplied
    expect_gt(abs(.lt[5]), 0)

    # A smoke band, NOT the pin: a wrong inner theta shifts the target without
    # crashing, so the estimates are the last thing to notice it.  Plain SAEM is
    # itself not converged at 90 iterations here (measured gap 0.135 on tcl).
    expect_lt(max(abs(.fe - fixef(.ss))), 0.3)
  })

  test_that("est='fsaem' fills an IOV model's phi0 theta from its own phi column", {
    skip_if_not_installed("nlmixr2data")
    # tv carries no eta, so it is a phi0 sitting BETWEEN two phi1 columns.  The
    # positional read walked mprior_phi1 instead, which skips phi0 columns, so
    # tv was filled from an IOV eta's prior mean (0) rather than from log(V).
    .iovPhi0 <- function() {
      ini({
        tka <- log(1.2); tcl <- log(0.25); tv <- log(5)
        eta.ka ~ 0.1; eta.cl ~ 0.09
        iov.cl ~ 0.05 | occ
        prop.sd <- 0.15
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl + iov.cl)
        v  <- exp(tv)
        d/dt(depot) <- -ka*depot
        d/dt(central) <- ka*depot - cl/v*central
        cp <- central/v
        cp ~ prop(prop.sd)
      })
    }
    .d <- .iovData()
    .fs <- suppressMessages(nlmixr2(.iovPhi0, .d, est = "fsaem", control = .iovCtl()))
    .ss <- suppressMessages(nlmixr2(.iovPhi0, .d, est = "saem", control = .iovCtl()))

    .lt <- .fs$fsaemDiag$lastTheta
    .fe <- fixef(.fs)
    expect_equal(length(.lt), 5L)
    # the phi0 slot, which is what the fallback got wrong (inner tv was 0)
    expect_lt(abs(.lt[3] - .fe[["tv"]]), 0.05)
    expect_lt(max(abs(.lt[1:4] - .fe[c("tka", "tcl", "tv", "prop.sd")])), 0.05)

    # a mis-filled inner theta shows up in the acceptance rate long before the
    # estimates: this model ran at 0.50 with tv = 0 and at 0.89 with tv = log(V)
    expect_gt(.fs$fsaemDiag$accRate, 0.7)
    expect_equal(.fs$fsaemDiag$nMapFail, 0)
    expect_lt(max(abs(.fe - fixef(.ss))), 0.1)
  })

  test_that("est='fsaem' degrades covariates + IOV instead of dying", {
    skip_if_not_installed("nlmixr2data")
    # The covariate inner absorbs each mu-ref intercept into a per-subject data
    # column and is therefore sized by muRefDataFrame; an IOV eta has no row
    # there, so the setup used to abort the whole fit with "etaMat must have the
    # same number of ETAs".  Degrade to standard SAEM instead.
    .iovCov <- function() {
      ini({
        tka <- log(1.2); tcl <- log(0.25); tv <- log(5); covwt <- 0.1
        eta.ka ~ 0.1; eta.cl ~ 0.09; eta.v ~ 0.05
        iov.cl ~ 0.05 | occ
        prop.sd <- 0.15
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl + covwt*log(WT/70) + iov.cl)
        v  <- exp(tv + eta.v)
        d/dt(depot) <- -ka*depot
        d/dt(central) <- ka*depot - cl/v*central
        cp <- central/v
        cp ~ prop(prop.sd)
      })
    }
    .d <- .iovData()
    .d$WT <- ifelse(.d$ID %% 2 == 0, 70, 90)
    .fs <- suppressMessages(nlmixr2(.iovCov, .d, est = "fsaem",
                                    control = saemControl(nBurn = 20, nEm = 10, nmc = 3,
                                                          seed = 42, print = 0L,
                                                          calcTables = FALSE)))
    expect_equal(.fs$fsaemDiag$nStep, 0)   # the fast kernel was not installed
    expect_true(all(is.finite(fixef(.fs))))
  })
})
