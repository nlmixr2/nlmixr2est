nmTest({
  # Issue #1010 -- foceiSetupEta_() sized the per-subject residual variance
  # block gVid as nall^2 when it is only ever indexed by sum(nobs_i^2).  nall
  # counts dose (and evid=2) records too, so the `nall > 65535` guard that
  # turned the resulting 32-bit overflow into a clean error rejected ordinary
  # population datasets needing well under a megabyte.

  .i1010Model <- function() {
    ini({
      tka <- log(1)
      tcl <- log(10)
      tv <- log(100)
      eta.cl ~ 0.1
      prop.err <- 0.1
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      linCmt() ~ prop(prop.err)
    })
  }

  test_that("focei fits a dataset whose expanded record count exceeds 65535 (#1010)", {
    # 20 subjects x (3400 daily doses + 5 samples) = 68100 expanded records,
    # but only 5 observations each: gVid needs 20*5^2 = 500 doubles, while the
    # old sizing asked for 68100^2 (37 GB) and the guard refused outright.
    nsub <- 20
    ndose <- 3400
    nobs <- 5
    d <- do.call(rbind, lapply(seq_len(nsub), function(id) {
      rbind(data.frame(ID = id, TIME = 0, DV = NA, AMT = 100, EVID = 1,
                       CMT = 1, II = 1, ADDL = ndose - 1),
            data.frame(ID = id, TIME = seq(0.5, ndose - 0.5, length.out = nobs),
                       DV = 5, AMT = 0, EVID = 0, CMT = 1, II = 0, ADDL = 0))
    }))
    expect_equal(nrow(d), nsub * (nobs + 1))

    fit <- suppressMessages(suppressWarnings(
      nlmixr(.i1010Model(), d, "focei",
             control = foceiControl(maxOuterIterations = 0, covMethod = "",
                                    print = 0))))
    # dataSav is the ADDL-expanded event table, i.e. exactly the `nall` the old
    # guard tested -- assert the fit really did cross the old threshold.
    expect_gt(nrow(fit$dataSav), 65535)
    expect_true(is.finite(fit$objf))
    expect_true(all(is.finite(fit$CWRES)))
    expect_equal(nrow(fit), nsub * nobs)
  })

  test_that("gVid holds sum(nobs_i^2) with unequal per-subject obs counts", {
    # gcHff/gcHfr/gcHrr are allocated immediately after gVid, so an undersized
    # gVid corrupts the censored inner-Hessian coefficients.  Unequal nobs_i is
    # what distinguishes sum(nobs_i^2) from any per-subject-constant sizing.
    nobsI <- c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11)
    d <- do.call(rbind, lapply(seq_along(nobsI), function(id) {
      rbind(data.frame(ID = id, TIME = 0, DV = NA, AMT = 100, EVID = 1,
                       CMT = 1, II = 24, ADDL = 6),
            data.frame(ID = id, TIME = seq(0.5, 24, length.out = nobsI[id]),
                       DV = 5, AMT = 0, EVID = 0, CMT = 1, II = 0, ADDL = 0))
    }))
    d$CENS <- 0
    # censor each subject's first observation so the M3 path (and with it
    # gcHff/gcHfr/gcHrr) is exercised
    .obs <- which(d$EVID == 0)
    d$CENS[.obs[!duplicated(d$ID[.obs])]] <- 1
    d$DV[d$CENS == 1] <- 1

    fit <- suppressMessages(suppressWarnings(
      nlmixr(.i1010Model(), d, "focei",
             control = foceiControl(maxOuterIterations = 0, covMethod = "",
                                    print = 0))))
    expect_equal(nrow(fit), sum(nobsI))
    expect_match(as.character(fit$censInformation), "^M3 censoring")
    expect_true(is.finite(fit$objf))
    expect_true(all(is.finite(fit$IPRED)))
    expect_true(all(is.finite(as.numeric(fit$eta.cl))))
  })

  test_that("a no-eta (neta == 0) focei fit gives thetaGrad its own block", {
    # foceiSetupNoEta_ sized the gthetaGrad block with gEtaGTransN, which is
    # neta * nsub -- and this path only runs when neta == 0, so the block was
    # empty and llikObsFull started at the same address.  The per-subject
    # thetaGrad writes (npars per subject) then ran past the allocation
    # whenever npars * nsub exceeded nall.
    d <- do.call(rbind, lapply(1:12, function(id) {
      rbind(data.frame(ID = id, TIME = 0, DV = NA, AMT = 320, EVID = 1, CMT = 1),
            data.frame(ID = id, TIME = 2 + id / 12, DV = 6 + id / 10, AMT = 0,
                       EVID = 0, CMT = 1))
    }))

    noEta <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl)
        v <- exp(tv)
        linCmt() ~ add(add.sd)
      })
    }

    fit <- suppressMessages(suppressWarnings(
      nlmixr(noEta(), d, "focei",
             control = foceiControl(covMethod = "", print = 0))))
    # the setup really is in the regime the old sizing overflowed: npars per
    # subject over 12 subjects is more than the whole event table
    expect_gt(length(fit$theta) * 12, nrow(fit$dataSav))
    expect_true(is.finite(fit$objf))
    # llikObsFull shared storage with gthetaGrad before the fix; the per-record
    # log-likelihoods must still add up to the objective
    expect_equal(-2 * sum(fit$llikObs, na.rm = TRUE), fit$objf)
  })

  test_that("llikObs comes back in record order, not reversed by subject", {
    # llikObsFull is the one per-subject block handed to R as a contiguous
    # array, so its per-subject offsets have to follow the subject order.  They
    # used to be a running total taken in the setup loop's BACKWARDS order, so
    # subject 1 was given the tail of the buffer and the exported per-record
    # log-likelihoods came back reversed by subject.
    d <- data.frame(ID = rep(1:3, each = 2), TIME = rep(c(0, 1), 3),
                    DV = c(NA, 1, NA, 100, NA, 20), AMT = c(320, 0, 320, 0, 320, 0),
                    EVID = c(1, 0, 1, 0, 1, 0), CMT = 1)

    # subject-order-dependent by construction: the three residuals are wildly
    # different, so a reversal is unmistakable
    .check <- function(fit) {
      ll <- fit$llikObs
      ll <- ll[!is.na(ll)]                     # drop the dose records
      sd <- fit$theta[["add.sd"]]
      expect_equal(ll, -log(sd) - 0.5 * ((fit$DV - fit$IPRED) / sd)^2)
    }

    # neta == 0 -> foceiSetupNoEta_
    noEta <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; add.sd <- 0.7 })
      model({ ka <- exp(tka); cl <- exp(tcl); v <- exp(tv)
              linCmt() ~ add(add.sd) })
    }
    .fitNoEta <- suppressMessages(suppressWarnings(
      nlmixr(noEta(), d, "focei",
             control = foceiControl(maxOuterIterations = 0, covMethod = "",
                                    print = 0))))
    .check(.fitNoEta)
    # with no etas the objective IS -2 * the per-record log-likelihood
    expect_equal(-2 * sum(.fitNoEta$llikObs, na.rm = TRUE), .fitNoEta$objf)

    # neta > 0 -> foceiSetupEta_
    wEta <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; eta.cl ~ 0.1; add.sd <- 0.7 })
      model({ ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv)
              linCmt() ~ add(add.sd) })
    }
    .check(suppressMessages(suppressWarnings(
      nlmixr(wEta(), d, "focei",
             control = foceiControl(maxOuterIterations = 0, covMethod = "",
                                    print = 0)))))
  })

  test_that("gVid scales with the mixture replicate count", {
    # gVid is sized (mixIdxN + 1) * sum(nobs_i^2) and the setup loop walks it
    # once per mixture replicate.  With mixIdxN == 0 the multiplier is 1, so a
    # dropped factor would be invisible; a mix() model with unequal
    # per-subject observation counts is what actually exercises it.
    d <- nlmixr2data::theo_sd
    # thin each subject's observations to a different count (theo_sd has 11
    # observations for every subject, which would hide a per-subject stride)
    .keep <- unlist(lapply(split(seq_len(nrow(d)), d$ID), function(w) {
      .obs <- w[d$EVID[w] == 0]
      c(w[d$EVID[w] != 0], utils::head(.obs, 2 + (as.integer(d$ID[w[1]]) %% 8)))
    }))
    d <- d[sort(.keep), ]
    nobsI <- tapply(d$EVID == 0, d$ID, sum)
    expect_gt(length(unique(as.integer(nobsI))), 1)

    mixMod <- function() {
      ini({
        tka <- 0.45
        tcl1 <- log(c(0, 2.7, 100))
        tcl2 <- log(c(0, 0.1, 120))
        tv <- 3.45
        p1 <- 0.3
        eta.cl ~ 0.3
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka)
        cl <- mix(exp(tcl1 + eta.cl), p1, exp(tcl2 + eta.cl))
        v <- exp(tv)
        linCmt() ~ add(add.sd)
      })
    }

    fit <- suppressMessages(suppressWarnings(
      nlmixr(mixMod(), d, "focei",
             control = foceiControl(maxOuterIterations = 0, covMethod = "",
                                    print = 0))))
    expect_equal(length(fit$mixList), 2L)
    expect_equal(nrow(fit), sum(as.integer(nobsI)))
    expect_true(is.finite(fit$objf))
    expect_true(all(is.finite(fit$CWRES)))
    expect_true(all(fit$mixNum$mixnum %in% c(1L, 2L)))
  })
})
