nmTest({
  # nlmixr2/nlmixr2est#455: combined IV + depot crossover encoded with two
  # overlapping-time evid=4 reset episodes per subject.  The SAEM kernel solves
  # each subject in the ODE solver's internal time-sorted order, which relocated
  # the reset record ahead of the first episode's observations and collapsed the
  # two episodes into one merged trajectory -- SAEM then reported an essentially
  # constant PRED and a grossly inflated residual.  .saemMonotonicResetTime()
  # offsets each reset episode as a block so the solve times increase within a
  # subject and the kernel's time-sort matches input order (predictions unchanged).

  .nm <- loadNamespace("nlmixr2est")

  test_that(".saemMonotonicResetTime offsets overlapping reset episodes", {
    # post-etTrans layout: ID, TIME, EVID, ... ; evid=3 marks the reset.  Two
    # arms share clock times 0..12; the second arm is introduced by a reset.
    dat <- data.frame(
      ID   = rep(1L, 15),
      TIME = c(0, 0.25, 0.5, 1, 2, 5, 12,   0,   0, 0.25, 0.5, 1, 2, 5, 12),
      EVID = c(101, 0, 0, 0, 0, 0, 0,       3, 201,    0,   0, 0, 0, 0, 0))
    out <- .nm$.saemMonotonicResetTime(dat)
    # first episode is untouched
    expect_equal(out$TIME[1:7], dat$TIME[1:7])
    # the reset and everything after it are shifted strictly past the running max
    expect_true(all(diff(out$TIME[7:15]) >= 0))
    expect_true(out$TIME[8] > out$TIME[7])
    # within-episode spacing is preserved (offset applied as a block)
    expect_equal(diff(out$TIME[8:15]), diff(dat$TIME[8:15]))
  })

  test_that(".saemMonotonicResetTime is a no-op for monotonic data", {
    # single-episode, monotonically increasing time -- unchanged
    dat1 <- data.frame(ID = rep(1L, 7),
                       TIME = c(0, 0.25, 0.5, 1, 2, 5, 12),
                       EVID = c(101, 0, 0, 0, 0, 0, 0))
    expect_equal(.nm$.saemMonotonicResetTime(dat1)$TIME, dat1$TIME)
    # a reset whose time is already monotonically later must not be shifted
    dat2 <- data.frame(ID = rep(1L, 4),
                       TIME = c(0, 12, 100, 112),
                       EVID = c(101, 0, 3, 0))
    expect_equal(.nm$.saemMonotonicResetTime(dat2)$TIME, dat2$TIME)
  })

  test_that(".saemMonotonicResetTime offsets each subject independently", {
    dat <- data.frame(
      ID   = c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L),
      TIME = c(0, 12, 0, 12,  0, 12, 0, 12),
      EVID = c(101, 0, 3, 0, 101, 0, 3, 0))
    out <- .nm$.saemMonotonicResetTime(dat)
    # each subject resets its own offset; times increase within each subject
    expect_true(all(diff(out$TIME[out$ID == 1L]) >= 0))
    expect_true(all(diff(out$TIME[out$ID == 2L]) >= 0))
    expect_equal(out$TIME[out$ID == 1L], out$TIME[out$ID == 2L])
  })

  test_that("SAEM recovers parameters for combined IV + depot crossover (evid=4)", {
    combined <- function() {
      ini({
        tcl1 <- log(2.1)
        tvc <- log(25)
        tka <- log(.8)
        flogit <- log(4)
        eta.cl1 ~ 0.3
        eta.vc ~ 0.3
        eta.ka ~ 0.3
        eta.f ~ 0.3
        add.sd <- c(0, 0.2, 5)
      })
      model({
        cl <- exp(tcl1 + eta.cl1)
        vc <- exp(tvc + eta.vc)
        ka <- exp(tka + eta.ka)
        ftemp <- exp(flogit + eta.f)
        f <- ftemp / (1 + ftemp)
        f(depot) <- f
        cp <- linCmt()
        lcp <- log(cp)
        lcp ~ add(add.sd)
      })
    }
    nID <- 12
    pkTimes <- c(0, .25, .5, 1, 2, 5, 12)
    mk <- function(i) {
      b <- data.frame(id = i, time = pkTimes)
      b$AMT <- ifelse(b$time == 0, 1000, NA)
      b$CMT <- 1 + is.na(b$AMT)
      b$DV <- NA
      b$evid <- ifelse(b$time == 0, 4, 0)  # depot arm
      b2 <- b
      b2$CMT <- 2                          # IV arm, same clock times, own reset
      rbind(b, b2)
    }
    template <- do.call(rbind, lapply(seq_len(nID), mk))
    set.seed(1)
    sim <- rxode2::rxSolve(combined, template)
    dat <- template
    dat$DV[is.na(dat$AMT)] <- sim$sim

    fit <- nlmixr2(combined, dat, "saem",
                   control = saemControl(covMethod = "", print = 0,
                                         nBurn = 100, nEm = 100, seed = 1))
    .f <- as.numeric(exp(fit$theta["flogit"]) / (1 + exp(fit$theta["flogit"])))
    # Before the fix SAEM collapsed both arms: add.sd inflated to ~0.78 and f
    # pushed to ~0.98 with a near-constant PRED.  With episodes kept separate the
    # data-generating values (add.sd = 0.2, f = 0.8) are recovered.
    expect_true(as.numeric(fit$theta["add.sd"]) < 0.4)
    expect_true(.f > 0.6 && .f < 0.9)
    # PRED is no longer collapsed to a single value
    expect_true(stats::sd(fit$PRED) > 0.3)

    # The low-level saemDopred* diagnostics rebuild the event table from the saved
    # data; they must apply the same reset offset as the fit, otherwise they solve a
    # merged trajectory.  With the offset the individual predictions track the
    # observations (cor ~0.97); a merged solve collapses this to ~0.5.
    .ip <- fit$saemDopredIpred
    .dvObs <- dat$DV[!is.na(dat$DV)]
    expect_equal(length(.ip), length(.dvObs))
    expect_true(stats::cor(.dvObs, .ip) > 0.85)
  })
})
