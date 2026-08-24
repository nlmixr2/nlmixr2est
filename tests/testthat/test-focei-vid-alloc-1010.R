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
})
