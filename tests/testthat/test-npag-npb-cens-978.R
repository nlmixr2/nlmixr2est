## #978: npag/npb's residual-error moment (npResidMoments(), inner.cpp) folded an
## M3/M4 censored row's recorded LOQ/limit into the additive/proportional
## sum-of-squares as if it had been a real measurement.  For the common
## allSimpleScale configuration (one add()/prop() SD per endpoint, no regressor
## theta), npOptimizeResid() installs that moment DIRECTLY with no further
## optimizer correction, so the biased value was the final residual-error
## estimate, not just a warm start. The model below deliberately uses proper
## mu-referencing (ke <- exp(tvK + bsvK), not tvK * exp(bsvK)) -- npag treats a
## non-mu-referenced theta as a regressor, which disqualifies allSimpleScale on
## its own and would let this test pass even with the bug still present.
##
## The regression check below is framework-independent: it does not compare
## against focei/saem's residual SD (npag/npb's ELS-at-fixed-eta moment is a
## different object from focei/saem's marginal-likelihood residual by design --
## see npag.cpp's "final support refinement" comment -- so those are not a valid
## ground truth here).  Instead it compares a fit WITHOUT any censored rows
## against the same fit with EXTRA M3-censored rows appended whose recorded DV
## is deliberately absurd (50, far outside the model's range). If those rows
## are correctly excluded from the moment, the residual estimate barely moves;
## if their absurd recorded DV leaks into the moment (the bug), the residual
## estimate is grossly distorted. A parallel check with M2 (not M3/M4) rows
## confirms those are correctly NOT excluded -- an inverted cens check would
## fail that one instead. Real fit -> weekly slow batch.

nmTest({
  .cens978Mod <- function() {
    ini({
      tvK <- log(0.5)
      bsvK ~ 0.04
      prop.sd <- sqrt(0.1)
    })
    model({
      ke <- exp(tvK + bsvK)
      v <- 1
      ipre <- 10 * exp(-ke * t)
      ipre ~ prop(prop.sd)
    })
  }

  .cens978Dat <- nlmixr2data::Wang2007
  .cens978Dat$DV <- .cens978Dat$Y
  .cens978Dat <- .cens978Dat[, names(.cens978Dat) != "Y"]
  .cens978Dat$cens <- 0

  .cens978DatCensored <- rbind(.cens978Dat,
                               data.frame(ID = 1:10, Time = 1.5, DV = 50, cens = 1))
  .cens978DatCensored <- .cens978DatCensored[order(.cens978DatCensored$ID,
                                                   .cens978DatCensored$Time), ]

  # M2: cens == 0 with a finite `limit` -- a real, defined observation that must
  # NOT be excluded from the moment (censEst.h's isM2()).  An absurd DV (50)
  # marked M2 instead of M3/M4 should distort the moment just like an ordinary
  # uncensored outlier would -- i.e. the OPPOSITE of the M3/M4 check above.
  .cens978DatM2 <- rbind(.cens978Dat,
                         data.frame(ID = 1:10, Time = 1.5, DV = 50, cens = 0))
  .cens978DatM2$limit <- 0
  .cens978DatM2 <- .cens978DatM2[order(.cens978DatM2$ID, .cens978DatM2$Time), ]

  test_that("est='npag' residual moment excludes M3/M4 rows from the moment (#978)", {
    .ctl <- npagControl(points = 128L, cycles = 20L, gammaOptimize = FALSE,
                        calcTables = FALSE, seed = 1L)
    f.base <- suppressMessages(suppressWarnings(
      nlmixr2(.cens978Mod, .cens978Dat, est = "npag", control = .ctl)
    ))
    f.cens <- suppressMessages(suppressWarnings(
      nlmixr2(.cens978Mod, .cens978DatCensored, est = "npag", control = .ctl)
    ))
    expect_match(as.character(f.cens$censInformation), "^M3 censoring")
    # An absurd recorded DV (50) leaking into the moment as if observed would
    # blow the proportional SD up by orders of magnitude (measured pre-fix:
    # 0.030 -> 5.34); properly excluded, it barely moves the estimate.
    expect_equal(as.numeric(f.cens$theta[["prop.sd"]]),
                 as.numeric(f.base$theta[["prop.sd"]]),
                 tolerance = 0.05)

    f.m2 <- suppressMessages(suppressWarnings(
      nlmixr2(.cens978Mod, .cens978DatM2, est = "npag", control = .ctl)
    ))
    expect_match(as.character(f.m2$censInformation), "^M2 censoring")
    # M2 keeps its real (here: absurd) DV, so this should NOT stay close to
    # baseline -- confirms cens==0 rows are not being wrongly excluded too.
    expect_false(isTRUE(all.equal(as.numeric(f.m2$theta[["prop.sd"]]),
                                  as.numeric(f.base$theta[["prop.sd"]]),
                                  tolerance = 0.05)))
  })

  test_that("est='npb' residual moment excludes M3/M4 rows from the moment (#978)", {
    .ctl <- npbControl(points = 128L, burnin = 200L, nsamp = 200L,
                       calcTables = FALSE, seed = 1L)
    f.base <- suppressMessages(suppressWarnings(
      nlmixr2(.cens978Mod, .cens978Dat, est = "npb", control = .ctl)
    ))
    f.cens <- suppressMessages(suppressWarnings(
      nlmixr2(.cens978Mod, .cens978DatCensored, est = "npb", control = .ctl)
    ))
    expect_match(as.character(f.cens$censInformation), "^M3 censoring")
    # pre-fix this collapsed the moment to exactly 0 (measured: 0.026 -> 0)
    expect_equal(as.numeric(f.cens$theta[["prop.sd"]]),
                 as.numeric(f.base$theta[["prop.sd"]]),
                 tolerance = 0.05)
    expect_true(as.numeric(f.cens$theta[["prop.sd"]]) > 0.001)
  })
})
