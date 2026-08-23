nmTest({
  # This file pins a "twin" pair -- one endpoint written as the closed-form
  # add() residual, one written with the exact same likelihood via a
  # general-likelihood family (dt() with a large, near-Gaussian degrees of
  # freedom) -- so any divergence between the twins is attributable to the
  # general-likelihood theta machinery, not to a genuinely different optimum.
  #
  # FIXED (nlmixr2/nlmixr2est#999): the t(df=30) twin used to diverge from its
  # add() twin.  Root cause was DV never reaching a general-likelihood SAEM
  # model's solve at all: rxUiGet.saemInParsAndMuRefCovariates
  # (R/saemRxUiGet.R) built inPars purely from ui$covariates, which never
  # includes DV (rxode2's etTran.cpp deliberately excludes any "dv"-named
  # column from covariate detection), so .configsaem's own DV-append step
  # (R/saem_fit.R, gated on "DV" already being listed in inPars) silently
  # never fired -- DV read as a constant 0 for the life of the fit, for EVERY
  # general-likelihood distribution (dnorm()/dt()/dcauchy()/literal ll()),
  # not just t().  Fixing inPars to include DV for a general-lik model fixed
  # this twin (and the SAEM general-likelihood MM-elimination corpus
  # divergence this was investigated alongside) without touching the
  # phi1Hessian/do_mcmc machinery at all -- confirming that machinery was
  # never the problem for this twin.
  mAdd <- function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45
      add.sd <- fixed(0.7)
      eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
    })
    model({
      ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }

  mT <- function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45
      add.sd <- fixed(0.7)
      nu <- fixed(30)
      eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
    })
    model({
      ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
      cp <- linCmt()
      cp ~ add(add.sd) + dt(nu)
    })
  }

  test_that("near-Gaussian t() twin matches its add() twin (SAEM general-lik phi1, #999)", {
    ctl <- saemControl(nBurn = 200, nEm = 300, seed = 42L, print = 0L,
                       covMethod = "", calcTables = FALSE)
    fA <- suppressWarnings(.nlmixr(mAdd, theo_sd, est = "saem", control = ctl))
    fT <- suppressWarnings(.nlmixr(mT, theo_sd, est = "saem", control = ctl))

    # the add() twin recovers the known-good estimate (matches #871's twin)
    expect_equal(unname(fixef(fA)[["tka"]]), 0.45, tolerance = 0.1)

    # the t(df=30) twin, mathematically equivalent to add() here, lands
    # within the same tolerance band -- confirms #999 stays fixed.
    .tkaT <- unname(fixef(fT)[["tka"]])
    .tkaA <- unname(fixef(fA)[["tka"]])
    expect_equal(.tkaT, .tkaA, tolerance = 0.1)
  })
})
