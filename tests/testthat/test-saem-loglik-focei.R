nmTest({
  # Phase 1 diagnostic (see the SAEM general-likelihood theta plan): SAEM's
  # general-likelihood theta step (do_mcmc's distribution==4 path, followed by
  # the plain stochastic-approximation SA-recursion for mu-referenced/phi1
  # thetas) does not minimize successfully, even on a small, well-identified
  # model.  This file pins the failure with a "twin" pair -- one endpoint
  # written as the closed-form add() residual, one written with the exact
  # same likelihood via a general-likelihood family (dt() with a large,
  # near-Gaussian degrees of freedom) -- so any divergence between the twins
  # is attributable to the general-likelihood theta machinery, not to a
  # genuinely different optimum.
  #
  # Root-cause checkpoint (measured, not assumed): a t(df=30) endpoint is
  # statistically indistinguishable from add() (df=30 is nearly Gaussian), so
  # if the twins land on different estimates, the estimation MACHINERY is at
  # fault, not the model's own likelihood shape.  A separate iteration-count
  # sweep (200/300 up to 1000/1500 burn-in/EM) at this file's authoring time
  # showed the t()-twin's tka estimate does NOT converge toward the add()
  # twin's value as iterations increase -- it stays in the wrong basin (and
  # drifts slightly further from truth), ruling out "just needs more
  # iterations" as the explanation.  The residual-family theta (add.sd) and
  # the degrees of freedom (nu) are both FIXED in the t() twin, isolating this
  # from the separate, already-fixed err-tagged-theta declaration/reporting
  # bug (see test-saem-loglik.R's "(#Phase0)" tests) -- what's left free here
  # is exactly the mu-referenced (phi1) structural thetas the wider plan's
  # Hessian-corrected theta step targets.
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

  test_that("KNOWN FAILURE: a near-Gaussian t() twin diverges from its add() twin (SAEM general-lik phi1)", {
    ctl <- saemControl(nBurn = 200, nEm = 300, seed = 42L, print = 0L,
                       covMethod = "", calcTables = FALSE)
    fA <- suppressWarnings(.nlmixr(mAdd, theo_sd, est = "saem", control = ctl))
    fT <- suppressWarnings(.nlmixr(mT, theo_sd, est = "saem", control = ctl))

    # the add() twin recovers the known-good estimate (matches #871's twin)
    expect_equal(unname(fixef(fA)[["tka"]]), 0.45, tolerance = 0.1)

    # KNOWN FAILURE (as of this writing): the t(df=30) twin, mathematically
    # equivalent to add() here, should land within the same tolerance band --
    # it instead lands nowhere near it.  This assertion is expected to FAIL
    # until the phi1 Hessian-corrected theta step (later plan phases) is
    # wired in; when it passes, tighten the tolerance to match fA's and
    # promote this from a diagnostic to a real regression guard.
    .tkaT <- unname(fixef(fT)[["tka"]])
    .tkaA <- unname(fixef(fA)[["tka"]])
    expect_equal(.tkaT, .tkaA, tolerance = 0.1)
  })
})
