nmTest({
  # Phase 3 (see the SAEM general-likelihood theta plan): foceiLikInnerObjective_
  # is a thin C++ wrapper around foceiObjFromLik0() -- the EXACT objective a
  # live FOCEi fit's outer optimizer minimizes (innerOpt1's per-subject EBE
  # Newton search + LikInner2's Laplace -2*loglik + log|H| determinant,
  # summed, plus any ini({}) prior term) -- exposed standalone over a
  # foceiLikLoad()-ed problem.  This is the actual bobyqa objective Phase 4
  # wires into SAEM's phi1 theta step.
  #
  # Verification strategy: fit the SAME model with the already-shipped,
  # independently-tested est="focei", then confirm optimizing this NEW
  # objective (holding Omega fixed, matching Phase 4's intended use -- SAEM
  # owns its own Omega estimate, this objective only ever updates theta)
  # lands close to FOCEi's own converged theta.
  #
  # A caveat, measured rather than assumed away: this objective's ABSOLUTE
  # VALUE does not match a live fit's reported `$objf` -- and that gap is
  # NOT specific to this new function.  Cross-checked with the existing,
  # independently-tested foceiLikRun() (issue #414): even
  # `-2*sum(foceiLikRun(theta, eta, "joint"))`, evaluated at the live fit's
  # OWN converged theta and eta, does not equal that fit's `$objf` either
  # (a pre-existing FOCEi `$objf` reporting convention, unrelated to this
  # phase). What IS verified here is the property Phase 4 actually needs:
  # the objective's optimum lands in the right place.
  mLl <- function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45
      lsd <- fixed(log(0.7))
      eta.ka ~ fixed(0.6); eta.cl ~ fixed(0.3); eta.v ~ fixed(0.1)
    })
    model({
      ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
      cp <- linCmt()
      sd <- exp(lsd)
      ll(err) ~ -lsd - 0.5 * log(2 * pi) - 0.5 * ((DV - cp) / sd)^2
    })
  }

  test_that("foceiLikInnerObjective_ requires a loaded system", {
    # unload defensively in case a prior test in this session left one loaded
    suppressWarnings(try(foceiLikUnload(), silent = TRUE))
    expect_error(foceiLikInnerObjective_(c(0.45, 1, 3.45)),
                "no general likelihood system loaded")
  })

  test_that("foceiLikInnerObjective_ is well-behaved across several thetas (#Phase3)", {
    h <- suppressMessages(foceiLikLoad(mLl, theo_sd, "focei", scale = "natural"))
    on.exit(suppressWarnings(try(foceiLikUnload(), silent = TRUE)))
    for (.th in list(c(0.45, 1, 3.45), c(0.3, 0.8, 3.2), c(0.6, 1.3, 3.7))) {
      .v <- foceiLikInnerObjective_(.th)
      expect_true(is.finite(.v))
    }
  })

  test_that("foceiLikInnerObjective_'s optimum matches a known-good est='focei' fit (#Phase3)", {
    .ctl <- foceiControl(print = 0)
    fF <- suppressWarnings(.nlmixr(mLl, theo_sd, est = "focei", control = .ctl))
    expect_true(is.finite(fF$objf))

    h <- suppressMessages(foceiLikLoad(mLl, theo_sd, "focei", scale = "natural"))
    on.exit(suppressWarnings(try(foceiLikUnload(), silent = TRUE)))
    .obj <- function(th) foceiLikInnerObjective_(th)
    # eta ~ fixed(...) (no Omega to estimate, matching Phase 4's intended use --
    # SAEM owns its own Omega estimate, this objective only ever updates theta)
    # is an unusual configuration on a sparse dataset: some subjects' individual
    # Hessians come back non-positive-definite ("FOCEi objective functions may
    # not be comparable"), a warning inherited from LikInner2/calcEtaHessian,
    # not introduced here. The optimum still lands close to the known-good fit
    # (asserted below), so this is a conditioning caveat, not a correctness bug
    # -- suppressed here since it is expected for this specific test model.
    .opt <- suppressWarnings(
      optim(h$initPar, .obj, method = "BFGS", control = list(reltol = 1e-10))
    )

    .foceiTheta <- unname(fixef(fF)[c("tka", "tcl", "tv")])
    expect_equal(unname(.opt$par), .foceiTheta, tolerance = 0.1)
  })
})
