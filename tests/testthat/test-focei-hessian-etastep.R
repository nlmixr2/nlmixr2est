nmTest({
  # Regression test for the innerOpt="trust" general-likelihood failure.
  #
  # A dnorm()/ll() endpoint sets needOptimHess, so the per-subject inner eta
  # Hessian has no Gauss-Newton shortcut and is finite-differenced from the
  # analytic eta gradient.  shi21's step search is told that function's noise
  # floor is rxControl(atolSens=), but the gradient comes out of a sensitivity
  # solve that rtolSens governs as well -- so the noise is understated, shiRC()'s
  # ratio is inflated past `ru`, and the search shrinks the step until the
  # difference is taken inside the noise.  n1qn1 absorbed that (it uses this
  # Hessian once, as a warm-start seed, and steers by exact gradients);
  # innerOpt="trust" rebuilds it at every trial point and adds its
  # log-determinant to the reported objective, so the noise steers the fit.
  #
  # foceiControl(hessEtaStepMin=) floors the step at a fraction of each eta's own
  # SD, which stops the runaway without paying for a tighter solve.
  #
  # The failure needs ORAL absorption (or two compartments) to show up: on the
  # 1-cmt IV bolus fixture in test-focei-hessian-method.R the unfloored search is
  # still fine, which is why that test could not catch this.  Measured on the CPT
  # benchmark corpus, 1-cmt bolus/infusion models (U001-U022) agreed with plain
  # focei to within 2.2% on Vc while every oral / 2-cmt model from U023 on
  # diverged -- U024 gave Vc 98.3 against 67.1.  On U023 itself (1-cmt oral, 120
  # subjects) the unfloored fit gave Vc 90.6 against a plain-focei 66.4, with eta
  # variance shrinkage 56/52/34% against 8/10/13% and the omegas pinned at their
  # starting values; the floor recovers Vc 65.5.
  .oralSim <- rxode2::rxode2({
    Cl <- exp(lCl + eta.Cl)
    Vc <- exp(lVc + eta.Vc)
    KA <- exp(lKA + eta.KA)
    d / dt(depot) <- -KA * depot
    d / dt(centr) <- KA * depot - (Cl / Vc) * centr
    cp <- centr / Vc
  })

  .oral <- function() {
    ini({
      lCl <- 1.8
      lVc <- 4.7
      lKA <- 0.2
      prop.err <- c(0, 0.3, 1)
      eta.Cl ~ 0.15
      eta.Vc ~ 0.15
      eta.KA ~ 0.15
    })
    model({
      Cl <- exp(lCl + eta.Cl)
      Vc <- exp(lVc + eta.Vc)
      KA <- exp(lKA + eta.KA)
      d / dt(depot) <- -KA * depot
      d / dt(centr) <- KA * depot - (Cl / Vc) * centr
      cp <- centr / Vc
      cp ~ prop(prop.err)
    })
  }

  .mkOralData <- function(nsub = 40L) {
    .th <- c(lCl = log(4), lVc = log(70), lKA = log(1))
    .ev <- rxode2::et(amt = 10000, cmt = "depot")
    .ev <- rxode2::et(.ev, c(0.25, 0.5, 1, 2, 4, 6, 8, 12, 16, 24))
    .ev <- rxode2::et(.ev, id = seq_len(nsub))
    .sim <- rxode2::rxSolve(.oralSim, .th, .ev,
      omega = lotri::lotri(eta.Cl ~ 0.09, eta.Vc ~ 0.09, eta.KA ~ 0.09),
      addDosing = FALSE, returnType = "data.frame"
    )
    .obs <- data.frame(
      ID = .sim$id, TIME = .sim$time,
      DV = .sim$cp * (1 + 0.2 * stats::rnorm(nrow(.sim))),
      AMT = NA_real_, EVID = 0
    )
    .dose <- data.frame(
      ID = seq_len(nsub), TIME = 0, DV = NA_real_,
      AMT = 10000, EVID = 1
    )
    rbind(.dose, .obs)
  }

  .shrinkVar <- function(f) {
    .s <- f$shrink
    as.numeric(unlist(.s["var shrinkage (%)", grep("^eta", colnames(.s))]))
  }

  test_that("a dnorm() oral fit under innerOpt='trust' tracks its plain-focei counterpart", {
    skip_on_cran()
    .testSeed(2308)
    .dat <- .mkOralData()
    .oralLL <- .oral |> model(cp ~ prop(prop.err) + dnorm())

    # Reference: the same statistical model with the Gauss-Newton inner Hessian
    # (a normal endpoint never reaches the needOptimHess branch).
    .fPlain <- suppressWarnings(nlmixr2(.oral, .dat,
      est = "focei",
      control = foceiControl(print = 0L, covMethod = "")
    ))
    expect_true(is.finite(.fPlain$objf))

    # Default control: the step floor is on (hessEtaStepMin = 0.05).
    .fLL <- suppressWarnings(nlmixr2(.oralLL, .dat,
      est = "focei",
      control = foceiControl(print = 0L, covMethod = "", innerOpt = "trust")
    ))
    expect_true(is.finite(.fLL$objf))
    expect_equal(unname(.fLL$theta), unname(.fPlain$theta), tolerance = 0.05)
    # The etas have to actually reach the conditional mode.  A stalled inner
    # solve shows up as shrinkage far above the reference fit's -- and as omegas
    # left near their starting values -- long before it shows in the objective,
    # so assert on those rather than on the objective alone.
    expect_true(max(.shrinkVar(.fLL)) < max(.shrinkVar(.fPlain)) + 20)
    expect_equal(unname(diag(.fLL$omega)), unname(diag(.fPlain$omega)),
      tolerance = 0.3
    )
  })

  test_that("hessEtaStepMin is what makes that work (the floor is load bearing)", {
    skip_on_cran()
    .testSeed(2308)
    .dat <- .mkOralData()
    .oralLL <- .oral |> model(cp ~ prop(prop.err) + dnorm())
    .fOn <- suppressWarnings(nlmixr2(.oralLL, .dat,
      est = "focei",
      control = foceiControl(print = 0L, covMethod = "", innerOpt = "trust")
    ))
    # hessEtaStepMin = 0 restores the plain absolute shi21hMin floor, i.e. the
    # behavior that produced the divergence.  Assert the floor actually changes
    # the answer, so a future change that quietly stops applying it fails here
    # rather than silently reintroducing the bias.
    .fOff <- suppressWarnings(nlmixr2(.oralLL, .dat,
      est = "focei",
      control = foceiControl(
        print = 0L, covMethod = "", innerOpt = "trust",
        hessEtaStepMin = 0
      )
    ))
    expect_true(is.finite(.fOff$objf))
    expect_false(isTRUE(all.equal(unname(.fOn$theta), unname(.fOff$theta),
      tolerance = 1e-6
    )))
  })
})
