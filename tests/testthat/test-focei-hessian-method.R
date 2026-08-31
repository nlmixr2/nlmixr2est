nmTest({
  test_that("foceiControl() hessianMethod validation", {
    expect_equal(foceiControl()$hessianMethod, 1L)
    expect_equal(foceiControl(hessianMethod = "fd")$hessianMethod, 1L)
    expect_equal(foceiControl(hessianMethod = "bfgs")$hessianMethod, 2L)
    expect_equal(foceiControl(hessianMethod = "sr1")$hessianMethod, 3L)
    expect_equal(foceiControl(hessianMethod = "bofill")$hessianMethod, 4L)
    expect_equal(foceiControl(hessianMethod = 3L)$hessianMethod, 3L)
    expect_error(foceiControl(hessianMethod = "not-a-method"))
  })

  # A non-normal (Poisson) endpoint: the per-subject inner Hessian has no
  # Gaussian Gauss-Newton shortcut, so calcEtaHessian()'s needOptimHess
  # branch (finite-difference by default, or hessianMethod=)'s quasi-Newton
  # update) is the ONLY inner-Hessian path exercised here.
  .poisMod <- function() {
    ini({
      te0 <- log(5)
      eta.e0 ~ 0.3
    })
    model({
      e0 <- exp(te0 + eta.e0)
      effect <- e0
      effect ~ dpois(effect)
    })
  }
  .poisData <- data.frame(ID = rep(1:20, each = 3), TIME = rep(c(0, 1, 2), 20))
  .poisData$DV <- c(5L, 4L, 6L, 3L, 5L, 4L, 6L, 5L, 7L, 4L, 5L, 5L,
                     6L, 6L, 5L, 4L, 3L, 5L, 5L, 6L, 4L, 7L, 5L, 6L,
                     4L, 5L, 5L, 3L, 4L, 6L, 6L, 5L, 4L, 5L, 6L, 5L,
                     4L, 4L, 6L, 5L, 7L, 5L, 6L, 5L, 4L, 3L, 5L, 6L,
                     5L, 5L, 6L, 4L, 6L, 5L, 5L, 4L, 5L, 6L, 4L, 5L)

  test_that("hessianMethod converges close to fd on a non-normal endpoint", {
    skip_on_cran()
    .fFd <- .nlmixr(.poisMod, .poisData, est = "focei",
                    control = foceiControl(print = 0L, hessianMethod = "fd"))
    expect_true(is.finite(.fFd$objf))
    for (.hm in c("bfgs", "sr1", "bofill")) {
      .fT <- .nlmixr(.poisMod, .poisData, est = "focei",
                     control = foceiControl(print = 0L, hessianMethod = .hm))
      expect_true(is.finite(.fT$objf), info = .hm)
      expect_equal(.fT$objf, .fFd$objf, tolerance = 1e-2, info = .hm)
      expect_equal(unname(.fT$theta), unname(.fFd$theta), tolerance = 1e-2, info = .hm)
      # Positive evidence the quasi-Newton mechanism actually ran (mirrors
      # test-focei-trust-inner.R's own .nTrustInner() convention).
      expect_true(.nHessianQN() > 0L, info = .hm)
    }
  })

  # Regression test for a real correctness bug (not just a Poisson-toy-model
  # difference): on a one-compartment IV bolus PK model fit as a general
  # ll()/dnorm() endpoint (needOptimHess), hessianMethod="bfgs"/"sr1"/"bofill"
  # all converged to the SAME wrong Vc (~90, vs a simulated 70 and a plain
  # (non-ll()) focei fit's ~67-68) with a WORSE reported objective than
  # "fd"'s correct answer. Unlike trustControl()'s outer-theta Hessian, this
  # inner Hessian's log-determinant is added directly into the reported
  # Laplace objective (LikInner2(), src/inner.cpp) -- a quasi-Newton estimate
  # accurate enough to guide the per-subject Newton step is not accurate
  # enough to serve as that objective term, and the outer optimizer chases
  # the resulting bias to a worse point. This is why "fd" stays the default
  # (see hessianMethod= roxygen); this test guards against a future default
  # flip silently reintroducing the bias.
  .ivBolus <- function() {
    ini({
      lCl <- log(4)
      lVc <- log(70)
      prop.err <- c(0, 0.1, 1)
      eta.Vc ~ 0.1
      eta.Cl ~ 0.1
    })
    model({
      Vc <- exp(lVc + eta.Vc)
      Cl <- exp(lCl + eta.Cl)
      d / dt(centr) <- -(Cl / Vc) * centr
      cp <- centr / Vc
      cp ~ prop(prop.err)
    })
  }

  test_that("default hessianMethod (fd) avoids the QN inner-Hessian log-det bias on a real PK ll() fit", {
    skip_on_cran()
    .testSeed(1104)
    .simMod <- rxode2::rxode2({
      Vc <- exp(lVc + eta.Vc)
      Cl <- exp(lCl + eta.Cl)
      d / dt(centr) <- -(Cl / Vc) * centr
      cp <- centr / Vc
    })
    .th <- c(lCl = log(4), lVc = log(70))
    .ev <- rxode2::et(amt = 60000, cmt = "centr")
    .ev <- rxode2::et(.ev, seq(0.25, 24, by = 2))
    .ev <- rxode2::et(.ev, id = seq_len(30))
    .sim <- rxode2::rxSolve(.simMod, .th, .ev,
                            omega = lotri::lotri(eta.Vc ~ 0.1, eta.Cl ~ 0.1),
                            addDosing = FALSE, returnType = "data.frame")
    .obs <- data.frame(ID = .sim$id, TIME = .sim$time,
                       DV = .sim$cp * (1 + 0.1 * stats::rnorm(nrow(.sim))),
                       AMT = NA_real_, EVID = 0)
    .dose <- data.frame(ID = seq_len(30), TIME = 0, DV = NA_real_,
                        AMT = 60000, EVID = 1)
    .dat <- rbind(.dose, .obs)
    .ivBolusLL <- .ivBolus |> model(cp ~ prop(prop.err) + dnorm())
    .fPlain <- suppressWarnings(nlmixr2(.ivBolus, .dat, est = "focei",
                                        control = foceiControl(print = 0L)))
    .fLLFd <- suppressWarnings(nlmixr2(.ivBolusLL, .dat, est = "focei",
                                       control = foceiControl(print = 0L)))
    expect_true(is.finite(.fPlain$objf))
    expect_true(is.finite(.fLLFd$objf))
    # fd (the default) should track the plain (Gauss-Newton) focei fit closely
    expect_equal(unname(.fLLFd$theta), unname(.fPlain$theta), tolerance = 0.05)
    for (.hm in c("bfgs", "sr1", "bofill")) {
      .fLLQn <- suppressWarnings(nlmixr2(.ivBolusLL, .dat, est = "focei",
                                         control = foceiControl(print = 0L, hessianMethod = .hm)))
      # A QN method drifting off is the known, not-yet-fixed bias -- assert
      # only that a caller who explicitly opts into it still gets a finite
      # result, not that it matches (it currently does not).
      expect_true(is.finite(.fLLQn$objf), info = .hm)
    }
  })

  .oneCmt <- function() {
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
      d / dt(depot) <- -ka * depot
      d / dt(centr) <- ka * depot - cl / v * centr
      cp <- centr / v
      cp ~ add(add.sd)
    })
  }

  test_that("hessianMethod has no effect on a normal-endpoint model (needOptimHess gates it off)", {
    skip_on_cran()
    .fN <- .nlmixr(.oneCmt, nlmixr2data::theo_sd, est = "focei",
                   control = foceiControl(print = 0L, hessianMethod = "sr1"))
    expect_true(is.finite(.fN$objf))
    # The Gauss-Newton inner Hessian is used unconditionally for a normal
    # endpoint -- the quasi-Newton mechanism should never engage here.
    expect_equal(.nHessianQN(), 0L)
  })

  test_that("hessianMethod is inert unless innerOpt='trust' (calcEtaHessian() is also reached from warmZm()/LikInner2() for other inner optimizers, where the per-attempt reset does not apply)", {
    skip_on_cran()
    .fN1qn1 <- .nlmixr(.poisMod, .poisData, est = "focei",
                       control = foceiControl(print = 0L, innerOpt = "n1qn1",
                                              hessianMethod = "sr1"))
    expect_true(is.finite(.fN1qn1$objf))
    # The quasi-Newton mechanism must never engage for a non-trust inner
    # optimizer, regardless of hessianMethod -- calcEtaHessian() is reached
    # from warmZm()'s one-time n1qn1 warm-start seed and from the final
    # objective recompute, neither of which is the repeated, reset-per-attempt
    # Newton-step loop the quasi-Newton state assumes.
    expect_equal(.nHessianQN(), 0L)
  })

  test_that("foceiControl() hessEtaStepMin validation and round-trip", {
    expect_equal(foceiControl()$hessEtaStepMin, 0.05)
    expect_equal(foceiControl(hessEtaStepMin = 0.1)$hessEtaStepMin, 0.1)
    expect_equal(foceiControl(hessEtaStepMin = 0)$hessEtaStepMin, 0)
    expect_error(foceiControl(hessEtaStepMin = -1))
    expect_error(foceiControl(hessEtaStepMin = c(0.1, 0.2)))
    expect_equal(do.call(foceiControl, foceiControl())$hessEtaStepMin, 0.05)
  })

  test_that("hessEtaStepMin does not disturb a normal endpoint", {
    skip_on_cran()
    # A normal endpoint uses the Gauss-Newton inner Hessian unconditionally, so
    # calcEtaHessian()'s needOptimHess finite-difference branch -- the only place
    # hessEtaStepMin is read -- is never reached.
    .fNorm0 <- .nlmixr(.oneCmt, nlmixr2data::theo_sd, est = "focei",
                       control = foceiControl(print = 0L, hessEtaStepMin = 0))
    .fNorm1 <- .nlmixr(.oneCmt, nlmixr2data::theo_sd, est = "focei",
                       control = foceiControl(print = 0L, hessEtaStepMin = 0.5))
    expect_equal(.fNorm1$objf, .fNorm0$objf, tolerance = 1e-8)
  })

  test_that("hessEtaStepMin applies to every inner optimizer, not just trust", {
    skip_on_cran()
    # The step a finite difference needs is a property of the problem, not of
    # who consumes the Hessian, so the floor must reach n1qn1 too.  This is why
    # the scale is read from curOmegaInv() (what lpInner() itself uses, always
    # current) rather than op_focei.omega, which innerOpt() refreshes only on
    # the trust branch.
    for (.io in c("trust", "n1qn1")) {
      .f0 <- .nlmixr(.poisMod, .poisData, est = "focei",
                     control = foceiControl(print = 0L, innerOpt = .io,
                                            hessEtaStepMin = 0))
      .f1 <- .nlmixr(.poisMod, .poisData, est = "focei",
                     control = foceiControl(print = 0L, innerOpt = .io,
                                            hessEtaStepMin = 0.5))
      expect_true(is.finite(.f0$objf), info = .io)
      expect_true(is.finite(.f1$objf), info = .io)
    }
    # Only trust is asserted to MOVE.  n1qn1 uses this Hessian as a warmZm seed
    # and then corrects it with its own quasi-Newton updates as it iterates, so
    # a changed seed can wash out by convergence -- that self-correction is
    # exactly why the unfloored step never broke n1qn1.  Trust re-derives the
    # Hessian as its trust-region model at every trial point, with nothing to
    # correct it, so there the floor must show.
    .fT0 <- .nlmixr(.poisMod, .poisData, est = "focei",
                    control = foceiControl(print = 0L, innerOpt = "trust",
                                           hessEtaStepMin = 0))
    .fT1 <- .nlmixr(.poisMod, .poisData, est = "focei",
                    control = foceiControl(print = 0L, innerOpt = "trust",
                                           hessEtaStepMin = 0.5))
    expect_false(isTRUE(all.equal(.fT1$objf, .fT0$objf, tolerance = 1e-8)))
  })
})
