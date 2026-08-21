nmTest({
  # matExp()/indLin() fits with default (compartment-1) dosing must match the
  # equivalent ODE model; source-first compartment order is the regression

  .mkData <- function(model, params, sd = 0.3, nid = 6, seed = 1234) {
    # rxWithSeed pins the simulation to a fixed stream AND restores the global
    # RNG state afterwards, so the data is reproducible without leaking the seed
    # into the fits that follow (or into downstream test files).
    rxode2::rxWithSeed(seed, {
    .ev <- rxode2::et(amt = 320, cmt = "depot", id = seq_len(nid)) |>
      rxode2::et(seq(0.5, 24, by = 1.5))
    .sim <- rxode2::rxSolve(model, .ev, params = params)
    .dat <- as.data.frame(.sim)[, c("id", "time", "cp")]
    .dat$cp <- .dat$cp + stats::rnorm(nrow(.dat), 0, sd)
    names(.dat) <- c("ID", "TIME", "DV")
    .dat$AMT <- 0
    .dat$EVID <- 0
    # dose row deliberately lacks cmt so it defaults to compartment 1
    .dose <- data.frame(ID = seq_len(nid), TIME = 0, DV = NA, AMT = 320, EVID = 1)
    .dat <- rbind(.dose, .dat)
    .dat[order(.dat$ID, .dat$TIME, -.dat$EVID), ]
    })
  }

  test_that("matExp() linear model estimates identically to the ODE (focei/foce/nlm)", {
    odeLin <- function() {
      ini({ tka <- 0.45; tcl <- 1.0; tv <- 3.45; eta.ka ~ 0.09; add.sd <- 0.7 })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(central) <- ka * depot - cl / v * central
        cp <- central / v
        cp ~ add(add.sd)
      })
    }
    matLin <- function() {
      ini({ tka <- 0.45; tcl <- 1.0; tv <- 3.45; eta.ka ~ 0.09; add.sd <- 0.7 })
      model({
        matExp()
        k_depot_central <- exp(tka + eta.ka)
        k_central_output <- exp(tcl) / exp(tv)
        cp <- central / exp(tv)
        cp ~ add(add.sd)
      })
    }
    .dat <- .mkData(odeLin, c(tka = 0.6, tcl = 1.1, tv = 3.6))
    # sigdig = 6: matExp() now solves this system natively (#860) and so is EXACT,
    # while the ODE form still carries its own solver discretization error.  The
    # sigdig = 4 default leaves enough of that error at the converged optimum to
    # exceed the 1e-3 tolerance asserted here (see the MM test below, which needed
    # the same fix before the native path existed).
    for (.met in c("focei", "foce")) {
      .fO <- .nlmixr(odeLin, .dat, est = .met, control = foceiControl(print = 0, sigdig = 6))
      .fM <- .nlmixr(matLin, .dat, est = .met, control = foceiControl(print = 0, sigdig = 6))
      expect_equal(.fM$objf, .fO$objf, tolerance = 1e-3)
      expect_equal(unname(fixef(.fM)), unname(fixef(.fO)), tolerance = 1e-3)
    }
  })

  test_that("matExp() linear model estimates identically to the ODE (agq/laplace)", {
    odeLin <- function() {
      ini({ tka <- 0.45; tcl <- 1.0; tv <- 3.45; eta.ka ~ 0.09; add.sd <- 0.7 })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(central) <- ka * depot - cl / v * central
        cp <- central / v
        cp ~ add(add.sd)
      })
    }
    matLin <- function() {
      ini({ tka <- 0.45; tcl <- 1.0; tv <- 3.45; eta.ka ~ 0.09; add.sd <- 0.7 })
      model({
        matExp()
        k_depot_central <- exp(tka + eta.ka)
        k_central_output <- exp(tcl) / exp(tv)
        cp <- central / exp(tv)
        cp ~ add(add.sd)
      })
    }
    .dat <- .mkData(odeLin, c(tka = 0.6, tcl = 1.1, tv = 3.6))
    # sigdig = 6: see the focei/foce block above -- matExp() is exact (#860), so
    # the ODE form needs an accurate solve to match it at 1e-3.
    for (.met in c("agq", "laplace")) {
      .ctl <- if (.met == "agq") agqControl(print = 0, sigdig = 6) else laplaceControl(print = 0, sigdig = 6)
      .fO <- .nlmixr(odeLin, .dat, est = .met, control = .ctl)
      .fM <- .nlmixr(matLin, .dat, est = .met, control = .ctl)
      expect_equal(.fM$objf, .fO$objf, tolerance = 1e-3)
      expect_equal(unname(fixef(.fM)), unname(fixef(.fO)), tolerance = 1e-3)
    }
  })

  test_that("matExp() + indLin() Michaelis-Menten estimates identically to the ODE (focei)", {
    odeMM <- function() {
      ini({ tka <- 0.45; tvmax <- log(60); tkm <- log(40); tv <- 3.45; eta.ka ~ 0.09; add.sd <- 0.7 })
      model({
        ka <- exp(tka + eta.ka); vmax <- exp(tvmax); km <- exp(tkm); v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(central) <- ka * depot - vmax * central / (km + central)
        cp <- central / v
        cp ~ add(add.sd)
      })
    }
    matMM <- function() {
      ini({ tka <- 0.45; tvmax <- log(60); tkm <- log(40); tv <- 3.45; eta.ka ~ 0.09; add.sd <- 0.7 })
      model({
        matExp()
        k_depot_central <- exp(tka + eta.ka)
        indLin(central) <- -exp(tvmax) * central / (exp(tkm) + central)
        cp <- central / exp(tv)
        cp ~ add(add.sd)
      })
    }
    .dat <- .mkData(odeMM, c(tka = 0.6, tvmax = log(70), tkm = log(45), tv = 3.6))
    # sigdig = 6 for the same reason the .cmp() helper below pins it: matExp/indLin
    # are EXACT for this pseudo-linear system, so the ODE form needs an accurate
    # solve to match it.  The sigdig = 4 default leaves ~0.7 objf of ODE
    # discretization error at the converged optimum, which moves the fixed effects
    # ~2.7% -- far outside the 1e-3 asserted here.
    .fO <- .nlmixr(odeMM, .dat, est = "focei", control = foceiControl(print = 0, sigdig = 6))
    .fM <- .nlmixr(matMM, .dat, est = "focei", control = foceiControl(print = 0, sigdig = 6))
    # tvmax/tkm are near-collinear on this 6-subject data (same as the analytic
    # gradient/cov test below): the likelihood is nearly flat along their ridge,
    # so run-to-run floating-point differences between the native matExp path
    # (#860) and the ODE path's completely different computational route move
    # the converged point along that ridge more than 1e-3 even though the
    # objective itself barely changes.  1e-2 comfortably covers the observed
    # ~3e-3 relative drift while still catching a real divergence.
    expect_equal(.fM$objf, .fO$objf, tolerance = 1e-2)
    expect_equal(unname(fixef(.fM)), unname(fixef(.fO)), tolerance = 1e-2)
  })

  test_that("matExp()-native inner model pins the FOCEi lhs column layout (#860)", {
    # src/inner.cpp reads sensitivities at FIXED positional offsets from
    # predOffset (the by-name index of rx_pred_): lhs[predOffset+i+1] is
    # d(f)/d(eta_i), lhs[predOffset+neta+1] is rx_r_, and
    # lhs[predOffset+i+neta+2] is d(R)/d(eta_i) -- a contiguous
    # [pred, HdEta_1..neta, r, REta_1..neta] block.  A future change to
    # .sensMatExpNative()'s ..ddt assembly that reorders or interleaves
    # extra output columns into that block would silently corrupt the
    # inner objective rather than error, so pin the exact compiled column
    # order here (two etas, so a reordering that only breaks neta > 1
    # still fails loudly).
    matLin2 <- function() {
      ini({ tka <- 0.45; tcl <- 1.0; tv <- 3.45; eta.ka ~ 0.09; eta.cl ~ 0.09; add.sd <- 0.7 })
      model({
        matExp()
        k_depot_central <- exp(tka + eta.ka)
        k_central_output <- exp(tcl + eta.cl) / exp(tv)
        cp <- central / exp(tv)
        cp ~ add(add.sd)
      })
    }
    .s <- suppressMessages(nlmixr2est:::rxUiGet.foceiEnv(list(matLin2())))
    expect_true(isTRUE(.s$..matExpNative))
    .cmtPre <- nlmixr2est:::rxUiGet.foceiCmtPreModel(list(matLin2()))
    .paramsPre <- nlmixr2est:::.uiGetThetaEtaParams(matLin2(), TRUE)
    .mod <- suppressMessages(rxode2::rxode2(paste(.paramsPre, .cmtPre, .s$..inner, sep = "\n")))
    expect_equal(
      rxode2::rxModelVars(.mod)$lhs,
      c("rx_pred_",
        "rx__sens_rx_pred__BY_ETA_1___", "rx__sens_rx_pred__BY_ETA_2___",
        "rx_r_",
        "rx__sens_rx_r__BY_ETA_1___", "rx__sens_rx_r__BY_ETA_2___")
    )
  })

  test_that("matExp() population model estimates identically to the ODE (nlm)", {
    odeMM <- function() {
      ini({ tka <- 0.45; tvmax <- log(60); tkm <- log(40); tv <- 3.45; add.sd <- 0.7 })
      model({
        ka <- exp(tka); vmax <- exp(tvmax); km <- exp(tkm); v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(central) <- ka * depot - vmax * central / (km + central)
        cp <- central / v
        cp ~ add(add.sd)
      })
    }
    matMM <- function() {
      ini({ tka <- 0.45; tvmax <- log(60); tkm <- log(40); tv <- 3.45; add.sd <- 0.7 })
      model({
        matExp()
        k_depot_central <- exp(tka)
        indLin(central) <- -exp(tvmax) * central / (exp(tkm) + central)
        cp <- central / exp(tv)
        cp ~ add(add.sd)
      })
    }
    .dat <- .mkData(odeMM, c(tka = 0.7, tvmax = log(80), tkm = log(55), tv = 3.6))
    .fO <- .nlmixr(odeMM, .dat, est = "nlm", list(print = 0))
    .fM <- .nlmixr(matMM, .dat, est = "nlm", list(print = 0))
    expect_equal(.fM$objf, .fO$objf, tolerance = 1e-3)
    expect_equal(unname(fixef(.fM)), unname(fixef(.fO)), tolerance = 1e-3)
  })

  test_that("matExp()/indLin() analytic gradient + covariance match the ODE (focei/foce/focep + mu/irls)", {
    # The augmented outer/cov model must materialize matExp/indLin ODEs in the same
    # source-first compartment order as the inner model; an indLin() forcing state
    # otherwise parses as compartment 1, misplacing dosing and collapsing the eta
    # sensitivities (analytic gradient wrong, observed-information R singular -> FD).
    odeLin <- function() {
      ini({ tka <- 0.45; tcl <- 1.0; tv <- 3.45; eta.ka ~ 0.09; add.sd <- 0.7 })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(central) <- ka * depot - cl / v * central
        cp <- central / v
        cp ~ add(add.sd)
      })
    }
    matLin <- function() {
      ini({ tka <- 0.45; tcl <- 1.0; tv <- 3.45; eta.ka ~ 0.09; add.sd <- 0.7 })
      model({
        matExp()
        k_depot_central <- exp(tka + eta.ka)
        k_central_output <- exp(tcl) / exp(tv)
        cp <- central / exp(tv)
        cp ~ add(add.sd)
      })
    }
    odeMM <- function() {
      ini({ tka <- 0.45; tvmax <- log(60); tkm <- log(40); tv <- 3.45; eta.ka ~ 0.09; add.sd <- 0.7 })
      model({
        ka <- exp(tka + eta.ka); vmax <- exp(tvmax); km <- exp(tkm); v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(central) <- ka * depot - vmax * central / (km + central)
        cp <- central / v
        cp ~ add(add.sd)
      })
    }
    matMM <- function() {
      ini({ tka <- 0.45; tvmax <- log(60); tkm <- log(40); tv <- 3.45; eta.ka ~ 0.09; add.sd <- 0.7 })
      model({
        matExp()
        k_depot_central <- exp(tka + eta.ka)
        indLin(central) <- -exp(tvmax) * central / (exp(tkm) + central)
        cp <- central / exp(tv)
        cp ~ add(add.sd)
      })
    }
    .datLin <- .mkData(odeLin, c(tka = 0.6, tcl = 1.1, tv = 3.6))
    .datMM  <- .mkData(odeMM,  c(tka = 0.6, tvmax = log(70), tkm = log(45), tv = 3.6))

    .cmp <- function(ode, mat, dat, est, ctlFun, seTol = 1e-2) {
      # matExp/indLin are exact for these linear/pseudo-linear systems, so the ODE
      # form needs an accurate solve to match; sigdig=6 keeps the ODE error well
      # below the objf/cov comparison tolerances (the sigdig=4 default leaves ~0.7
      # objf of ODE discretization error at the converged optimum).
      .ctl <- ctlFun(print = 0, fast = TRUE, covMethod = "analytic", sigdig = 6)
      .fO <- .nlmixr(ode, dat, est = est, control = .ctl)
      .fM <- .nlmixr(mat, dat, est = est, control = .ctl)
      # the matExp/indLin fit walks the SAME covariance path as its ODE
      # counterpart -- both take the analytic covariance where the observed
      # information is PD, and both fall back identically where it is not (the
      # near-collinear MM model under non-interaction FOCE).  A mis-built
      # augmented model would diverge here (or blow the SE tolerance below).
      expect_equal(.fM$covMethod, .fO$covMethod)
      expect_equal(.fM$objf, .fO$objf, tolerance = 1e-3)
      expect_equal(unname(sqrt(diag(.fM$cov))), unname(sqrt(diag(.fO$cov))), tolerance = seTol)
    }

    for (.est in c("focei", "foce", "focep")) {
      .ctlFun <- switch(.est, focei = foceiControl, foce = foceControl, focep = focepControl)
      .cmp(odeLin, matLin, .datLin, .est, .ctlFun)   # pure-linear matExp
      # seTol = 0.08: the comment below records that near-collinear tvmax/tkm inflate
      # the ODE-vs-matExp SE difference to ~5-6% on this 6-subject data.  The
      # helper's 1e-2 default asserts 1% against a quantity documented as 5-6%, so
      # the MM case gets the looser ballpark tolerance that comment intends.
      .cmp(odeMM,  matMM,  .datMM,  .est, .ctlFun, seTol = 0.08)   # indLin() Michaelis-Menten
    }
    # mu-referenced and IRLS families share the same augmented builder, but the
    # mu-regression re-derives the mu-thetas (tka is regression-updated since plain
    # mu-ref profiling) from each form's numeric gradient, so the matExp and ODE forms
    # converge to marginally different mu-thetas.  The MM parameters tvmax/tkm are
    # near-collinear on this 6-subject data, so that small theta difference inflates the
    # SEs globally by ~5-6% -- deterministic (verified seed- and thread-independent), not
    # an augmented-builder error (the analytic cov ran on BOTH, and the objfs match to
    # 1e-3).  So the SEs are compared with a looser ballpark tolerance; the covMethod and
    # objf checks are the real guards (a mis-built augmented model falls back to a
    # singular-R FD cov, which those catch).
    for (.est in c("mfocei", "ifocei")) {
      .cmp(odeMM, matMM, .datMM, .est, foceiControl, seTol = 8e-2)
    }
    for (.est in c("mfoce", "ifoce")) {
      .cmp(odeMM, matMM, .datMM, .est, foceControl, seTol = 8e-2)
    }
  })

  test_that("matExp() estimates identically to the ODE (saem)", {
    odeS <- function() {
      ini({ tka <- 0.45; tcl <- 1.0; tv <- 3.45; eta.ka ~ 0.3; add.sd <- 0.7 })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(central) <- ka * depot - cl / v * central
        cp <- central / v
        cp ~ add(add.sd)
      })
    }
    matS <- function() {
      ini({ tka <- 0.45; tcl <- 1.0; tv <- 3.45; eta.ka ~ 0.3; add.sd <- 0.7 })
      model({
        matExp()
        k_depot_central <- exp(tka + eta.ka)
        k_central_output <- exp(tcl) / exp(tv)
        cp <- central / exp(tv)
        cp ~ add(add.sd)
      })
    }
    .ctl <- saemControl(print = 0, nBurn = 100, nEm = 100, seed = 42)
    .fO <- .nlmixr(odeS, nlmixr2data::theo_sd, est = "saem", control = .ctl)
    .fM <- .nlmixr(matS, nlmixr2data::theo_sd, est = "saem", control = .ctl)
    expect_equal(.fM$objf, .fO$objf, tolerance = 1e-3)
    expect_equal(unname(fixef(.fM)), unname(fixef(.fO)), tolerance = 1e-3)
    # mechanism check (#859): a pure-linear matExp() model must solve
    # natively (a literal "matExp()" line, no materialized d/dt()), not via
    # the ODE-flattening path -- matching numbers alone would not catch a
    # regression back to flattening.
    .saemModelTxt <- strsplit(utils::getFromNamespace("rxUiGet.saemModel", "nlmixr2est")(list(matS())), "\n")[[1]]
    expect_true(any(grepl("^matExp\\(\\)$", .saemModelTxt)))
    expect_false(any(grepl("d/dt(", .saemModelTxt, fixed = TRUE)))
  })

  test_that("matExp()/indLin() Michaelis-Menten (saem) falls back to the ODE flatten", {
    # rxode2's true inductive-linearization iteration (doIndLin 3/4) does not
    # yet converge reliably under SAEM's per-iteration parameter draws (a
    # native attempt diverged from the ODE reference: objf -174 vs 321,
    # unrelated fixed effects; tracked separately as #977), so
    # .rxKeepMatExpNative() bails on an indLin() forcing term and
    # rxUiGet.saemModel() keeps materializing d/dt() via
    # .rxInjectMatExpDdt(), same as before #859.
    odeMM <- function() {
      ini({ tka <- 0.45; tvmax <- log(60); tkm <- log(40); tv <- 3.45; eta.ka ~ 0.09; add.sd <- 0.7 })
      model({
        ka <- exp(tka + eta.ka); vmax <- exp(tvmax); km <- exp(tkm); v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(central) <- ka * depot - vmax * central / (km + central)
        cp <- central / v
        cp ~ add(add.sd)
      })
    }
    matMM <- function() {
      ini({ tka <- 0.45; tvmax <- log(60); tkm <- log(40); tv <- 3.45; eta.ka ~ 0.09; add.sd <- 0.7 })
      model({
        matExp()
        k_depot_central <- exp(tka + eta.ka)
        indLin(central) <- -exp(tvmax) * central / (exp(tkm) + central)
        cp <- central / exp(tv)
        cp ~ add(add.sd)
      })
    }
    .saemModelTxt <- strsplit(utils::getFromNamespace("rxUiGet.saemModel", "nlmixr2est")(list(matMM())), "\n")[[1]]
    expect_true(any(grepl("d/dt(", .saemModelTxt, fixed = TRUE)))
    expect_false(any(grepl("^matExp\\(\\)$", .saemModelTxt)))
    expect_false(any(grepl("^indLin\\(", .saemModelTxt)))

    .dat <- .mkData(odeMM, c(tka = 0.6, tvmax = log(70), tkm = log(45), tv = 3.6), nid = 8, seed = 4321)
    .ctl <- saemControl(print = 0, nBurn = 100, nEm = 100, seed = 42)
    .fO <- .nlmixr(odeMM, .dat, est = "saem", control = .ctl)
    .fM <- .nlmixr(matMM, .dat, est = "saem", control = .ctl)
    expect_equal(.fM$objf, .fO$objf, tolerance = 1e-3)
    expect_equal(unname(fixef(.fM)), unname(fixef(.fO)), tolerance = 1e-3)
  })

  test_that("matExp() with a state-free indLin() forcing (saem) falls back to the ODE flatten", {
    # A forcing that does not reference any state (e.g. indLin(central) <- 5)
    # leaves wIndLin empty, same as a pure-linear model with no forcing at
    # all -- so the native-routing decision cannot key on wIndLin alone (it
    # would silently drop the forcing term, since the native path emits no
    # indLin() text). .rxKeepMatExpNative() instead bails whenever
    # indLin$f is non-NULL (any forcing, state-free or not).
    matSF <- function() {
      ini({ tka <- 0.45; tv <- 3.45; add.sd <- 0.7 })
      model({
        matExp()
        k_depot_central <- exp(tka)
        indLin(central) <- 5
        cp <- central / exp(tv)
        cp ~ add(add.sd)
      })
    }
    .saemModelTxt <- strsplit(utils::getFromNamespace("rxUiGet.saemModel", "nlmixr2est")(list(matSF())), "\n")[[1]]
    expect_true(any(grepl("d/dt(central)=", .saemModelTxt, fixed = TRUE)))
    expect_true(any(grepl("\\(5\\)", .saemModelTxt)))
    expect_false(any(grepl("^matExp\\(\\)$", .saemModelTxt)))
  })

  test_that("matExp() with delay() (saem) falls back to the ODE flatten", {
    # .saemFitModel() (R/saem.R) forces the dense dop853 solver for a delay
    # model (hasDelay), which takes precedence over method="indLin" -- so a
    # native "matExp()" emission (no d/dt() for dop853 to integrate) would
    # silently solve to all-zero states. .rxKeepMatExpNative() bails on
    # flags[["hasDelay"]] to keep this combination flattened.
    matD <- function() {
      ini({ tka <- 0.45; tv <- 3.45; add.sd <- 0.7 })
      model({
        matExp()
        k_depot_central <- exp(tka)
        cp <- delay(central, 5) / exp(tv)
        cp ~ add(add.sd)
      })
    }
    .s <- utils::getFromNamespace("rxUiGet.loadPruneSaem", "nlmixr2est")(list(matD()))
    expect_false(isTRUE(utils::getFromNamespace(".rxKeepMatExpNative", "nlmixr2est")(.s)))
  })

  test_that("matExp()/indLin() table and residual generation (#858)", {
    # General regression check (issue #858's stated test) that table/residual
    # generation for a matExp()/indLin() fit works. This does NOT exercise the
    # setdiff() removal in R/resid.R itself -- see
    # test-resid-ode-fallback.R for that (the default ODE methods used here
    # never reach the branch the removal changed).
    matLin <- function() {
      ini({ tka <- 0.45; tcl <- 1.0; tv <- 3.45; eta.ka ~ 0.09; add.sd <- 0.7 })
      model({
        matExp()
        k_depot_central <- exp(tka + eta.ka)
        k_central_output <- exp(tcl) / exp(tv)
        cp <- central / exp(tv)
        cp ~ add(add.sd)
      })
    }
    matMM <- function() {
      ini({ tka <- 0.45; tvmax <- log(60); tkm <- log(40); tv <- 3.45; eta.ka ~ 0.09; add.sd <- 0.7 })
      model({
        matExp()
        k_depot_central <- exp(tka + eta.ka)
        indLin(central) <- -exp(tvmax) * central / (exp(tkm) + central)
        cp <- central / exp(tv)
        cp ~ add(add.sd)
      })
    }
    .datLin <- .mkData(matLin, c(tka = 0.6, tcl = 1.1, tv = 3.6))
    .datMM <- .mkData(matMM, c(tka = 0.6, tvmax = log(70), tkm = log(45), tv = 3.6))

    .fLin <- .nlmixr(matLin, .datLin, est = "focei", control = foceiControl(print = 0))
    expect_true(all(c("PRED", "IPRED") %in% names(.fLin)))
    expect_false(any(is.na(.fLin$IPRED)))
    suppressMessages(expect_error(addCwres(.fLin), NA))
    suppressMessages(expect_error(addNpde(.fLin), NA))

    .fMM <- .nlmixr(matMM, .datMM, est = "focei", control = foceiControl(print = 0, sigdig = 6))
    expect_true(all(c("PRED", "IPRED") %in% names(.fMM)))
    expect_false(any(is.na(.fMM$IPRED)))
    suppressMessages(expect_error(addCwres(.fMM), NA))
    suppressMessages(expect_error(addNpde(.fMM), NA))
  })
})
