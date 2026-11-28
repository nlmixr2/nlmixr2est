nmTest({
  # The event-sensitivity ("jump") shape is a PROCESS GLOBAL: rxode2 holds exactly
  # one model's shape installed at a time.  The SAEM prediction model carries no
  # jump sensitivities and must never be solved while another model's shape is
  # live, which is why saem_fit / saem_do_pred / fsaemRestoreSaemSolve all call
  # odeSwapEsOff().
  #
  # A fit-level test cannot pin that guard, and the fsaem suite's existing
  # modeled-lag fit does not.  Measured on the same model: removing all four
  # halves of the guard (the three odeSwapEsOff() calls and
  # OdeSwapEsBatch(odeSlotInner)) leaves an fsaem fit BIT-IDENTICAL, because a
  # fit's FIRST load of the SAEM model rebinds rxode2's event globals itself --
  # a shape installed beforehand is already gone by the time saem_fit runs
  # (measured: esOffN moves, esDroppedN stays 0).
  #
  # What the guard really protects is RE-ENTRY into an already-loaded SAEM
  # solve while another component's shape is live -- focei's fit-wide
  # rxEventSensLoadModel, or the nlm theta-gradient path.  There rxLoad is a
  # no-op, so nothing else takes the shape down.  The test below constructs
  # exactly that, and it is falsifiable: with odeSwapEsOff() removed from
  # saem_do_pred, esLive stays 1 and esDroppedN does not move.
  .mkInnerCtl <- function(cov = "jacobian") {
    list(rxControl = rxode2::rxControl(), fastCov = cov, fastLik = "focei",
         fastInnerIt = 100L, sumProd = FALSE, optExpression = TRUE, literalFix = FALSE,
         addProp = "combined2", eventSens = "jump", indTolRelax = TRUE,
         maxOdeRecalc = 5L, odeRecalcFactor = 10^0.5)
  }
  # modeled dosing (lag) is the ONLY thing that gives a FOCEi inner event
  # sensitivities; a bolus linCmt() model has eventSensInfo NULL
  .lagm <- function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45; tlag <- -0.7
      eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
      lag(depot) <- exp(tlag)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }
  .bolus <- function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45
      eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }

  test_that("only a modeled-dosing model gives the fsaem inner a jump shape", {
    skip_if_not_installed("nlmixr2data")
    # The premise every other assertion in this file rests on, and the reason the
    # rest of the fsaem suite is blind to this path: it all uses bolus models.
    .data <- nlmixr2data::theo_sd
    .N <- length(unique(.data$ID))
    .esOf <- function(fn) {
      .ui <- rxode2::rxUiDecompress(rxode2::rxode2(fn))
      .fsaemInnerSetup(.ui, .data, matrix(0, .N, 3L), .mkInnerCtl())
      on.exit(.fsaemInnerFree(), add = TRUE)
      .m <- .odeSwapInfo()$models
      .m$esActive[match("inner", .m$name)]
    }
    expect_equal(.esOf(.lagm), 1L)
    expect_equal(.esOf(.bolus), 0L)
  })

  test_that("saem_do_pred drops a jump shape it inherited", {
    skip_if_not_installed("nlmixr2data")
    .data <- nlmixr2data::theo_sd
    .N <- length(unique(.data$ID))
    .ctl <- saemControl(nBurn = 20, nEm = 5, nmc = 3, seed = 11, print = 0L,
                        calcTables = FALSE)
    .fs <- suppressMessages(nlmixr2(.lagm, .data, est = "fsaem", control = .ctl))
    # the fast kernel really ran, and the fit left nothing installed
    expect_gt(.fs$fsaemDiag$nStep, 0)
    expect_equal(.odeSwapInfo()$esLive, 0L)

    # fit$saemDopredIpred re-enters saem_do_pred at the fit's own phi, so it is
    # an exact before/after of the guarded SAEM solve -- no stochastic tolerance.
    .ref <- .fs$saemDopredIpred

    # construct the hazard: install the FOCEi inner's jump shape, as focei.R does
    # fit-wide, and go back through the SAEM prediction entry
    .ui <- rxode2::rxUiDecompress(rxode2::rxode2(.lagm))
    .env <- .fsaemInnerSetup(.ui, .data, matrix(0, .N, 3L), .mkInnerCtl())
    on.exit(.fsaemInnerFree(), add = TRUE)
    expect_true(rxode2::rxEventSensLoadModel(.env$model$inner))
    .b <- .odeSwapInfo()
    expect_equal(.b$esLive, 1L)                       # premise: a shape IS live
    expect_gt(.b$esLiveNstate, 0L)

    .got <- .fs$saemDopredIpred
    .a <- .odeSwapInfo()
    expect_equal(.a$esLive, 0L)                       # the guard took it down
    expect_equal(.a$esDroppedN - .b$esDroppedN, 1)    # exactly once, and it counted
    expect_equal(.got, .ref)                          # and the solve is unchanged
  })

  test_that("an fsaem fit never leaves a jump shape installed", {
    skip_if_not_installed("nlmixr2data")
    # Entering a fit with a shape live is the OTHER half of the invariant.  The
    # first load of the SAEM model rebinds the globals, so this does not exercise
    # saem_fit's odeSwapEsOff() (esDroppedN stays 0 here, deliberately asserted);
    # what it pins is that the whole fsaem path -- inner setup, the MAP+IMH
    # batch, fsaemRestoreSaemSolve -- is shape-neutral end to end and leaves
    # nothing behind for the next fit in the session to inherit.
    .data <- nlmixr2data::theo_sd
    .N <- length(unique(.data$ID))
    .ui <- rxode2::rxUiDecompress(rxode2::rxode2(.lagm))
    .env <- .fsaemInnerSetup(.ui, .data, matrix(0, .N, 3L), .mkInnerCtl())
    on.exit(.fsaemInnerFree(), add = TRUE)
    expect_true(rxode2::rxEventSensLoadModel(.env$model$inner))
    .b <- .odeSwapInfo()
    expect_equal(.b$esLive, 1L)

    .ctl <- saemControl(nBurn = 20, nEm = 5, nmc = 3, seed = 11, print = 0L,
                        calcTables = FALSE)
    .fs <- suppressMessages(nlmixr2(.lagm, .data, est = "fsaem", control = .ctl))
    .a <- .odeSwapInfo()
    expect_gt(.fs$fsaemDiag$nStep, 0)
    expect_gt(.a$esOffN - .b$esOffN, 0)   # the guard ran on the fit path
    expect_equal(.a$esLive, 0L)           # and nothing is left installed
  })
})
