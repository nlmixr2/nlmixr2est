nmTest({
  # f-SAEM proposal core: reuse the FOCEi inner to get the per-subject conditional
  # MAP and the proposal precision H = Gamma_i^-1.  Validate that (a) the returned
  # eta is a stationary point of the individual objective and (b) H matches the
  # independent paper Eq-17 information J' J / sigma^2 + Omega^-1.
  test_that("fsaemInnerMap: MAP + Jacobian (Eq 17) proposal covariance", {
    skip_if_not_installed("numDeriv")
    skip_if_not_installed("nlmixr2data")

    one.cmt <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v  <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    .data <- nlmixr2data::theo_sd
    .ui <- rxode2::rxUiDecompress(rxode2::rxode2(one.cmt))
    .neta <- 3L
    .N <- length(unique(.data$ID))

    .ctl <- list(rxControl = rxode2::rxControl(), fastCov = "jacobian", fastLik = "focei",
                 fastInnerIt = 100L, sumProd = FALSE, optExpression = TRUE, literalFix = FALSE,
                 addProp = "combined2", eventSens = "jump", indTolRelax = TRUE,
                 maxOdeRecalc = 5L, odeRecalcFactor = 10^0.5)

    .env <- .fsaemInnerSetup(.ui, .data, matrix(0, .N, .neta), .ctl)
    on.exit(.fsaemInnerFree(), add = TRUE)
    .map <- .fsaemInnerMap(.ctl, .neta)

    expect_equal(nrow(.map$eta), .N)
    expect_true(all(.map$ok == 1L))

    .sigma2 <- 0.7^2
    .OmegaInv <- solve(diag(c(0.6, 0.3, 0.1)))

    # prediction f_i(eta) via the inner driver (perturb only subject i's row)
    .predOf <- function(i, e, base) {
      .M <- base; .M[i, ] <- e
      as.numeric(vaeInnerLik(.M, 1L, FALSE, TRUE)$f[[i]])
    }
    # individual objective for subject i (for the MAP stationarity check)
    .objOf <- function(i, e) {
      .M <- .map$eta; .M[i, ] <- e
      vaeInnerLik(.M, 1L, FALSE, FALSE)$obj[i]
    }

    .maxRel <- 0
    .maxGrad <- 0
    for (i in seq_len(.N)) {
      .J <- numDeriv::jacobian(function(e) .predOf(i, e, .map$eta), .map$eta[i, ])
      .Hindep <- t(.J) %*% .J / .sigma2 + .OmegaInv
      .Hinner <- .map$hess[[i]]
      .maxRel <- max(.maxRel, max(abs(.Hinner - .Hindep) / pmax(abs(.Hindep), 1)))
      .maxGrad <- max(.maxGrad, max(abs(numDeriv::grad(function(e) .objOf(i, e), .map$eta[i, ]))))
    }
    # MAP is a stationary point of the individual objective
    expect_lt(.maxGrad, 1e-3)
    # proposal precision matches paper Eq 17 (relative; entries reach ~700)
    expect_lt(.maxRel, 1e-3)

    # proposal covariance is the inverse of the precision and is positive definite
    for (i in seq_len(.N)) {
      expect_true(all(is.finite(.map$gamma[[i]])))
      expect_true(all(eigen(.map$gamma[[i]], only.values = TRUE)$values > 0))
    }
  })

  test_that("fsaemInnerUpdate bridge: re-parameterizing theta/omega moves the MAP", {
    skip_if_not_installed("nlmixr2data")
    # The SAEM loop feeds the inner its current estimate as theta = [structural
    # fixed effects, residual] (THETA order) and omega = diag(Gamma2).  Verify
    # that update path re-parameterizes the inner correctly: a tighter prior must
    # shrink every subject's MAP toward 0.
    one.cmt <- function() {
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
    .data <- nlmixr2data::theo_sd
    .ui <- rxode2::rxUiDecompress(rxode2::rxode2(one.cmt))
    .N <- length(unique(.data$ID))
    .ctl <- list(rxControl = rxode2::rxControl(), fastCov = "jacobian", fastLik = "focei",
                 fastInnerIt = 100L, sumProd = FALSE, optExpression = TRUE, literalFix = FALSE,
                 addProp = "combined2", eventSens = "jump", indTolRelax = TRUE,
                 maxOdeRecalc = 5L, odeRecalcFactor = 10^0.5)
    .env <- .fsaemInnerSetup(.ui, .data, matrix(0, .N, 3L), .ctl)
    on.exit(.fsaemInnerFree(), add = TRUE)
    # thetaIni layout is [tka, tcl, tv, add.sd]
    expect_equal(unname(.env$thetaIni), c(0.45, 1, 3.45, 0.7))

    .wide <- .fsaemInnerMap(.ctl, 3L)$eta
    .fsaemInnerUpdate(.env, theta = c(0.45, 1, 3.45, 0.7), omega = c(0.02, 0.02, 0.02),
                      matrix(0, .N, 3L))
    .tight <- .fsaemInnerMap(.ctl, 3L)$eta
    # tighter prior -> MAP shrinks toward 0 in aggregate (the update re-solved the
    # inner at the new omega; individual components balance prior vs likelihood)
    expect_lt(sum(abs(.tight)), sum(abs(.wide)))
    expect_gt(mean(abs(.tight) <= abs(.wide) + 1e-8), 0.6)
  })

  # issue #874: the SAEM E-step's residual SD is unconditionally combined1
  # (ares + bres*|f|), so the f-SAEM inner has to be built combined1 too or the
  # IMH proposal is preconditioned against a variance the chain never uses.  A
  # wrong proposal costs acceptance, not the estimate, so pin the inner itself.
  test_that("fsaem inner is combined1 whatever addProp the fit asks for", {
    skip_if_not_installed("nlmixr2data")

    comb <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
        add.sd <- 0.3; prop.sd <- 0.1
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd) + prop(prop.sd)
      })
    }
    .data <- nlmixr2data::theo_sd
    .ui <- rxode2::rxUiDecompress(rxode2::rxode2(comb))
    # combined2 is saemControl()'s default, i.e. what the inner used to inherit
    .ui$control <- saemControl(addProp = "combined2", print = 0L, calcTables = FALSE)
    .neta <- 3L
    .N <- length(unique(.data$ID))

    .cfg <- .fsaemInstallStep(.ui, .data, rxode2::rxControl(), list())
    expect_false(is.null(.cfg$fsaemStep))
    # the mechanism: the inner's own foceiControl, not the fit's addProp
    expect_equal(.cfg$fsaemInnerEnv$control$addProp, "combined1")
    .hessInstalled <- .fsaemInnerMap(list(rxControl = rxode2::rxControl()), .neta)$hess
    .fsaemInnerFree()

    # ...and the proposal precision it produces is the one a combined1 inner
    # gives, not a combined2 one.  Compared against the inner built each way
    # rather than an analytic Eq-17 -- the control has to survive the focei
    # model cache, which keys on it, and reach the residual the MAP is scored
    # against.
    .innerHess <- function(addProp) {
      .ctl <- list(rxControl = rxode2::rxControl(), fastCov = "jacobian", fastLik = "focei",
                   fastInnerIt = 100L, sumProd = FALSE, optExpression = TRUE, literalFix = FALSE,
                   addProp = addProp, eventSens = "jump", indTolRelax = TRUE,
                   maxOdeRecalc = 5L, odeRecalcFactor = 10^0.5)
      .fsaemInnerSetup(.ui, .data, matrix(0, .N, .neta), .ctl)
      on.exit(.fsaemInnerFree(), add = TRUE)
      .fsaemInnerMap(.ctl, .neta)$hess
    }
    .h1 <- .innerHess("combined1")
    .h2 <- .innerHess("combined2")
    expect_equal(.hessInstalled, .h1)
    # the two residual forms really are distinguishable here (guards the guard)
    expect_false(isTRUE(all.equal(.h1, .h2)))
  })

  test_that("a model DECLARING combined2 degrades instead of mis-targeting", {
    # the control cannot override a model-level declaration, so such an endpoint
    # falls back to standard SAEM rather than building a combined2 proposal
    default <- function() {
      ini({tka <- 0.45; tcl <- 1; tv <- 3.45; eta.ka ~ 0.6; add.sd <- 0.3; prop.sd <- 0.1})
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv)
        linCmt() ~ add(add.sd) + prop(prop.sd)
      })
    }
    declared1 <- function() {
      ini({tka <- 0.45; tcl <- 1; tv <- 3.45; eta.ka ~ 0.6; add.sd <- 0.3; prop.sd <- 0.1})
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv)
        linCmt() ~ add(add.sd) + prop(prop.sd) + combined1()
      })
    }
    declared2 <- function() {
      ini({tka <- 0.45; tcl <- 1; tv <- 3.45; eta.ka ~ 0.6; add.sd <- 0.3; prop.sd <- 0.1})
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv)
        linCmt() ~ add(add.sd) + prop(prop.sd) + combined2()
      })
    }
    .sup <- function(m) .fsaemSupported(rxode2::rxUiDecompress(rxode2::rxode2(m)))
    expect_true(.sup(default))     # default -> the control makes it combined1
    expect_true(.sup(declared1))   # declared, and it agrees with the E-step
    expect_false(.sup(declared2))  # declared, and it does not

    # a lone add() endpoint is unaffected: the two forms coincide there
    addOnly <- function() {
      ini({tka <- 0.45; tcl <- 1; tv <- 3.45; eta.ka ~ 0.6; add.sd <- 0.3})
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv)
        linCmt() ~ add(add.sd) + combined2()
      })
    }
    expect_true(.fsaemSupported(rxode2::rxUiDecompress(rxode2::rxode2(addOnly))))
  })
})
