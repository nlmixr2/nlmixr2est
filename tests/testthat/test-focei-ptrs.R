# FOCEi conditional-likelihood C API (#937): the _nlmixr2est_foceiPtrs table
# and its entry points, exercised through the R shims (an external pointer
# cannot be invoked from R).  Same small analytic model as test-foceiLik.R so
# it compiles fast and the conditional density is hand-computable.

.foceiPtrMod <- function() {
  ini({
    tcl <- 1
    tv <- 3
    add.sd <- 0.5
    eta.cl ~ 0.1
  })
  model({
    cl <- exp(tcl + eta.cl)
    v <- exp(tv)
    cp <- 100 / v * exp(-cl / v * time)
    cp ~ add(add.sd)
  })
}

.foceiPtrData <- function() {
  .testSeed(42)
  do.call(rbind, lapply(1:4, function(id) {
    tt <- c(0.5, 1, 2, 4, 8)
    data.frame(ID = id, TIME = tt,
               DV = 5 * exp(-0.05 * tt) + stats::rnorm(length(tt), 0, 0.5),
               AMT = 0, EVID = 0)
  }))
}

test_that("the foceiPtrs table has the documented shape (#937 + #955)", {
  .p <- .nlmixr2estFoceiPtrs()
  expect_length(.p, 10L)
  expect_equal(names(.p), c("apiVersion", "dims", "setTheta", "condBatch",
                            "setOmegaInv", "thetaSensIdx", "condThetaGrad",
                            "nMix", "iterPrintRow", "condBatchThetaGrad"))
  for (.i in seq_along(.p)) expect_true(inherits(.p[[.i]], "externalptr"))
})

test_that("entries report 'not loaded' by return code, not error (#937)", {
  skip_on_cran()
  foceiLikUnload()  # no-op if nothing is loaded
  expect_equal(foceiLikDims_()$status, -1L)
  expect_equal(foceiLikDims_()$apiVersion, 1L)
  expect_equal(foceiLikSetThetaC_(c(1, 2, 3, 4)), -1L)
  expect_equal(foceiLikSetOmegaInvC_(diag(1)), -1L)
  expect_error(foceiLikCondGrad_(matrix(0, 4, 1), 1L), "status -1")
  expect_error(foceiLikThetaSensIdxC_(), "no general likelihood system loaded")
})

test_that("dims + flags describe the loaded problem (#937)", {
  skip_on_cran()
  d <- .foceiPtrData()
  h <- foceiLikLoad(.foceiPtrMod, d, "focei")
  on.exit(foceiLikUnload(), add = TRUE)
  .d <- foceiLikDims_()
  expect_equal(.d$status, 0L)
  expect_equal(.d$nid, 4L)
  expect_equal(.d$neta, 1L)
  expect_equal(.d$ntheta, 3L)
  expect_equal(.d$npars, 4L)
  # focei => interaction bit on; not focep/fo/mixture
  expect_equal(bitwAnd(.d$flags, 0x01), 0x01)
  expect_equal(bitwAnd(.d$flags, 0x02), 0L)
  expect_equal(bitwAnd(.d$flags, 0x04), 0L)
  expect_equal(bitwAnd(.d$flags, 0x10), 0L)
  # no theta-sensitivity model on a plain load
  expect_equal(bitwAnd(.d$flags, 0x40), 0L)
  expect_equal(foceiLikThetaSensIdxC_(), integer(0))
  foceiLikUnload()
  # the focep hazard flag is set so a gradient-based caller can refuse
  h2 <- foceiLikLoad(.foceiPtrMod, d, "focep")
  expect_equal(bitwAnd(foceiLikDims_()$flags, 0x02), 0x02)
})

test_that("condBatch matches foceiLikRun(type='cond') exactly (#937)", {
  skip_on_cran()
  d <- .foceiPtrData()
  h <- foceiLikLoad(.foceiPtrMod, d, "focei")
  on.exit(foceiLikUnload(), add = TRUE)
  .testSeed(7)
  eta <- matrix(stats::rnorm(h$nid * h$neta, 0, 0.2), h$nid, h$neta)
  ref <- foceiLikRun(h$initPar, eta, type = "cond")   # also sets theta
  got <- foceiLikCondGrad_(eta, 1L)
  expect_equal(got$nBad, 0L)
  expect_equal(as.numeric(got$value), as.numeric(ref), tolerance = 1e-12)
})

test_that("condBatch gradient matches central differences of the value (#937)", {
  # This is the test that catches the sign convention in the NEW assembly
  # step (existing focei consumers use lp against likInner0, the objective it
  # is the gradient of, so they are unaffected): fInd->lp is stored as
  # -(dlogp/deta) + Omega^-1 eta, so the conditional gradient must be
  # Omega^-1 eta - lp -- PLUS, not minus.  The flipped assembly is wrong by
  # 2*Omega^-1*eta and raises no error anywhere: an MH-corrected sampler
  # degrades into divergences/zero ESS, and gradient-only consumers
  # (optimization, VI, autodiff composition) converge to wrong answers.
  skip_on_cran()
  d <- .foceiPtrData()
  h <- foceiLikLoad(.foceiPtrMod, d, "focei")
  on.exit(foceiLikUnload(), add = TRUE)
  .testSeed(11)
  eta <- matrix(stats::rnorm(h$nid * h$neta, 0, 0.25), h$nid, h$neta)
  expect_equal(foceiLikSetThetaC_(h$initPar), 0L)
  got <- foceiLikCondGrad_(eta, 1L)
  .h <- 1e-5
  fd <- matrix(0, h$nid, h$neta)
  for (k in seq_len(h$neta)) {
    up <- eta; up[, k] <- up[, k] + .h
    dn <- eta; dn[, k] <- dn[, k] - .h
    vUp <- foceiLikCondGrad_(up, 1L)$value
    vDn <- foceiLikCondGrad_(dn, 1L)$value
    fd[, k] <- (vUp - vDn) / (2 * .h)
  }
  expect_equal(as.numeric(got$grad), as.numeric(fd), tolerance = 1e-4)
  # and the wrong sign would NOT pass: the gradient is not symmetric in eta
  expect_false(isTRUE(all.equal(as.numeric(got$grad), as.numeric(-fd),
                                tolerance = 1e-2)))
})

test_that("condBatch is deterministic and history-independent (#937)", {
  skip_on_cran()
  d <- .foceiPtrData()
  h <- foceiLikLoad(.foceiPtrMod, d, "focei")
  on.exit(foceiLikUnload(), add = TRUE)
  expect_equal(foceiLikSetThetaC_(h$initPar), 0L)
  etaA <- matrix(c(-0.2, 0.05, 0.3, -0.1), h$nid, h$neta)
  etaB <- matrix(c(0.4, -0.3, 0.1, 0.2), h$nid, h$neta)
  r1 <- foceiLikCondGrad_(etaA, 1L)
  # interleave a different point, then return: bitwise identical
  invisible(foceiLikCondGrad_(etaB, 1L))
  r2 <- foceiLikCondGrad_(etaA, 1L)
  expect_identical(r1$value, r2$value)
  expect_identical(r1$grad, r2$grad)
  # and thread-count invariant
  r4 <- foceiLikCondGrad_(etaA, 4L)
  expect_equal(r1$value, r4$value, tolerance = 1e-12)
  expect_equal(r1$grad, r4$grad, tolerance = 1e-12)
})

test_that("the conditional value is Omega-free and its gradient Omega-invariant (#937)", {
  # setOmegaInv exists for numerical conditioning only: the conditional
  # log-likelihood does not depend on Omega, and because the gradient is
  # assembled as Omega^-1 eta - lp with the SAME Omega^-1 likInner0 used,
  # the assembled gradient must not change either.
  skip_on_cran()
  d <- .foceiPtrData()
  h <- foceiLikLoad(.foceiPtrMod, d, "focei")
  on.exit(foceiLikUnload(), add = TRUE)
  expect_equal(foceiLikSetThetaC_(h$initPar), 0L)
  eta <- matrix(c(-0.2, 0.05, 0.3, -0.1), h$nid, h$neta)
  r0 <- foceiLikCondGrad_(eta, 1L)
  for (.c in c(0.01, 1, 100)) {
    expect_equal(foceiLikSetOmegaInvC_(diag(.c, h$neta)), 0L)
    r <- foceiLikCondGrad_(eta, 1L)
    expect_equal(r$value, r0$value, tolerance = 1e-10)
    expect_equal(r$grad, r0$grad, tolerance = 1e-7)
  }
  # bad inputs by return code
  expect_equal(foceiLikSetOmegaInvC_(diag(1, h$neta + 1L)), -2L)
  expect_equal(foceiLikSetOmegaInvC_(matrix(-1, 1, 1)), -3L)
})

test_that("setTheta validates and applies by return code (#937)", {
  skip_on_cran()
  d <- .foceiPtrData()
  h <- foceiLikLoad(.foceiPtrMod, d, "focei")
  on.exit(foceiLikUnload(), add = TRUE)
  eta0 <- matrix(0, h$nid, h$neta)
  expect_equal(foceiLikSetThetaC_(h$initPar[-1]), -2L)
  expect_equal(foceiLikSetThetaC_(h$initPar), 0L)
  v0 <- foceiLikCondGrad_(eta0, 1L)$value
  th2 <- h$initPar
  th2[3] <- th2[3] * 2
  expect_equal(foceiLikSetThetaC_(th2), 0L)
  v2 <- foceiLikCondGrad_(eta0, 1L)$value
  expect_false(isTRUE(all.equal(v0, v2)))
  # round-trip: the theta really is re-applied, not memoized
  expect_equal(foceiLikSetThetaC_(h$initPar), 0L)
  expect_equal(foceiLikCondGrad_(eta0, 1L)$value, v0, tolerance = 1e-12)
})

test_that("condThetaGrad requires the wired sensitivity model (#937)", {
  skip_on_cran()
  d <- .foceiPtrData()
  h <- foceiLikLoad(.foceiPtrMod, d, "focei")
  on.exit(foceiLikUnload(), add = TRUE)
  # plain load: flag 0x40 clear -> -4 by contract
  expect_equal(bitwAnd(foceiLikDims_()$flags, 0x40), 0L)
  expect_error(foceiLikCondThetaGrad_(matrix(0, h$nid, h$neta), 1L),
               "status -4")
})

test_that("condThetaGrad matches central differences when wired (#937 + #939)", {
  skip_on_cran()
  # foceiLikLoad() grows thetaSens=/scale= in #939; until that merges this
  # functional check cannot wire the sensitivity model and skips
  skip_if_not(all(c("thetaSens", "scale") %in% names(formals(foceiLikLoad))),
              "foceiLikLoad() without thetaSens/scale (#939 not merged)")
  d <- .foceiPtrData()
  # scale="natural" so the theta the FD perturbs is the same natural-scale
  # theta impThetaScore differentiates (its forward sensitivities are w.r.t.
  # the model's THETA directly)
  h <- foceiLikLoad(.foceiPtrMod, d, "focei", scale = "natural",
                    thetaSens = TRUE)
  on.exit(foceiLikUnload(), add = TRUE)
  skip_if_not(isTRUE(h$thetaSens), "theta-sensitivity model not built")
  expect_equal(bitwAnd(foceiLikDims_()$flags, 0x40), 0x40)
  # tcl is mu-referenced -> its column stays 0; tv and add.sd carry scores
  sensIdx <- foceiLikThetaSensIdxC_() + 1L
  expect_equal(sensIdx, h$thetaSensIdx)
  .testSeed(5)
  eta <- matrix(stats::rnorm(h$nid * h$neta, 0, 0.2), h$nid, h$neta)
  th <- h$initPar
  expect_equal(foceiLikSetThetaC_(th), 0L)
  got <- foceiLikCondThetaGrad_(eta, 1L)
  expect_equal(got$nBad, 0L)
  .h <- 1e-5
  for (t in sensIdx) {
    up <- th; up[t] <- up[t] + .h
    dn <- th; dn[t] <- dn[t] - .h
    expect_equal(foceiLikSetThetaC_(up), 0L)
    vUp <- foceiLikCondGrad_(eta, 1L)$value
    expect_equal(foceiLikSetThetaC_(dn), 0L)
    vDn <- foceiLikCondGrad_(eta, 1L)$value
    fd <- (vUp - vDn) / (2 * .h)
    expect_equal(as.numeric(got$dTheta[, t]), as.numeric(fd), tolerance = 1e-3)
  }
  # mu-referenced theta columns are the caller's (zero here)
  muCols <- setdiff(seq_len(foceiLikDims_()$ntheta), sensIdx)
  for (t in muCols) expect_true(all(got$dTheta[, t] == 0))
})

test_that("condBatch value/gradient FD-agree across error families (#937)", {
  # Direct gradient validation, not convergence side-effects: the NONMEM
  # comparison suite pins the gradient only through one proportional run's
  # per-subject outputs (the transform-family expectations are nlmixr2's own
  # golden values), and an optimizer's line search on the true objective can
  # rescue a wrong gradient without moving the converged objective past
  # tolerance.  Here every family's exposed (value, gradient) pair must be
  # internally consistent by central differences -- this fails on its first
  # assertion for a sign or term error, for error models NONMEM cannot run.
  skip_on_cran()
  d <- .foceiPtrData()
  .base <- rxode2::rxode2(.foceiPtrMod)
  .fam <- list(
    add      = list(mod = function(f) f,
                    lik = "focei"),
    prop     = list(mod = function(f) f |>
                      rxode2::model(cp ~ prop(prop.sd)) |>
                      rxode2::ini(prop.sd = 0.1),
                    lik = "focei"),   # R depends on eta: exercises dR/deta
    propFoce = list(mod = function(f) f |>
                      rxode2::model(cp ~ prop(prop.sd)) |>
                      rxode2::ini(prop.sd = 0.1),
                    lik = "foce"),    # R frozen at eta=0: consistent pair
    propT    = list(mod = function(f) f |>
                      rxode2::model(cp ~ propT(prop.sd)) |>
                      rxode2::ini(prop.sd = 0.1),
                    lik = "focei"),
    pow      = list(mod = function(f) f |>
                      rxode2::model(cp ~ pow(pow.sd, pw)) |>
                      rxode2::ini(pow.sd = 0.1, pw = 0.5),
                    lik = "focei"),
    boxCox   = list(mod = function(f) f |>
                      rxode2::model(cp ~ add(add.sd) + boxCox(lambda)) |>
                      rxode2::ini(lambda = 0.5),
                    lik = "focei"),
    yeoJohnson = list(mod = function(f) f |>
                      rxode2::model(cp ~ add(add.sd) + yeoJohnson(lambda)) |>
                      rxode2::ini(lambda = 0.5),
                    lik = "focei"),
    lnorm    = list(mod = function(f) f |>
                      rxode2::model(cp ~ lnorm(lnorm.sd)) |>
                      rxode2::ini(lnorm.sd = 0.1),
                    lik = "focei"))
  .testSeed(13)
  for (.n in names(.fam)) {
    .s <- .fam[[.n]]
    h <- foceiLikLoad(.s$mod(.base), d, .s$lik)
    eta <- matrix(stats::rnorm(h$nid * h$neta, 0, 0.2), h$nid, h$neta)
    expect_equal(foceiLikSetThetaC_(h$initPar), 0L, info = .n)
    got <- foceiLikCondGrad_(eta, 1L)
    .h <- 1e-5
    fd <- matrix(0, h$nid, h$neta)
    for (k in seq_len(h$neta)) {
      up <- eta; up[, k] <- up[, k] + .h
      dn <- eta; dn[, k] <- dn[, k] - .h
      fd[, k] <- (foceiLikCondGrad_(up, 1L)$value -
                    foceiLikCondGrad_(dn, 1L)$value) / (2 * .h)
    }
    expect_equal(as.numeric(got$grad), as.numeric(fd), tolerance = 1e-4,
                 info = .n)
    foceiLikUnload()
  }

  # focep ("foce+") keeps the live conditional R in the value while lp omits
  # its dR/deta term: value and gradient are gradients of DIFFERENT functions.
  # That inconsistency is exactly why dims flags it (0x02) for refusal by
  # gradient-based callers -- lock the reason in, not just the flag.
  h <- foceiLikLoad(.base |> rxode2::model(cp ~ prop(prop.sd)) |>
                      rxode2::ini(prop.sd = 0.1), d, "focep")
  on.exit(foceiLikUnload(), add = TRUE)
  expect_equal(bitwAnd(foceiLikDims_()$flags, 0x02), 0x02)
  eta <- matrix(c(-0.2, 0.1, 0.3, -0.15), h$nid, h$neta)
  expect_equal(foceiLikSetThetaC_(h$initPar), 0L)
  got <- foceiLikCondGrad_(eta, 1L)
  .h <- 1e-5
  fd <- matrix(0, h$nid, h$neta)
  for (k in seq_len(h$neta)) {
    up <- eta; up[, k] <- up[, k] + .h
    dn <- eta; dn[, k] <- dn[, k] - .h
    fd[, k] <- (foceiLikCondGrad_(up, 1L)$value -
                  foceiLikCondGrad_(dn, 1L)$value) / (2 * .h)
  }
  expect_false(isTRUE(all.equal(as.numeric(got$grad), as.numeric(fd),
                                tolerance = 1e-3)))
})

test_that("dose-handling theta sensitivities carry the event jump (#946)", {
  # An estimated alag theta's derivative needs a jump condition at the dose
  # event.  The theta-sensitivity model is now compiled with rxode2's
  # analytic event ("jump") sensitivities and solved under its own event
  # shape (OdeSwapEsBatch), so the column is real rather than silently zero.
  skip_on_cran()
  .lagMod <- function() {
    ini({
      tcl <- 1
      tv <- 3
      tlag <- -1
      add.sd <- 0.5
      eta.cl ~ 0.1
    })
    model({
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      d / dt(central) <- -cl / v * central
      alag(central) <- exp(tlag)
      cp <- central / v
      cp ~ add(add.sd)
    })
  }
  .testSeed(42)
  d <- do.call(rbind, lapply(1:4, function(id) {
    rbind(data.frame(ID = id, TIME = 0, DV = NA_real_, AMT = 100, EVID = 1),
          data.frame(ID = id, TIME = c(0.5, 1, 2, 4, 8),
                     DV = 5 * exp(-0.05 * c(0.5, 1, 2, 4, 8)) +
                       stats::rnorm(5, 0, 0.5),
                     AMT = 0, EVID = 0))
  }))
  h <- foceiLikLoad(.lagMod, d, "focei", scale = "natural", thetaSens = TRUE)
  on.exit(foceiLikUnload(), add = TRUE)
  # tlag (ntheta 3) is a non-mu structural theta and carries a sensitivity
  expect_true(3L %in% h$thetaSensIdx)
  expect_equal(bitwAnd(foceiLikDims_()$flags, 0x40), 0x40)
  th <- c(1, 3, -1, 0.5, h$initPar[5])
  expect_equal(foceiLikSetThetaC_(th), 0L)
  eta <- matrix(c(-0.1, 0.05, 0.2, -0.15), 4, 1)
  got <- foceiLikCondThetaGrad_(eta, 1L)
  expect_equal(got$nBad, 0L)
  # the derivative through the event is real (nonzero) ...
  expect_true(all(abs(got$dTheta[, 3]) > 1e-3))
  # ... and every sensitivity column (tv, tlag, add.sd) FD-agrees
  .h <- 1e-5
  for (t in h$thetaSensIdx) {
    up <- th
    up[t] <- up[t] + .h
    dn <- th
    dn[t] <- dn[t] - .h
    expect_equal(foceiLikSetThetaC_(up), 0L)
    vUp <- foceiLikCondGrad_(eta, 1L)$value
    expect_equal(foceiLikSetThetaC_(dn), 0L)
    vDn <- foceiLikCondGrad_(eta, 1L)$value
    fd <- (vUp - vDn) / (2 * .h)
    expect_equal(as.numeric(got$dTheta[, t]), as.numeric(fd),
                 tolerance = 1e-3, info = paste0("theta ", t))
  }
})

test_that("an eta entering dose handling gets its jump in the eta gradient (#946)", {
  # The inner batch installs the INNER model's event shape, so d/d(eta) of the
  # conditional through an eta-in-alag dose event is real, not silently zero.
  skip_on_cran()
  .lagEtaMod <- function() {
    ini({
      tcl <- 1
      tv <- 3
      tlag <- -1
      add.sd <- 0.5
      eta.cl ~ 0.1
      eta.lag ~ 0.05
    })
    model({
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      d / dt(central) <- -cl / v * central
      alag(central) <- exp(tlag + eta.lag)
      cp <- central / v
      cp ~ add(add.sd)
    })
  }
  .testSeed(42)
  d <- do.call(rbind, lapply(1:4, function(id) {
    rbind(data.frame(ID = id, TIME = 0, DV = NA_real_, AMT = 100, EVID = 1),
          data.frame(ID = id, TIME = c(0.5, 1, 2, 4, 8),
                     DV = 5 * exp(-0.05 * c(0.5, 1, 2, 4, 8)) +
                       stats::rnorm(5, 0, 0.5),
                     AMT = 0, EVID = 0))
  }))
  h <- foceiLikLoad(.lagEtaMod, d, "focei", scale = "natural")
  on.exit(foceiLikUnload(), add = TRUE)
  # full par vector: 4 thetas + the 2 omega parameters at their initials
  th <- h$initPar
  th[1:4] <- c(1, 3, -1, 0.5)
  expect_equal(foceiLikSetThetaC_(th), 0L)
  eta <- matrix(c(-0.1, 0.05, 0.2, -0.15,
                  0.08, -0.04, 0.1, -0.06), 4, 2)
  got <- foceiLikCondGrad_(eta, 1L)
  expect_equal(got$nBad, 0L)
  # the eta.lag column runs only through the dose event: it must be real
  expect_true(all(abs(got$grad[, 2]) > 1e-3))
  # and both eta columns FD-agree
  .h <- 1e-5
  for (k in 1:2) {
    up <- eta
    up[, k] <- up[, k] + .h
    dn <- eta
    dn[, k] <- dn[, k] - .h
    fd <- (foceiLikCondGrad_(up, 1L)$value -
             foceiLikCondGrad_(dn, 1L)$value) / (2 * .h)
    expect_equal(as.numeric(got$grad[, k]), as.numeric(fd),
                 tolerance = 1e-3, info = paste0("eta ", k))
  }
})

test_that("an estimated transform-both-sides lambda carries its column (#949)", {
  # lambda is a residual-error theta, so the sensitivity model hard-coded
  # d(f)/d(lambda) = 0 -- but the conditional depends on lambda through BOTH
  # sides of err = h(y; lambda) - h(f; lambda).  The f side now comes from the
  # model (rxTBSdL) and the DV side from d(lambda)/d(theta) times the analytic
  # d(h(y))/d(lambda), applied where the DV transform is (in C++).
  skip_on_cran()
  d <- .foceiPtrData()
  .base <- rxode2::rxode2(.foceiPtrMod)
  .fam <- list(
    boxCox = .base |> rxode2::model(cp ~ add(add.sd) + boxCox(lambda)) |>
      rxode2::ini(lambda = 0.5),
    yeoJohnson = .base |> rxode2::model(cp ~ add(add.sd) + yeoJohnson(lambda)) |>
      rxode2::ini(lambda = 0.5),
    boxCoxProp = .base |> rxode2::model(cp ~ prop(prop.sd) + boxCox(lambda)) |>
      rxode2::ini(prop.sd = 0.1, lambda = 0.5),
    # composed transforms: lambda applies to the logit/probit-mapped scale
    logitYj = .base |>
      rxode2::model(cp ~ logitNorm(lg.sd, 0, 10) + yeoJohnson(lambda)) |>
      rxode2::ini(lg.sd = 0.4, lambda = 0.5),
    probitYj = .base |>
      rxode2::model(cp ~ probitNorm(pb.sd, 0, 10) + yeoJohnson(lambda)) |>
      rxode2::ini(pb.sd = 0.4, lambda = 0.5))
  eta <- matrix(c(-0.1, 0.05, 0.2, -0.15), 4, 1)
  for (.n in names(.fam)) {
    h <- foceiLikLoad(.fam[[.n]], d, "focei", scale = "natural", thetaSens = TRUE)
    expect_true(isTRUE(h$thetaSens), info = .n)
    # lambda is ntheta 4 (tcl, tv, <resid sd>, lambda); tcl is mu-referenced
    expect_equal(h$thetaSensIdx, c(2L, 3L, 4L), info = .n)
    th <- h$initPar
    th[1:4] <- c(1, 3, switch(.n, boxCoxProp = 0.1, logitYj = 0.4,
                              probitYj = 0.4, 0.5), 0.5)
    expect_equal(foceiLikSetThetaC_(th), 0L, info = .n)
    got <- foceiLikCondThetaGrad_(eta, 1L)
    expect_equal(got$nBad, 0L, info = .n)
    # the historical zero is gone
    expect_true(all(abs(got$dTheta[, 4]) > 1e-3), info = .n)
    .h <- 1e-5
    for (t in h$thetaSensIdx) {
      up <- th
      up[t] <- up[t] + .h
      dn <- th
      dn[t] <- dn[t] - .h
      expect_equal(foceiLikSetThetaC_(up), 0L)
      vUp <- foceiLikCondGrad_(eta, 1L)$value
      expect_equal(foceiLikSetThetaC_(dn), 0L)
      vDn <- foceiLikCondGrad_(eta, 1L)$value
      fd <- (vUp - vDn) / (2 * .h)
      expect_equal(as.numeric(got$dTheta[, t]), as.numeric(fd),
                   tolerance = 1e-4, info = paste0(.n, " theta ", t))
    }
    foceiLikUnload()
  }
})

test_that("the lambda column is right for censored and multi-endpoint data (#949)", {
  # The DV-side term also moves the censoring bounds (M3/M4 read the TRANSFORMED
  # DV and LIMIT), and with several endpoints rx_lambda_ is a CMT switch so only
  # the transformed endpoint's rows may pick the term up.
  skip_on_cran()
  .testSeed(42)
  .censData <- function(limit = FALSE) {
    do.call(rbind, lapply(1:4, function(id) {
      tt <- c(0.5, 1, 2, 4, 8)
      dv <- 5 * exp(-0.05 * tt) + stats::rnorm(5, 0, 0.5)
      d <- data.frame(ID = id, TIME = tt, DV = dv, AMT = 0, EVID = 0,
                      CENS = c(0, 0, 0, -1, -1))
      d$DV[d$CENS == -1] <- 4.6
      if (limit) d$LIMIT <- ifelse(d$CENS == -1, 3.5, NA_real_)
      d
    }))
  }
  .tbsMod <- rxode2::rxode2(.foceiPtrMod) |>
    rxode2::model(cp ~ add(add.sd) + boxCox(lambda)) |>
    rxode2::ini(lambda = 0.5)
  .tbsYj <- rxode2::rxode2(.foceiPtrMod) |>
    rxode2::model(cp ~ add(add.sd) + yeoJohnson(lambda)) |>
    rxode2::ini(lambda = 0.5)
  eta <- matrix(c(-0.1, 0.05, 0.2, -0.15), 4, 1)
  .fdCheck <- function(h, th, info) {
    expect_equal(foceiLikSetThetaC_(th), 0L, info = info)
    got <- foceiLikCondThetaGrad_(eta, 1L)
    expect_equal(got$nBad, 0L, info = info)
    .h <- 1e-5
    for (t in h$thetaSensIdx) {
      up <- th
      up[t] <- up[t] + .h
      dn <- th
      dn[t] <- dn[t] - .h
      expect_equal(foceiLikSetThetaC_(up), 0L)
      vUp <- foceiLikCondGrad_(eta, 1L)$value
      expect_equal(foceiLikSetThetaC_(dn), 0L)
      vDn <- foceiLikCondGrad_(eta, 1L)$value
      expect_equal(as.numeric(got$dTheta[, t]),
                   as.numeric((vUp - vDn) / (2 * .h)),
                   tolerance = 1e-4, info = paste0(info, " theta ", t))
    }
    got
  }
  # M3 (CENS only), M4 (CENS + LIMIT) and M2 (finite LIMIT, CENS 0): each takes a
  # different branch of the dv/limit partials
  .cases <- list(M3 = .censData(FALSE), M4 = .censData(TRUE))
  .cases$M2 <- .censData(TRUE)
  .cases$M2$DV <- .cases$M3$DV
  .cases$M2$CENS <- NULL
  # both transform families: the DV/LIMIT partials run through _powerDLambda, whose
  # branch differs per transform
  for (.tr in c("boxCox", "yeoJohnson")) {
    for (.nm in names(.cases)) {
      h <- foceiLikLoad(if (.tr == "boxCox") .tbsMod else .tbsYj, .cases[[.nm]],
                        "focei", scale = "natural", thetaSens = TRUE)
      th <- h$initPar
      th[1:4] <- c(1, 3, 0.5, 0.5)
      got <- .fdCheck(h, th, paste0(.tr, " ", .nm))
      expect_true(all(abs(got$dTheta[, 4]) > 1e-3), info = paste0(.tr, " ", .nm))
      foceiLikUnload()
    }
  }
  # two endpoints, only the first transformed
  .multiMod <- function() {
    ini({
      tcl <- 1
      tv <- 3
      add.sd <- 0.5
      lambda <- c(-2, 0.5, 2)
      ef.sd <- 0.7
      eta.cl ~ 0.1
    })
    model({
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      cp <- 100 / v * exp(-cl / v * time)
      ef <- 10 - cp
      cp ~ add(add.sd) + boxCox(lambda)
      ef ~ add(ef.sd)
    })
  }
  .testSeed(42)
  dm <- do.call(rbind, lapply(1:4, function(id) {
    tt <- c(0.5, 1, 2, 4, 8)
    rbind(data.frame(ID = id, TIME = tt, DVID = "cp",
                     DV = 5 * exp(-0.05 * tt) + stats::rnorm(5, 0, 0.5)),
          data.frame(ID = id, TIME = tt, DVID = "ef",
                     DV = 5 + 0.05 * tt + stats::rnorm(5, 0, 0.7)))
  }))
  dm$AMT <- 0
  dm$EVID <- 0
  h <- foceiLikLoad(.multiMod, dm, "focei", scale = "natural", thetaSens = TRUE)
  on.exit(foceiLikUnload(), add = TRUE)
  expect_equal(h$thetaSensIdx, c(2L, 3L, 4L, 5L))
  th <- h$initPar
  th[1:5] <- c(1, 3, 0.5, 0.5, 0.7)
  got <- .fdCheck(h, th, "multi")
  expect_true(all(abs(got$dTheta[, 4]) > 1e-3))
  # and thread-count invariant with the extra per-row term
  expect_equal(foceiLikSetThetaC_(th), 0L)
  expect_equal(foceiLikCondThetaGrad_(eta, 4L)$dTheta, got$dTheta,
               tolerance = 1e-12)
})

test_that("mixture models: the component-major batch layout is blessed (#955)", {
  # nid = nsub * nMix with id = component*nsub + subject; each row is the
  # COMPONENT-conditional log-likelihood (p-free -- the mixing probability
  # enters only the caller's mixture weights, e.g. Stan's log_sum_exp)
  skip_on_cran()
  .mixMod <- function() {
    ini({
      tcl1 <- 1
      tcl2 <- 2
      tv <- 3
      p1 <- 0.3
      add.sd <- 0.5
      eta.cl ~ 0.1
    })
    model({
      cl <- mix(exp(tcl1 + eta.cl), p1, exp(tcl2 + eta.cl))
      v <- exp(tv)
      cp <- 100 / v * exp(-cl / v * time)
      cp ~ add(add.sd)
    })
  }
  .testSeed(42)
  .tt <- c(0.5, 1, 2, 4, 8)
  d <- do.call(rbind, lapply(1:4, function(id) {
    data.frame(ID = id, TIME = .tt,
               DV = 5 * exp(-0.05 * .tt) + stats::rnorm(5, 0, 0.5),
               AMT = 0, EVID = 0)
  }))
  h <- foceiLikLoad(.mixMod, d, "focei", scale = "natural", thetaSens = TRUE)
  on.exit(foceiLikUnload(), add = TRUE)
  expect_equal(foceiLikNMix_(), 2L)
  # flags 0x10 stays set as the LAYOUT marker (no longer a refusal)
  expect_equal(bitwAnd(foceiLikDims_()$flags, 0x10), 0x10)
  .th <- h$initPar
  expect_equal(foceiLikSetThetaC_(.th), 0L)
  set.seed(7)
  eta <- matrix(stats::rnorm(8, 0, 0.2), 8, 1) # 2 components x 4 subjects
  got <- foceiLikCondGrad_(eta, 1L)
  expect_equal(got$nBad, 0L)
  # each row is the COMPONENT-conditional density: hand-computable with that
  # component's cl, up to the adjLik constant (+nobs * 0.5*log(2pi))
  .hand <- vapply(1:8, function(r) {
    .i <- ((r - 1) %% 4) + 1
    .m <- ((r - 1) %/% 4) + 1
    .di <- d[d$ID == .i, ]
    .cl <- exp(c(1, 2)[.m] + eta[r, 1])
    .f <- 100 / exp(3) * exp(-.cl / exp(3) * .di$TIME)
    sum(stats::dnorm(.di$DV, .f, 0.5, log = TRUE))
  }, numeric(1))
  expect_equal(got$value, .hand + 5 * 0.5 * log(2 * pi), tolerance = 1e-8)
  # eta gradient FD-agrees per expanded row
  .h <- 1e-5
  fd <- (foceiLikCondGrad_(eta + .h, 1L)$value -
           foceiLikCondGrad_(eta - .h, 1L)$value) / (2 * .h)
  expect_equal(as.numeric(got$grad), as.numeric(fd), tolerance = 1e-4)
  # theta gradients: component 1 rows respond to tcl1 only, component 2 to
  # tcl2 only, and NO row responds to the mixing probability (the
  # conditional is p-free by construction)
  gt <- foceiLikCondThetaGrad_(eta, 1L)
  expect_equal(gt$nBad, 0L)
  for (t in h$thetaSensIdx) {
    up <- .th
    up[t] <- up[t] + .h
    dn <- .th
    dn[t] <- dn[t] - .h
    expect_equal(foceiLikSetThetaC_(up), 0L)
    vUp <- foceiLikCondGrad_(eta, 1L)$value
    expect_equal(foceiLikSetThetaC_(dn), 0L)
    vDn <- foceiLikCondGrad_(eta, 1L)$value
    expect_equal(as.numeric(gt$dTheta[, t]),
                 as.numeric((vUp - vDn) / (2 * .h)),
                 tolerance = 1e-3, info = paste0("theta ", t))
  }
  expect_equal(foceiLikSetThetaC_(.th), 0L)
  # determinism in the expanded layout
  expect_identical(foceiLikCondGrad_(eta, 1L), foceiLikCondGrad_(eta, 1L))
  # and the WRONG shape still refuses
  expect_error(foceiLikCondGrad_(eta[1:4, , drop = FALSE], 1L), "-2")
})

test_that("sampler iteration print: scale.h rows + parHistData over the residency", {
  skip_on_cran()
  .mod <- function() {
    ini({
      tcl <- 1
      tv <- 3
      add.sd <- 0.5
      eta.cl ~ 0.1
    })
    model({
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      cp <- 100 / v * exp(-cl / v * time)
      cp ~ add(add.sd)
    })
  }
  .testSeed(42)
  d <- do.call(rbind, lapply(1:4, function(id) {
    data.frame(ID = id, TIME = c(0.5, 1, 2, 4, 8),
               DV = 5 * exp(-0.05 * c(0.5, 1, 2, 4, 8)) +
                 stats::rnorm(5, 0, 0.5),
               AMT = 0, EVID = 0)
  }))
  h <- foceiLikLoad(.mod, d, "focei", scale = "natural")
  on.exit(foceiLikUnload(), add = TRUE)
  # display vector: natural-scale thetas + the sampler's current ACTUAL
  # omega entries (not the internal chol(Omega^-1) diagXform tail), named
  # om.<eta> for variances / cov.<eta1>.<eta2> for covariances
  .par <- c(h$initPar[seq_len(h$ntheta)], 0.1)
  .nm <- c(h$thetaNames, "om.eta.cl")
  .out <- utils::capture.output(type = "message", {
    foceiLikIterPrintStart_(2L, .par, .nm)
    # 5 ticks at cadence 2: the FIRST tick prints, then every 2nd after ->
    # rows at ticks 1, 3 and 5
    for (i in 1:5) {
      expect_equal(foceiLikRowTick_(.par, -100 - i), 0L)
    }
    .ph <- foceiLikIterPrintEnd_()
  })
  expect_true(is.data.frame(.ph))
  # two printed rows, each recorded in the history (plus back-transform rows)
  expect_equal(sum(.ph$type == "Scaled"), 3L)
  expect_true(all(c("tcl", "tv", "add.sd", "om.eta.cl") %in% names(.ph)))
  # inactive after End
  expect_equal(foceiLikRowTick_(h$initPar, -1), -1L)
})

test_that("combined eta+theta sensitivity build (#958): layout, parity, fused entry", {
  skip_on_cran()
  .mod <- function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45
      kout <- 0.2
      eta.ka ~ 0.6; eta.cl ~ 0.3
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center - kout * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }
  .d <- nlmixr2data::theo_sd
  .testSeed(11)
  .eta <- matrix(stats::rnorm(24, 0, 0.2), 12, 2)
  # reference: the two-model path
  .h1 <- foceiLikLoad(.mod, .d, "focei", scale = "natural", thetaSens = TRUE)
  .v1 <- foceiLikCondGrad_(.eta, 1L)
  .t1 <- foceiLikCondThetaGrad_(.eta, 1L)
  foceiLikUnload()
  # combined build
  .h2 <- foceiLikLoad(.mod, .d, "focei", scale = "natural", combSens = TRUE)
  on.exit(foceiLikUnload(), add = TRUE)
  .d2 <- foceiLikDims_()
  expect_true(bitwAnd(.d2[["flags"]], 0x40L) != 0L) # theta sens wired
  expect_true(bitwAnd(.d2[["flags"]], 0x80L) != 0L) # combined build
  .v2 <- foceiLikCondGrad_(.eta, 1L)
  .t2 <- foceiLikCondThetaGrad_(.eta, 1L)
  # value/eta grad agree at solver tolerance (the theta-sens states ride the
  # same adaptive integration, so bit-parity is impossible by construction)
  expect_equal(.v1$value, .v2$value, tolerance = 1e-6)
  expect_equal(.v1$grad, .v2$grad, tolerance = 1e-5)
  # theta gradient from the combined solve == two-model path (solver tol)
  expect_equal(.t1$dTheta, .t2$dTheta, tolerance = 1e-5)
  # fused entry: one solve per subject, all three outputs; value/etaGrad
  # BITWISE equal to the separate combined-path calls, dTheta bitwise too
  # (the reuse pass reads the very solve the value pass produced)
  .f <- foceiLikCondBatchThetaGrad_(.eta, 1L)
  expect_identical(.f$value, .v2$value)
  expect_identical(.f$gradEta, .v2$grad)
  expect_identical(.f$dTheta, .t2$dTheta)
  expect_equal(.f$nBad, 0L)
})

test_that("combined build: sigma-only theta set adds columns, no states (#958)", {
  skip_on_cran()
  # every structural theta mu-referenced; only add.sd (sigma) remains -- the
  # paste0-recycling guard keeps a malformed THETA__ variable out of the build
  .mod <- function() {
    ini({
      tcl <- 1
      eta.cl ~ 0.1
      add.sd <- 0.5
    })
    model({
      cl <- exp(tcl + eta.cl)
      cp <- 100 * exp(-cl * time)
      cp ~ add(add.sd)
    })
  }
  .d <- data.frame(ID = rep(1:4, each = 3), TIME = rep(c(1, 2, 4), 4),
                   DV = 3, AMT = 0, EVID = 0)
  .h <- foceiLikLoad(.mod, .d, "focei", scale = "natural", combSens = TRUE)
  on.exit(foceiLikUnload(), add = TRUE)
  .f <- foceiLikCondBatchThetaGrad_(matrix(0, 4, 1), 1L)
  expect_equal(.f$nBad, 0L)
  expect_true(all(is.finite(.f$value)))
})

test_that("IOV magnitude theta sensitivity column FD-agrees (#952)", {
  skip_on_cran()
  # the sd/var transforms used abs(theta), whose symengine rewrite carries
  # an indicator with a tanh-SMOOTHED derivative -- the column came out
  # wrong by 1 + 10*theta*(1 - tanh(10*theta)^2) (x1.42 at theta=0.1).
  # |theta| is now emitted as sqrt(theta^2) (identical value, exact
  # derivative), and this locks the FD agreement in for every transform.
  .mod <- function() {
    ini({
      tcl <- 1.05; tv <- 2.95; add.sd <- 0.55
      eta.cl ~ 0.1
      iov.cl ~ 0.02 | OCC
    })
    model({
      cl <- exp(tcl + eta.cl + iov.cl)
      v <- exp(tv)
      cp <- 100 / v * exp(-cl / v * time)
      cp ~ add(add.sd)
    })
  }
  .testSeed(7)
  .d <- expand.grid(ID = 1:4, OCC = 1:2, TIME = c(1, 2, 4))
  .d <- .d[order(.d$ID, .d$OCC, .d$TIME), ]
  .d$DV <- 3 + stats::rnorm(nrow(.d), 0, 0.5)
  .d$AMT <- 0; .d$EVID <- 0
  for (.xf in c("sd", "logsd", "var", "logvar")) {
    .h <- foceiLikLoad(.mod, .d, "focei", scale = "natural",
                       thetaSens = TRUE, iovXform = .xf)
    .th <- .h$initPar
    .iov <- grep("iov", .h$thetaNames)
    .testSeed(3)
    .eta <- matrix(stats::rnorm(.h$nid * .h$neta, 0, 0.3), .h$nid, .h$neta)
    .th2 <- .th
    .th2[.iov] <- if (.xf %in% c("logsd", "logvar")) -1.5 else 0.15
    foceiLikSetTheta_(.th2)
    .an <- foceiLikCondThetaGrad_(.eta, 1L)$dTheta[, .iov]
    .hs <- 1e-5
    .tp <- .th2; .tp[.iov] <- .tp[.iov] + .hs
    .tm <- .th2; .tm[.iov] <- .tm[.iov] - .hs
    foceiLikSetTheta_(.tp); .vp <- foceiLikCondGrad_(.eta, 1L)$value
    foceiLikSetTheta_(.tm); .vm <- foceiLikCondGrad_(.eta, 1L)$value
    .fd <- (.vp - .vm) / (2 * .hs)
    expect_lt(max(abs(.an / .fd - 1)), 1e-6,
              label = paste0("iovXform=", .xf, " analytic/FD"))
    foceiLikUnload()
  }
})
