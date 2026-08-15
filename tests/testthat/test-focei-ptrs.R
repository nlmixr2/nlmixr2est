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

test_that("the foceiPtrs table has the documented shape (#937)", {
  .p <- .nlmixr2estFoceiPtrs()
  expect_length(.p, 7L)
  expect_equal(names(.p), c("apiVersion", "dims", "setTheta", "condBatch",
                            "setOmegaInv", "thetaSensIdx", "condThetaGrad"))
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
  # This is the test that catches the sign convention: fInd->lp is stored as
  # -(dlogp/deta) + Omega^-1 eta, so the conditional gradient must be
  # Omega^-1 eta - lp -- PLUS, not minus.  A sign error here is wrong by
  # 2*Omega^-1*eta and never errors, it just samples the wrong distribution.
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
