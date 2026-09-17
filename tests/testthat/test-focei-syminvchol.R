## #1079 / rxode2#1365: rxSymInvCholCreate() refuses a POSITIVE-DEFINITE omega
## whose zero pattern is not block-decomposable, so the old chol() guard let it
## through and the setup died with "theta has to have N elements".  Cheap (one
## small fit), so this stays in the essential push/PR subset.

nmTest({
  .blockZeroMod <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      ## the (eta.ka, eta.v) covariance is declared at exactly 0 INSIDE the
      ## block -- positive definite, but not block-decomposable
      eta.ka + eta.cl + eta.v ~ c(0.1,
                                  0.01, 0.1,
                                  0.00, 0.01, 0.1)
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }

  test_that("the declared omega really is the rxode2#1365 shape", {
    .ui <- rxode2::rxUiDecompress(rxode2::rxode2(.blockZeroMod))
    .om <- .ui$omega
    expect_true(all(eigen(.om)$values > 0))
    expect_equal(nrow(.omegaBlockZeros(.om)), 1L)
    ## the raw call is what used to abort the fit
    expect_error(rxode2::rxSymInvCholCreate(mat=.om, diag.xform="sqrt"),
                 "theta has to have")
  })

  test_that(".foceiSymInvCholCreate fills a block-internal zero", {
    .ui <- rxode2::rxUiDecompress(rxode2::rxode2(.blockZeroMod))
    .om <- .ui$omega
    expect_warning(.r <- .foceiSymInvCholCreate(.om, "sqrt", NULL),
                   "omega block zero cov is estimated")
    ## the mechanism: the returned matrix is the FILLED one and the inverse
    ## carries the full dense-block parameter count (3 diag + 3 off-diag)
    expect_equal(length(.r$rxInv$theta), 6L)
    expect_equal(nrow(.omegaBlockZeros(.r$mat)), 0L)
    ## it names which random effects were involved
    .w <- tryCatch(.foceiSymInvCholCreate(.om, "sqrt", NULL),
                   warning=function(w) conditionMessage(w))
    expect_true(grepl("eta.ka", .w, fixed=TRUE))
    expect_true(grepl("eta.v", .w, fixed=TRUE))
  })

  test_that(".foceiSymInvCholCreate is a no-op on an acceptable omega", {
    .om <- matrix(c(0.1, 0.01, 0.01, 0.1), 2, 2,
                  dimnames=list(c("eta.ka", "eta.cl"), c("eta.ka", "eta.cl")))
    expect_warning(.r <- .foceiSymInvCholCreate(.om, "sqrt", NULL), NA)
    expect_equal(.r$mat, .om)
    expect_equal(length(.r$rxInv$theta), 3L)
  })

  test_that(".foceiOptEnvSetupBounds carries the fill into the bounds", {
    .ui <- rxode2::rxUiDecompress(rxode2::rxode2(.blockZeroMod))
    assign("control", foceiControl(), envir=.ui)
    .env <- new.env(parent=emptyenv())
    .env$etaNames <- c("eta.ka", "eta.cl", "eta.v")
    expect_warning(.foceiOptEnvSetupBounds(.ui, .env),
                   "omega block zero cov is estimated")
    expect_false(is.null(.env$rxInv))
    ## nomega and the bound vectors agree with the (filled) parameter count --
    ## the a/b bound matrices are built from the SAME repaired omega, so their
    ## theta vectors line up row for row with env$rxInv's
    expect_equal(rxode2::rxGetControl(.ui, "nomega", -1L), 6L)
    expect_equal(length(.env$lower), length(.env$upper))
    expect_equal(length(.env$lower),
                 rxode2::rxGetControl(.ui, "ntheta", 0L) + 6L)
  })

  test_that("a FOCEi fit with a block-internal omega zero runs (#1079)", {
    .fit <- nlmixr2(.blockZeroMod, nlmixr2data::theo_sd, est="focei",
                    control=foceiControl(maxOuterIterations=0,
                                         maxInnerIterations=10,
                                         covMethod="", calcTables=FALSE,
                                         print=0))
    expect_true(inherits(.fit, "nlmixr2FitCore"))
    ## the fit estimated the previously-structural zero rather than erroring
    expect_equal(dim(.fit$omega), c(3L, 3L))
    ## and it said so, naming the random effects
    expect_true(any(grepl("omega block zero cov is estimated", .fit$runInfo)))
  })

  test_that("est='vae' survives the same omega (#1079)", {
    .fit <- suppressMessages(
      nlmixr2(.blockZeroMod, nlmixr2data::theo_sd, est="vae",
              control=vaeControl(iters=3, itersBurnIn=2, calcTables=FALSE)))
    expect_true(inherits(.fit, "nlmixr2FitCore"))
    expect_equal(dim(.fit$omega), c(3L, 3L))
  })


  test_that("a non-contiguous correlated block is filled, not flattened", {
    ## eta1 correlates with eta3 and eta2 sits between them: every component is
    ## dense, so a component-only rule calls this fine -- but the call refuses
    ## it, and the repair ladder would then have dropped the 0.5 covariance for
    ## a floored diagonal.
    .nm <- c("eta.a", "eta.b", "eta.c")
    .om <- matrix(c(0.1, 0, 0.05, 0, 0.1, 0, 0.05, 0, 0.1), 3, 3,
                  dimnames=list(.nm, .nm))
    expect_error(rxode2::rxSymInvCholCreate(mat=.om, diag.xform="sqrt"))
    expect_warning(.r <- .foceiSymInvCholCreate(.om, "sqrt", NULL),
                   "omega block zero cov is estimated")
    ## the covariance SURVIVED -- this is the check a floored-diagonal fallback
    ## would fail
    expect_equal(.r$mat[1, 3], 0.05)
    expect_equal(diag(.r$mat), diag(.om))
    expect_equal(length(.r$rxInv$theta), 6L)
  })


  test_that("a 2x2 block declaring a 0 covariance is left alone", {
    ## Two etas whose covariance is declared at exactly 0 are structurally
    ## uncorrelated -- the call accepts that, and the repair must not decide it
    ## knows better and start estimating a covariance the model did not ask for.
    .nm <- c("eta.ka", "eta.cl")
    .om <- matrix(c(0.1, 0, 0, 0.1), 2, 2, dimnames=list(.nm, .nm))
    expect_equal(nrow(.omegaBlockZeros(.om)), 0L)
    expect_warning(.r <- .foceiSymInvCholCreate(.om, "sqrt", NULL), NA)
    expect_equal(.r$mat, .om)
    expect_equal(length(.r$rxInv$theta), 2L)
  })

  test_that("the vae omega position list follows the FILLED omega (#1079)", {
    ## vaeOmegaSel is the 0-based position list the C++ fast path packs
    ## chol(Omega^-1) into, and it has to match rxSymInvCholCreate's parameter
    ## order.  Building it from the unrepaired omega would give 5 positions for
    ## a 6-parameter inverse.
    .theoZ <- function() {
      ini({
        lka <- log(1.8)
        lke <- log(0.086)
        lV <- log(32)
        eta.ka + eta.ke + eta.V ~ c(0.3,
                                    0.00, 0.03,
                                    0.02, 0.005, 0.03)
        add.err <- 0.7
      })
      model({
        ka <- exp(lka + eta.ka)
        ke <- exp(lke + eta.ke)
        V <- exp(lV + eta.V)
        d/dt(depot) = -ka * depot
        d/dt(central) = ka * depot - ke * central
        cp <- central / V
        cp ~ add(add.err)
      })
    }
    .ui <- rxode2::assertRxUi(.theoZ)
    ## the declared omega really needs the repair
    expect_equal(nrow(.omegaBlockZeros(.ui$omega)), 1L)
    .ctl <- vaeControl()
    .n <- length(unique(nlmixr2data::theo_sd$ID))
    set.seed(3)
    .etaMat <- matrix(rnorm(.n * 3, 0, 0.1), .n, 3)
    .prep <- .vaeDataPrep(.ui, nlmixr2data::theo_sd)
    .env <- .vaeInnerSetup(.ui, nlmixr2data::theo_sd, .etaMat, .ctl)
    on.exit(.vaeInnerFree(), add=TRUE)
    expect_equal(nrow(.env$vaeOmegaSel), length(.env$rxInv$theta))
    expect_equal(length(.env$rxInv$theta), 6L)
    ## the C++ fast path packs the OUTER (unrepaired) omega; it must still land
    ## where the repaired full re-setup does
    vaeInnerUpdatePar_(as.numeric(.prep$th), .prep$omegaMat)
    .fast <- .vaeInnerEval(.etaMat, .ctl, grad=TRUE)
    .vaeInnerUpdate(.env, .prep$th, .prep$omegaMat, .etaMat)
    .ref <- .vaeInnerEval(.etaMat, .ctl, grad=TRUE)
    expect_lt(max(abs(.fast$obj - .ref$obj)), 1e-8)
    expect_lt(max(abs(.fast$lp - .ref$lp)), 1e-8)
    ## and the 1e-10 fill is inert: declaring that covariance explicitly gives
    ## the same objective as leaving it at 0
    .omTiny <- .prep$omegaMat
    .omTiny[1, 3] <- .omTiny[3, 1] <-
      1e-10 * sqrt(.omTiny[1, 1] * .omTiny[3, 3])
    .vaeInnerUpdate(.env, .prep$th, .omTiny, .etaMat)
    .tiny <- .vaeInnerEval(.etaMat, .ctl, grad=TRUE)
    expect_equal(.ref$obj, .tiny$obj, tolerance=1e-10)
  })


  test_that("the last rungs: a floored diagonal, then an actionable error", {
    ## A structurally acceptable but NON positive-definite omega: there is
    ## nothing to fill, so the ladder has to run out to the floored diagonal.
    .nm <- c("eta1", "eta2")
    .om <- matrix(c(1, 2, 2, 1), 2, 2, dimnames=list(.nm, .nm))
    expect_true(any(eigen(.om)$values < 0))
    expect_equal(nrow(.omegaBlockZeros(.om)), 0L)
    expect_null(.omegaFillBlockZeros(.om))
    expect_warning(.r <- .foceiSymInvCholCreate(.om, "sqrt", NULL),
                   "floored diagonal")
    ## the diagonal is kept (it is already above the floor), the covariance goes
    expect_equal(diag(.r$mat), c(eta1=1, eta2=1))
    expect_equal(.r$mat[1, 2], 0)
    ## a non-finite omega floors to the minimum instead of erroring
    .nan <- matrix(NaN, 2, 2, dimnames=list(.nm, .nm))
    expect_warning(.rn <- .foceiSymInvCholCreate(.nan, "sqrt", NULL),
                   "floored diagonal")
    expect_equal(unname(diag(.rn$mat)), c(1e-6, 1e-6))
    ## with fallback=FALSE (the per-step vae rebuild) only the fill may run, so
    ## the same omega errors -- and the message names the random effects
    expect_error(.foceiSymInvCholCreate(.om, "sqrt", NULL, fallback=FALSE),
                 "eta1, eta2")
  })


  test_that("a repeated same() block keeps its sharing through the fill", {
    ## The fill is the SAME rule in every block (cor * sqrt(d_i d_j)), so
    ## identical repeated blocks stay identical and the same() map survives --
    ## which it must, or the bound matrices below it no longer mirror the
    ## right rows.
    .d <- nlmixr2data::theo_sd
    .d$occ <- 1 + (.d$TIME >= 5)
    .f <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        add.sd <- 0.7
        eta.ka ~ 0.6
        iov.ka + iov.cl + iov.v ~ c(0.1,
                                    0.02, 0.2,
                                    0.01, 0.03, 0.3) | occ
      })
      model({
        ka <- exp(tka + eta.ka + iov.ka)
        cl <- exp(tcl + iov.cl)
        v <- exp(tv + iov.v)
        linCmt() ~ add(add.sd)
      })
    }
    .ui <- .uiApplyIov(rxode2::rxode2(.f()), "focei", .d,
                       foceiControl(iovMethod="omega"))$ui
    .om <- .ui$omega
    .sm <- .ui$omegaSameMap
    ## eta.ka plus two repeated 3x3 blocks, sharing 7 parameters
    expect_equal(dim(.om), c(7L, 7L))
    expect_equal(.sm, c(0L, 0L, 0L, 0L, 2L, 3L, 4L))
    .nTheta <- length(rxode2::rxSymInvCholCreate(.om, "sqrt", same=.sm)$theta)
    expect_equal(.nTheta, 7L)
    ## zero the SAME within-block cell in both repeats
    .z <- .om
    .z[2, 4] <- .z[4, 2] <- 0
    .z[5, 7] <- .z[7, 5] <- 0
    expect_equal(nrow(.omegaBlockZeros(.z)), 2L)
    expect_warning(.r <- .foceiSymInvCholCreate(.z, "sqrt", .sm),
                   "omega block zero cov is estimated")
    ## the sharing survived: same map kept, parameter count unchanged, and the
    ## two blocks are still identical after the fill
    expect_equal(.r$same, .sm)
    expect_equal(length(.r$rxInv$theta), .nTheta)
    expect_equal(unname(.r$mat[2:4, 2:4]), unname(.r$mat[5:7, 5:7]))
    ## with same(), only the MASTER block's pattern is parameterized, so a zero
    ## in a REPEAT alone needs no repair at all
    .one <- .om
    .one[5, 7] <- .one[7, 5] <- 0
    expect_warning(.r2 <- .foceiSymInvCholCreate(.one, "sqrt", .sm), NA)
    expect_equal(.r2$same, .sm)
    expect_equal(length(.r2$rxInv$theta), .nTheta)
  })

  test_that("dropping same() sharing is never silent", {
    ## The ladder may estimate the repeated blocks independently to get an
    ## inverse at all.  That changes what is estimated, so it has to be said --
    ## a non-PD omega runs out to the floored diagonal, which cannot share.
    .nm <- paste0("eta", 1:4)
    .om <- matrix(0, 4, 4, dimnames=list(.nm, .nm))
    diag(.om) <- 1
    .om[1, 2] <- .om[2, 1] <- 2
    .om[3, 4] <- .om[4, 3] <- 2
    .sm <- c(0L, 0L, 1L, 2L)
    .w <- NULL
    withCallingHandlers(.foceiSymInvCholCreate(.om, "sqrt", .sm),
                        warning=function(w) {
                          .w <<- c(.w, conditionMessage(w))
                          invokeRestart("muffleWarning")
                        })
    expect_true(any(grepl("floored diagonal", .w, fixed=TRUE)))
    expect_true(any(grepl("same() sharing dropped", .w, fixed=TRUE)))
    ## and every note stays on one $runInfo line
    expect_true(all(nchar(.w) < 75L))
  })

})
