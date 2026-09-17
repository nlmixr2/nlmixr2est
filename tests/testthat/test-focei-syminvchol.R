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

})
