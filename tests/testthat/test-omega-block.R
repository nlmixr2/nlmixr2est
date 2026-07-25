## Unit tests for the shared omega-block helpers (R/omegaBlock.R) and for the
## omega parameter ORDER the vae inner fast path relies on.  Fast (no fits), so
## these stay in the essential push/PR subset -- the fit-based cross-method
## checks live in test-omega-offdiag.R (weekly batch).

nmTest({
  test_that(".omegaBlockFromIniDf builds the declared block + fix mask", {
    .idf <- data.frame(
      name = c("tka", "eta.cl", "(eta.cl,eta.v)", "eta.v", "eta.q"),
      ntheta = c(1, NA, NA, NA, NA),
      neta1 = c(NA, 1, 2, 2, 3), neta2 = c(NA, 1, 1, 2, 3),
      est = c(0.4, 0.1, 0.01, 0.2, 0.3),
      fix = c(FALSE, FALSE, FALSE, FALSE, TRUE),
      stringsAsFactors = FALSE)
    .b <- .omegaBlockFromIniDf(.idf, c("eta.cl", "eta.v", "eta.q"))
    expect_equal(unname(.b$mat[1L, 2L]), 0.01)
    expect_equal(unname(.b$mat[2L, 1L]), 0.01)          # symmetric
    expect_equal(unname(diag(.b$mat)), c(0.1, 0.2, 0.3))
    expect_equal(unname(.b$mat[1L, 3L]), 0)             # undeclared stays 0
    expect_true(.b$fixMat[3L, 3L])
    expect_false(.b$fixMat[1L, 2L])
    expect_true(.omegaHasOffDiag(.b$mat))
    # an eta absent from etaNames is skipped, not an error
    .b2 <- .omegaBlockFromIniDf(.idf, c("eta.cl", "eta.v"))
    expect_equal(dim(.b2$mat), c(2L, 2L))
    expect_equal(unname(.b2$mat[1L, 2L]), 0.01)
    # a diagonal-only model reports no off-diagonal
    expect_false(.omegaHasOffDiag(diag(c(0.1, 0.2))))
  })

  test_that("the vae omega position list matches rxSymInvCholCreate's order", {
    ## vaeInnerUpdatePar_ writes chol(Omega^-1) entries into the reduced par
    ## vector positionally, using which(upper.tri & nonzero, arr.ind=TRUE).
    ## That order MUST equal rxSymInvCholCreate(diag.xform="sqrt")'s theta
    ## order (sqrt on the diagonal), or the inner prior is silently wrong.
    .check <- function(om) {
      .sel <- upper.tri(om, diag = TRUE) & om != 0
      diag(.sel) <- TRUE
      .pos <- which(.sel, arr.ind = TRUE)                # column-major
      .u <- chol(solve(om))
      .expect <- vapply(seq_len(nrow(.pos)), function(k) {
        .i <- .pos[k, 1L]; .j <- .pos[k, 2L]
        if (.i == .j) sqrt(.u[.i, .i]) else .u[.i, .j]
      }, numeric(1))
      expect_equal(rxode2::rxSymInvCholCreate(mat = om, diag.xform = "sqrt")$theta,
                   .expect, tolerance = 1e-10)
    }
    ## full 3x3 block
    .check(matrix(c(0.1, 0.01, 0.02,
                    0.01, 0.2, 0.03,
                    0.02, 0.03, 0.3), 3, 3))
    ## partial block: etas 1-2 correlated, eta 3 independent
    .check(matrix(c(0.1, 0.01, 0,
                    0.01, 0.2, 0,
                    0, 0, 0.3), 3, 3))
    ## diagonal: the parameters are omega_kk^(-1/4), the historic closed form
    .om <- diag(c(0.1, 0.2))
    expect_equal(rxode2::rxSymInvCholCreate(mat = .om, diag.xform = "sqrt")$theta,
                 diag(.om)^(-0.25), tolerance = 1e-10)
  })

  test_that(".omegaWriteIni writes blocks and singletons back into a model", {
    .m <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45
        eta.ka ~ 0.2
        eta.cl + eta.v ~ c(0.1,
                           0.01, 0.1)
        add.sd <- 0.7 })
      model({ ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd) })
    }
    .u <- rxode2::rxUiDecompress(rxode2::assertRxUi(.m))
    .om <- matrix(0, 3, 3, dimnames = list(c("eta.ka", "eta.cl", "eta.v"),
                                           c("eta.ka", "eta.cl", "eta.v")))
    diag(.om) <- c(0.33, 0.44, 0.55)
    .om["eta.cl", "eta.v"] <- .om["eta.v", "eta.cl"] <- 0.066
    .u2 <- suppressMessages(.omegaWriteIni(.u, .om))
    expect_equal(.u2$omega[c("eta.ka", "eta.cl", "eta.v"),
                           c("eta.ka", "eta.cl", "eta.v")], .om)
  })
})
