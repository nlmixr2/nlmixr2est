## Unit tests for the shared omega-block helpers (R/omegaBlock.R) and for the
## omega parameter ORDER the vae inner fast path relies on.  Fast (no fits), so
## these stay in the essential push/PR subset -- the fit-based cross-method
## checks live in test-omega-offdiag.R (weekly batch).

nmTest({
  test_that(".omegaBlockFromIniDf builds the declared block + fix mask", {
    .idf <- data.frame(
      name = c("tka", "eta.cl", "(eta.cl,eta.v)", "eta.v", "eta.q"),
      ntheta = c(1, NA, NA, NA, NA),
      neta1 = c(NA, 1, 2, 2, 3),
      neta2 = c(NA, 1, 1, 2, 3),
      est = c(0.4, 0.1, 0.01, 0.2, 0.3),
      fix = c(FALSE, FALSE, FALSE, FALSE, TRUE),
      stringsAsFactors = FALSE
    )
    .b <- .omegaBlockFromIniDf(.idf, c("eta.cl", "eta.v", "eta.q"))
    expect_equal(unname(.b$mat[1L, 2L]), 0.01)
    expect_equal(unname(.b$mat[2L, 1L]), 0.01) # symmetric
    expect_equal(unname(diag(.b$mat)), c(0.1, 0.2, 0.3))
    expect_equal(unname(.b$mat[1L, 3L]), 0) # undeclared stays 0
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
      .pos <- which(.sel, arr.ind = TRUE) # column-major
      .u <- chol(solve(om))
      .expect <- vapply(
        seq_len(nrow(.pos)),
        function(k) {
          .i <- .pos[k, 1L]
          .j <- .pos[k, 2L]
          if (.i == .j) sqrt(.u[.i, .i]) else .u[.i, .j]
        },
        numeric(1)
      )
      expect_equal(rxode2::rxSymInvCholCreate(mat = om, diag.xform = "sqrt")$theta, .expect, tolerance = 1e-10)
    }
    ## full 3x3 block
    .check(matrix(c(0.1, 0.01, 0.02, 0.01, 0.2, 0.03, 0.02, 0.03, 0.3), 3, 3))
    ## partial block: etas 1-2 correlated, eta 3 independent
    .check(matrix(c(0.1, 0.01, 0, 0.01, 0.2, 0, 0, 0, 0.3), 3, 3))
    ## diagonal: the parameters are omega_kk^(-1/4), the historic closed form
    .om <- diag(c(0.1, 0.2))
    expect_equal(rxode2::rxSymInvCholCreate(mat = .om, diag.xform = "sqrt")$theta, diag(.om)^(-0.25), tolerance = 1e-10)
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
    .om <- matrix(0, 3, 3, dimnames = list(c("eta.ka", "eta.cl", "eta.v"), c("eta.ka", "eta.cl", "eta.v")))
    diag(.om) <- c(0.33, 0.44, 0.55)
    .om["eta.cl", "eta.v"] <- .om["eta.v", "eta.cl"] <- 0.066
    .u2 <- suppressMessages(.omegaWriteIni(.u, .om))
    expect_equal(.u2$omega[c("eta.ka", "eta.cl", "eta.v"), c("eta.ka", "eta.cl", "eta.v")], .om)
  })

  test_that(".omegaBlockIds finds the connected components", {
    ## fully diagonal -> one block per eta
    expect_equal(.omegaBlockIds(diag(3)), c(1L, 2L, 3L))
    ## 4x4 as two 2x2 blocks
    .m <- matrix(0, 4, 4)
    diag(.m) <- 1
    .m[1, 2] <- .m[2, 1] <- 0.1
    .m[3, 4] <- .m[4, 3] <- 0.1
    expect_equal(.omegaBlockIds(.m), c(1L, 1L, 2L, 2L))
    ## a chain 1-2-3 is ONE block even though (1,3) is zero
    .c <- matrix(c(1, .1, 0, .1, 1, .1, 0, .1, 1), 3, 3)
    expect_equal(.omegaBlockIds(.c), c(1L, 1L, 1L))
  })

  test_that(".omegaBlockZeros names the zeros rxSymInvCholCreate cannot hold", {
    ## acceptable patterns have none
    expect_equal(nrow(.omegaBlockZeros(diag(3))), 0L)
    .two <- matrix(0, 4, 4)
    diag(.two) <- 1
    .two[1, 2] <- .two[2, 1] <- 0.1
    .two[3, 4] <- .two[4, 3] <- 0.1
    expect_equal(nrow(.omegaBlockZeros(.two)), 0L)
    ## the rxode2#1365 matrix: (2,3) is zero INSIDE the 1-2-3 block
    .bad <- matrix(c(1, .1, .1, .1, 1, 0, .1, 0, 1), 3, 3)
    expect_equal(unname(.omegaBlockZeros(.bad)), matrix(c(2L, 3L), 1, 2))
    ## and that is exactly the matrix rxSymInvCholCreate cannot hold as given
    expect_true(.rxInvPatternMismatch(.bad))
    ## a NON-CONTIGUOUS component cannot be held either, even though every
    ## component is dense: eta1 correlates with eta3 and eta2 sits between
    ## them.  The whole 1..3 span has to be filled, not just the component.
    .gap <- matrix(c(1, 0, .5, 0, 1, 0, .5, 0, 1), 3, 3)
    expect_true(.rxInvPatternMismatch(.gap))
    expect_equal(nrow(.omegaBlockZeros(.gap)), 2L)
  })

  test_that(".omegaBlockZeros matches rxSymInvCholCreate on EVERY 4x4 pattern", {
    ## The predicate is the whole fix: it decides whether an omega needs
    ## filling before the call, so it must agree with the call itself rather
    ## than with a plausible story about it.  Enumerate every off-diagonal
    ## zero pattern on 4 etas and check both directions, plus that the fill
    ## gives a matrix the call parameterizes by its own pattern.
    .pairs <- which(upper.tri(diag(4)), arr.ind = TRUE)
    .nClean <- 0L
    .nFill <- 0L
    for (.b in 0:63) {
      .bits <- as.integer(intToBits(.b))[1:6]
      .m <- diag(4)
      for (.k in which(.bits == 1)) {
        .i <- .pairs[.k, 1]
        .j <- .pairs[.k, 2]
        .m[.i, .j] <- .m[.j, .i] <- 0.15
      }
      .ok <- !.rxInvPatternMismatch(.m)
      expect_equal(nrow(.omegaBlockZeros(.m)) == 0L, .ok, info = paste("pattern", .b))
      if (.ok) {
        .nClean <- .nClean + 1L
      } else {
        .nFill <- .nFill + 1L
        .f <- .omegaFillBlockZeros(.m)
        expect_false(is.null(.f), info = paste("pattern", .b))
        expect_false(.rxInvPatternMismatch(.f), info = paste("pattern", .b))
      }
    }
    ## the sweep really covered both outcomes
    expect_equal(.nClean, 8L)
    expect_equal(.nFill, 56L)
  })

  test_that(".omegaFillBlockZeros makes the pattern acceptable", {
    .bad <- matrix(c(1, .1, .1, .1, 1, 0, .1, 0, 1), 3, 3, dimnames = list(c("a", "b", "c"), c("a", "b", "c")))
    .fill <- .omegaFillBlockZeros(.bad)
    expect_false(is.null(.fill))
    ## only the offending cell moved, and it moved by a negligible amount
    expect_true(.fill[2, 3] > 0)
    expect_true(.fill[2, 3] < 1e-9)
    ## every other cell is untouched
    .chk <- .fill
    .chk[2, 3] <- .chk[3, 2] <- 0
    expect_equal(.chk, .bad)
    expect_equal(nrow(.omegaBlockZeros(.fill)), 0L)
    ## the mechanism: the filled matrix is parameterized by its own pattern,
    ## the full dense-block parameter count
    .r <- rxode2::rxSymInvCholCreate(mat = .fill, diag.xform = "sqrt")
    expect_equal(length(.r$theta), 6L)
    ## nothing to fill -> NULL, so the caller knows this rung does not apply
    expect_null(.omegaFillBlockZeros(diag(3)))
    ## a non-positive diagonal cannot be scaled into a covariance
    .zero <- .bad
    .zero[2, 2] <- 0
    expect_null(.omegaFillBlockZeros(.zero))
  })

  test_that(".omegaBlockZeroNames names the random effects, truncated", {
    .bad <- matrix(
      c(1, .1, .1, .1, 1, 0, .1, 0, 1),
      3,
      3,
      dimnames = list(c("eta.cl", "eta.v", "eta.ka"), c("eta.cl", "eta.v", "eta.ka"))
    )
    expect_equal(.omegaBlockZeroNames(.bad, .omegaBlockZeros(.bad)), "eta.v, eta.ka")
    expect_true(nchar(.omegaBlockZeroNames(.bad, .omegaBlockZeros(.bad), width = 8L)) <= 8L)
  })
})
