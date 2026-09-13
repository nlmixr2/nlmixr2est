## The copula correlation update, as a formula rather than as a fit outcome.
##
## The M-steps used to disagree about this: saem estimated the correlation from
## the standardized sufficient statistic of the RAW latents, while imp and the
## FOCEi family took a product-moment correlation of the COMBINED ones.  Those
## are not two implementations of one estimator -- the combined latent
##
##     w_k = rho*z_j + sqrt(1 - rho^2)*z_k
##
## is built FROM the current rho, so correlating it against w_j returns the
## value it was handed whenever the latent second moments happen to be equal.
## A fixed point at the current estimate rather than at the data's, which is
## how a declared fit could report its copula still sitting on its ini() value
## after the M-step had "run" a dozen times.
##
## `rxEtaDistCorFromRz()` is now the single copy, shared by all three.

nmTest({

  test_that("rz = 0 leaves the correlation exactly where it was", {
    ## THE property that makes this an M-step rather than a random walk: S_z == I
    ## means the current rho already explains the draws, so there is nothing to
    ## move.  A route without this fixed point drifts every time it fires.
    for (.r in c(-0.9, -0.5, 0, 0.3, 0.5, 0.9)) {
      expect_equal(rxEtaDistCorFromRz_(.r, 0), .r, tolerance = 1e-12)
    }
  })

  test_that("residual latent correlation moves rho in its own direction", {
    ## rz > 0 says the latents are still correlated after the current rho was
    ## taken out, so the correlation was too small; rz < 0 the reverse.
    expect_gt(rxEtaDistCorFromRz_(0.3, 0.2), 0.3)
    expect_lt(rxEtaDistCorFromRz_(0.3, -0.2), 0.3)
    expect_gt(rxEtaDistCorFromRz_(-0.3, 0.2), -0.3)
    ## and it is monotone in rz
    .v <- vapply(seq(-0.9, 0.9, by = 0.1),
                 function(z) rxEtaDistCorFromRz_(0.4, z), numeric(1))
    expect_true(all(diff(.v) > 0))
  })

  test_that("the update matches the S_n = L S_z L' definition", {
    ## Independent of the implementation: build L from rho, build S_z as a
    ## correlation matrix with off-diagonal rz, and read the normalized
    ## off-diagonal of L S_z L'.
    .ref <- function(rho, rz) {
      .L <- matrix(c(1, rho, 0, sqrt(1 - rho^2)), 2, 2)
      .Sz <- matrix(c(1, rz, rz, 1), 2, 2)
      .Sn <- .L %*% .Sz %*% t(.L)
      .Sn[2, 1] / sqrt(.Sn[1, 1] * .Sn[2, 2])
    }
    for (.r in c(-0.7, -0.2, 0, 0.25, 0.6, 0.85)) {
      for (.z in c(-0.6, -0.1, 0, 0.15, 0.5)) {
        expect_equal(rxEtaDistCorFromRz_(.r, .z), .ref(.r, .z),
                     tolerance = 1e-10,
                     label = sprintf("rho=%.2f rz=%.2f", .r, .z))
      }
    }
  })

  test_that("the result is clamped where the copula stays identified", {
    ## 0.99, not 0.999: at the boundary the partner's latent becomes numerically
    ## its partner's, the family M-step's own spread guard then fails, and every
    ## family theta freezes at its starting value.  Recorded on Bauer's g1 at
    ## MARE 41.2% against a baseline of 18.3%.
    expect_lte(rxEtaDistCorFromRz_(0.98, 0.99), 0.99)
    expect_gte(rxEtaDistCorFromRz_(-0.98, -0.99), -0.99)
    expect_true(is.finite(rxEtaDistCorFromRz_(0.999999, 0.999999)))
    ## a non-finite statistic is treated as "no information", i.e. no move
    expect_equal(rxEtaDistCorFromRz_(0.4, NaN), 0.4, tolerance = 1e-12)
    expect_equal(rxEtaDistCorFromRz_(0.4, NA_real_), 0.4, tolerance = 1e-12)
  })

  test_that("it is NOT the product-moment correlation of the combined latents", {
    ## The whole point.  Build the combined latent the way the model does and
    ## correlate it the way the old route did: with EQUAL latent spreads the two
    ## agree, and with unequal ones -- which is the normal state, an
    ## over-dispersed latent being signal here rather than pathology -- they do
    ## not, and the old one is the biased one.
    set.seed(42)
    .n <- 20000L
    .rho <- 0.3
    .zj <- rnorm(.n)
    .zk <- rnorm(.n)
    .w <- function(zj, zk) .rho * zj + sqrt(1 - .rho^2) * zk
    ## equal spreads, no residual correlation: both routes return rho
    expect_equal(rxEtaDistCorTest_(.zj, .w(.zj, .zk)), .rho, tolerance = 0.02)
    expect_equal(rxEtaDistCorFromRz_(.rho, 0), .rho, tolerance = 1e-12)
    ## now inflate ONE latent, still with no residual correlation between them.
    ## Nothing about the copula changed, so a correct update must not move.
    .zjBig <- 3 * .zj
    .rzBig <- sum(.zjBig * .zk) / sqrt(sum(.zjBig^2) * sum(.zk^2))
    expect_equal(rxEtaDistCorFromRz_(.rho, .rzBig), .rho, tolerance = 0.02)
    ## the product-moment route on the combined pair does move, and that drift
    ## is what accumulates into the clamp over an entire fit
    .pm <- rxEtaDistCorTest_(.zjBig, .w(.zjBig, .zk))
    expect_gt(abs(.pm - .rho), 0.1)
  })

  test_that("the two controls default DIFFERENTLY, on measurement", {
    ## focei TRUE: better or equal on all four of Bauer's declared-gamma arms,
    ## and widest exactly where the product-moment route is predicted to fail
    ## (g1, smallest latent spread: rho +0.691 vs +0.540 against truth 0.50).
    expect_true(foceiControl()$etaDistCorSuff)
    expect_false(foceiControl(etaDistCorSuff = FALSE)$etaDistCorSuff)
    ## imp FALSE: the same route leaves imp's correlation accurate but sends the
    ## LOCATION away (g2 CL 1823, g3 CL 306204 against truth 5.10).  The
    ## asymmetry is deliberate; this pins it so neither drifts to match the
    ## other by accident.
    expect_false(impmapControl()$etaDistCorSuff)
    expect_true(impmapControl(etaDistCorSuff = TRUE)$etaDistCorSuff)
    expect_error(impmapControl(etaDistCorSuff = "yes"))
    ## and they round-trip through do.call(), which is how the control is rebuilt
    .c <- impmapControl(etaDistCorSuff = TRUE)
    expect_true(do.call(impmapControl,
                        .c[names(.c) %in% names(formals(impmapControl))])$etaDistCorSuff)
  })

  test_that("focei's default cannot leak into imp", {
    ## impmapControl() builds itself with do.call(foceiControl, .) and THEN
    ## overwrites this field from its own argument (R/impmap.R).  If that
    ## overwrite were ever removed, imp would silently inherit focei's TRUE --
    ## which is the setting measured to send imp's location parameters to
    ## CL 306204.  Assert the overwrite, not just the default.
    expect_false(impmapControl()$etaDistCorSuff)
    expect_false(impmapControl(etaDistCorSuff = FALSE)$etaDistCorSuff)
    ## the C++ reads this off the control list by name (src/inner.cpp ~8488),
    ## so the name must actually be PRESENT and logical, not absent-and-defaulted
    .c <- impmapControl()
    expect_true("etaDistCorSuff" %in% names(.c))
    expect_type(.c$etaDistCorSuff, "logical")
  })

  test_that("the copula write has its own counter", {
    ## _foceiEtaDistN counts an ATTEMPT where anything moved, so a run with a
    ## dozen family fits and zero copula writes is indistinguishable from one
    ## with both -- and those have very different answers.  Measured on Bauer's
    ## g1: 13 M-steps, copula still on its ini() value at the end.
    expect_true(is.numeric(foceiEtaDistCorN_()))
    expect_gte(foceiEtaDistCorN_(), 0)
  })

})
