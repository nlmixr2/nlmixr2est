## Colinearity clusters for the VAE covariate M-step (R/vaeCovColinear.R).
## A cluster NAMES interchangeable covariates; it never constrains what may be
## selected.  These tests pin the clustering rule, the group/block exclusion and
## the gate that decides whether the vector is sent to C++ at all.  Pure linear
## algebra, no ODE solve: essential.

nmTest({

  ## three groups: cols 1-2 are two shapes of one covariate, cols 3 and 4 are
  ## separate covariates
  .design <- function(N = 80L, seed = 3L) {
    .testSeed(seed)
    X <- matrix(rnorm(N * 4L), N, 4L)
    X[, 2] <- X[, 1] * 0.99 + rnorm(N, sd = 0.05)   # shape mate of col 1
    colnames(X) <- c("A_power", "A_lin", "B", "C")
    X
  }
  .grp <- c(1L, 1L, 2L, 3L)

  test_that("an uncorrelated design leaves every group its own cluster", {
    X <- .design()
    cl <- nlmixr2est:::.vaeCovCluster(X, .grp)
    ## the identity coarsening: cluster and group induce the same partition
    expect_identical(cl, match(.grp, unique(.grp)))
    expect_false(nlmixr2est:::.vaeClusterBinds(cl, .grp))
  })

  test_that("two shapes of one covariate never cluster, however correlated", {
    X <- .design()
    ## cols 1 and 2 correlate at ~0.999 and must still not bind: they share a
    ## group, and the mutual-exclusion machinery already arbitrates them
    expect_gt(abs(stats::cor(X[, 1], X[, 2])), 0.98)
    expect_false(nlmixr2est:::.vaeClusterBinds(
      nlmixr2est:::.vaeCovCluster(X, .grp), .grp))
  })

  test_that("a colinear pair of DIFFERENT covariates merges their groups", {
    X <- .design()
    X[, 4] <- X[, 3] * 0.97 + rnorm(nrow(X), sd = 0.1)  # B ~ C
    cl <- nlmixr2est:::.vaeCovCluster(X, .grp)
    ## groups 2 and 3 land in one cluster; group 1 stays alone
    expect_identical(length(unique(cl)), 2L)
    expect_identical(cl[3], cl[4])
    expect_false(cl[1] == cl[3])
    expect_true(nlmixr2est:::.vaeClusterBinds(cl, .grp))
  })

  test_that("the cluster is always a coarsening of the group", {
    X <- .design()
    X[, 4] <- X[, 3] * 0.97 + rnorm(nrow(X), sd = 0.1)
    cl <- nlmixr2est:::.vaeCovCluster(X, .grp)
    ## every group lies entirely inside one cluster -- the interface guarantee
    ## the C++ side relies on
    expect_true(all(vapply(split(cl, .grp),
                           function(z) length(unique(z)) == 1L, logical(1))))
  })

  test_that("single linkage chains a~b~c into one component", {
    .testSeed(4L)
    N <- 200L
    X <- matrix(rnorm(N * 3L), N, 3L)
    X[, 2] <- X[, 1] * 0.80 + rnorm(N, sd = 0.60)
    X[, 3] <- X[, 2] * 0.80 + rnorm(N, sd = 0.60)
    ## a~b and b~c clear the cut while a~c does not (the correlation compounds
    ## to ~0.64 across the chain): the documented single-linkage behavior,
    ## harmless because a cluster is a label rather than a constraint
    cut <- 0.70
    expect_gt(abs(stats::cor(X[, 1], X[, 2])), cut)
    expect_gt(abs(stats::cor(X[, 2], X[, 3])), cut)
    expect_lt(abs(stats::cor(X[, 1], X[, 3])), cut)
    expect_identical(length(unique(nlmixr2est:::.vaeCovCluster(X, 1:3, cut))), 1L)
  })

  test_that("clustering is invariant to shifting and scaling a column", {
    X <- .design()
    X[, 4] <- X[, 3] * 0.97 + rnorm(nrow(X), sd = 0.1)
    ref <- nlmixr2est:::.vaeCovCluster(X, .grp)
    Y <- X
    Y[, 3] <- Y[, 3] * 1000 + 17      # affine, as a re-centering would be
    Y[, 4] <- Y[, 4] - 4
    expect_identical(nlmixr2est:::.vaeCovCluster(Y, .grp), ref)
  })

  test_that("degenerate designs fall back to the groups without warning", {
    X <- .design()
    id <- match(.grp, unique(.grp))
    expect_identical(nlmixr2est:::.vaeCovCluster(X[1:2, ], .grp), id)   # nrow < 3
    expect_identical(nlmixr2est:::.vaeCovCluster(matrix(0, 80L, 0L), integer(0)),
                     integer(0))
    ## a constant column correlates with nothing, and must not reach stats::cor
    Z <- X; Z[, 4] <- 5
    expect_silent(zc <- nlmixr2est:::.vaeCovCluster(Z, .grp))
    expect_identical(zc, id)
    ## a non-finite entry is a fallback, not an error
    W <- X; W[1, 1] <- NA_real_
    expect_identical(nlmixr2est:::.vaeCovCluster(W, .grp), id)
  })

  test_that("a mismatched group length is an error, not silent recycling", {
    expect_error(nlmixr2est:::.vaeCovCluster(.design(), 1:3),
                 "one entry per covariate column")
  })

  test_that("clusterBinds is not anyDuplicated -- the theo_sd regression", {
    skip_if_not_installed("nlmixr2data")
    res <- vaeCovariates(nlmixr2data::theo_sd)
    ## theo_sd offers four shape columns of ONE covariate, so cluster ids repeat
    expect_gt(anyDuplicated(res$cluster), 0L)
    ## ...but nothing merged two groups, so the vector must NOT be shipped.
    ## Gating on anyDuplicated() would switch the new mechanisms on here.
    expect_false(nlmixr2est:::.vaeClusterBinds(res$cluster, res$group))
    expect_identical(length(unique(res$group)), 1L)
  })

  test_that("vaeCovariates reports the cluster and honors colinearCut", {
    ## both covariates must be constant WITHIN subject or they are excluded as
    ## time-varying and never reach the search at all
    .wt <- c(60, 70, 80, 90, 65, 75, 85, 95)
    .lbm <- .wt * 0.8 + c(1, -1, 1, -1, 1, -1, 1, -1) * 0.05
    d <- data.frame(id = rep(1:8, each = 2), time = rep(0:1, 8),
                    dv = rnorm(16),
                    wt = rep(.wt, each = 2), lbm = rep(.lbm, each = 2))
    res <- vaeCovariates(d, warn = FALSE)
    expect_true("cluster" %in% names(res))
    expect_true(nlmixr2est:::.vaeClusterBinds(res$cluster, res$group))
    ## a cut of 1 admits only exact duplicates, so nothing binds
    res1 <- vaeCovariates(d, warn = FALSE, colinearCut = 1)
    expect_false(nlmixr2est:::.vaeClusterBinds(res1$cluster, res1$group))
    expect_error(vaeCovariates(d, warn = FALSE, colinearCut = 1.5))
  })

  test_that("the near-tie note fits the runInfo one-line budget", {
    expect_identical(nlmixr2est:::.vaeColinearMsg(0L), character(0))
    expect_identical(nlmixr2est:::.vaeColinearMsg(NA_integer_), character(0))
    msg <- nlmixr2est:::.vaeColinearMsg(2L)
    expect_gt(length(msg), 0L)
    expect_true(all(nchar(msg) <= 75L))
    ## no method-name prefix: the fit already reports which method ran
    expect_false(any(grepl("^vae", msg)))
  })

  test_that("alternates pair each mate with the covariate it stands in for", {
    ## C++ marks only the MATE; the covariate it is an alternative to is the
    ## selected column of the same cluster on that dimension
    etaNames <- c("eta.cl", "eta.v")
    covNames <- c("WT_power", "LBM_power", "AGE_power")
    cluster <- c(1L, 1L, 2L)
    selected <- matrix(c(TRUE, FALSE, FALSE,
                         FALSE, FALSE, TRUE), nrow = 2L, byrow = TRUE)
    tie <- matrix(0L, 2L, 3L); tie[1L, 2L] <- 1L      # LBM near-ties WT on eta.cl
    gap <- matrix(0, 2L, 3L);  gap[1L, 2L] <- 0.02
    got <- nlmixr2est:::.vaeColinearAlt(list(altTie = tie, altGap = gap),
                                        selected, cluster, etaNames, covNames)
    expect_identical(nrow(got), 1L)
    expect_identical(got$param, "eta.cl")
    expect_identical(got$covariate, "WT_power")
    expect_identical(got$alternate, "LBM_power")
    expect_equal(got$delta, 0.02)
  })

  test_that("an empty alternates table is a zero-row frame, never NULL", {
    e <- nlmixr2est:::.vaeColinearAlt(NULL, matrix(FALSE, 1L, 1L), 1L, "eta", "WT")
    expect_s3_class(e, "data.frame")
    expect_identical(nrow(e), 0L)
    expect_identical(names(e), c("param", "covariate", "alternate", "delta"))
    ## a marked mate with no selected cluster-mate is dropped rather than
    ## producing an NA row
    tie <- matrix(1L, 1L, 1L)
    z <- nlmixr2est:::.vaeColinearAlt(list(altTie = tie, altGap = matrix(0, 1L, 1L)),
                                      matrix(FALSE, 1L, 1L), 1L, "eta", "WT")
    expect_identical(nrow(z), 0L)
  })

  test_that("vaeControl exposes the colinearity knobs with measured defaults", {
    expect_true(vaeControl()$covSelectColinear)
    expect_equal(vaeControl()$covSelectColinearCut, 0.9)
    ## the accessor default must not drift from the control default
    expect_equal(formals(vaeCovariates)$colinearCut,
                 formals(vaeControl)$covSelectColinearCut)
    expect_false(vaeControl(covSelectColinear = FALSE)$covSelectColinear)
    expect_error(vaeControl(covSelectColinearCut = 1.5))
    expect_error(vaeControl(covSelectColinearCut = -0.1))
    expect_error(vaeControl(covSelectColinear = "yes"))
    expect_equal(vaeControl()$covSelectHysteresis, 0.25)
    expect_equal(vaeControl()$covSelectAltTol, 0.1)
    expect_error(vaeControl(covSelectHysteresis = -1))
    expect_error(vaeControl(covSelectAltTol = -1))
  })
})
