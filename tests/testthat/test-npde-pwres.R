# PWRES (#1390): Sigma_sim^(-1/2) (y - EPRED) from npdeCalc()'s simulations
nmTest({
  # rxode2's per-thread engines exist once something has been solved
  invisible(rxode2::rxSolve(rxode2::rxode2("y = t"), rxode2::et(1), returnType = "data.frame"))

  # K replicates of nid subjects with n observations each, drawn from N(mu, S)
  .pwresSims <- function(mu, S, nid, K) {
    .n <- nrow(S)
    .c <- chol(S)
    .sim <- do.call(
      rbind,
      lapply(seq_len(K), function(k) {
        .z <- matrix(rnorm(nid * .n), nid) %*% .c
        .z <- sweep(.z, 2, mu, "+")
        data.frame(simId = k, id = rep(seq_len(nid), each = .n), sim = as.vector(t(.z)))
      })
    )
    .sim$rxLambda <- 1
    .sim$rxYj <- 2
    .sim$rxLow <- 0
    .sim$rxHi <- 1
    .sim
  }

  .pwres <- function(sim, dv, cens = rep(0L, length(dv)), opt = list()) {
    .r <- .Call(
      `_nlmixr2est_npdeCalc`,
      sim,
      dv,
      rep(0L, length(dv)),
      cens,
      NULL,
      c(list(seed = 7L), opt)
    )
    .r[[2]]$PWRES
  }

  # symmetric inverse square root
  .isqrt <- function(S) {
    .e <- eigen(S, symmetric = TRUE)
    .e$vectors %*% diag(1 / sqrt(.e$values)) %*% t(.e$vectors)
  }

  .mu <- c(1, 2, 3)
  .S <- matrix(c(1, 0.6, 0.3, 0.6, 1.5, 0.5, 0.3, 0.5, 2), 3)
  .y <- c(2.5, 0.5, 4)

  test_that("PWRES approaches the analytic Sigma^(-1/2)(y - mu) as K grows", {
    .true <- as.vector(.isqrt(.S) %*% (.y - .mu))
    set.seed(42)
    .err <- vapply(
      c(100L, 20000L),
      function(K) {
        .p <- .pwres(.pwresSims(.mu, .S, 1L, K), .y)
        expect_false(anyNA(.p))
        max(abs(.p - .true))
      },
      numeric(1)
    )
    expect_lt(.err[2], .err[1])
    expect_lt(.err[2], 0.05)
  })

  test_that("PWRES is the whitened residual of the simulations' own moments", {
    set.seed(3)
    .sim <- .pwresSims(.mu, .S, 1L, 50L)
    .m <- matrix(.sim$sim, nrow = 3)
    .p <- .pwres(.sim, .y)
    expect_equal(.p, as.vector(.isqrt(cov(t(.m))) %*% (.y - rowMeans(.m))), tolerance = 1e-8)
  })

  test_that("a censored subject or K <= n_obs gives NA PWRES", {
    set.seed(11)
    .sim <- .pwresSims(.mu, .S, 2L, 300L)
    .p <- .pwres(.sim, c(.y, .y), cens = c(0L, 0L, 0L, 0L, 1L, 0L))
    expect_false(anyNA(.p[1:3]))
    expect_true(all(is.na(.p[4:6])))
    # the uncensored subject is not changed by the other's censoring
    expect_equal(.p[1:3], .pwres(.sim, c(.y, .y))[1:3])
    # EPRED is still given for the censored subject
    .r <- .Call(
      `_nlmixr2est_npdeCalc`,
      .sim,
      c(.y, .y),
      rep(0L, 6),
      c(0L, 0L, 0L, 0L, 1L, 0L),
      NULL,
      list(seed = 7L)
    )
    expect_false(anyNA(.r[[2]]$EPRED))

    .p <- .pwres(.pwresSims(.mu, .S, 1L, 3L), .y)
    expect_true(all(is.na(.p)))
    .p <- .pwres(.pwresSims(.mu, .S, 1L, 4L), .y)
    expect_false(anyNA(.p))
  })

  test_that("eigen and Cholesky PWRES have the same per-subject sum of squares", {
    set.seed(5)
    .sim <- .pwresSims(.mu, .S, 3L, 200L)
    .dv <- c(.y, .y - 1, .y + 0.5)
    .e <- .pwres(.sim, .dv)
    .c <- .pwres(.sim, .dv, opt = list(pwresMethod = 2L))
    expect_false(isTRUE(all.equal(.e, .c)))
    .ss <- function(p) tapply(p^2, rep(1:3, each = 3), sum)
    expect_equal(.ss(.e), .ss(.c), tolerance = 1e-10)
    # the lower Cholesky factor whitens the first residual by its own SD
    .m <- matrix(.sim$sim[.sim$id == 1], nrow = 3)
    expect_equal(.c[1], (.y[1] - mean(.m[1, ])) / sd(.m[1, ]), tolerance = 1e-10)
    expect_identical(.pwres(.sim, .dv, opt = list(pwresMethod = 1L)), .e)
  })

  test_that("tableControl(pwresMethod=) maps to the C option", {
    expect_identical(tableControl()$pwresMethod, 1L)
    expect_identical(tableControl(pwresMethod = "cholesky")$pwresMethod, 2L)
    expect_error(tableControl(pwresMethod = "svd"))
    # a stored control is passed back through tableControl()
    expect_identical(do.call(tableControl, tableControl(pwresMethod = "cholesky"))$pwresMethod, 2L)
  })

  test_that("fit$EWRES is the PWRES column", {
    .fit <- data.frame(PWRES = c(0.1, -0.2))
    expect_identical(nmObjGetData.EWRES(list(.fit, FALSE)), c(0.1, -0.2))
    expect_null(nmObjGetData.EWRES(list(data.frame(NPDE = 1), FALSE)))
  })
})
