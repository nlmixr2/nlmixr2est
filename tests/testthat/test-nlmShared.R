# The nlm-family covariance from the optimizer's Hessian (R/nlmShared.R).

test_that("a positive-definite Hessian is inverted as is", {
  .h <- matrix(c(4, 1, 1, 3), 2)
  .r <- .nlmCovFromHessian(.h)
  expect_identical(.r$type, "r")
  expect_identical(.r$r, .h)
  expect_null(.r$warning)
  # the factor inverted is the Hessian's own: nothing added
  expect_equal(crossprod(.r$u), .h, tolerance = 1e-14)
})

test_that("a positive-definite but nearly singular Hessian is repaired as \"r+\", with a warning (issue 1140)", {
  # eigenvalues 2 - 1e-7 and 1e-7: positive definite, but Schnabel-Eskow's
  # modified Cholesky adds to the diagonal (as FOCEi's cholSE0 does, "r+")
  .h <- matrix(c(1, 1 - 1e-7, 1 - 1e-7, 1), 2)
  expect_gt(min(eigen(.h, symmetric = TRUE, only.values = TRUE)$values), 0)
  .r <- .nlmCovFromHessian(.h)
  expect_identical(.r$type, "r+")
  expect_identical(.r$warning, "R matrix is nearly singular; corrected as \"r+\"")
  .e <- diag(crossprod(.r$u)) - diag(.h)
  expect_true(all(.e > 0) && all(.e <= foceiControl()$cholAccept))
})

test_that("a Hessian that is not positive definite is repaired as FOCEi repairs R, under its labels (issue 1140)", {
  # nearly positive definite: Schnabel-Eskow's modified Cholesky adds at most
  # cholAccept (2.2e-5 here) to the diagonal, "r+" as in foceiCovUsable()
  .h <- diag(c(2, 1, -1e-5))
  .r <- .nlmCovFromHessian(.h)
  expect_identical(.r$type, "r+")
  expect_identical(.r$r, .h)
  expect_identical(.r$warning, "R matrix is not positive definite; corrected as \"r+\"")
  expect_true(all(diag(crossprod(.r$u)) - diag(.h) <= foceiControl()$cholAccept))
  # otherwise sqrtm(R %*% R), "|r|"
  .h <- matrix(c(1, 2, 2, 1), 2) # eigenvalues 3, -1
  .r <- .nlmCovFromHessian(.h)
  expect_identical(.r$type, "|r|")
  expect_equal(.r$r, sqrtm(.h %*% .h))
  expect_equal(eigen(.r$r, symmetric = TRUE, only.values = TRUE)$values, c(3, 1))
  expect_equal(crossprod(.r$u), .r$r, tolerance = 1e-12)
  expect_identical(.r$warning, "R matrix is not positive definite; corrected as \"|r|\"")
})

test_that("a numerically singular Hessian is repaired as FOCEi's R is", {
  # the shared rule accepts "r+" on a rank-deficient R when cholAccept allows it;
  # otherwise there is no covariance
  .fail <- "R matrix is not positive definite; covariance step failed"
  for (.h in list(diag(c(2, 0, -1)), tcrossprod(1:4))) {
    .r <- .nlmCovFromHessian(.h)
    expect_identical(.r$type, "failed")
    expect_null(.r$r)
    expect_identical(.r$warning, .fail)
  }
  .r <- .nlmCovFromHessian(diag(c(2, -1e-17, 1)))
  expect_identical(.r$type, "r+")
  expect_identical(.r$warning, "R matrix is not positive definite; corrected as \"r+\"")
  # a well-conditioned indefinite Hessian is still repaired
  expect_identical(.nlmCovFromHessian(diag(c(2, -1, 1)))$type, "|r|")
})

test_that("a Hessian that cannot be repaired gives no covariance", {
  .r <- .nlmCovFromHessian(matrix(c(1, NA, NA, 1), 2))
  expect_identical(.r$type, "failed")
  expect_null(.r$r)
  expect_identical(.r$warning, "R matrix is not finite; covariance step failed")
  .r <- .nlmCovFromHessian(matrix(c(1, Inf, Inf, 1), 2))
  expect_identical(.r$warning, "R matrix is not finite; covariance step failed")
  # the zero-filled Hessian of a failed solve
  .r <- .nlmCovFromHessian(matrix(0, 2, 2))
  expect_identical(.r$type, "failed")
  expect_identical(.r$warning, "R matrix is not positive definite; covariance step failed")
})

test_that("every nlm-family control takes covMethod = \"\" (issue 1140)", {
  # "" skips the covariance step; match.arg() cannot match it
  for (.f in c(
    "nlmControl",
    "nlminbControl",
    "optimControl",
    "bobyqaControl",
    "newuoaControl",
    "uobyqaControl",
    "n1qn1Control",
    "lbfgsb3cControl",
    "trustControl"
  )) {
    expect_identical(get(.f)(covMethod = "")$covMethod, "", info = .f)
    expect_identical(get(.f)(covMethod = "r")$covMethod, "r", info = .f)
  }
  expect_identical(nlmControl()$covMethod, "nlm")
  expect_identical(nlmControl(solveType = "fun")$covMethod, "r")
  expect_identical(bobyqaControl()$covMethod, "r")
  expect_error(nlmControl(covMethod = "s"))
})

# sensMethod of the foceiControl an nlm-family control finalizes with
.nlmFinalSensMethod <- function(m, ...) {
  .env <- new.env(parent = emptyenv())
  assign(paste0(m, "Control"), do.call(paste0(m, "Control"), list(...)), envir = .env)
  get(paste0(".", m, "ControlToFoceiControl"))(.env, assign = FALSE)$sensMethod
}

test_that("the finalization foceiControl keeps the sensMethod of every nlm-family control", {
  for (.m in c("nlm", "nlminb", "optim", "lbfgsb3c", "n1qn1")) {
    expect_identical(.nlmFinalSensMethod(.m, sensMethod = "forward"), "forward", info = .m)
    expect_identical(.nlmFinalSensMethod(.m), "default", info = .m)
  }
  # controls without a sensMethod finalize with the foceiControl() default
  for (.m in c("bobyqa", "newuoa", "uobyqa", "trust", "nls")) {
    expect_identical(.nlmFinalSensMethod(.m), "default", info = .m)
  }
})

nmTest({
  .pk <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl)
      v <- exp(tv)
      d / dt(depot) <- -ka * depot
      d / dt(centr) <- ka * depot - cl / v * centr
      cp <- centr / v
      cp ~ add(add.sd)
    })
  }

  test_that("a derivative-free fit with an indefinite Hessian gets a repaired positive-definite covariance", {
    # at the covariance tolerances this model's Hessian is positive definite, so
    # its smallest eigenvalue is flipped
    .hess <- nlmixr2Hess
    local_mocked_bindings(nlmixr2Hess = function(...) {
      .e <- eigen(.hess(...), symmetric = TRUE)
      .v <- .e$values
      .v[length(.v)] <- -abs(.v[length(.v)])
      .e$vectors %*% diag(.v) %*% t(.e$vectors)
    })
    .fit <- suppressMessages(nlmixr2(.pk, nlmixr2data::theo_sd, est = "bobyqa", control = bobyqaControl(print = 0L)))
    .h <- .fit$env$bobyqa$r
    expect_lt(min(eigen(.h, symmetric = TRUE, only.values = TRUE)$values), 0)
    expect_identical(.fit$covMethod, "|r|")
    expect_gt(min(eigen(.fit$cov, symmetric = TRUE, only.values = TRUE)$values), 0)
    expect_equal(unname(.fit$env$bobyqa$cov.scaled), unname(solve(sqrtm(.h %*% .h))), tolerance = 1e-8)
    expect_true("R matrix is not positive definite; corrected as \"|r|\"" %in% .fit$runInfo)
    expect_null(.fit$env$bobyqa$covWarning)
  })

  test_that("an nlm-family fit with covMethod = \"\" computes no covariance (issue 1140)", {
    .fit <- .nlmixr(.pk, nlmixr2data::theo_sd, est = "n1qn1", control = n1qn1Control(print = 0L, covMethod = ""))
    expect_null(.fit$cov)
    expect_null(.fit$env$n1qn1$r)
  })

  test_that("the warnings of every nlm-family run reach $runInfo (issue 1140)", {
    # censored observations are finite-differenced, which nlmWarnings()
    # reports during the run
    .d <- nlmixr2data::theo_sd
    .d$CENS <- ifelse(.d$DV < 2 & .d$EVID == 0, 1L, 0L)
    .d$DV[.d$CENS == 1] <- 2
    for (.est in c("nlm", "nlminb", "n1qn1")) {
      .fit <- .nlmixr(.pk, .d, est = .est, control = list(print = 0L))
      expect_true(
        "NaN symbolic gradients were resolved with finite differences" %in% .fit$runInfo,
        info = .est
      )
    }
  })

  test_that("the covariance is mapped to the natural scale by the diagonal scaling Jacobian (issue 1140)", {
    # scaleType = "mult" estimates x = u * scaleTo / init, so du/dx = init / scaleTo
    # and the covariance of u is J Cov(x) J with J = diag(init / scaleTo).  A
    # characterization test: J's zero off-diagonal is explicit, and Armadillo
    # (>= 10.5) also zero-fills, so it cannot tell the two apart.
    .x <- nlmObjectiveSetup(
      .pk,
      nlmixr2data::theo_sd,
      control = nlmControl(print = 0L, scaleType = "mult", scaleTo = 2)
    )
    on.exit(.nlmFreeEnv())
    .init <- c(0.45, 1, 3.45, 0.7)
    .n <- c("tka", "tcl", "tv", "add.sd")
    .cov <- matrix(c(4, 1, 0.5, 0.2, 1, 3, 0.1, 0.3, 0.5, 0.1, 2, 0.4, 0.2, 0.3, 0.4, 1), 4, dimnames = list(.n, .n))
    expect_equal(.nlmAdjustCov(.cov, .x), .cov * tcrossprod(.init / 2), tolerance = 1e-14)
  })
})

test_that("covAccept_() takes every branch of the covariance acceptance rule", {
  # the one rule FOCEi and the nlm family use: "" as it is, "+" within cholAccept
  # (a rank-deficient matrix included), "|" for sqrtm(A %*% A)
  .tol <- .Machine$double.eps^(1 / 3)
  .cases <- list(
    pd = list(diag(2), ""),
    one = list(matrix(2), ""),
    zero1 = list(matrix(0), "failed"),
    neg1 = list(matrix(-1), "|"),
    nearPd = list(matrix(c(1, 1, 1, 1 - 1e-8), 2), "+"),
    indefinite = list(matrix(c(1, 0, 0, -1), 2), "|"),
    rank1 = list(matrix(c(1, 1, 1, 1), 2), "+"),
    zero2 = list(matrix(0, 2, 2), "failed"),
    nan = list(matrix(c(1, NaN, NaN, 1), 2), "failed"),
    inf = list(matrix(c(Inf, 0, 0, 1), 2), "failed")
  )
  for (.n in names(.cases)) {
    .c <- .cases[[.n]]
    expect_identical(covAccept_(.c[[1]], .tol, 1e-3)$type, .c[[2]], label = .n)
  }
  # "|" factors sqrtm(A %*% A) and returns it; "+" keeps A and factors A + diag(E)
  .a <- covAccept_(matrix(c(1, 0, 0, -1), 2), .tol, 1e-3)
  expect_equal(.a$M, sqrtm(matrix(c(1, 0, 0, -1), 2) %*% matrix(c(1, 0, 0, -1), 2)))
  expect_equal(crossprod(.a$U), .a$M)
  .p <- covAccept_(matrix(c(1, 1, 1, 1 - 1e-8), 2), .tol, 1e-3)
  expect_identical(.p$M, matrix(c(1, 1, 1, 1 - 1e-8), 2))
  # cholAccept bounds the "+" correction; past it the matrix is "|"
  .small <- diag(c(1, -1e-6))
  expect_identical(covAccept_(.small, .tol, 1e-3)$type, "+")
  expect_identical(covAccept_(.small, .tol, 1e-9)$type, "|")
})
