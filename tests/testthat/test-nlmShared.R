# The nlm-family covariance from the optimizer's Hessian (R/nlmShared.R).

test_that("a positive-definite Hessian is inverted as is", {
  .h <- matrix(c(4, 1, 1, 3), 2)
  .r <- .nlmCovFromHessian(.h)
  expect_identical(.r$type, "r")
  expect_identical(.r$r, .h)
  expect_null(.r$warning)
})

test_that("an indefinite Hessian is repaired as |r|, else as the nearest positive-definite matrix", {
  .h <- matrix(c(1, 2, 2, 1), 2) # eigenvalues 3, -1
  .r <- .nlmCovFromHessian(.h)
  expect_identical(.r$type, "|r|")
  expect_equal(.r$r, sqrtm(.h %*% .h))
  expect_equal(eigen(.r$r, symmetric = TRUE, only.values = TRUE)$values, c(3, 1))
  expect_identical(.r$warning, "R matrix is not positive definite; corrected as \"|r|\"")
  # a singular |R| cannot be inverted either
  .h <- diag(c(2, 0, -1))
  .r <- .nlmCovFromHessian(.h)
  expect_identical(.r$type, "r+")
  expect_true(min(eigen(.r$r, symmetric = TRUE, only.values = TRUE)$values) > 0)
  expect_equal(.r$r, nmNearPD(.h))
  expect_identical(.r$warning, "R matrix is not positive definite; corrected as \"r+\"")
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
  # match.arg() cannot match "": every one of them stopped
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
    # reports; only est = "nlm" passed the warnings of its run on
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
    # and the covariance of u is J Cov(x) J with J = diag(init / scaleTo)
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
