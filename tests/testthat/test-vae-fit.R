## Phase 6: est="vae" result -- the updated (final) covariate model with exact
## centered expressions, the original model, the linearization -2LL, the
## linearization-Hessian covariance, and encoder EBEs. (The standard
## nlmixr2FitData wrapper is the remaining Phase 6 wiring.) Slow; skipped on CRAN.

test_that("est=\"vae\" returns the updated covariate model + linearization -2LL/cov", {
  skip_on_cran()
  theo <- function() {
    ini({
      lka <- log(1.8); lke <- log(0.086); lV <- log(32)
      eta.ka ~ 0.3; eta.ke ~ 0.03; eta.V ~ 0.03
      add.err <- 0.7
    })
    model({
      ka <- exp(lka + eta.ka); ke <- exp(lke + eta.ke); V <- exp(lV + eta.V)
      d/dt(depot) = -ka * depot
      d/dt(central) = ka * depot - ke * central
      cp <- central / V
      cp ~ add(add.err)
    })
  }
  f <- nlmixr2(theo, nlmixr2data::theo_sd, est = "vae",
               control = vaeControl(itersBurnIn = 80L, klWarmup = 40L, gammaIter = 120L,
                                    iters = 160L, hiddenDim = 25L, seed = 1L))
  expect_s3_class(f, "vaeFit")

  ## final model carries the exact centered covariate expression (WT on ka and V)
  .lines <- vapply(f$finalUi$lstExpr, function(e) deparse1(e), character(1))
  expect_true(any(grepl("beta_lka_WT \\* log\\(WT/", .lines)))
  expect_true(any(grepl("beta_lV_WT \\* log\\(WT/", .lines)))
  ## ke unmodified (WT not selected on ke)
  expect_true(any(grepl("^ke <- exp\\(lke \\+ eta.ke\\)$", .lines)))
  ## original model preserved and structurally different (no covariates)
  .oLines <- vapply(f$uiIni$lstExpr, function(e) deparse1(e), character(1))
  expect_false(any(grepl("beta_lka_WT", .oLines)))

  ## linearization -2LL near the paper (~331; loose band for the short schedule)
  expect_true(is.finite(f$objective))
  expect_lt(abs(f$objective - 331) / 331, 0.10)
  ## linearization-Hessian covariance: finite positive SEs, or NULL if the
  ## short-schedule Hessian is not positive definite (hardening)
  if (!is.null(f$cov)) expect_true(all(is.finite(sqrt(diag(f$cov))) & diag(f$cov) >= 0))
  ## encoder EBEs present
  expect_true(all(c("eta.ka", "eta.ke", "eta.V") %in% names(f$eta)))
})
