nmTest({
  # The trust inner solve stopped up to sqrt(trustFterm) short of each EBE, and the
  # FOCEi log|H| term moved the objective with it, so the outer objective depended on
  # the warm start.  bobyqa then stopped ~0.03 short on pheno and the full FD
  # sandwich came out with an indefinite R (#1152).
  .pheno <- function() {
    ini({
      tcl <- log(0.008)
      tv <- log(0.6)
      eta.cl ~ 0.1
      eta.v ~ 0.1
      prop.sd <- 0.1
    })
    model({
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      d / dt(central) <- -cl / v * central
      cp <- central / v
      cp ~ prop(prop.sd)
    })
  }

  test_that("trustPolish is opt-in", {
    expect_false(foceiControl()$trustPolish)
    expect_true(foceiControl(trustPolish = TRUE)$trustPolish)
    expect_error(foceiControl(trustPolish = NA))
    .ctl <- foceiControl(trustPolish = TRUE)
    .dep <- rxUiDeparse(.ctl, "ctl")
    expect_true(any(grepl("trustPolish = TRUE", deparse(.dep), fixed = TRUE)))
    expect_true(eval(.dep[[3]])$trustPolish)
  })

  test_that("trustPolish restarts the bobyqa outer search once, within its budget", {
    .fn <- function(x) sum((x - c(1, 2))^2)
    .ctl <- list(npt = 5, rhobeg = 0.2, rhoend = 1e-6, iprint = 0L, maxfun = 1000)
    # a search that stopped short is improved and its evaluations counted
    .stop <- list(par = c(0.5, 1.5), fval = .fn(c(0.5, 1.5)), feval = 40L, ierr = 0L)
    .r <- .bobyqaRestart(.fn, c(-5, -5), c(5, 5), .ctl, .stop)
    expect_lt(.r$fval, 1e-8)
    expect_gt(.r$feval, 40L)
    # no restart after a non-normal exit or with no budget left
    expect_identical(
      .bobyqaRestart(.fn, c(-5, -5), c(5, 5), .ctl, replace(.stop, "ierr", 1L)),
      replace(.stop, "ierr", 1L)
    )
    expect_identical(.bobyqaRestart(.fn, c(-5, -5), c(5, 5), replace(.ctl, "maxfun", 45), .stop), .stop)
    # a restart that does not improve keeps the first result
    .at <- list(par = c(1, 2), fval = 0, feval = 40L, ierr = 0L)
    expect_identical(.bobyqaRestart(.fn, c(-5, -5), c(5, 5), .ctl, .at)$par, c(1, 2))
  })

  test_that("trustPolish=TRUE reaches the pheno minimum with a sane full sandwich (#1152)", {
    skip_on_cran()
    .fit <- suppressMessages(suppressWarnings(nlmixr2(
      .pheno,
      nlmixr2data::pheno_sd,
      est = "focei",
      control = foceiControl(print = 0, calcTables = FALSE, trustPolish = TRUE)
    )))
    .cnt <- .fit$env$nTrustInner
    expect_gt(.cnt[["calls"]], 0L)
    expect_gt(.cnt[["polish"]], 0L)
    # innerOpt="n1qn1" reaches 730.866 here; the default (unpolished) trust fit
    # stops at 730.891 with omega^2 CL 0.122-0.126
    expect_lt(.fit$objf, 730.875)
    expect_lt(.fit$omega["eta.cl", "eta.cl"], 0.117)
    expect_equal(.fit$covMethod, "r,s (full)")
    .R <- solve(.fit$env$.fdFullCov)
    expect_gt(min(eigen(.R, symmetric = TRUE, only.values = TRUE)$values), 0)
    # NONMEM's sandwich SE for omega^2 CL is 0.141; the stopped fit gave 1.16
    expect_lt(sqrt(.fit$cov["om.eta.cl", "om.eta.cl"]), 0.2)
    .ev <- Re(eigen(.fit$env$.fdFullCov %*% .fit$env$.fdFullS, only.values = TRUE)$values)
    expect_lt(max(.ev), 6)
  })

  test_that("the polish stays off by default", {
    skip_on_cran()
    .fit <- suppressMessages(suppressWarnings(nlmixr2(
      .pheno,
      nlmixr2data::pheno_sd,
      est = "focei",
      control = foceiControl(print = 0, calcTables = FALSE, covMethod = "", maxOuterIterations = 2L)
    )))
    expect_gt(.fit$env$nTrustInner[["calls"]], 0L)
    expect_equal(.fit$env$nTrustInner[["polish"]], 0L)
  })
})
