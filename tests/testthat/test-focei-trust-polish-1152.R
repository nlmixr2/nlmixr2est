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

  test_that("default FOCEi reaches the pheno minimum with a sane full sandwich (#1152)", {
    skip_on_cran()
    .fit <- suppressMessages(suppressWarnings(nlmixr2(
      .pheno,
      nlmixr2data::pheno_sd,
      est = "focei",
      control = foceiControl(print = 0, calcTables = FALSE)
    )))
    .cnt <- .fit$env$nTrustInner
    expect_gt(.cnt[["calls"]], 0L)
    expect_gt(.cnt[["polish"]], 0L)
    # innerOpt="n1qn1" (unaffected by the polish) reaches 730.866 here; the
    # unpolished trust fit stopped at 730.891 with omega^2 CL 0.122-0.126
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
})
