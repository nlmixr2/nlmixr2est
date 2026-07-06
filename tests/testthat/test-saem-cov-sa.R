nmTest({
  test_that("SAEM covMethod='sa' (stochastic-approximation Louis FIM) full covariance", {
    # The dedicated fixed-theta cov phase must (a) leave the estimate unperturbed,
    # (b) produce a PD full theta + Omega + residual covariance, and (c) agree with
    # the linearized FIM.

    one.cmt <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45; add.sd <- 0.7
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v  <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }

    ctlL <- saemControl(nBurn = 200, nEm = 300, print = 0, seed = 1L, covMethod = "linFim")
    ctlS <- saemControl(nBurn = 200, nEm = 300, print = 0, seed = 1L, covMethod = "sa",
                        nSaCov = 1000)

    fL <- .nlmixr(one.cmt, theo_sd, est = "saem", control = ctlL)
    fS <- .nlmixr(one.cmt, theo_sd, est = "saem", control = ctlS)

    # (a) the cov phase is frozen at theta_hat -- the reported estimate is unchanged
    expect_equal(unname(fS$theta), unname(fL$theta), tolerance = 1e-6)

    # (b) full PD covariance, method label preserved
    expect_equal(fS$covMethod, "sa")
    expect_true(is.matrix(fS$cov) && nrow(fS$cov) >= 4L)
    expect_true(all(is.finite(fS$cov)))
    expect_true(min(eigen(fS$cov, symmetric = TRUE, only.values = TRUE)$values) > 0)

    # Omega diagonal and residual rows are present and named on the eta
    expect_true(all(c("om.eta.ka", "om.eta.cl", "om.eta.v", "add.sd") %in% rownames(fS$cov)))

    # (c) SA and linFim SEs agree (different estimators -> allow Monte-Carlo tolerance)
    .cmn <- intersect(rownames(fS$cov), rownames(fL$cov))
    .seS <- sqrt(diag(fS$cov))[.cmn]
    .seL <- sqrt(diag(fL$cov))[.cmn]
    expect_equal(unname(.seS), unname(.seL), tolerance = 0.25)

    # residual SE surfaced in the parameter table
    expect_true(is.finite(fS$parFixedDf["add.sd", "SE"]))
    expect_gt(fS$parFixedDf["add.sd", "SE"], 0)
  })

  test_that("SAEM covMethod='fim' inverts the (mu-block-corrected) estimation-phase FIM", {
    # Regression: the shared per-iteration integrand omits the deterministic mu-block
    # complete Hessian, so before the correction Ha's theta block was indefinite and
    # covMethod='fim' produced NaN SEs.  With the correction fim is a valid PD full
    # covariance agreeing with the linearized FIM.

    one.cmt <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45; add.sd <- 0.7
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v  <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }

    ctlL <- saemControl(nBurn = 200, nEm = 300, print = 0, seed = 1L, covMethod = "linFim")
    ctlF <- saemControl(nBurn = 200, nEm = 300, print = 0, seed = 1L, covMethod = "fim")

    fL <- .nlmixr(one.cmt, theo_sd, est = "saem", control = ctlL)
    fF <- .nlmixr(one.cmt, theo_sd, est = "saem", control = ctlF)

    expect_equal(fF$covMethod, "fim")
    expect_true(is.matrix(fF$cov) && nrow(fF$cov) >= 4L)
    expect_true(all(is.finite(fF$cov)))
    expect_true(min(eigen(fF$cov, symmetric = TRUE, only.values = TRUE)$values) > 0)
    expect_true(all(c("om.eta.ka", "om.eta.cl", "om.eta.v", "add.sd") %in% rownames(fF$cov)))

    .cmn <- intersect(rownames(fF$cov), rownames(fL$cov))
    expect_equal(unname(sqrt(diag(fF$cov))[.cmn]),
                 unname(sqrt(diag(fL$cov))[.cmn]), tolerance = 0.25)

    expect_true(is.finite(fF$parFixedDf["add.sd", "SE"]))
    expect_gt(fF$parFixedDf["add.sd", "SE"], 0)
  })

  test_that("SAEM iteration history records off-diagonal Omega covariances", {
    # block-Omega model: the off-diagonal covariance trajectory (cov.<eta>.<eta>) is
    # recorded in parHistData, named consistently with the covariance matrix, and its
    # final value matches the fitted Omega off-diagonal.
    blk <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45; add.sd <- 0.3; prop.sd <- 0.1
        eta.ka ~ 0.6
        eta.cl + eta.v ~ c(0.3, 0.05, 0.1)
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd) + prop(prop.sd)
      })
    }
    f <- .nlmixr(blk, theo_sd, est = "saem",
                 control = saemControl(nBurn = 150, nEm = 200, print = 0, seed = 1L,
                                       covMethod = "linFim"))
    ph <- f$parHistData
    expect_true("cov.eta.v.eta.cl" %in% names(ph))
    expect_equal(ph[["cov.eta.v.eta.cl"]][nrow(ph)],
                 unname(f$omega["eta.cl", "eta.v"]), tolerance = 0.05)

    # a diagonal-Omega model must not gain any covariance columns (no regression)
    diagM <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; add.sd <- 0.7
            eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1 })
      model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd) })
    }
    fd <- .nlmixr(diagM, theo_sd, est = "saem",
                  control = saemControl(nBurn = 100, nEm = 120, print = 0, seed = 1L))
    expect_false(any(grepl("^cov\\.", names(fd$parHistData))))
  })
})
