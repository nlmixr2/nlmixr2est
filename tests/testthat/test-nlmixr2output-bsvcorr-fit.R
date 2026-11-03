nmTest({
  test_that("a fit puts its non-mu-referenced eta's row before the BSV correlation rows", {
    skip_on_cran()
    .ctl <- foceiControl(print = 0, maxOuterIterations = 0, calcTables = FALSE, covMethod = "")
    fit <- .nlmixr(.bsvCorMixedFit, theo_sd, "focei", control = .ctl)
    .cor <- "bsvcorr(eta.ka,tcl)"
    expect_identical(rownames(fit$parFixedDf), c("tka", "tcl", "tv", "add.sd", "eta.ka", .cor))
    # the model has labels, so the label column describes the row
    .desc <- "BSV correlation between eta.ka and tcl"
    expect_identical(unname(fit$parFixedDf[.cor, "Parameter"]), .desc)
    expect_equal(
      fit$parFixedDf[.cor, "BSV(CV%)"],
      stats::cov2cor(fit$omega)[["eta.ka", "eta.cl"]],
      tolerance = 1e-12
    )
    expect_identical(fit$parFixed[.cor, "Parameter"], .desc)
    expect_identical(fit$parFixed[["Parameter"]][1:4], c("Ka", "", "", ""))
    expect_identical(fit$parFixed[.cor, "BSV(CV%)"], "0.0236")
  })

  test_that("an IOV fit reports only its subject-level BSV correlations", {
    skip_on_cran()
    .d <- theo_sd
    .d$occ <- ifelse(.d$TIME < 6, 1L, 2L)
    .ctl <- foceiControl(print = 0, maxOuterIterations = 0, calcTables = FALSE, covMethod = "", iovMethod = "omega")
    fit <- .nlmixr(.bsvCorIovFit, .d, "focei", control = .ctl)
    .cor <- "bsvcorr(tcl,tv)"
    expect_identical(rownames(fit$parFixedDf), c("tka", "tcl", "tv", "add.sd", "iov.cl", "iov.v", .cor))
    expect_identical(rownames(fit$parFixed), rownames(fit$parFixedDf))
    expect_equal(
      fit$parFixedDf[.cor, "BSV(CV% or SD)"],
      stats::cov2cor(fit$omega$id)[["eta.cl", "eta.v"]],
      tolerance = 1e-12
    )
    expect_identical(fit$parFixed[.cor, "BSV(CV% or SD)"], "0.289")
  })

  test_that("a fit with a theta mu-referenced by two correlated etas reports every pair", {
    skip_on_cran()
    .ctl <- foceiControl(print = 0, maxOuterIterations = 0, calcTables = FALSE, covMethod = "")
    fit <- .nlmixr(.bsvCorSharedFit, theo_sd, "focei", control = .ctl)
    .rows <- c("bsvcorr(eta.cl,eta.cl2)", "bsvcorr(eta.cl,tv)", "bsvcorr(eta.cl2,tv)")
    expect_identical(rownames(fit$parFixedDf), c("tcl", "tv", "add.sd", "eta.cl", "eta.cl2", .rows))
    expect_identical(rownames(fit$parFixed), rownames(fit$parFixedDf))
    .cor <- stats::cov2cor(fit$omega)
    expect_equal(
      unname(fit$parFixedDf[.rows, "BSV(CV%)"]),
      c(.cor[["eta.cl", "eta.cl2"]], .cor[["eta.cl", "eta.v"]], .cor[["eta.cl2", "eta.v"]]),
      tolerance = 1e-12
    )
    expect_identical(unname(fit$parFixed[.rows, "BSV(CV%)"]), c("0.0816", "0.173", "0.283"))
  })

  test_that("a fit reports the copy of a repeated (same()) block as its own correlation", {
    skip_on_cran()
    .ctl <- foceiControl(print = 0, maxOuterIterations = 0, calcTables = FALSE, covMethod = "")
    fit <- .nlmixr(.bsvCorSameFit, theo_sd, "focei", control = .ctl)
    .rows <- c("bsvcorr(tcl,tv)", "bsvcorr(tq,tvp)")
    expect_identical(rownames(fit$parFixedDf), c("tka", "tcl", "tv", "tq", "tvp", "add.sd", .rows))
    .cor <- stats::cov2cor(fit$omega)
    expect_equal(
      unname(fit$parFixedDf[.rows, "BSV(CV%)"]),
      c(.cor[["eta.cl", "eta.v"]], .cor[["eta.q", "eta.vp"]]),
      tolerance = 1e-12
    )
    expect_identical(unname(fit$parFixed[.rows, "BSV(CV%)"]), c("0.289", "0.289"))
  })

  test_that("a saem fit reports its BSV correlation, even one estimated at zero", {
    skip_on_cran()
    .ctl <- saemControl(print = 0, nBurn = 1, nEm = 1, nmc = 1, nu = c(1, 1, 1), calcTables = FALSE)
    fit <- .nlmixr(.bsvCorBlockFit, theo_sd, "saem", control = .ctl)
    .row <- "bsvcorr(tcl,tv)"
    # two iterations leave the covariance at exactly zero; saem's ui, rebuilt
    # from that omega, has no row for it, so the row comes from iniDf0
    expect_identical(fit$omega[["eta.cl", "eta.v"]], 0)
    expect_false("(eta.cl,eta.v)" %in% fit$ui$iniDf$name)
    expect_identical(rownames(fit$parFixedDf), c("tka", "tcl", "tv", "add.sd", .row))
    expect_identical(rownames(fit$parFixed), rownames(fit$parFixedDf))
    .value <- stats::cov2cor(fit$omega)[["eta.cl", "eta.v"]]
    expect_equal(fit$parFixedDf[.row, "BSV(CV%)"], .value, tolerance = 1e-12)
    expect_identical(fit$parFixed[.row, "BSV(CV%)"], formatMinWidth(.value, digits = fit$control$sigdigTable))
    expect_identical(fit$parFixed[.row, "Est."], "")
  })

  test_that("the BSV correlation is printed to the fit's sigdigTable", {
    skip_on_cran()
    .ctl <- foceiControl(print = 0, maxOuterIterations = 0, calcTables = FALSE, covMethod = "", sigdigTable = 5)
    fit <- .nlmixr(.bsvCorBlockFit, theo_sd, "focei", control = .ctl)
    .row <- "bsvcorr(tcl,tv)"
    expect_equal(fit$control$sigdigTable, 5)
    expect_identical(fit$parFixed[.row, "BSV(CV%)"], "0.28868")
    expect_identical(fit$parFixed["tcl", "BSV(CV%)"], "59.149")
  })

  test_that("a saem two-level IOV fit reports only its subject-level BSV correlation", {
    skip_on_cran()
    .d <- theo_sd
    .d$occ <- ifelse(.d$TIME < 6, 1L, 2L)
    .ctl <- saemControl(print = 0, nBurn = 3, nEm = 3, calcTables = FALSE, iovMethod = "twoLevel")
    fit <- .nlmixr(.bsvCorIovSaemFit, .d, "saem", control = .ctl)
    .cor <- "bsvcorr(tcl,tv)"
    .rn <- rownames(fit$parFixedDf)
    expect_identical(.rn[.parFixedIsBsvCor(.rn)], .cor)
    expect_identical(.rn[length(.rn)], .cor)
    expect_identical(rownames(fit$parFixed), .rn)
    .om <- if (is.list(fit$omega)) fit$omega$id else fit$omega
    expect_equal(
      fit$parFixedDf[.cor, "BSV(CV% or SD)"],
      stats::cov2cor(.om)[["eta.cl", "eta.v"]],
      tolerance = 1e-12
    )
  })
})
