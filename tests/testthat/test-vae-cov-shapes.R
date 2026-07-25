nmTest({
  test_that("shapes collapse to the two searchable families", {
    expect_equal(.vaeShapeFamily(c("power", "log")), c("log", "log"))
    expect_equal(.vaeShapeFamily(c("lin", "identity", "center")),
                 c("lin", "lin", "lin"))
    expect_equal(.vaeShapeFamily("cat"), "cat")
    expect_error(.vaeShapeFamily("hockey"), "unknown covariate shape")
  })

  test_that("shape expressions match the nlmixr2scm vocabulary", {
    expect_equal(.vaeShapeExpr("power", "WT", 70.5), "log(WT/70.5)")
    expect_equal(.vaeShapeExpr("log", "WT", 70.5), "log(WT)")
    expect_equal(.vaeShapeExpr("lin", "WT", 70.5), "(WT - 70.5)")
    expect_equal(.vaeShapeExpr("identity", "WT", 70.5), "WT")
    expect_equal(.vaeShapeExpr("center", "WT", 70.5), "(WT/70.5)")
    expect_error(.vaeShapeExpr("hockey", "WT", 1), "unknown covariate shape")
  })

  test_that("categorical expressions compare against the level", {
    expect_equal(.vaeShapeExpr("cat", "SEX", level = "F"), "(SEX == \"F\")")
    expect_equal(.vaeShapeExpr("cat", "GRP", level = 2), "(GRP == 2)")
    ## a level stored as the character "2" still compares numerically
    expect_equal(.vaeShapeExpr("cat", "GRP", level = "2"), "(GRP == 2)")
    ## a bare 0/1 indicator keeps its natural parameterization
    expect_equal(.vaeShapeExpr("cat", "SEXF", raw = TRUE), "SEXF")
  })

  test_that("shape re-parameterization leaves the prediction unchanged", {
    ## fitted family model: mu = ic + b * x, x the family column
    cov <- c(50, 70, 90, 110)
    ctr <- 70
    ic <- 1.25
    b <- 0.8

    ## log family: x = log(cov/ctr)
    x <- log(cov / ctr)
    mu <- ic + b * x
    for (sh in c("power", "log")) {
      r <- .vaeShapeBeta(sh, ctr, b)
      e <- switch(sh, power = log(cov / ctr), log = log(cov))
      expect_equal(ic + r$interceptAdj + r$beta * e, mu, info = sh)
    }

    ## linear family: x = cov - ctr
    x <- cov - ctr
    mu <- ic + b * x
    for (sh in c("lin", "identity", "center")) {
      r <- .vaeShapeBeta(sh, ctr, b)
      e <- switch(sh, lin = cov - ctr, identity = cov, center = cov / ctr)
      expect_equal(ic + r$interceptAdj + r$beta * e, mu, info = sh)
    }
  })

  test_that("a character shapes= vector is one global rule", {
    r <- .vaeResolveShapes(c("power", "lin"))
    expect_equal(nrow(r), 1L)
    expect_true(is.na(r$var) && is.na(r$cov))
    expect_equal(.vaeShapesFor(r, "cl", "WT"), c("power", "lin"))
    ## NULL means every continuous shape
    expect_equal(.vaeShapesFor(.vaeResolveShapes(NULL), "cl", "WT"),
                 .vaeContShapes)
  })

  test_that("a covariate-named list restricts that covariate only", {
    r <- .vaeResolveShapes(list(wt = "power"))
    ## matching is case-insensitive because the VAE upper-cases data columns
    expect_equal(.vaeShapesFor(r, "cl", "WT"), "power")
    ## an unmatched covariate stays free
    expect_equal(.vaeShapesFor(r, "cl", "AGE"), .vaeContShapes)
  })

  test_that("a pairsVec list restricts one parameter/covariate pair", {
    r <- .vaeResolveShapes(list(
      list(var = "cl", covar = "wt", shapes = "power"),
      list(var = "v", covar = "wt", shapes = c("lin", "identity"))))
    expect_equal(.vaeShapesFor(r, "cl", "WT"), "power")
    expect_equal(.vaeShapesFor(r, "v", "WT"), c("lin", "identity"))
    ## shapes= restricts parameterizations only -- a pair no rule mentions is
    ## still searched, with every shape available
    expect_equal(.vaeShapesFor(r, "ka", "WT"), .vaeContShapes)
    expect_equal(.vaeShapesFor(r, "cl", "AGE"), .vaeContShapes)
  })

  test_that("the most specific rule wins", {
    r <- .vaeResolveShapes(list(
      list(covar = "wt", shapes = "lin"),
      list(var = "cl", covar = "wt", shapes = "power")))
    expect_equal(.vaeShapesFor(r, "cl", "WT"), "power")
    expect_equal(.vaeShapesFor(r, "v", "WT"), "lin")
  })

  test_that("a parameter matches any of its aliases", {
    r <- .vaeResolveShapes(list(list(var = "tcl", covar = "wt", shapes = "lin")))
    expect_equal(.vaeShapesFor(r, c("eta.cl", "tcl", "cl"), "WT"), "lin")
    ## a different parameter is unaffected by the rule
    expect_equal(.vaeShapesFor(r, c("eta.v", "tv", "v"), "WT"), .vaeContShapes)
  })

  test_that("bad shape specifications are rejected", {
    expect_error(.vaeResolveShapes("hockey"), "unknown covariate shape")
    expect_error(.vaeResolveShapes(c("power", "power")), "duplicate")
    expect_error(.vaeResolveShapes(list("power")), "named by covariate")
    expect_error(.vaeResolveShapes(1L), "character vector or a list")
    ## "cat" is never user-selectable
    expect_error(.vaeResolveShapes("cat"), "unknown covariate shape")
  })
})

## Regressions for three bugs found in independent review of this feature.
nmTest({
  .pinLog <- function() {
    ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; wt.cl <- 0.1; add.err <- 0.7
      eta.cl ~ 0.1 })
    model({ ka <- exp(tka); cl <- exp(tcl + wt.cl * log(WT / 70) + eta.cl)
      v <- exp(tv)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v; cp ~ add(add.err) })
  }

  test_that("a per-parameter shapes= rule is enforced, not just widened", {
    ## WT is power-only on cl and lin-only on v.  Both columns must exist (the
    ## design is shared), so the restriction has to land in the covAllow mask --
    ## otherwise cl could be given the linear column it forbids.
    ctl <- vaeControl(shapes = list(list(var = "cl", covar = "wt", shapes = "power"),
                                    list(var = "v", covar = "wt", shapes = "lin")),
                      muRefCovAlg = FALSE)
    m <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; add.err <- 0.7
        eta.cl ~ 0.1; eta.v ~ 0.1 })
      model({ ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v; cp ~ add(add.err) })
    }
    p <- suppressWarnings(.vaeDataPrep(rxode2::assertRxUi(m), nlmixr2data::theo_sd, ctl))
    expect_setequal(p$covShape, c("power", "lin"))
    expect_false(is.null(p$covAllow))
    .k <- match("eta.cl", p$etaNames)
    .kv <- match("eta.v", p$etaNames)
    .jp <- match("power", p$covShape)
    .jl <- match("lin", p$covShape)
    ## cl may only take the log family, v may only take the linear one
    expect_equal(p$covAllow[.k, .jp], 1L)
    expect_equal(p$covAllow[.k, .jl], 0L)
    expect_equal(p$covAllow[.kv, .jp], 0L)
    expect_equal(p$covAllow[.kv, .jl], 1L)
  })

  test_that("a pinned column that is not its group's first still reaches the encoder", {
    ## shapes= lists lin first, so WT_lin is the group's leading column while the
    ## pinned log(WT/70) pair takes WT_power -- deduping the encoder input against
    ## a fixed canonical column would drop WT from the encoder entirely
    ctl <- vaeControl(shapes = c("lin", "power"), muRefCovAlg = FALSE)
    p <- suppressWarnings(.vaeDataPrep(rxode2::assertRxUi(.pinLog()),
                                       nlmixr2data::theo_sd, ctl))
    expect_equal(p$covNames, c("WT_lin", "WT_power"))
    expect_equal(which(colSums(p$covAllow) > 0L), 2L)
    expect_equal(ncol(p$covIn), 1L)
    expect_equal(unname(p$covIn[, 1]), unname(p$covMat[, 2]))
  })

  test_that("a shape that cannot be written at its center is not used", {
    ## "center" would emit beta*(COV/0) and rescale the coefficient to exactly 0
    expect_false(.vaeShapeUsable("center", 0))
    expect_false(.vaeShapeUsable("log", 0))
    expect_true(.vaeShapeUsable("lin", 0))
    expect_true(.vaeShapeUsable("identity", 0))
    ## a covariate whose center is 0 falls back to a writable sibling shape
    d <- data.frame(id = rep(1:4, each = 3), time = rep(0:2, 4), dv = 1:12,
                    score = rep(c(-3, -1, 1, 3), each = 3))
    res <- vaeCovariates(d, shapes = c("center", "lin"))
    expect_equal(unique(res$center), 0)
    expect_false(any(res$shape == "center"))
    expect_true(all(res$shape == "lin"))
    ## when NOTHING requested is writable the plain form of the family is
    ## substituted, so the covariate stays searchable rather than vanishing
    res0 <- vaeCovariates(d, shapes = "center", covCenter = c(score = 0))
    expect_equal(res0$shape, "lin")
    expect_equal(nrow(res0), 1L)
  })
})
