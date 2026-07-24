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
