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
    ## a CHARACTER level is quoted even when it looks numeric -- comparing a
    ## character covariate against a bare number would test the wrong thing and
    ## would not match its own column when the model is piped back
    expect_equal(.vaeShapeExpr("cat", "GRP", level = "2"), "(GRP == \"2\")")
    expect_equal(.vaeShapeExpr("cat", "GRP", level = "01"), "(GRP == \"01\")")
    ## the numeric/character decision comes from the COVARIATE's type
    expect_identical(.vaeCovLevelValue(c(1, 2, 3), "2"), 2)
    expect_identical(.vaeCovLevelValue(c("A", "01"), "01"), "01")
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

## Piping: a model the VAE wrote must feed back into another fit with every
## covariate form recognized and pinned to the shape it was written in.
nmTest({
  .mk <- function(term) {
    .txt <- paste0(
      "function() {\n",
      "  ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; b1 <- 0.1; add.err <- 0.7\n",
      "        eta.cl ~ 0.1 })\n",
      "  model({ ka <- exp(tka)\n",
      "    cl <- exp(tcl + b1 * ", term, " + eta.cl)\n",
      "    v <- exp(tv)\n",
      "    d/dt(depot) <- -ka * depot\n",
      "    d/dt(center) <- ka * depot - cl / v * center\n",
      "    cp <- center / v; cp ~ add(add.err) })\n}")
    eval(parse(text = .txt))
  }

  test_that("every written shape is detected with its center", {
    expect_equal(.vaeDetectShape(quote(log(WT/70)), "WT")[c("shape", "center")],
                 list(shape = "power", center = 70))
    expect_equal(.vaeDetectShape(quote(log(WT)), "WT")[c("shape", "center")],
                 list(shape = "log", center = 1))
    expect_equal(.vaeDetectShape(quote((WT - 70)), "WT")[c("shape", "center")],
                 list(shape = "lin", center = 70))
    expect_equal(.vaeDetectShape(quote(WT), "WT")[c("shape", "center")],
                 list(shape = "identity", center = 0))
    expect_equal(.vaeDetectShape(quote((WT/70)), "WT")[c("shape", "center")],
                 list(shape = "center", center = 70))
    expect_equal(.vaeDetectShape(quote((SEX == "F")), "SEX")$shape, "cat")
    ## an untransferable form stays unrecognized
    expect_true(is.na(.vaeDetectShape(quote(log(WT/AGE)), "WT")$shape))
    expect_true(is.na(.vaeDetectShape(quote(sqrt(WT)), "WT")$shape))
  })

  test_that("the coefficient's own factor is what is read", {
    .l <- quote(cl <- exp(tcl + b1 * log(WT/70) + b2 * (AGE - 40) + eta.cl))
    expect_equal(.vaeDetectShape(.vaeCoefFactor(.l, "b1"), "WT")$shape, "power")
    expect_equal(.vaeDetectShape(.vaeCoefFactor(.l, "b2"), "AGE")$shape, "lin")
  })

  test_that("each written shape pins to its own column and transfers directly", {
    d <- as.data.frame(nlmixr2data::theo_sd)
    .wt <- vapply(unique(d$ID), function(i) d$WT[d$ID == i][1], numeric(1))
    .cases <- list(power = list("log(WT/70)", function(v) log(v / 70)),
                   log = list("log(WT)", function(v) log(v)),
                   lin = list("(WT - 70)", function(v) v - 70),
                   identity = list("WT", function(v) v),
                   center = list("(WT/70)", function(v) v / 70))
    for (.s in names(.cases)) {
      .term <- .cases[[.s]][[1]]
      .val <- .cases[[.s]][[2]]
      p <- suppressWarnings(.vaeDataPrep(rxode2::assertRxUi(.mk(.term)), d,
                                         vaeControl(muRefCovAlg = FALSE)))
      expect_true(p$pinActive, info = .s)
      ## detected as the shape it was written in, and pinned
      expect_equal(p$pinPairs$shape, .s, info = .s)
      expect_true(all(p$pinPairs$inPool), info = .s)
      ## exactly one cell allowed, on the column of that shape's family
      expect_equal(sum(p$covAllow), 1L, info = .s)
      .j <- which(colSums(p$covAllow) > 0L)
      expect_equal(p$covFamily[.j], .vaeShapeFamily(.s), info = .s)
      ## the column IS the model's own expression, so the slope transfers with
      ## no correction (this is the pinned contract)
      expect_equal(unname(p$covMat[, .j]), unname(.val(.wt)), info = .s)
      expect_equal(p$covPop[.j], 0, info = .s)
      ## and the coefficient is searched, not routed to the regress M-step
      expect_false("b1" %in% p$regressNames, info = .s)
    }
  })

  test_that("an unrecognized form still falls back to the regress M-step", {
    d <- as.data.frame(nlmixr2data::theo_sd)
    p <- suppressWarnings(.vaeDataPrep(rxode2::assertRxUi(.mk("sqrt(WT)")), d,
                                       vaeControl(muRefCovAlg = FALSE)))
    expect_false(any(p$pinPairs$inPool))
    expect_true("b1" %in% p$regressNames)
  })
})

nmTest({
  test_that("a shapes= rule never silently deletes a declared covariate effect", {
    ## The model declares a LINEAR WT effect on cl while shapes= allows lin only
    ## on v.  Intersecting the search restriction with the pin mask would empty
    ## cl's row, leaving the declared coefficient neither searched nor regressed
    ## -- it would then be written back as exactly 0.  Pinning must win.
    d <- as.data.frame(nlmixr2data::theo_sd)
    m <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; cl.wt <- 0.1; add.err <- 0.7
        eta.cl ~ 0.1; eta.v ~ 0.1 })
      model({ ka <- exp(tka); cl <- exp(tcl + cl.wt * (WT - 70) + eta.cl)
        v <- exp(tv + eta.v)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v; cp ~ add(add.err) })
    }
    p <- suppressWarnings(.vaeDataPrep(
      rxode2::assertRxUi(m), d,
      vaeControl(shapes = list(list(var = "cl", covar = "wt", shapes = "power"),
                               list(var = "v", covar = "wt", shapes = "lin")),
                 muRefCovAlg = FALSE)))
    expect_true(p$pinActive)
    expect_equal(p$pinPairs$shape, "lin")
    ## the declared effect is still reachable: searched on its own column ...
    expect_true(all(p$pinPairs$inPool))
    expect_equal(sum(p$covAllow), 1L)
    .j <- which(colSums(p$covAllow) > 0L)
    expect_equal(p$covFamily[.j], "lin")
    ## ... or, failing that, estimated in place -- never both dropped
    expect_true(sum(p$covAllow) > 0L || "cl.wt" %in% p$regressNames)
  })

  test_that("a declared shape with no column to pin to falls back to regress", {
    ## shapes="power" builds only the log column, so a declared linear effect has
    ## nothing to pin to and must be estimated in place rather than vanish
    d <- as.data.frame(nlmixr2data::theo_sd)
    m <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; cl.wt <- 0.1; add.err <- 0.7
        eta.cl ~ 0.1 })
      model({ ka <- exp(tka); cl <- exp(tcl + cl.wt * (WT - 70) + eta.cl)
        v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v; cp ~ add(add.err) })
    }
    p <- suppressWarnings(.vaeDataPrep(rxode2::assertRxUi(m), d,
                                       vaeControl(shapes = "power",
                                                  muRefCovAlg = FALSE)))
    expect_equal(p$covNames, "WT_power")
    expect_false(any(p$pinPairs$inPool))
    expect_true("cl.wt" %in% p$regressNames)
  })
})

## Regressions for the second independent-review pass.
nmTest({
  .d <- as.data.frame(nlmixr2data::theo_sd)

  test_that("a coefficient must multiply a standalone additive term", {
    ## searching anywhere in the line returns WT for both of these, which would
    ## pin a slope that does not transfer and silently replace the user's term
    .sq <- quote(cl <- exp(tcl + sqrt(b1 * WT) + eta.cl))
    expect_null(.vaeCoefFactor(.sq, "b1"))
    .ix <- quote(cl <- exp(tcl + b1 * WT * AGE + eta.cl))
    expect_null(.vaeCoefFactor(.ix, "b1"))
    ## the legitimate additive forms still resolve, in either factor order
    expect_equal(.vaeCoefFactor(quote(cl <- exp(tcl + b1 * log(WT/70) + eta.cl)), "b1"),
                 quote(log(WT/70)))
    expect_equal(.vaeCoefFactor(quote(cl <- exp(tcl + log(WT/70) * b1 + eta.cl)), "b1"),
                 quote(log(WT/70)))
    expect_equal(.vaeCoefFactor(quote(cl <- exp(tcl + b1 * (WT - 70) + eta.cl)), "b1"),
                 quote((WT - 70)))
  })

  test_that("an interaction or wrapped term regresses instead of pinning", {
    mk <- function(term) {
      eval(parse(text = paste0(
        "function() {\n",
        "  ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; cl.wt <- 0.1; add.err <- 0.7\n",
        "        eta.cl ~ 0.1 })\n",
        "  model({ ka <- exp(tka)\n",
        "    cl <- exp(tcl + ", term, " + eta.cl)\n",
        "    v <- exp(tv)\n",
        "    d/dt(depot) <- -ka * depot\n",
        "    d/dt(center) <- ka * depot - cl / v * center\n",
        "    cp <- center / v; cp ~ add(add.err) })\n}")))
    }
    for (.t in c("sqrt(cl.wt * WT)", "cl.wt * WT * AGE")) {
      d2 <- .d
      d2$AGE <- 40 + (as.integer(d2$ID) %% 5)
      p <- suppressWarnings(.vaeDataPrep(rxode2::assertRxUi(mk(.t)), d2,
                                         vaeControl(muRefCovAlg = FALSE)))
      expect_false(any(p$pinPairs$inPool), info = .t)
      expect_true("cl.wt" %in% p$regressNames, info = .t)
    }
  })

  test_that("a written factor level pins to that level's indicator", {
    d2 <- .d
    ids <- unique(d2$ID)
    lv <- c("W", "W", "W", "W", "W", "B", "B", "B", "U", "U", "U", "U")
    d2$RACE <- lv[match(d2$ID, ids)]
    mk <- function(term) {
      eval(parse(text = paste0(
        "function() {\n",
        "  ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; r1 <- 0.1; add.err <- 0.7\n",
        "        eta.cl ~ 0.1 })\n",
        "  model({ ka <- exp(tka)\n",
        "    cl <- exp(tcl + r1 * ", term, " + eta.cl)\n",
        "    v <- exp(tv)\n",
        "    d/dt(depot) <- -ka * depot\n",
        "    d/dt(center) <- ka * depot - cl / v * center\n",
        "    cp <- center / v; cp ~ add(add.err) })\n}")))
    }
    for (.l in c("B", "U")) {
      p <- suppressWarnings(.vaeDataPrep(
        rxode2::assertRxUi(mk(paste0('(RACE == "', .l, '")'))), d2,
        vaeControl(muRefCovAlg = FALSE)))
      ## the indicator round-trips: recognized, pinned, and searched
      expect_equal(p$pinPairs$shape, "cat", info = .l)
      expect_true(all(p$pinPairs$inPool), info = .l)
      expect_false("r1" %in% p$regressNames, info = .l)
      ## and pinned to the column for THAT level, not the covariate's first one
      .j <- which(colSums(p$covAllow) > 0L)
      expect_equal(length(.j), 1L, info = .l)
      expect_equal(p$covLevel[.j], .l, info = .l)
    }
  })

  test_that("a fallback column is not masked away by shapes=", {
    ## WT has non-positive values so no log-family column can be built; the
    ## linear fallback keeps it searchable and must not then be masked out
    d2 <- data.frame(id = rep(1:6, each = 3), time = rep(0:2, 6), dv = 1:18,
                     wt = rep(c(-10, 0, 10, 20, 30, 40), each = 3))
    res <- suppressWarnings(vaeCovariates(d2, shapes = "power"))
    expect_equal(nrow(res), 1L)
    expect_equal(res$shape, "lin")
    m <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; add.err <- 0.7; eta.cl ~ 0.1 })
      model({ ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v; cp ~ add(add.err) })
    }
    names(d2) <- toupper(names(d2))
    d2$AMT <- 0
    p <- suppressWarnings(.vaeDataPrep(rxode2::assertRxUi(m), d2,
                                       vaeControl(shapes = "power")))
    ## the fallback column stays selectable rather than being zeroed everywhere
    expect_true(is.null(p$covAllow) || sum(p$covAllow) > 0L)
  })
})

## Regressions for the third independent-review pass.
nmTest({
  test_that("a subtracted coefficient is not pinned with a flipped sign", {
    ## `theta - beta*cov` fits beta for +cov, and writing that beta back into the
    ## model applies it as -beta*cov: the covariate effect inverts
    expect_null(.vaeCoefFactor(quote(cl <- exp(tcl - b1 * (WT - 70) + eta.cl)), "b1"))
    expect_null(.vaeCoefFactor(quote(cl <- exp(tcl + -(b1 * (WT - 70)) + eta.cl)), "b1"))
    ## a positive term is still found when something else is subtracted
    expect_equal(.vaeCoefFactor(quote(cl <- exp(tcl + b1 * (WT - 70) - eta.cl)), "b1"),
                 quote((WT - 70)))
  })

  test_that("a coefficient used more than once is not pinned", {
    ## fitting on the first term and writing the estimate back to both corrupts
    ## the prediction, so this must regress instead
    expect_null(.vaeCoefFactor(
      quote(cl <- exp(tcl + b1 * (WT - 70) + b1 * (AGE - 40) + eta.cl)), "b1"))
    ## and a coefficient reused inside another call escapes the additive walk
    expect_null(.vaeCoefFactor(
      quote(cl <- exp(tcl + b1 * (WT - 70) + sqrt(b1) + eta.cl)), "b1"))
  })

  test_that("a subtracted or duplicated coefficient regresses end-to-end", {
    d <- as.data.frame(nlmixr2data::theo_sd)
    m <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; b1 <- 0.1; add.err <- 0.7
        eta.cl ~ 0.1 })
      model({ ka <- exp(tka); cl <- exp(tcl - b1 * (WT - 70) + eta.cl)
        v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v; cp ~ add(add.err) })
    }
    p <- suppressWarnings(.vaeDataPrep(rxode2::assertRxUi(m), d,
                                       vaeControl(muRefCovAlg = FALSE)))
    expect_false(any(p$pinPairs$inPool))
    expect_true("b1" %in% p$regressNames)
    expect_equal(sum(p$covAllow), 0L)
  })

  test_that("a written level with no column regresses instead of mis-pinning", {
    ## "Asian" is held by 1/12 subjects, below catCutoff, so it is lumped into the
    ## reference and has no indicator column.  Pinning it to some OTHER level's
    ## column would fit the wrong indicator entirely.
    d <- as.data.frame(nlmixr2data::theo_sd)
    ids <- unique(d$ID)
    d$RACE <- c(rep("W", 10), "B", "Asian")[match(d$ID, ids)]
    m <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; r1 <- 0.1; add.err <- 0.7
        eta.cl ~ 0.1 })
      model({ ka <- exp(tka); cl <- exp(tcl + r1 * (RACE == "Asian") + eta.cl)
        v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v; cp ~ add(add.err) })
    }
    p <- suppressWarnings(.vaeDataPrep(rxode2::assertRxUi(m), d,
                                       vaeControl(muRefCovAlg = FALSE,
                                                  catCutoff = 0.10)))
    ## no RACE_Asian column exists
    expect_false("Asian" %in% p$covLevel)
    ## so the declared effect is estimated in place, not pinned to RACE_B
    expect_false(any(p$pinPairs$inPool))
    expect_true("r1" %in% p$regressNames)
    expect_equal(sum(p$covAllow), 0L)
  })
})

nmTest({
  test_that("a numeric-looking character level round-trips", {
    ## level "01" must be written as == "01" and pin back to its own column
    d <- as.data.frame(nlmixr2data::theo_sd)
    ids <- unique(d$ID)
    ## "A" must be clearly modal, else a frequency tie makes "01" the reference
    d$GRP <- c(rep("A", 6), rep("01", 3), rep("02", 3))[match(d$ID, ids)]
    res <- vaeCovariates(d)
    expect_true("01" %in% res$level)
    m <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; g1 <- 0.1; add.err <- 0.7
        eta.cl ~ 0.1 })
      model({ ka <- exp(tka); cl <- exp(tcl + g1 * (GRP == "01") + eta.cl)
        v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v; cp ~ add(add.err) })
    }
    p <- suppressWarnings(.vaeDataPrep(rxode2::assertRxUi(m), d,
                                       vaeControl(muRefCovAlg = FALSE)))
    expect_true(all(p$pinPairs$inPool))
    .j <- which(colSums(p$covAllow) > 0L)
    expect_equal(p$covLevel[.j], "01")
    expect_false("g1" %in% p$regressNames)
  })
})

## Regressions for the fourth independent-review pass.
nmTest({
  test_that("a bounded mu transform is still walked", {
    ## expit/logit/probit take limit arguments, so requiring a one-argument call
    ## demoted an ordinary bounded parameter to the regress M-step
    p <- function(s) parse(text = s)[[1]]
    expect_equal(.vaeCoefFactor(p("f1 <- expit(tf1 + b1 * (WT - 70) + eta.f1, 0, 1)"), "b1"),
                 quote((WT - 70)))
    expect_equal(.vaeCoefFactor(p("f1 <- probitInv(tf1 + b1 * (WT - 70) + eta.f1, 0, 1)"), "b1"),
                 quote((WT - 70)))
    expect_equal(.vaeCoefFactor(p("f1 <- logit(tf1 + b1 * log(WT/70) + eta.f1, 0, 10)"), "b1"),
                 quote(log(WT/70)))
    ## the unbounded spellings keep working
    expect_equal(.vaeCoefFactor(p("f1 <- expit(tf1 + b1 * (WT - 70) + eta.f1)"), "b1"),
                 quote((WT - 70)))
  })

  test_that("a level containing quotes or backslashes emits parseable text", {
    for (.l in c("A\"B", "A\\B", "A B", "A'B")) {
      .txt <- .vaeShapeExpr("cat", "GRP", level = .l)
      .e <- tryCatch(parse(text = .txt)[[1]], error = function(e) NULL)
      expect_false(is.null(.e), info = .l)
      ## and it round-trips: the level reads back exactly as written
      expect_equal(.vaeDetectShape(.e, "GRP")$level, .l, info = .l)
    }
  })
})
