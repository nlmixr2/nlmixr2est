nmTest({
  test_that("shapes collapse to the searchable families", {
    expect_equal(.vaeShapeFamily(c("power", "log")), c("log", "log"))
    expect_equal(.vaeShapeFamily(c("lin", "identity", "center")),
                 c("lin", "lin", "lin"))
    expect_equal(.vaeShapeFamily("cat"), "cat")
    expect_error(.vaeShapeFamily("notAShape"), "unknown covariate shape")
  })

  test_that("hockey is its own family, reachable by request or by arm", {
    ## a family of its own: its span CONTAINS lin's, so it is not another
    ## parameterization of lin and must not be swapped for one
    expect_equal(.vaeShapeFamily("hockey"), "hockey")
    expect_equal(.vaeShapeFamily(.vaeHockeyArms), c("hockey", "hockey"))
    expect_identical(.vaeHockeyArms, c("hockeyLow", "hockeyHi"))
  })

  test_that("hockey is nameable but the arms are internal", {
    ## `shapes=` takes the RELATIONSHIP; the arms are how it is built, and
    ## naming one on its own would ask for half a hockey stick
    expect_true("hockey" %in% .vaeContShapes)
    expect_silent(.vaeAssertContShapes("hockey"))
    expect_error(.vaeAssertContShapes("hockeyLow"), "unknown covariate shape")
    expect_error(.vaeAssertContShapes("hockeyHi"), "unknown covariate shape")
    ## naming a shape and defaulting to it stay separate decisions, even though
    ## hockey now happens to be both
    expect_true(all(.vaeDefaultShapes %in% .vaeContShapes))
  })

  test_that("the hockey arms are disjoint and sum to the lin column", {
    ## the property the atomic block rests on, and the reason a subject sitting
    ## exactly ON the knot must belong to one arm only
    cov <- c(50, 70, 70.5, 90, 110)
    ctr <- 70.5
    lo <- .vaeShapeValue("hockeyLow", cov, ctr)
    hi <- .vaeShapeValue("hockeyHi", cov, ctr)
    expect_equal(lo + hi, .vaeShapeValue("lin", cov, ctr))
    ## no subject is in both arms, and the subject at the knot is in the high one
    expect_true(all(lo == 0 | hi == 0))
    expect_equal(lo[cov == ctr], 0)
    expect_equal(hi[cov == ctr], 0)   # both vanish AT the knot
    ## below the knot only the low arm is non-zero, above only the high arm
    expect_true(all(hi[cov < ctr] == 0))
    expect_true(all(lo[cov > ctr] == 0))
    expect_true(all(lo[cov < ctr] < 0))
    expect_true(all(hi[cov > ctr] > 0))
  })

  test_that("hockey arm text evaluates to the arm's design column", {
    ## .vaeShapeExpr and .vaeShapeValue must agree EXACTLY, or the written model
    ## would not reproduce the fit -- including at the knot, where the `<` / `>=`
    ## boundary decides which arm a subject lands in
    WT <- c(50, 70, 70.5, 90, 110)
    ctr <- 70.5
    for (arm in .vaeHockeyArms) {
      txt <- .vaeShapeExpr(arm, "WT", ctr)
      expect_equal(eval(str2lang(txt)), .vaeShapeValue(arm, WT, ctr), info = arm)
    }
    expect_equal(.vaeShapeExpr("hockeyLow", "WT", 70.5),
                 "(WT < 70.5)*(WT - 70.5)")
    expect_equal(.vaeShapeExpr("hockeyHi", "WT", 70.5),
                 "(WT >= 70.5)*(WT - 70.5)")
  })

  test_that("a hockey arm needs no intercept correction", {
    ## the arm's written expression IS its design column, so unlike "identity" or
    ## "center" no part of the effect moves into the structural theta
    for (arm in .vaeHockeyArms) {
      r <- .vaeShapeBeta(arm, 70.5, -0.4)
      expect_equal(r$beta, -0.4, info = arm)
      expect_equal(r$interceptAdj, 0, info = arm)
    }
  })

  test_that("hockey is usable at any finite knot", {
    ## unlike log/power (need a positive center) and center (needs a non-zero
    ## one), bending is well defined wherever the knot is
    expect_true(all(.vaeShapeUsable(c("hockey", .vaeHockeyArms), 0)))
    expect_true(all(.vaeShapeUsable(c("hockey", .vaeHockeyArms), -12)))
    expect_false(any(.vaeShapeUsable(c("hockey", .vaeHockeyArms), NA_real_)))
    ## and it does not disturb the existing rules
    expect_false(.vaeShapeUsable("center", 0))
    expect_false(.vaeShapeUsable("log", -1))
  })

  test_that("arm coefficients are named hockey.low and hockey.hi", {
    expect_equal(.vaeShapeCoefTag("hockeyLow"), "hockey.low")
    expect_equal(.vaeShapeCoefTag("hockeyHi"), "hockey.hi")
    ## every other shape names itself
    for (sh in c("power", "lin", "log", "identity", "center", "cat")) {
      expect_equal(.vaeShapeCoefTag(sh), sh, info = sh)
    }
  })

  test_that("shape expressions match the nlmixr2scm vocabulary", {
    expect_equal(.vaeShapeExpr("power", "WT", 70.5), "log(WT/70.5)")
    expect_equal(.vaeShapeExpr("log", "WT", 70.5), "log(WT)")
    expect_equal(.vaeShapeExpr("lin", "WT", 70.5), "(WT - 70.5)")
    expect_equal(.vaeShapeExpr("identity", "WT", 70.5), "WT")
    expect_equal(.vaeShapeExpr("center", "WT", 70.5), "(WT/70.5)")
    expect_error(.vaeShapeExpr("notAShape", "WT", 1), "unknown covariate shape")
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
    r <- .vaeResolveShapes(c("power", "lin"))$rules
    expect_equal(nrow(r), 1L)
    expect_true(is.na(r$var) && is.na(r$cov))
    expect_equal(.vaeShapesFor(r, "cl", "WT"), c("power", "lin"))
    ## NULL means the default shapes
    expect_equal(.vaeShapesFor(.vaeResolveShapes(NULL)$rules, "cl", "WT"),
                 .vaeDefaultShapes)
  })

  test_that("a covariate-named list restricts that covariate only", {
    r <- .vaeResolveShapes(list(wt = "power"))$rules
    ## matching is case-insensitive because the VAE upper-cases data columns
    expect_equal(.vaeShapesFor(r, "cl", "WT"), "power")
    ## an unmatched covariate stays free
    expect_equal(.vaeShapesFor(r, "cl", "AGE"), .vaeDefaultShapes)
  })

  test_that("a pairsVec list restricts one parameter/covariate pair", {
    r <- .vaeResolveShapes(list(
      list(var = "cl", covar = "wt", shapes = "power"),
      list(var = "v", covar = "wt", shapes = c("lin", "identity"))))$rules
    expect_equal(.vaeShapesFor(r, "cl", "WT"), "power")
    expect_equal(.vaeShapesFor(r, "v", "WT"), c("lin", "identity"))
    ## shapes= restricts parameterizations only -- a pair no rule mentions is
    ## still searched, with the default shapes available
    expect_equal(.vaeShapesFor(r, "ka", "WT"), .vaeDefaultShapes)
    expect_equal(.vaeShapesFor(r, "cl", "AGE"), .vaeDefaultShapes)
  })

  test_that("the most specific rule wins", {
    r <- .vaeResolveShapes(list(
      list(covar = "wt", shapes = "lin"),
      list(var = "cl", covar = "wt", shapes = "power")))$rules
    expect_equal(.vaeShapesFor(r, "cl", "WT"), "power")
    expect_equal(.vaeShapesFor(r, "v", "WT"), "lin")
  })

  test_that("a parameter matches any of its aliases", {
    r <- .vaeResolveShapes(list(list(var = "tcl", covar = "wt", shapes = "lin")))$rules
    expect_equal(.vaeShapesFor(r, c("eta.cl", "tcl", "cl"), "WT"), "lin")
    ## a different parameter is unaffected by the rule
    expect_equal(.vaeShapesFor(r, c("eta.v", "tv", "v"), "WT"), .vaeDefaultShapes)
  })

  test_that("bad shape specifications are rejected", {
    expect_error(.vaeResolveShapes("notAShape"), "unknown covariate shape")
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
    ## "B" must clear catCutoff so a RACE_B column EXISTS -- otherwise the pool
    ## is empty and the test passes trivially instead of showing that "Asian"
    ## refuses to pin to some other level's column
    d$RACE <- c(rep("W", 8), rep("B", 3), "Asian")[match(d$ID, ids)]
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
    ## RACE_B exists but RACE_Asian does not, so there IS a wrong column to
    ## mis-pin to
    expect_true("B" %in% p$covLevel)
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

## fixCov: naming a covariate in shapes= is the statement that it belongs in the
## search.  Eligibility is asserted on the matrix rather than on a fit outcome --
## exact, and it does not need a converged model to be meaningful.
nmTest({
  .etas <- c("eta.cl", "eta.v", "eta.ka")
  .thetas <- c("tcl", "tv", "tka")
  .raw <- c("WT", "WT", "AGE", "AGE", "SEX")
  .el <- function(spec) {
    r <- .vaeResolveShapes(spec)
    .vaeEligible(r$rules, r$fixCov, .etas, .thetas, .raw)
  }

  test_that("fixCov=TRUE is the default and fixes the searched set", {
    m <- .el(list(wt = "power"))
    expect_true(all(m[, "WT"]))
    expect_false(any(m[, "AGE"]))
    expect_false(any(m[, "SEX"]))
    ## the default only applies to the list form -- a character vector names no
    ## covariate, so there is nothing to fix
    expect_true(all(.el(c("power", "lin"))))
    expect_true(all(.el(NULL)))
  })

  test_that("fixCov=FALSE reproduces the parameterization-only behavior", {
    expect_true(all(.el(list(wt = "power", fixCov = FALSE))))
    ## and it leaves the shape restriction itself intact
    r <- .vaeResolveShapes(list(wt = "power", fixCov = FALSE))
    expect_false(r$fixCov)
    expect_equal(.vaeShapesFor(r$rules, "cl", "WT"), "power")
  })

  test_that("a pair rule makes only that pair eligible", {
    m <- .el(list(list(var = "cl", covar = "wt", shapes = "power")))
    expect_true(m[1L, "WT"])
    expect_false(any(m[-1L, "WT"]))
    ## covar-only: that covariate on every parameter
    m2 <- .el(list(list(covar = "wt", shapes = "power")))
    expect_true(all(m2[, "WT"]))
    expect_false(any(m2[, "AGE"]))
    ## var-only: every covariate, on that parameter alone
    m3 <- .el(list(list(var = "v", shapes = "power")))
    expect_true(all(m3[2L, ]))
    expect_false(any(m3[-2L, ]))
  })

  test_that("named and pair entries mix in one list", {
    r <- .vaeResolveShapes(list(list(var = "cl", covar = "wt", shapes = "power"),
                                sex = TRUE))
    expect_equal(nrow(r$rules), 2L)
    ## a named entry is exactly the covar-only rule it is shorthand for
    expect_equal(.vaeShapesFor(r$rules, "cl", "WT"), "power")
    expect_equal(.vaeShapesFor(r$rules, "ka", "SEX"), .vaeDefaultShapes)
    m <- .vaeEligible(r$rules, r$fixCov, .etas, .thetas, .raw)
    expect_true(m[1L, "WT"])
    expect_false(any(m[-1L, "WT"]))
    expect_true(all(m[, "SEX"]))
    expect_false(any(m[, "AGE"]))
  })

  test_that("fixCov survives the pair form it is mixed into", {
    ## the regression this ordering exists to prevent: a bare logical element is
    ## not a list, so leaving it in would flip the per-element dispatch
    r <- .vaeResolveShapes(list(list(var = "cl", covar = "wt", shapes = "power"),
                                fixCov = FALSE))
    expect_false(r$fixCov)
    expect_equal(nrow(r$rules), 1L)
    expect_true(all(.el(list(list(var = "cl", covar = "wt", shapes = "power"),
                             fixCov = FALSE))))
  })

  test_that("TRUE means eligible with the default shapes", {
    r <- .vaeResolveShapes(list(wt = TRUE))
    expect_equal(.vaeShapesFor(r$rules, "cl", "WT"), .vaeDefaultShapes)
    expect_true(all(.el(list(wt = TRUE))[, "WT"]))
  })

  test_that("contradictory fixCov specifications are rejected", {
    ## a rule naming neither var nor covar makes everything eligible
    expect_error(.el(list(list(shapes = "lin"), wt = "power")),
                 "contradicts fixCov")
    expect_error(.vaeResolveShapes(list(wt = FALSE)),
                 "TRUE or a shape vector")
    expect_error(.vaeResolveShapes(list(wt = "power", fixCov = TRUE, fixCov = FALSE)),
                 "more than once")
    expect_error(.vaeResolveShapes(list(wt = "power", fixCov = "yes")), "fixCov")
    ## a data column colliding with the flag name
    expect_error(.vaeEligible(.vaeResolveShapes(list(wt = "power"))$rules, TRUE,
                              .etas, .thetas, c("WT", "FIXCOV")),
                 "collides with the eligibility flag")
  })

  test_that("an unnamed non-list element is still an error", {
    expect_error(.vaeResolveShapes(list("power")), "named by covariate")
    expect_error(.vaeResolveShapes(list(list(covar = "wt", shapes = "lin"), "power")),
                 "named by covariate")
  })
})

## fixCov end to end: the mask .vaeShapeAllowMask actually produces, and the
## three .vaeDataPrep paths (exclusion note, nothing-searchable error, and the
## declaring model that overrides the flag).
nmTest({
  .d <- nlmixr2data::theo_sd
  .d$WT <- rep(c(60, 70, 80), length.out = length(unique(.d$ID)))[
    match(.d$ID, unique(.d$ID))]
  .d$AGE <- rep(c(30, 40, 50), length.out = length(unique(.d$ID)))[
    match(.d$ID, unique(.d$ID))]
  names(.d) <- toupper(names(.d))

  ## built as text like the pinning tests above: model()/ini() piping inside
  ## test_that() resolves the covariate symbol in the wrong environment
  .mk <- function(term = NULL) {
    .ini <- if (is.null(term)) "" else " b1 <- 0.1;"
    .cl <- if (is.null(term)) "exp(tcl + eta.cl)"
           else paste0("exp(tcl + b1 * ", term, " + eta.cl)")
    eval(parse(text = paste0(
      "function() {\n",
      "  ini({ tka <- 0.45; tcl <- 1; tv <- 3.45;", .ini, " add.sd <- 0.7\n",
      "        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1 })\n",
      "  model({ ka <- exp(tka + eta.ka)\n",
      "    cl <- ", .cl, "\n",
      "    v <- exp(tv + eta.v)\n",
      "    d/dt(depot) <- -ka * depot\n",
      "    d/dt(center) <- ka * depot - cl / v * center\n",
      "    cp <- center / v; cp ~ add(add.sd) })\n}")))
  }

  test_that("fixCov masks whole covariates in .vaeShapeAllowMask", {
    r <- .vaeResolveShapes(list(wt = "power"))
    cov <- .vaeCovariateSearch(.d, unique(.d$ID), r$rules)
    m <- .vaeShapeAllowMask(cov, r, c("eta.ka", "eta.cl"), c("tka", "tcl"))
    ## every WT column stays selectable, every AGE column is masked out
    expect_true(all(m[, cov$covRaw == "WT"] == 1L))
    expect_true(all(m[, cov$covRaw == "AGE"] == 0L))
    ## fixCov=FALSE leaves AGE alone
    r2 <- .vaeResolveShapes(list(wt = "power", fixCov = FALSE))
    m2 <- .vaeShapeAllowMask(cov, r2, c("eta.ka", "eta.cl"), c("tka", "tcl"))
    expect_true(all(m2[, cov$covRaw == "AGE"] == 1L))
  })

  ## .vaeDataPrep emits several warnings and their order is not part of the
  ## contract, so collect them all rather than relying on which one comes first
  .warns <- function(expr) {
    .w <- character(0)
    withCallingHandlers(invisible(force(expr)),
                        warning = function(cond) {
                          .w <<- c(.w, conditionMessage(cond))
                          invokeRestart("muffleWarning")
                        })
    .w
  }

  test_that("a covariate restricted to one shape is not reported as excluded", {
    ## the false positive an independent review predicted: WT's non-power
    ## columns are masked, but WT itself is still searched through WT_power
    .w <- .warns(.vaeDataPrep(rxode2::assertRxUi(.mk()), .d,
                              vaeControl(shapes = list(wt = "power"),
                                         muRefCovAlg = FALSE)))
    .fx <- grep("fixCov=TRUE, covariate", .w, value = TRUE)
    expect_length(.fx, 1L)
    expect_match(.fx, "AGE")
    expect_false(grepl("WT", .fx))
  })

  test_that("fixCov=TRUE with nothing searchable is an error", {
    expect_error(
      .vaeDataPrep(rxode2::assertRxUi(.mk()), .d,
                   vaeControl(shapes = list(noSuchCov = "power"),
                              muRefCovAlg = FALSE)),
      "leaves no covariate searchable")
  })

  test_that("covariateSelection=FALSE has nothing for fixCov to narrow", {
    ## the error above told the user to set covariateSelection=FALSE, so it must
    ## not fire at someone who already has -- there is no search to restrict
    expect_error(
      suppressWarnings(
        .vaeDataPrep(rxode2::assertRxUi(.mk()), .d,
                     vaeControl(covariateSelection = FALSE,
                                shapes = list(noSuchCov = "power"),
                                muRefCovAlg = FALSE))),
      NA)
    ## nor may it blame fixCov for a search that was off regardless
    .w <- .warns(.vaeDataPrep(rxode2::assertRxUi(.mk()), .d,
                              vaeControl(covariateSelection = FALSE,
                                         shapes = list(wt = "power"),
                                         muRefCovAlg = FALSE)))
    expect_false(any(grepl("fixCov=TRUE, covariate", .w)))
  })

  test_that("a declaring model overrides fixCov", {
    .w <- .warns(.vaeDataPrep(rxode2::assertRxUi(.mk("log(WT/70)")), .d,
                              vaeControl(shapes = list(age = "power"),
                                         muRefCovAlg = FALSE)))
    expect_true(any(grepl("fixCov=TRUE ignored", .w)))
    ## the declaration is what actually restricts the search
    expect_true(any(grepl("pinned to model-specified covariates", .w)))
  })

  test_that("fixCov=TRUE naming no covariate is rejected outright", {
    expect_error(.vaeResolveShapes(list(fixCov = TRUE)), "no covariate is named")
    ## ... while an explicit FALSE with no rules is simply unrestricted
    expect_false(.vaeResolveShapes(list(fixCov = FALSE))$fixCov)
  })
})
