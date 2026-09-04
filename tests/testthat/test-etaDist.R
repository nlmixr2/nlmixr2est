nmTest({
  # Declared non-Gaussian random effect distributions.
  #
  # rxode2 does the model rewrite (`rxEtaDistExpand()`); what is tested
  # here is nlmixr2est's two jobs -- running that rewrite before
  # estimation, and refusing the methods for which the rewritten model is
  # not something they can fit.

  .declMod <- function() {
    .f <- function() {
      ini({
        tka <- 0.45
        lclm  <- log(1.1)
        lclrv <- log(0.09)
        tv <- 3.45
        eta.ka ~ 0.6
        eta.cl ~ 1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- eta.cl
        v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl/v * center
        cp <- center/v
        cp ~ add(add.sd)
      })
    }
    .b <- body(.f)
    .i <- .b[[2]]
    .i[[2]][[length(.i[[2]]) + 1L]] <-
      str2lang("dist(eta.cl) ~ dgamma(shape=1/exp(lclrv), rate=1/(exp(lclrv)*exp(lclm)))")
    .b[[2]] <- .i
    body(.f) <- .b
    eval(parse(text=paste(deparse(.f), collapse="\n")))
  }

  test_that("the pre-processing hook expands a declaration", {
    .u <- .declMod()
    expect_equal(nrow(rxode2::rxUiEtaDists(.u)), 1L)
    .r <- nlmixr2est:::.preProcessEtaDist(.u, "focei", NULL, NULL)
    expect_true(is.list(.r))
    expect_equal(nrow(rxode2::rxUiEtaDists(.r$ui)), 0L)
    expect_true("rxz.eta.cl" %in% .r$ui$iniDf$name)
    ## a model with no declaration is left alone entirely
    .p <- function() {
      ini({
        tka <- 0.45
        tv <- 3.45
        eta.ka ~ 0.6
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - center/v
        cp <- center/v
        cp ~ add(add.sd)
      })
    }
    expect_null(nlmixr2est:::.preProcessEtaDist(.p(), "focei", NULL, NULL))
  })

  test_that("the expansion hook runs before the bounded transform", {
    ## the bounded transform wraps whatever thetas are there by then,
    ## including the `rxCor.*` ones the expansion adds
    .n <- preProcessHooks()
    if (".preProcessEtaDist" %in% .n && ".preProcessBoundedTransform" %in% .n) {
      expect_lt(which(.n == ".preProcessEtaDist"),
                which(.n == ".preProcessBoundedTransform"))
      expect_equal(.n[1], ".preProcessEtaDist")
    }
  })

  test_that("methods declare whether they support a declaration", {
    ## the FOCEi family, SAEM and simulation do
    for (.m in c("focei", "foce", "fo", "laplace", "agq", "imp", "impmap",
                 "saem", "posthoc", "rxSolve", "simulate")) {
      expect_true(nlmixr2est:::.isEtaDistMethod(.m), info=.m)
    }
    ## a nonparametric random effect distribution contradicts a declared
    ## one; nlme/nls are Gaussian by construction; the variational methods
    ## hardcode the normal family in their ELBO
    for (.m in c("npag", "npb", "nlme", "nls", "vae", "emvi", "fbvi")) {
      expect_false(nlmixr2est:::.isEtaDistMethod(.m), info=.m)
    }
    expect_false(nlmixr2est:::.isEtaDistMethod("notAMethod"))
  })

  test_that("an unsupported method refuses rather than ignoring", {
    .env <- new.env(parent=emptyenv())
    assign("ui", .declMod(), envir=.env)
    assign("control", NULL, envir=.env)
    class(.env) <- c("npag", "nlmixr2Est")
    expect_error(nlmixr2est:::.nlmixr2AssertEtaDist(.env), "eta.cl")
    expect_error(nlmixr2est:::.nlmixr2AssertEtaDist(.env), "npag")
    class(.env) <- c("focei", "nlmixr2Est")
    expect_silent(nlmixr2est:::.nlmixr2AssertEtaDist(.env))
  })

  test_that("the copula correlation is recovered from the rxCor thetas", {
    ## round trip: a correlation matrix -> the unconstrained parameters
    ## rxode2 writes -> back again
    .R <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
    .y <- rxode2:::.rxEtaDistCorToY(.R)
    .back <- nlmixr2est:::.etaDistCorFromY(c("a", "b"), list(rxCor.b.a=.y[2, 1]))
    expect_equal(unname(.back), .R, tolerance=1e-10)
    ## and a 3x3, where the parameters are partial correlations
    .R3 <- matrix(c(1, 0.5, 0.2,
                    0.5, 1, -0.3,
                    0.2, -0.3, 1), 3, 3)
    .y3 <- rxode2:::.rxEtaDistCorToY(.R3)
    .back3 <- nlmixr2est:::.etaDistCorFromY(c("a", "b", "c"),
                               list(rxCor.b.a=.y3[2, 1],
                                    rxCor.c.a=.y3[3, 1],
                                    rxCor.c.b=.y3[3, 2]))
    expect_equal(unname(.back3), .R3, tolerance=1e-10)
  })

  test_that("a declared normal fits exactly like the model it replaces", {
    skip_on_cran()
    ## the identity test: the phiU/inverse-CDF/copula pipeline must not
    ## move the objective function at all when the declared family IS the
    ## normal it replaces
    .dec <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        sdcl <- 0.3
        eta.cl ~ 1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl/v * center
        cp <- center/v
        cp ~ add(add.sd)
      })
    }
    .b <- body(.dec); .i <- .b[[2]]
    .i[[2]][[length(.i[[2]]) + 1L]] <- str2lang("dist(eta.cl) ~ dnorm(0, sdcl)")
    .b[[2]] <- .i; body(.dec) <- .b
    .dec <- eval(parse(text=paste(deparse(.dec), collapse="\n")))
    .plain <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        eta.cl ~ 0.09
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl/v * center
        cp <- center/v
        cp ~ add(add.sd)
      })
    }
    .ctl <- foceiControl(maxOuterIterations=0, maxInnerIterations=100,
                         covMethod="", print=0)
    .a <- suppressWarnings(nlmixr2(.dec, nlmixr2data::theo_sd, est="focei",
                                   control=.ctl))
    .b2 <- suppressWarnings(nlmixr2(.plain, nlmixr2data::theo_sd, est="focei",
                                    control=.ctl))
    expect_equal(.a$objf, .b2$objf, tolerance=1e-4)
  })
})
