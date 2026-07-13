nmTest({
  # Plain (covariate-free) mu-referenced theta+eta pairs form intercept-only
  # regression groups (.muRefClassify/.muRefGroups, R/muRefClassify.R) so the
  # mu-FOCEI family profiles them out of the outer optimizer like the
  # covariate groups. UI-only unit tests; the end-to-end fits live in
  # test-mu-plain-fit.R (weekly batch).

  .ocmt <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }

  test_that("plain mu-ref pairs are classified and grouped (covariate-free model)", {
    ui <- rxode2::rxode2(.ocmt)
    cls <- nlmixr2est:::.muRefClassify(ui)
    expect_equal(cls$muPlainThetas, c("tka", "tcl", "tv"))
    expect_equal(cls$muPlainEtas, c("eta.ka", "eta.cl", "eta.v"))
    expect_equal(cls$muCovThetas, character(0))

    # plain=FALSE (the default, and what the imp family uses) is unchanged:
    # no groups without a covariate relationship
    expect_length(nlmixr2est:::.muRefGroups(ui), 0L)
    s0 <- nlmixr2est:::.muRefCppGroupSetup(ui)
    expect_equal(s0$muGroupTheta, integer(0))

    g <- nlmixr2est:::.muRefGroups(ui, plain = TRUE)
    expect_length(g, 3L)
    expect_equal(vapply(g, function(x) x$theta, character(1)),
                 c("tka", "tcl", "tv"))
    expect_true(all(vapply(g, function(x) nrow(x$covariates) == 0L, logical(1))))

    s <- nlmixr2est:::.muRefCppGroupSetup(ui, plain = TRUE)
    expect_equal(s$muGroupTheta, 0:2)
    expect_equal(s$muGroupEta, 0:2)
    expect_equal(s$muGroupCovCount, c(0L, 0L, 0L))
    expect_equal(s$muGroupCovTheta, integer(0))
    expect_equal(s$muGroupCovNames, character(0))
    # unbounded model: clamp bounds default to +-Inf
    expect_equal(s$muGroupThetaLower, rep(-Inf, 3))
    expect_equal(s$muGroupThetaUpper, rep(Inf, 3))
    expect_equal(s$muGroupCovLower, numeric(0))
    expect_equal(s$muGroupCovUpper, numeric(0))
  })

  test_that("a shared eta excludes its theta from the plain groups", {
    shared <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v + 0.5 * eta.cl)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }
    cls <- nlmixr2est:::.muRefClassify(rxode2::rxode2(shared))
    # eta.cl enters the model in two positions, so tcl+eta.cl cannot be
    # rewritten by the regression; tka/tv are unaffected
    expect_equal(cls$muPlainThetas, c("tka", "tv"))
    expect_false("tcl" %in% cls$muPlainThetas)
  })

  test_that("bounded plain mu theta is grouped with clamp bounds; fixed excluded silently", {
    bnd <- function() {
      ini({
        tka <- c(-2, 0.45, 2)
        tcl <- 1
        tv <- 3.45
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }
    ui <- rxode2::rxode2(bnd)
    expect_no_warning(g <- suppressMessages(nlmixr2est:::.muRefGroups(ui, plain = TRUE)))
    # bounded tka is regression-updated too (update clamped to [-2, 2])
    expect_equal(vapply(g, function(x) x$theta, character(1)),
                 c("tka", "tcl", "tv"))
    s <- nlmixr2est:::.muRefCppGroupSetup(ui, plain = TRUE)
    expect_equal(s$muGroupThetaLower, c(-2, -Inf, -Inf))
    expect_equal(s$muGroupThetaUpper, c(2, Inf, Inf))
    expect_null(s$muGroupCovBounded)

    fx <- function() {
      ini({
        tka <- fix(0.45)
        tcl <- 1
        tv <- 3.45
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }
    g2 <- nlmixr2est:::.muRefGroups(rxode2::rxode2(fx), plain = TRUE)
    expect_equal(vapply(g2, function(x) x$theta, character(1)), c("tcl", "tv"))
  })

  test_that("mixed covariate + plain groups flatten with aligned covariate arrays", {
    mixed <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        allo.cl <- 0.75
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl + allo.cl * logWT)
        v <- exp(tv + eta.v)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }
    ui <- rxode2::rxode2(mixed)
    s <- nlmixr2est:::.muRefCppGroupSetup(ui, plain = TRUE)
    thNames <- ui$iniDf$name[!is.na(ui$iniDf$ntheta)]
    # covariate group (tcl) first, then the plain groups (tka, tv)
    expect_equal(thNames[s$muGroupTheta + 1L], c("tcl", "tka", "tv"))
    expect_equal(s$muGroupCovCount, c(1L, 0L, 0L))
    expect_equal(s$muGroupCovStart, c(0L, 1L, 1L))
    expect_equal(thNames[s$muGroupCovTheta + 1L], "allo.cl")
    expect_equal(s$muGroupCovNames, "logWT")
  })

  test_that("bounds on a pop theta and a coefficient flatten aligned", {
    mixedBnd <- function() {
      ini({
        tka <- 0.45
        tcl <- c(0, 1, 5)
        tv <- 3.45
        allo.cl <- c(-1, 0.75, 2)
        allo.cl2 <- 0.1
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl + allo.cl * logWT + allo.cl2 * sexf)
        v <- exp(tv + eta.v)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }
    ui <- rxode2::rxode2(mixedBnd)
    expect_no_warning(s <- suppressMessages(
      nlmixr2est:::.muRefCppGroupSetup(ui, plain = TRUE)))
    thNames <- ui$iniDf$name[!is.na(ui$iniDf$ntheta)]
    expect_equal(thNames[s$muGroupTheta + 1L], c("tcl", "tka", "tv"))
    # bounded tcl group carries its clamp bounds; plain groups are infinite
    expect_equal(s$muGroupThetaLower, c(0, -Inf, -Inf))
    expect_equal(s$muGroupThetaUpper, c(5, Inf, Inf))
    # both coefficients (bounded and not) are in the flattened arrays; the
    # bounds stay aligned to their coefficient whatever the flattening order
    .covNm <- thNames[s$muGroupCovTheta + 1L]
    expect_setequal(.covNm, c("allo.cl", "allo.cl2"))
    .w <- match("allo.cl", .covNm)
    expect_equal(s$muGroupCovLower[.w], -1)
    expect_equal(s$muGroupCovUpper[.w], 2)
    .w2 <- match("allo.cl2", .covNm)
    expect_equal(s$muGroupCovLower[.w2], -Inf)
    expect_equal(s$muGroupCovUpper[.w2], Inf)
    expect_equal(s$muGroupCovCount, c(2L, 0L, 0L))
  })

  test_that("muModelClampRetries is validated", {
    expect_equal(foceiControl()$muModelClampRetries, 10L)
    expect_error(foceiControl(muModelClampRetries = 0L))
    expect_equal(irlsfoceiControl(muModelClampRetries = 3L)$muModelClampRetries, 3L)
  })

  test_that("plain groups are dropped when nothing would be left for the outer optimizer", {
    degen <- function() {
      ini({
        tka <- 0.45
        eta.ka ~ fix(0.6)
        add.sd <- fix(0.7)
      })
      model({
        ka <- exp(tka + eta.ka)
        d/dt(depot) <- -ka * depot
        depot ~ add(add.sd)
      })
    }
    g <- suppressMessages(nlmixr2est:::.muRefGroups(rxode2::rxode2(degen), plain = TRUE))
    expect_length(g, 0L)
  })
})
