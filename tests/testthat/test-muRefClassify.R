test_that(".muRefClassify buckets theta/eta pairs correctly", {

  mod <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      allo.cl <- 0.75      # covariate coefficient WITH an eta (cl has eta.cl)
      tf1 <- 0             # theta WITH covariate, NO eta
      allo.f1 <- 0.1
      tstruct <- 2         # theta by itself: no eta, no covariate
      eta.ka ~ 0.6         # eta by itself: no covariate
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl + allo.cl * logWT)
      v <- exp(tv + eta.v)
      f1 <- expit(tf1 + allo.f1 * logWT)
      dummy <- tstruct
      linCmt() ~ add(add.sd)
    })
  }

  ui <- rxode2::rxode2(mod)
  .cls <- .muRefClassify(ui)

  # theta+eta+covariate -> mu-ref covariate eligible
  expect_true("tcl" %in% .cls$muCovThetas)
  expect_true("eta.cl" %in% .cls$muCovEtas)
  expect_true("allo.cl" %in% .cls$muCovCovariateParams)

  # theta+covariate, no eta -> also mu-ref covariate eligible (no eta produced)
  expect_true("tf1" %in% .cls$muCovThetas)
  expect_true("allo.f1" %in% .cls$muCovCovariateParams)

  # eta by itself (no covariate) -> standard, not touched
  expect_true("eta.ka" %in% .cls$standardEtas)
  expect_true("eta.v" %in% .cls$standardEtas)
  expect_false("eta.ka" %in% .cls$muCovEtas)
  expect_false("eta.v" %in% .cls$muCovEtas)

  # theta by itself (no eta, no covariate) -> standard
  expect_true("tstruct" %in% .cls$standardThetas)

  # plain thetas with an eta but no covariate ("eta by itself" from the
  # theta side) -> standard, not mu-ref covariate eligible
  expect_true("tka" %in% .cls$standardThetas)
  expect_true("tv" %in% .cls$standardThetas)
  expect_false("tka" %in% .cls$muCovThetas)
  expect_false("tv" %in% .cls$muCovThetas)

  # sigma/residual-error parameters are never classified as thetas at all
  expect_false("add.sd" %in% .cls$standardThetas)
  expect_false("add.sd" %in% .cls$muCovThetas)

  # covariate-coefficient thetas ("slopes") must never leak into
  # standardThetas -- they are accounted for separately
  expect_false("allo.cl" %in% .cls$standardThetas)
  expect_false("allo.f1" %in% .cls$standardThetas)

  # no double counting: a theta cannot be in both muCovThetas and standardThetas
  expect_length(intersect(.cls$muCovThetas, .cls$standardThetas), 0L)
  expect_length(intersect(.cls$muCovEtas, .cls$standardEtas), 0L)
})

test_that(".muRefClassify handles a model with no mu-ref covariates at all", {

  mod <- function() {
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
      linCmt() ~ add(add.sd)
    })
  }

  ui <- rxode2::rxode2(mod)
  .cls <- .muRefClassify(ui)

  expect_length(.cls$muCovThetas, 0L)
  expect_length(.cls$muCovEtas, 0L)
  expect_length(.cls$muCovCovariateParams, 0L)
  expect_setequal(.cls$standardEtas, c("eta.ka", "eta.cl", "eta.v"))
  expect_true(all(c("tka", "tcl", "tv") %in% .cls$standardThetas))
})

test_that(".muRefClassify handles a model that is entirely mu-ref covariate", {

  mod <- function() {
    ini({
      tcl <- 1
      allo.cl <- 0.75
      eta.cl ~ 0.3
      add.sd <- 0.7
    })
    model({
      cl <- exp(tcl + eta.cl + allo.cl * logWT)
      v <- 10
      linCmt() ~ add(add.sd)
    })
  }

  ui <- rxode2::rxode2(mod)
  .cls <- .muRefClassify(ui)

  expect_setequal(.cls$muCovThetas, "tcl")
  expect_setequal(.cls$muCovEtas, "eta.cl")
  expect_setequal(.cls$muCovCovariateParams, "allo.cl")
  expect_length(.cls$standardEtas, 0L)
})

test_that(".muRefGroups keeps the group but marks a bounded covariate coefficient, with a warning", {
  # The closed-form/IRLS regression is an unconstrained solve -- it has no
  # way to respect a box constraint on a covariate coefficient -- so a
  # bounded covariate is treated as if it were time-varying and excluded
  # from the regression (marked bounded=TRUE, handled by the ordinary,
  # bounded outer optimizer instead), while the population theta (the
  # group's intercept) still benefits from the mu-ref speed-up. This is
  # narrower than a bounded *population* theta, which drops the whole
  # group (see the next test) -- a bound on a slope doesn't prevent
  # solving for the intercept.
  mod <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      allo.cl <- c(0, 0.75, 2) # bounded covariate coefficient
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl + allo.cl * logWT)
      v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }
  ui <- rxode2::rxode2(mod)
  expect_warning(.groups <- .muRefGroups(ui), "allo\\.cl.*boundar")
  expect_length(.groups, 1L)
  expect_equal(.groups[[1]]$theta, "tcl")
  expect_true(.groups[[1]]$covariates$bounded[.groups[[1]]$covariates$covariateParameter == "allo.cl"])
})

test_that(".muRefGroups excludes a group and warns when its population theta has boundaries", {
  mod <- function() {
    ini({
      tka <- 0.45
      tcl <- c(0, 1, 5) # bounded population theta
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
      linCmt() ~ add(add.sd)
    })
  }
  ui <- rxode2::rxode2(mod)
  expect_warning(.groups <- .muRefGroups(ui), "tcl.*boundar")
  expect_length(.groups, 0L)
})

test_that(".muRefGroups does not warn or exclude anything when there are no boundaries", {
  mod <- function() {
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
      linCmt() ~ add(add.sd)
    })
  }
  ui <- rxode2::rxode2(mod)
  expect_warning(.groups <- .muRefGroups(ui), NA)
  expect_length(.groups, 1L)
  expect_equal(.groups[[1]]$theta, "tcl")
})

test_that(".muRefGroups only marks the bounded covariate when a group has several, mixed covariates", {
  mod <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      allo.cl <- 0.75 # unbounded
      allo.cl2 <- c(0, 0.1, 1) # bounded
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl + allo.cl * logWT + allo.cl2 * sexf)
      v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }
  ui <- rxode2::rxode2(mod)
  expect_warning(.groups <- .muRefGroups(ui), "allo\\.cl2.*boundar")
  # the group (population theta + the *unbounded* covariate) stays
  # mu-ref eligible; only the bounded covariate is carved out
  expect_length(.groups, 1L)
  .cv <- .groups[[1]]$covariates
  expect_false(.cv$bounded[.cv$covariateParameter == "allo.cl"])
  expect_true(.cv$bounded[.cv$covariateParameter == "allo.cl2"])
})

test_that(".muRefCppCovData evaluates covariate expressions that are not data columns (#711)", {
  dataSav <- data.frame(
    ID = c(2L, 2L, 1L, 1L),
    WT = c(80, 80, 60, 60),
    logWT = log(c(80, 80, 60, 60) / 70)
  )
  # bare columns still subset directly
  .m <- .muRefCppCovData(c("logWT", "WT"), dataSav)
  expect_equal(dim(.m), c(2L, 2L))
  expect_equal(.m[, 1], log(c(60, 80) / 70))
  expect_equal(.m[, 2], c(60, 80))
  # expression names (e.g. from >=2 mu-ref covariate expressions) are
  # evaluated against the baseline rows instead of erroring
  .m2 <- .muRefCppCovData(c("log(WT/70)", "log(WT/70)", "WT"), dataSav)
  expect_equal(dim(.m2), c(2L, 3L))
  expect_equal(.m2[, 1], log(c(60, 80) / 70))
  expect_equal(.m2[, 2], .m2[, 1])
  expect_equal(.m2[, 3], c(60, 80))
})
