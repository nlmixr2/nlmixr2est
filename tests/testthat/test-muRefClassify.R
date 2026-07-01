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
