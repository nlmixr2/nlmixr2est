test_that(".covStoreKey() needs the hand-off and keeps every key field", {
  .env <- new.env(parent = emptyenv())
  expect_null(.covStoreKey(.env, list(hessEps = 1e-4)))
  .env$covHandoff <- list(theta = c(1, 2), omega = 0.5)
  .k <- .covStoreKey(.env, list(hessEps = 1e-4, gillKcov = 10L, unrelated = 3))
  expect_identical(.k$handoff, .env$covHandoff)
  expect_named(.k$versions, c("nlmixr2est", "rxode2"))
  expect_identical(names(.k$settings), c(.covStoreKeyFields, "covEtaLegs"))
  expect_identical(.k$settings$hessEps, 1e-4)
  expect_identical(.k$settings$gillKcov, 10L)
  expect_null(.k$settings$covSolveTol)
  expect_false("unrelated" %in% names(.k$settings))
  # the ETA policy of the legs: a fit's own budget, a held-ETA fit's refit budget, or fixed
  expect_identical(.k$settings$covEtaLegs, 0L)
  expect_identical(.covStoreKey(.env, list(maxInnerIterations = 1000))$settings$covEtaLegs, 1000L)
  expect_identical(
    .covStoreKey(.env, list(maxInnerIterations = 0L, covMaxInnerIterations = 500L))$settings$covEtaLegs,
    500L
  )
  expect_identical(.covStoreKey(.env, list(maxInnerIterations = 0L))$settings$covEtaLegs, 0L)
})

test_that(".covStoreRecord() merges what each covariance step computed under its key", {
  .fit <- new.env(parent = emptyenv())
  .fit$covHandoff <- list(theta = 1, omega = 2)
  .key <- .covStoreKey(.fit, list(hessEps = 1e-4))
  .R <- matrix(c(2, 0.1, 0.1, 3), 2, dimnames = list(c("a", "b"), c("a", "b")))
  .steps <- list(theta = 1, f0 = 10, aEps = 0.1, rEps = 0.1, aEpsC = 0.2, rEpsC = 0.2)
  # nothing computed: nothing written
  .none <- new.env(parent = emptyenv())
  expect_false(.covStoreRecord(.fit, .key, .none))
  expect_null(.fit$covStore)
  # an "r" step: the full R and steps, and the theta-only R and steps
  .r <- new.env(parent = emptyenv())
  .r$.fdFullR <- .R
  .r$.fdFullH <- c(0.01, 0.02)
  .r$.fdFullX0 <- c(1, 2)
  .r$covSteps <- .steps
  .r$R.0 <- matrix(2)
  expect_true(.covStoreRecord(.fit, .key, .r))
  expect_length(.fit$covStore, 1L)
  .e <- .covStoreGet(.fit, .key)
  expect_identical(.e$full$R, .R)
  expect_null(.e$full$S)
  expect_identical(.e$theta$R0, matrix(2))
  expect_null(.e$theta$S0)
  # an "s" step at the same steps adds S and keeps R
  .s <- new.env(parent = emptyenv())
  .s$.fdFullR <- .R
  .s$.fdFullH <- c(0.01, 0.02)
  .s$.fdFullX0 <- c(1, 2)
  .s$.fdFullS <- diag(2)
  .s$covSteps <- .steps
  .s$S0 <- matrix(4)
  .s$Sper <- 1
  .s$SHasZero <- TRUE
  expect_true(.covStoreRecord(.fit, .key, .s))
  expect_length(.fit$covStore, 1L)
  .e <- .covStoreGet(.fit, .key)
  expect_identical(.e$full$S, diag(2))
  expect_identical(.e$theta$R0, matrix(2))
  expect_identical(.e$theta$S0, matrix(4))
  expect_identical(.e$theta$Sper, 1)
  expect_true(.e$theta$SHasZero)
  expect_identical(.e$full$x0, c(1, 2))
  # a later step without S keeps the stored S when its R, steps and point are the same
  expect_true(.covStoreRecord(.fit, .key, .r))
  expect_identical(.covStoreGet(.fit, .key)$full$S, diag(2))
  .moved <- new.env(parent = emptyenv())
  for (.n in ls(.r, all.names = TRUE)) assign(.n, get(.n, envir = .r), envir = .moved)
  .moved$.fdFullX0 <- c(1, 2.5)
  expect_true(.covStoreRecord(.fit, .key, .moved))
  expect_null(.covStoreGet(.fit, .key)$full$S)
  expect_true(.covStoreRecord(.fit, .key, .s))
  # a full stage without its point is not stored
  .nox <- new.env(parent = emptyenv())
  .nox$.fdFullR <- .R
  .nox$.fdFullH <- c(0.01, 0.02)
  .fitNox <- new.env(parent = emptyenv())
  .fitNox$covHandoff <- .fit$covHandoff
  expect_false(.covStoreRecord(.fitNox, .key, .nox))
  # an analytic theta-only R is not stored as a finite-difference one
  .fit2 <- new.env(parent = emptyenv())
  .fit2$covHandoff <- .fit$covHandoff
  expect_true(.covStoreRecord(.fit2, .key, .r, fd = FALSE))
  expect_null(.covStoreGet(.fit2, .key)$theta$R0)
  # another key is another entry
  .key2 <- .covStoreKey(.fit, list(hessEps = 1e-5))
  expect_null(.covStoreGet(.fit, .key2))
  expect_true(.covStoreRecord(.fit, .key2, .s))
  expect_length(.fit$covStore, 2L)
  expect_identical(.covStoreIndex(.fit$covStore, .key2), 2L)
  expect_identical(.covStoreIndex(.fit$covStore, NULL), 0L)
  expect_null(.covStoreGet(.fit, NULL))
})

test_that("a refit with settings outside the key does not use the store", {
  expect_true(.covStoreRefitOk(list()))
  expect_true(.covStoreRefitOk(list(covMethod = "r,s", covFull = TRUE, covSmall = 1e-5, hessEps = 1e-4)))
  expect_false(.covStoreRefitOk(list(covMethod = "r", rxControl = list(atol = 1e-10))))
  expect_false(.covStoreRefitOk(list(interaction = 0L)))
})

test_that("every fit starts with an empty covariance store and only setCov()'s inputs", {
  .env <- new.env(parent = emptyenv())
  # left from an earlier fit in the same environment
  .env$covStore <- list(list(key = 1))
  .env$covHandoff <- list(theta = 1, omega = 2)
  .env$.fdFullStore <- list(R = diag(2))
  .env$covThetaStore <- list(R0 = diag(1))
  .covStoreFitStart(.env)
  expect_true(exists("covStore", envir = .env, inherits = FALSE))
  expect_null(.env$covStore)
  for (.n in c("covHandoff", ".fdFullStore", "covThetaStore")) {
    expect_false(exists(.n, envir = .env, inherits = FALSE))
  }
  # a setCov() refit hands its fit's in explicitly
  .env$covStore <- list(list(key = 1))
  .env$.covRefitInputs <- list(covHandoff = list(theta = 3, omega = 4), .fdFullStore = list(R = diag(3)), covThetaStore = NULL)
  .covStoreFitStart(.env)
  expect_null(.env$covStore)
  expect_identical(.env$covHandoff, list(theta = 3, omega = 4))
  expect_identical(.env$.fdFullStore, list(R = diag(3)))
  expect_false(exists("covThetaStore", envir = .env, inherits = FALSE))
  expect_false(exists(".covRefitInputs", envir = .env, inherits = FALSE))
})

nmTest({
  .storeOneCmt <- function() {
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
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }
  .storeRefit <- function(fit, reuse, ...) {
    withr::local_envvar(NLMIXR2EST_COV_NO_REUSE = if (reuse) "" else "1")
    suppressMessages(.setCovRefit(fit, ...))
  }

  test_that("a theta-only covariance request reads what the fit's covariance step computed", {
    skip_on_cran()
    .f <- .nlmixr(
      .storeOneCmt,
      theo_sd,
      "focei",
      foceiControl(print = 0, calcTables = FALSE, covMethod = "r", covFull = FALSE)
    )
    expect_length(.f$env$covStore, 1L)
    expect_false(is.null(.f$env$covStore[[1]]$theta$R0))
    # "r,s" computes only S: no step search and no R stencil
    .on <- .storeRefit(.f, TRUE, covMethod = "r,s", covFull = FALSE)
    .off <- .storeRefit(.f, FALSE, covMethod = "r,s", covFull = FALSE)
    expect_identical(.on$covMethod, "r,s")
    expect_identical(.on$cov, .off$cov)
    expect_identical(unname(.on$env$covEvals[c("gill", "r")]), c(0L, 0L))
    expect_gt(.off$env$covEvals[["r"]], 0L)
    # the S it computed is now stored too, so "s" computes nothing beyond the centre
    expect_false(is.null(.f$env$covStore[[1]]$theta$S0))
    .s <- .storeRefit(.f, TRUE, covMethod = "s", covFull = FALSE)
    expect_identical(.s$cov, .storeRefit(.f, FALSE, covMethod = "s", covFull = FALSE)$cov)
    expect_identical(unname(.s$env$covEvals[c("gill", "r", "s")]), c(0L, 0L, 0L))
  })

  test_that("a full-shape request reads the stored full stage and gives the fit's own result", {
    skip_on_cran()
    .f <- .nlmixr(
      .storeOneCmt,
      theo_sd,
      "focei",
      foceiControl(print = 0, calcTables = FALSE, covMethod = "r", covFull = TRUE)
    )
    expect_false(is.null(.f$env$covStore[[1]]$full$R))
    expect_null(.f$env$covStore[[1]]$full$S)
    .on <- .storeRefit(.f, TRUE, covMethod = "r,s", covFull = TRUE)
    # the same, bit for bit, as a fit that asked for "r,s" itself
    .nat <- .nlmixr(
      .storeOneCmt,
      theo_sd,
      "focei",
      foceiControl(print = 0, calcTables = FALSE, covMethod = "r,s", covFull = TRUE)
    )
    expect_identical(.on$covMethod, "r,s (full)")
    expect_identical(.on$cov, .nat$cov)
    # only the full S legs: 2 per parameter (4 thetas, 3 omegas), no step search or stencil
    expect_identical(unname(.on$env$covEvals[c("fullGill", "fullR", "fullS")]), c(0L, 0L, 14L))
    expect_identical(.f$env$covStore[[1]]$full$S, .nat$env$.fdFullS)
    # a changed step setting is another key: everything is computed again
    .k <- .storeRefit(.f, TRUE, covMethod = "r,s", covFull = TRUE, gillKcov = 5L)
    expect_gt(.k$env$covEvals[["fullR"]], 0L)
    expect_length(.f$env$covStore, 2L)
    # a setting outside the key neither reads nor writes the store
    .o <- .storeRefit(.f, TRUE, covMethod = "r,s", covFull = TRUE, epsilon = 1e-9)
    expect_gt(.o$env$covEvals[["fullR"]], 0L)
    expect_length(.f$env$covStore, 2L)
    # an entry taken about another point is not read
    .st <- .f$env$covStore
    .st[[1]]$full$x0[1] <- .st[[1]]$full$x0[1] * (1 + 1e-9)
    .f$env$covStore <- .st
    .x <- .storeRefit(.f, TRUE, covMethod = "r", covFull = TRUE)
    expect_gt(.x$env$covEvals[["fullR"]], 0L)
  })
})
