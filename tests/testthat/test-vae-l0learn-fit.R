## est="vae" covariate selection through vaeControl(covSelectMethod=) end to end:
## does the L0Learn path actually engage, does it report itself, does it agree
## with the exact branch-and-bound, and is it reproducible.  Runs real (short)
## fits and a wider kernel parity sweep, so this file is weekly, not essential.

nmTest({

  ## 48 subjects (theo_sd replicated) with 30 subject-constant nuisance
  ## covariates -- enough columns to cross the covSelectMaxExact threshold and
  ## enough subjects for the regression to be full rank.
  wideData <- function(nCov = 30L) {
    d <- do.call(rbind, lapply(0:3, function(b) {
      x <- nlmixr2data::theo_sd
      x$ID <- x$ID + b * 12L
      x
    }))
    .testSeed(42L)
    ids <- unique(d$ID)
    cv <- matrix(rnorm(length(ids) * nCov), length(ids), nCov)
    for (j in seq_len(nCov)) d[[paste0("C", j)]] <- cv[match(d$ID, ids), j]
    d
  }

  wideModel <- function() {
    ini({tka <- 0.45; tcl <- 1; tv <- 3.45
         eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1; add.sd <- 0.7})
    model({ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
           linCmt() ~ add(add.sd)})
  }

  shortCtl <- function(...) {
    vaeControl(itersBurnIn = 5L, iters = 8L, klWarmup = 4L, gammaIter = 6L,
               nGradStep = 2L, print = 0L, returnVae = TRUE, ...)
  }

  test_that("auto engages L0Learn past the threshold and agrees with the exact search", {
    skip_on_cran()
    skip_if_not_installed("L0Learn")
    skip_if_not_installed("nlmixr2data")
    d <- wideData(30L)

    fL <- suppressWarnings(suppressMessages(
      nlmixr2(wideModel, d, est = "vae", control = shortCtl())))
    ## mechanism, not just the result: the approximate path really ran
    expect_identical(fL$covSelectMethodUsed, "l0learn")

    fB <- suppressWarnings(suppressMessages(
      nlmixr2(wideModel, d, est = "vae", control = shortCtl(covSelectMethod = "bnb"))))
    expect_identical(fB$covSelectMethodUsed, "bnb")

    ## the two searches optimize the same objective, so they must agree on the
    ## selected covariate set
    expect_identical(fL$selected, fB$selected)
  })

  test_that("below the threshold nothing switches and nothing is said", {
    skip_on_cran()
    skip_if_not_installed("L0Learn")
    skip_if_not_installed("nlmixr2data")
    d <- wideData(6L)
    f <- suppressMessages(nlmixr2(wideModel, d, est = "vae", control = shortCtl()))
    expect_identical(f$covSelectMethodUsed, "bnb")
  })

  test_that("an explicit covSelectMethod overrides the threshold both ways", {
    skip_on_cran()
    skip_if_not_installed("L0Learn")
    skip_if_not_installed("nlmixr2data")
    d <- wideData(6L)
    ## forced on below the threshold
    f <- suppressWarnings(suppressMessages(
      nlmixr2(wideModel, d, est = "vae", control = shortCtl(covSelectMethod = "l0learn"))))
    expect_identical(f$covSelectMethodUsed, "l0learn")
    ## raising the threshold keeps a wide problem exact
    d30 <- wideData(30L)
    f <- suppressMessages(
      nlmixr2(wideModel, d30, est = "vae", control = shortCtl(covSelectMaxExact = 100L)))
    expect_identical(f$covSelectMethodUsed, "bnb")
    ## covSelectMaxExact=Inf forces the exact search everywhere
    f <- suppressMessages(
      nlmixr2(wideModel, d30, est = "vae", control = shortCtl(covSelectMaxExact = Inf)))
    expect_identical(f$covSelectMethodUsed, "bnb")
  })

  test_that("the approximate search reports itself in runInfo", {
    skip_on_cran()
    skip_if_not_installed("L0Learn")
    skip_if_not_installed("nlmixr2data")
    d <- wideData(30L)
    f <- suppressMessages(nlmixr2(wideModel, d, est = "vae",
                                  control = vaeControl(itersBurnIn = 4L, iters = 5L,
                                                       klWarmup = 3L, gammaIter = 4L,
                                                       nGradStep = 2L, print = 0L,
                                                       calcTables = FALSE)))
    expect_identical(f$vae$covSelectMethodUsed, "l0learn")
    ## a non-exact selection must never arrive silently
    expect_match(f$runInfo, "covariate search used L0Learn", all = FALSE)
  })

  test_that("the L0Learn path is reproducible", {
    skip_on_cran()
    skip_if_not_installed("L0Learn")
    skip_if_not_installed("nlmixr2data")
    d <- wideData(30L)
    a <- suppressWarnings(suppressMessages(
      nlmixr2(wideModel, d, est = "vae", control = shortCtl())))
    b <- suppressWarnings(suppressMessages(
      nlmixr2(wideModel, d, est = "vae", control = shortCtl())))
    expect_identical(a$selected, b$selected)
    expect_identical(a$zPop, b$zPop)
    expect_identical(a$omega, b$omega)
  })

  test_that("L0Learn candidates plus polish reproduce the exact optimum at scale", {
    skip_on_cran()
    skip_if_not_installed("L0Learn")
    ## wider version of the essential-suite parity check: more replicates, more
    ## covariates, and both independent and rho=0.7 correlated designs
    .testSeed(2026L)
    N <- 150L
    nRep <- 50L
    for (rep in seq_len(nRep)) {
      nCov <- sample(8:16, 1L)
      X <- matrix(rnorm(N * nCov), N, nCov)
      if (rep %% 2L == 0L) {
        f <- rnorm(N)
        X <- sqrt(0.7) * matrix(f, N, nCov) + sqrt(0.3) * X
      }
      k <- sample(0:5, 1L)
      sel <- if (k > 0) sort(sample.int(nCov, k)) else integer(0)
      y <- as.numeric(0.5 + (if (k > 0) X[, sel, drop = FALSE] %*%
                               (runif(k, 1, 3) * sample(c(-1, 1), k, TRUE)) else 0) +
                        rnorm(N, sd = runif(1, 0.4, 1.5)))
      omega <- runif(1, 0.2, 1)
      penalty <- log(N)
      got <- vaeScoreSupports_(y, X, omega, penalty,
                               nlmixr2est:::.vaeL0Supports(X, y), polish = TRUE)
      ref <- vaeBestSubset_(matrix(y, ncol = 1), X, omega, FALSE, penalty)
      expect_identical(which(got$selected == 1L), which(ref$selected[1, ] == 1L),
                       info = paste0("rep ", rep, " nCov ", nCov))
    }
  })
})
