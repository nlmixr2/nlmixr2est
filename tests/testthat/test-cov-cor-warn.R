nmTest({
  .corCov <- function(r = 0.99) {
    .s <- c(cl = 0.1, v = 0.2, ka = 0.3)
    .r <- diag(3)
    .r[1, 2] <- .r[2, 1] <- r
    .c <- .r * outer(.s, .s)
    dimnames(.c) <- list(names(.s), names(.s))
    .c
  }

  test_that("a 0.99 pair is named, strongest first, under 75 characters", {
    .p <- .covCorPairs(.corCov(), 0.95)
    expect_equal(nrow(.p), 1L)
    expect_equal(c(.p$p1, .p$p2), c("cl", "v"))
    expect_equal(.covCorMsg(.p), "highly correlated estimates: cl/v (0.99)")
    expect_equal(.covCorMsg(.covCorPairs(.corCov(-0.97), 0.95)), "highly correlated estimates: cl/v (-0.97)")
    expect_null(.covCorMsg(.covCorPairs(.corCov(0.9), 0.95)))
  })

  test_that("NA and Inf turn the warning off", {
    expect_equal(nrow(.covCorPairs(.corCov(), NA_real_)), 0L)
    expect_equal(nrow(.covCorPairs(.corCov(), Inf)), 0L)
    expect_error(foceiControl(corWarn = -1))
    expect_error(foceiControl(corWarn = c(0.9, 0.95)))
    expect_equal(foceiControl(corWarn = NA)$corWarn, NA)
  })

  test_that("zero-variance, NA-variance and fixed parameters are skipped", {
    .c <- .corCov()
    .c["cl", ] <- .c[, "cl"] <- 0
    expect_equal(nrow(.covCorPairs(.c, 0.95)), 0L)
    .c <- .corCov()
    .c["v", "v"] <- NA_real_
    expect_equal(nrow(.covCorPairs(.c, 0.95)), 0L)
    expect_equal(nrow(.covCorPairs(.corCov(), 0.95, skip = "v")), 0L)
  })

  test_that("a fixed omega is skipped under its full-covariance name", {
    .c <- .corCov()
    dimnames(.c) <- list(c("tcl", "om.eta.v", "tka"), c("tcl", "om.eta.v", "tka"))
    .env <- new.env(parent = emptyenv())
    .env$cov <- .c
    .fit <- list(
      env = .env,
      control = list(corWarn = 0.95),
      iniDf = data.frame(name = c("tcl", "eta.v"), fix = c(FALSE, TRUE))
    )
    expect_no_warning(expect_null(.covCorWarn(.fit)))
    .fit$iniDf$fix[2] <- FALSE
    expect_warning(.covCorWarn(.fit), "tcl/om.eta.v \\(0.99\\)")
  })

  test_that("a missing covariance gives no warning", {
    .env <- new.env(parent = emptyenv())
    expect_no_warning(expect_null(.covCorWarn(.env)))
    expect_equal(nrow(.covCorPairs(NULL, 0.95)), 0L)
  })

  test_that("the listed pairs are capped", {
    .s <- setNames(rep(0.1, 4), c("a", "b", "c", "d"))
    .c <- matrix(0.99, 4, 4, dimnames = list(names(.s), names(.s)))
    diag(.c) <- 1
    .c <- .c * 0.01
    .m <- .covCorMsg(.covCorPairs(.c, 0.95))
    expect_match(.m, ", \\.\\.\\.$")
    expect_lt(nchar(.m), 75L)
  })

  test_that("the warning reaches $runInfo and corWarn = NA suppresses it", {
    skip_on_cran()
    one.compartment <- function() {
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
    # corWarn = 0 names every pair, so the test does not hinge on the estimates
    fit <- .nlmixr(one.compartment, nlmixr2data::theo_sd, "focei", foceiControl(print = 0, corWarn = 0))
    expect_true(is.matrix(fit$cov))
    expect_true(any(grepl("^highly correlated estimates: ", fit$runInfo)))
    fit2 <- .nlmixr(one.compartment, nlmixr2data::theo_sd, "focei", foceiControl(print = 0, corWarn = NA))
    expect_false(any(grepl("^highly correlated estimates", fit2$runInfo)))
  })
})
