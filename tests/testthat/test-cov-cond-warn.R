nmTest({
  # covariance whose correlation matrix has condition number (1 + r) / (1 - r)
  .condCov <- function(r) {
    .s <- c(tcl = 0.1, tv = 0.2)
    .c <- matrix(c(1, r, r, 1), 2) * outer(.s, .s)
    dimnames(.c) <- list(names(.s), names(.s))
    .c
  }
  .condFit <- function(cov, condWarn = 1000, cond = NULL) {
    .env <- new.env(parent = emptyenv())
    .env$cov <- cov
    if (!is.null(cond)) {
      .env$conditionNumberCor <- cond
    }
    list(env = .env, control = list(condWarn = condWarn))
  }

  test_that("the condition number is that of the correlation matrix", {
    expect_equal(.covCondCor(.condCov(0.5)), 3)
    expect_equal(.covCondCor(.condCov(0.999)), 1999, tolerance = 1e-8)
    .c <- .condCov(0.999)
    .c["tv", "tv"] <- NA_real_
    expect_equal(.covCondCor(.c), 1)
    expect_true(is.na(.covCondCor(NULL)))
  })

  test_that("only a condition number above the threshold warns, under 75 characters", {
    expect_null(.covCondMsg(3, 1000))
    expect_null(.covCondMsg(1000, 1000))
    expect_null(.covCondMsg(NA_real_, 1000))
    expect_equal(
      .covCondMsg(1397.7, 1000),
      "covariance condition number is large (1397.7 > 1000)"
    )
    expect_lt(nchar(.covCondMsg(1.234567e15, 1000)), 75L)
    expect_null(.covCondMsg(1e6, NA_real_))
    expect_null(.covCondMsg(1e6, Inf))
  })

  test_that("the fit's covariance warns above the threshold and not below it", {
    .w <- .covSite(".covCondWarn", capture_warnings(.covCondWarn(.condFit(.condCov(0.999)))))
    expect_identical(.w, "covariance condition number is large (1999 > 1000)")
    expect_no_warning(expect_null(.covCondWarn(.condFit(.condCov(0.5)))))
    expect_no_warning(expect_null(.covCondWarn(.condFit(.condCov(0.999), NA_real_))))
    expect_no_warning(expect_null(.covCondWarn(.condFit(.condCov(0.999), Inf))))
    expect_no_warning(expect_null(.covCondWarn(new.env(parent = emptyenv()))))
  })

  test_that("the installed covariance is used, not a stale stored number", {
    .w <- .covSite(
      ".covCondWarn",
      capture_warnings(.covCondWarn(.condFit(.condCov(0.999), cond = 3)))
    )
    expect_identical(.w, "covariance condition number is large (1999 > 1000)")
    expect_no_warning(expect_null(.covCondWarn(.condFit(.condCov(0.5), cond = 1397.7))))
  })

  test_that("the printed number exceeds the printed threshold", {
    expect_equal(
      .covCondMsg(1000.01, 1000),
      "covariance condition number is large (1000.01 > 1000)"
    )
    expect_match(.covCondMsg(1000 + 1e-9, 1000), "\\(1000\\.0+1 > 1000\\)$")
    expect_null(.covCondMsg(1000 + 2^-43, 1000))
    withr::with_options(list(OutDec = ","), {
      expect_equal(
        .covCondMsg(1397.7, 1000),
        "covariance condition number is large (1397,7 > 1000)"
      )
    })
  })

  test_that("condWarn is validated by foceiControl", {
    expect_equal(foceiControl()$condWarn, 1000)
    expect_equal(foceiControl(condWarn = NA)$condWarn, NA)
    expect_equal(foceiControl(condWarn = Inf)$condWarn, Inf)
    expect_error(foceiControl(condWarn = 0.5))
    expect_error(foceiControl(condWarn = c(10, 100)))
    expect_error(foceiControl(condWarn = "a"))
  })

  test_that("the warning reaches $runInfo and condWarn = NA suppresses it", {
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
    # condWarn = 1 flags any covariance, so the test does not hinge on the estimates
    fit <- .nlmixr(one.compartment, nlmixr2data::theo_sd, "focei", foceiControl(print = 0, condWarn = 1))
    expect_true(is.matrix(fit$cov))
    .ri <- grep("^covariance condition number is large", fit$runInfo, value = TRUE)
    expect_length(.ri, 1L)
    expect_lt(nchar(.ri), 75L)
    fit2 <- .nlmixr(one.compartment, nlmixr2data::theo_sd, "focei", foceiControl(print = 0, condWarn = NA))
    expect_false(any(grepl("^covariance condition number", fit2$runInfo)))
  })
})
