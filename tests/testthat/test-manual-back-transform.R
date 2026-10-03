nmTest({
  test_that("manual back-transform", {
    # Throughout these tests, don't use .nlmixr() for quiet testing so that
    # `t100()` is in the environment.
    t100 <- function(x) {
      x * 100
    }

    one.cmt <- function() {
      ini({
        tka <- fix(0.45); backTransform("t100")
        tcl <- log(c(0, 2.7, 100)); backTransform("t100")
        tv <- 3.45; backTransform("none")
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

    fit <-
      suppressMessages(nlmixr(
        one.cmt,
        theo_sd,
        est = "saem",
        control = saemControlFast
      ))

    expect_equal(
      setNames(fit$parFixedDf["tka", "Estimate"] * 100, NULL),
      setNames(fit$parFixedDf["tka", "Back-transformed"], NULL)
    )

    expect_equal(
      fit$parFixed["tka", "Back-transformed(95%CI)"],
      formatMinWidth(fit$parFixedDf["tka", "Estimate"] * 100)
    )

    expect_equal(
      setNames(fit$parFixedDf["tcl", "Estimate"] * 100, NULL),
      setNames(fit$parFixedDf["tcl", "Back-transformed"], NULL)
    )

    qn <- qnorm(1.0 - (1 - 0.95) / 2)

    expect_equal(
      setNames(t100(fit$parFixedDf["tcl", "Estimate"] + qn * fit$parFixedDf["tcl", "SE"]), NULL),
      setNames(fit$parFixedDf["tcl", "CI Upper"], NULL)
    )

    expect_equal(
      setNames(t100(fit$parFixedDf["tcl", "Estimate"] - qn * fit$parFixedDf["tcl", "SE"]), NULL),
      setNames(fit$parFixedDf["tcl", "CI Lower"], NULL)
    )

    expect_equal(
      fit$parFixed["tcl", "Back-transformed(95%CI)"],
      sprintf(
        "%s (%s, %s)",
        formatMinWidth(t100(fit$parFixedDf["tcl", "Estimate"])),
        formatMinWidth(t100(fit$parFixedDf["tcl", "Estimate"] - qn * fit$parFixedDf["tcl", "SE"])),
        formatMinWidth(t100(fit$parFixedDf["tcl", "Estimate"] + qn * fit$parFixedDf["tcl", "SE"]))
      )
    )

    expect_equal(
      setNames(exp(fit$parFixedDf["tv", "Estimate"] + qn * fit$parFixedDf["tv", "SE"]), NULL),
      setNames(fit$parFixedDf["tv", "CI Upper"], NULL)
    )

    expect_equal(
      setNames(exp(fit$parFixedDf["tv", "Estimate"] - qn * fit$parFixedDf["tv", "SE"]), NULL),
      setNames(fit$parFixedDf["tv", "CI Lower"], NULL)
    )

    qn <- qnorm(1.0 - (1 - 0.80) / 2)

    # Test difference confidence intervals
    fit <-
      suppressMessages(nlmixr(
        one.cmt,
        theo_sd,
        est = "saem",
        control = saemControl(print = 0, nBurn = 1, nEm = 1, ci = 0.8)
      ))

    expect_equal(
      setNames(fit$parFixedDf["tka", "Estimate"] * 100, NULL),
      setNames(fit$parFixedDf["tka", "Back-transformed"], NULL)
    )

    expect_equal(
      fit$parFixed["tka", "Back-transformed(80%CI)"],
      formatMinWidth(fit$parFixedDf["tka", "Estimate"] * 100)
    )

    expect_equal(NA_real_, setNames(fit$parFixedDf["tka", "CI Lower"], NULL))
    expect_equal(NA_real_, setNames(fit$parFixedDf["tka", "CI Upper"], NULL))

    expect_equal(
      setNames(fit$parFixedDf["tcl", "Estimate"] * 100, NULL),
      setNames(fit$parFixedDf["tcl", "Back-transformed"], NULL)
    )

    expect_equal(
      setNames(t100(fit$parFixedDf["tcl", "Estimate"] + qn * fit$parFixedDf["tcl", "SE"]), NULL),
      setNames(fit$parFixedDf["tcl", "CI Upper"], NULL)
    )

    expect_equal(
      setNames(t100(fit$parFixedDf["tcl", "Estimate"] - qn * fit$parFixedDf["tcl", "SE"]), NULL),
      setNames(fit$parFixedDf["tcl", "CI Lower"], NULL)
    )

    expect_equal(
      fit$parFixed["tcl", "Back-transformed(80%CI)"],
      sprintf(
        "%s (%s, %s)",
        formatMinWidth(t100(fit$parFixedDf["tcl", "Estimate"])),
        formatMinWidth(t100(fit$parFixedDf["tcl", "Estimate"] - qn * fit$parFixedDf["tcl", "SE"])),
        formatMinWidth(t100(fit$parFixedDf["tcl", "Estimate"] + qn * fit$parFixedDf["tcl", "SE"]))
      )
    )

    expect_equal(
      setNames(exp(fit$parFixedDf["tv", "Estimate"] + qn * fit$parFixedDf["tv", "SE"]), NULL),
      setNames(fit$parFixedDf["tv", "CI Upper"], NULL)
    )

    expect_equal(
      setNames(exp(fit$parFixedDf["tv", "Estimate"] - qn * fit$parFixedDf["tv", "SE"]), NULL),
      setNames(fit$parFixedDf["tv", "CI Lower"], NULL)
    )

    expect_error(
      suppressMessages(
        nlmixr(
          one.cmt,
          theo_sd,
          est = "focei",
          control = foceiControl(print = 0, maxOuterIterations = 0, maxInnerIterations = 0, covMethod = "")
        )
      ),
      NA
    )
  })

  test_that("a refreshed covariance refreshes a manually back-transformed CI (issue 1140)", {
    t100 <- function(x) {
      x * 100
    }
    one.cmt <- function() {
      ini({
        tka <- 0.45
        tcl <- 1; backTransform("t100")
        tv <- 3.45; backTransform("none")
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
    fit <- suppressWarnings(suppressMessages(nlmixr(one.cmt, theo_sd, est = "saem", control = saemControlFast)))
    .th <- c("tka", "tcl", "tv")
    .se0 <- fit$parFixedDf[.th, "SE"]
    expect_true(all(is.finite(.se0)))
    # install a covariance with twice the standard errors
    .cov <- diag((2 * .se0)^2)
    dimnames(.cov) <- list(.th, .th)
    .updateParFixedRefreshSeFromCov(fit$env, .cov)
    .pf <- fit$parFixedDf
    .e <- .pf[.th, "Estimate"]
    .s <- .pf[.th, "SE"]
    expect_equal(.s, 2 * .se0)
    qn <- qnorm(0.975)
    # backTransform("t100"): the interval is t100() of the new one
    expect_equal(.pf["tcl", "CI Lower"], t100(.e[[2]] - qn * .s[[2]]))
    expect_equal(.pf["tcl", "CI Upper"], t100(.e[[2]] + qn * .s[[2]]))
    expect_equal(
      fit$parFixed["tcl", "Back-transformed(95%CI)"],
      sprintf(
        "%s (%s, %s)",
        formatMinWidth(t100(.e[[2]])),
        formatMinWidth(t100(.e[[2]] - qn * .s[[2]])),
        formatMinWidth(t100(.e[[2]] + qn * .s[[2]]))
      )
    )
    # backTransform("none") names no function: the default exp() stays
    expect_equal(.pf["tv", "CI Lower"], exp(.e[[3]] - qn * .s[[3]]))
    expect_equal(.pf["tv", "CI Upper"], exp(.e[[3]] + qn * .s[[3]]))
    expect_equal(.pf["tka", "CI Upper"], exp(.e[[1]] + qn * .s[[1]]))
  })
})
