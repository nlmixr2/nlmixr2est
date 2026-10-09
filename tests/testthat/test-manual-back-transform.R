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

  test_that("a refreshed covariance refreshes the CI of a theta with an extra mu term (issue 1140)", {
    # tv is mu-referenced inside exp() with an extra constant: the parameter
    # table back-transforms it with exp(), the literal-fix rule leaves it alone
    one.cmt <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 1.45
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v + 2)
        linCmt() ~ add(add.sd)
      })
    }
    fit <- suppressWarnings(suppressMessages(nlmixr(
      one.cmt,
      theo_sd,
      est = "focei",
      control = foceiControl(print = 0, maxOuterIterations = 0, covMethod = "r", covFull = FALSE, calcTables = FALSE)
    )))
    .th <- c("tka", "tcl", "tv")
    .se0 <- fit$parFixedDf[.th, "SE"]
    expect_true(all(is.finite(.se0)))
    .i <- match("tv", rownames(fit$parFixedDf))
    expect_equal(unname(fit$parFixedDf[["Back-transformed"]][.i]), exp(unname(fit$parFixedDf[["Estimate"]][.i])))
    .cov <- diag((2 * .se0)^2)
    dimnames(.cov) <- list(.th, .th)
    expect_no_warning(.updateParFixedRefreshSeFromCov(fit$env, .cov))
    .pf <- fit$parFixedDf
    qn <- qnorm(0.975)
    .e <- unname(.pf[["Estimate"]][.i])
    .s <- unname(.pf[["SE"]][.i])
    expect_equal(.s, 2 * unname(.se0[[3]]))
    expect_equal(unname(.pf[["CI Lower"]][.i]), exp(.e - qn * .s))
    expect_equal(unname(.pf[["CI Upper"]][.i]), exp(.e + qn * .s))
  })

  test_that("a refreshed covariance back-transforms each CI bound the way the table does", {
    # the table applies a backTransform() function one value at a time, so the
    # refresh does too: a function that only takes a scalar keeps its CI.  A
    # row whose back-transform cannot be reproduced loses its CI with a warning
    # rather than keeping the previous covariance's interval.
    assign("btScalar", function(x) if (x > 0) exp(x) else x, envir = globalenv())
    withr::defer(rm("btScalar", envir = globalenv()))
    btMod <- function() {
      ini({
        tka <- 0.45
        backTransform("btScalar")
        tcl <- 1
        tq <- 0.2
        backTransform("btMissing")
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl)
        cp <- ka * cl * tq
        cp ~ add(add.sd)
      })
    }
    .mkEnv <- function() {
      .e <- new.env(parent = emptyenv())
      .e$ui <- rxode2::rxode2(btMod)
      .est <- c(tka = 0.45, tcl = 1, tq = 0.2, add.sd = 0.7)
      .e$parFixedDf <- data.frame(
        Estimate = .est,
        SE = rep(0.1, 4),
        `%RSE` = 10,
        `Back-transformed` = c(exp(0.45), exp(1), 99, 0.7),
        `CI Lower` = c(1, 2, 3, 0.5),
        `CI Upper` = c(2, 3, 4, 0.9),
        check.names = FALSE,
        row.names = names(.est)
      )
      .e
    }
    .cov <- diag(c(0.2, 0.3, 0.4, 0.05)^2)
    dimnames(.cov) <- list(c("tka", "tcl", "tq", "add.sd"), c("tka", "tcl", "tq", "add.sd"))
    qn <- qnorm(0.975)
    .e <- .mkEnv()
    expect_warning(
      .updateParFixedRefreshSeFromCov(.e, .cov),
      "the confidence interval of 'tq' was dropped: its back-transform could not be reproduced",
      fixed = TRUE
    )
    .pf <- .e$parFixedDf
    expect_equal(unname(unlist(.pf["tka", c("CI Lower", "CI Upper")])), exp(0.45 + c(-1, 1) * qn * 0.2))
    expect_equal(unname(unlist(.pf["tcl", c("CI Lower", "CI Upper")])), exp(1 + c(-1, 1) * qn * 0.3))
    expect_equal(unname(unlist(.pf["tq", c("CI Lower", "CI Upper")])), c(NA_real_, NA_real_))
    expect_equal(unname(unlist(.pf["add.sd", c("CI Lower", "CI Upper")])), 0.7 + c(-1, 1) * qn * 0.05)
    expect_equal(unname(.pf[["SE"]]), c(0.2, 0.3, 0.4, 0.05))
    # ciIdentity (the variational covariance's rule): only a row reported
    # untransformed gets Estimate +/- z SE; the others keep their interval
    .e <- .mkEnv()
    expect_no_warning(.updateParFixedRefreshSeFromCov(.e, .cov, ciIdentity = TRUE))
    .pf <- .e$parFixedDf
    expect_equal(unname(unlist(.pf["add.sd", c("CI Lower", "CI Upper")])), 0.7 + c(-1, 1) * qn * 0.05)
    expect_equal(unname(unlist(.pf["tka", c("CI Lower", "CI Upper")])), c(1, 2))
    expect_equal(unname(unlist(.pf["tq", c("CI Lower", "CI Upper")])), c(3, 4))
  })
})
