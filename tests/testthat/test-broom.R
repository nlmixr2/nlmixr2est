nmTest({
  tol <- 1e-5
  ## Per-subject ran_coef are exp(theta + eta_i) at the converged EBEs, so they
  ## inherit the ODE solve accuracy.  posthocControl() defaults to sigdig=3
  ## (rtol=1e-3), so pinning these golden values to 1e-5 asserts far more
  ## precision than the solve is asked to deliver -- they drift ~2e-4 run to run
  ## and tighten to ~3e-6 only at sigdig=8.  Compare fit-to-fit results at `tol`;
  ## compare frozen golden numbers at the accuracy the default solve provides.
  tolEbe <- 1e-3
  ## From https://raw.githubusercontent.com/bbolker/broom.mixed/master/tests/testthat/helper-checkers.R

  ##' test the basics of tidy/augment/glance output: is a data frame, no row names
  check_tidiness <- function(o) {
    expect_s3_class(o, "tbl_df")
    expect_equal(rownames(o), as.character(seq_len(nrow(o))))
  }

  #' check the output of a tidy function
  check_tidy <- function(o, exp.row = NULL, exp.col = NULL, exp.names = NULL) {
    check_tidiness(o)
    if (!is.null(exp.row)) {
      expect_equal(nrow(o), exp.row)
    }
    if (!is.null(exp.col)) {
      expect_equal(ncol(o), exp.col)
    }
    if (!is.null(exp.names)) {
      expect_true(all(exp.names %in% colnames(o)))
    }
  }
  tmpDir <- tempdir()
  if (!dir.exists(tmpDir)) {
    dir.create(tmpDir)
  }
  options(
    nlmixr.save = TRUE,
    nlmixr.save.dir = tmpDir
  )

  # Use centralized fit from helper-fits.R
  fitS <- one.compartment.fit.saem

  test_that("tidy works on nlmixr fit SAEM fits", {
    td <- broom.mixed::tidy(fitS, exponentiate = NA)
    check_tidy(td, 7, 7, c("effect", "group", "term", "estimate", "std.error", "statistic", "p.value"))
    expect_equal(
      td$term,
      c(
        "tka",
        "tcl",
        "tv",
        "sd__eta.ka",
        "sd__eta.cl",
        "sd__eta.v",
        "add.sd"
      )
    )
    td <- broom.mixed::tidy(fitS, conf.level = 0.9, exponentiate = NA)
    check_tidy(
      td,
      7,
      9,
      c(
        "effect",
        "group",
        "term",
        "estimate",
        "std.error",
        "statistic",
        "p.value",
        "conf.low",
        "conf.high"
      )
    )
    expect_equal(
      td$term,
      c(
        "tka",
        "tcl",
        "tv",
        "sd__eta.ka",
        "sd__eta.cl",
        "sd__eta.v",
        "add.sd"
      )
    )
    .est <- td$estimate
    .stdErr <- td$std.error
    .confLow <- td$conf.low

    td <- broom.mixed::tidy(fitS, conf.level = 0.9, exponentiate = FALSE)
    check_tidy(td)
    expect_equal(
      td$estimate[1:3],
      log(.est)[1:3],
      tolerance = tol
    )
    expect_equal(
      td$estimate[-(1:3)],
      .est[-(1:3)],
      tolerance = tol
    )
    ## exp(.df$model.est[.exp])*.df$std.error[.exp]
    expect_equal(
      setNames(exp(td$estimate[1:3]) * td$std.error[1:3], NULL),
      .stdErr[1:3],
      tolerance = tol
    )
    expect_equal(exp(td$conf.low), .confLow, tolerance = tol)

    for (ef in c("ran_vals", "random")) {
      td <- broom.mixed::tidy(fitS, effects = ef, exponentiate = NA)
      td1 <- td$estimate
      check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

      td <- broom.mixed::tidy(fitS, effects = ef, exponentiate = FALSE)
      td2 <- td$estimate
      check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

      td <- broom.mixed::tidy(fitS, effects = ef, exponentiate = TRUE)
      td3 <- td$estimate
      check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

      expect_equal(td1, td2)
      expect_equal(td2, td3)
    }

    td <- broom.mixed::tidy(fitS, effects = "ran_coef", exponentiate = NA)
    td1 <- td$estimate
    check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

    .td1 <- td1

    td <- broom.mixed::tidy(fitS, effects = "ran_coef", exponentiate = FALSE)
    td2 <- td$estimate
    check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

    td <- broom.mixed::tidy(fitS, effects = "ran_coef", exponentiate = TRUE)
    td3 <- td$estimate
    check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

    expect_equal(log(td1), td2)
    expect_equal(td2, log(td3))

    td <- broom.mixed::tidy(fitS, effects = "ran_pars", exponentiate = NA)
    td1 <- td$estimate
    check_tidy(td, 4, 7, c("effect", "group", "term", "estimate", "std.error", "statistic", "p.value"))

    td <- broom.mixed::tidy(fitS, effects = "ran_pars", exponentiate = FALSE)
    td2 <- td$estimate
    check_tidy(td, 4, 7, c("effect", "group", "term", "estimate", "std.error", "statistic", "p.value"))

    td <- broom.mixed::tidy(fitS, effects = "ran_pars", exponentiate = TRUE)
    td3 <- td$estimate
    check_tidy(td, 4, 7, c("effect", "group", "term", "estimate", "std.error", "statistic", "p.value"))

    expect_equal(td1, td2)
    expect_equal(td2, td3)
  })

  skip_on_cran()

  for (f in c("focei", "foce")) {
    fitF <- switch(
      f,
      focei = one.compartment.fit.focei.fast,
      foce = one.compartment.fit.foce
    )

    test_that(sprintf("tidy works on nlmixr fit %s fits", f), {
      td <- broom.mixed::tidy(fitF, exponentiate = NA)
      check_tidy(td, 7, 7, c("effect", "group", "term", "estimate", "std.error", "statistic", "p.value"))
      expect_equal(
        td$term,
        c(
          "tka",
          "tcl",
          "tv",
          "sd__eta.ka",
          "sd__eta.cl",
          "sd__eta.v",
          "add.sd"
        )
      )
      td <- broom.mixed::tidy(fitF, conf.level = 0.9, exponentiate = NA)
      check_tidy(
        td,
        7,
        9,
        c(
          "effect",
          "group",
          "term",
          "estimate",
          "std.error",
          "statistic",
          "p.value",
          "conf.low",
          "conf.high"
        )
      )
      expect_equal(
        td$term,
        c(
          "tka",
          "tcl",
          "tv",
          "sd__eta.ka",
          "sd__eta.cl",
          "sd__eta.v",
          "add.sd"
        )
      )
      .est <- td$estimate
      .stdErr <- td$std.error
      .confLow <- td$conf.low

      td <- broom.mixed::tidy(fitF, conf.level = 0.9, exponentiate = FALSE)
      check_tidy(td)
      expect_equal(td$estimate[1:3], log(.est)[1:3], tolerance = tol)
      expect_equal(td$estimate[-(1:3)], .est[-(1:3)], tolerance = tol)
      expect_equal(setNames(exp(td$estimate[1:3]) * td$std.error[1:3], NULL), .stdErr[1:3], tolerance = tol)
      expect_equal(exp(td$conf.low), .confLow, tolerance = tol)

      for (ef in c("ran_vals", "random")) {
        td <- broom.mixed::tidy(fitF, effects = ef, exponentiate = NA)
        td1 <- td$estimate
        check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

        td <- broom.mixed::tidy(fitF, effects = ef, exponentiate = FALSE)
        td2 <- td$estimate
        check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

        td <- broom.mixed::tidy(fitF, effects = ef, exponentiate = TRUE)
        td3 <- td$estimate
        check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

        expect_equal(td1, td2)
        expect_equal(td2, td3)
      }

      td <- broom.mixed::tidy(fitF, effects = "ran_coef", exponentiate = NA)
      td1 <- td$estimate
      check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

      td <- broom.mixed::tidy(fitF, effects = "ran_coef", exponentiate = FALSE)
      td2 <- td$estimate
      check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

      td <- broom.mixed::tidy(fitF, effects = "ran_coef", exponentiate = TRUE)
      td3 <- td$estimate
      check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

      expect_equal(log(td1), td2)
      expect_equal(td2, log(td3))

      td <- broom.mixed::tidy(fitF, effects = "ran_pars", exponentiate = NA)
      td1 <- td$estimate
      check_tidy(td, 4, 7, c("effect", "group", "term", "estimate", "std.error", "statistic", "p.value"))

      td <- broom.mixed::tidy(fitF, effects = "ran_pars", exponentiate = FALSE)
      td2 <- td$estimate
      check_tidy(td, 4, 7, c("effect", "group", "term", "estimate", "std.error", "statistic", "p.value"))

      td <- broom.mixed::tidy(fitF, effects = "ran_pars", exponentiate = TRUE)
      td3 <- td$estimate
      check_tidy(td, 4, 7, c("effect", "group", "term", "estimate", "std.error", "statistic", "p.value"))
      expect_equal(td1, td2, tolerance = tol)
      expect_equal(td2, td3, tolerance = tol)
    })
  }

  for (f in c("foi", "fo")) {
    fitF <- switch(
      f,
      foi = one.compartment.fit.foi,
      fo = one.compartment.fit.fo
    )
    test_that(sprintf("tidy works on nlmixr fit %s fits", f), {
      td <- broom.mixed::tidy(fitF, exponentiate = NA)
      check_tidy(td, 7, 7, c("effect", "group", "term", "estimate", "std.error", "statistic", "p.value"))
      expect_equal(
        td$term,
        c(
          "tka",
          "tcl",
          "tv",
          "sd__eta.ka",
          "sd__eta.cl",
          "sd__eta.v",
          "add.sd"
        )
      )
      td <- broom.mixed::tidy(fitF, conf.level = 0.9, exponentiate = NA)
      check_tidy(
        td,
        7,
        9,
        c(
          "effect",
          "group",
          "term",
          "estimate",
          "std.error",
          "statistic",
          "p.value",
          "conf.low",
          "conf.high"
        )
      )
      expect_equal(
        td$term,
        c(
          "tka",
          "tcl",
          "tv",
          "sd__eta.ka",
          "sd__eta.cl",
          "sd__eta.v",
          "add.sd"
        )
      )
      .est <- td$estimate
      .stdErr <- td$std.error
      .confLow <- td$conf.low

      td <- broom.mixed::tidy(fitF, conf.level = 0.9, exponentiate = FALSE)
      check_tidy(td)
      expect_equal(td$estimate[1:3], log(.est)[1:3], tolerance = tol)
      expect_equal(td$estimate[-(1:3)], .est[-(1:3)], tolerance = tol)
      ## exp(.df$model.est[.exp])*.df$std.error[.exp]
      expect_equal(setNames(exp(td$estimate[1:3]) * td$std.error[1:3], NULL), .stdErr[1:3], tolerance = tol)
      expect_equal(exp(td$conf.low), .confLow, tolerance = tol)
      for (ef in c("ran_vals", "random")) {
        td <- broom.mixed::tidy(fitF, effects = ef, exponentiate = NA)
        td1 <- td$estimate
        check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

        td <- broom.mixed::tidy(fitF, effects = ef, exponentiate = FALSE)
        td2 <- td$estimate
        check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

        td <- broom.mixed::tidy(fitF, effects = ef, exponentiate = TRUE)
        td3 <- td$estimate
        check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

        expect_equal(td1, td2, tolerance = tol)
        expect_equal(td2, td3, tolerance = tol)
      }

      td <- broom.mixed::tidy(fitF, effects = "ran_coef", exponentiate = NA)
      td1 <- td$estimate
      check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

      td <- broom.mixed::tidy(fitF, effects = "ran_coef", exponentiate = FALSE)
      td2 <- td$estimate
      check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

      td <- broom.mixed::tidy(fitF, effects = "ran_coef", exponentiate = TRUE)
      td3 <- td$estimate
      check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

      expect_equal(log(td1), td2, tolerance = tol)
      expect_equal(td2, log(td3), tolerance = tol)

      td <- broom.mixed::tidy(fitF, effects = "ran_pars", exponentiate = NA)
      td1 <- td$estimate
      check_tidy(td, 4, 7, c("effect", "group", "term", "estimate", "std.error", "statistic", "p.value"))

      td <- broom.mixed::tidy(fitF, effects = "ran_pars", exponentiate = FALSE)
      td2 <- td$estimate
      check_tidy(td, 4, 7, c("effect", "group", "term", "estimate", "std.error", "statistic", "p.value"))

      td <- broom.mixed::tidy(fitF, effects = "ran_pars", exponentiate = TRUE)
      td3 <- td$estimate
      check_tidy(td, 4, 7, c("effect", "group", "term", "estimate", "std.error", "statistic", "p.value"))
      expect_equal(td1, td2, tolerance = tol)
      expect_equal(td2, td3, tolerance = tol)
    })
  }

  fitP <- one.compartment.fit.posthoc

  test_that("tidy works on posthoc fit fits", {
    td <- broom.mixed::tidy(fitP, exponentiate = NA)
    check_tidy(td, 7, 7, c("effect", "group", "term", "estimate", "std.error", "statistic", "p.value"))
    expect_equal(
      td$term,
      c(
        "tka",
        "tcl",
        "tv",
        "sd__eta.ka",
        "sd__eta.cl",
        "sd__eta.v",
        "add.sd"
      )
    )
    td <- broom.mixed::tidy(fitP, conf.level = 0.9, exponentiate = NA)
    check_tidy(
      td,
      7,
      9,
      c(
        "effect",
        "group",
        "term",
        "estimate",
        "std.error",
        "statistic",
        "p.value",
        "conf.low",
        "conf.high"
      )
    )
    expect_equal(
      td$term,
      c(
        "tka",
        "tcl",
        "tv",
        "sd__eta.ka",
        "sd__eta.cl",
        "sd__eta.v",
        "add.sd"
      )
    )
    expect_equal(
      td$estimate,
      c(
        1.56831218549017,
        2.71828182845905,
        31.5003923087479,
        0.774596669241483,
        0.547722557505166,
        0.316227766016838,
        0.7
      ),
      tolerance = tol
    )
    expect_equal(
      td$std.error,
      c(
        NA_real_,
        NA_real_,
        NA_real_,
        NA_real_,
        NA_real_,
        NA_real_,
        NA_real_
      )
    )
    expect_equal(
      td$conf.low,
      c(
        NA_real_,
        NA_real_,
        NA_real_,
        NA_real_,
        NA_real_,
        NA_real_,
        NA_real_
      )
    )

    td <- broom.mixed::tidy(fitP, conf.level = 0.9, exponentiate = FALSE)
    check_tidy(td)
    expect_equal(
      td$estimate,
      c(
        0.45,
        1,
        3.45,
        0.774596669241483,
        0.547722557505166,
        0.316227766016838,
        0.7
      ),
      tolerance = tol
    )

    expect_equal(
      setNames(td$std.error, NULL),
      c(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_),
      tolerance = tol
    )

    expect_equal(
      td$conf.low,
      c(
        NA_real_,
        NA_real_,
        NA_real_,
        NA_real_,
        NA_real_,
        NA_real_,
        NA_real_
      ),
      tolerance = tol
    )

    for (ef in c("ran_vals", "random")) {
      td <- broom.mixed::tidy(fitP, effects = ef, exponentiate = NA)
      td1 <- td$estimate
      check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

      td <- broom.mixed::tidy(fitP, effects = ef, exponentiate = FALSE)
      td2 <- td$estimate
      check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

      td <- broom.mixed::tidy(fitP, effects = ef, exponentiate = TRUE)
      td3 <- td$estimate
      check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

      expect_equal(td1, td2, tolerance = tol)
      expect_equal(td2, td3, tolerance = tol)
    }

    td <- broom.mixed::tidy(fitP, effects = "ran_coef", exponentiate = NA)
    td1 <- td$estimate
    check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

    expect_equal(
      td1,
      c(
        1.75611262206464,
        1.92886734474446,
        2.36872361414073,
        1.18776870402381,
        1.48421770349779,
        1.1603466406079,
        0.728548560764397,
        1.37496750947703,
        6.63504345280614,
        0.728113370661147,
        3.56192283422217,
        0.87667797314966,
        1.62338691901516,
        3.22874000491057,
        2.81088025117325,
        2.7068918775575,
        2.37701456492688,
        4.04231470515951,
        3.23451222726721,
        3.26106788182658,
        2.88069603771061,
        1.86803407932834,
        3.73237924876953,
        2.49388025078765,
        29.2305913871308,
        31.8379070004274,
        33.9077584685152,
        31.2454399174589,
        27.0633487903421,
        40.6859508443866,
        33.6487106457298,
        35.5118803877647,
        31.940988239621,
        26.0728299589223,
        37.2139995031868,
        24.7515885649664
      ),
      tolerance = tolEbe
    )

    td <- broom.mixed::tidy(fitP, effects = "ran_coef", exponentiate = FALSE)
    td2 <- td$estimate
    check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

    td <- broom.mixed::tidy(fitP, effects = "ran_coef", exponentiate = TRUE)
    td3 <- td$estimate
    check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

    expect_equal(log(td1), td2, tolerance = tol)
    expect_equal(td2, log(td3), tolerance = tol)

    td <- broom.mixed::tidy(fitP, effects = "ran_pars", exponentiate = NA)
    td1 <- td$estimate
    check_tidy(td, 4, 4, c("effect", "group", "term", "estimate"))
    expect_equal(td1, c(0.774596669241483, 0.547722557505166, 0.316227766016838, 0.7), tolerance = tol)

    td <- broom.mixed::tidy(fitP, effects = "ran_pars", exponentiate = FALSE)
    td2 <- td$estimate
    check_tidy(td, 4, 4, c("effect", "group", "term", "estimate"))

    td <- broom.mixed::tidy(fitP, effects = "ran_pars", exponentiate = TRUE)
    td3 <- td$estimate
    check_tidy(td, 4, 4, c("effect", "group", "term", "estimate"))

    expect_equal(td1, td2, tolerance = tol)
    expect_equal(td2, td3, tolerance = tol)
  })
})

# eta.ka is not mu-referenced, so $parFixedDf carries a row for it after the
# thetas; tcl is fixed and two thetas are labelled
.broomEtaRowMod <- function() {
  ini({
    tka <- 0.45
    label("Ka")
    tcl <- fix(1)
    label("Cl")
    tv <- 3.45
    eta.ka ~ 0.6
    eta.cl ~ 0.3
    eta.v ~ 0.1
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka) * exp(eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    linCmt() ~ add(add.sd)
  })
}

nmTest({
  test_that("tidy() and confint() report the thetas only, not a non-mu-referenced eta's row", {
    skip_if_not_installed("broom.mixed")
    .ctl <- foceiControl(print = 0, maxOuterIterations = 0, calcTables = FALSE)
    fit <- .nlmixr(.broomEtaRowMod, theo_sd, est = "focei", control = .ctl)
    .pf <- fit$parFixedDf
    expect_identical(rownames(.pf), c("tka", "tcl", "tv", "add.sd", "eta.ka"))
    expect_identical(.pf$Parameter[1:4], c("Ka", "Cl", "", ""))
    .th <- c("tka", "tcl", "tv")
    .fx <- as.data.frame(broom.mixed::tidy(fit, effects = "fixed"))
    expect_identical(.fx$term, .th)
    expect_equal(.fx$estimate, unname(.pf[.th, "Estimate"]))
    expect_equal(.fx$std.error, unname(.pf[.th, "SE"]))
    .rp <- as.data.frame(broom.mixed::tidy(fit, effects = "ran_pars"))
    expect_identical(.rp$term, c("sd__eta.ka", "sd__eta.cl", "sd__eta.v", "add.sd"))
    expect_identical(.rp$group, c("ID", "ID", "ID", "Residual(add)"))
    expect_equal(.rp$estimate[4], unname(.pf["add.sd", "Estimate"]))
    .all <- c(.th, "add.sd")
    .z <- stats::qnorm(0.975)
    .lo <- unname(.pf[.all, "Estimate"] - .z * .pf[.all, "SE"])
    .hi <- unname(.pf[.all, "Estimate"] + .z * .pf[.all, "SE"])
    .ci <- confint(fit)
    expect_identical(rownames(.ci), .all)
    expect_equal(.ci[["2.5 %"]], .lo)
    expect_equal(.ci[["97.5 %"]], .hi)
    # exponentiate = NA back-transforms the exp() thetas; the eta row's NA used
    # to stop it ("NAs are not allowed in subscripted assignments")
    .isExp <- c(TRUE, TRUE, TRUE, FALSE)
    .ciE <- confint(fit, exponentiate = NA)
    expect_identical(rownames(.ciE), .all)
    expect_equal(.ciE[["2.5 %"]], ifelse(.isExp, exp(.lo), .lo))
    expect_equal(.ciE[["97.5 %"]], ifelse(.isExp, exp(.hi), .hi))
    .fxE <- as.data.frame(broom.mixed::tidy(fit, effects = "fixed", exponentiate = NA))
    expect_identical(.fxE$term, .th)
    expect_equal(.fxE$std.error, unname(exp(.pf[.th, "Estimate"]) * .pf[.th, "SE"]))
    .fxC <- as.data.frame(broom.mixed::tidy(fit, effects = "fixed", conf.int = TRUE))
    expect_identical(.fxC$term, .th)
    expect_equal(.fxC$conf.low, .lo[1:3])
    expect_equal(.fxC$conf.high, .hi[1:3])
  })

  test_that("tidy() and confint() leave out an IOV fit's occasion variability row", {
    skip_if_not_installed("broom.mixed")
    fit <- .iovFit("focei", control = list(maxOuterIterations = 0, calcTables = FALSE))
    expect_identical(rownames(fit$parFixedDf), c("tka", "tcl", "tv", "add.sd", "iov.cl"))
    .fx <- as.data.frame(broom.mixed::tidy(fit, effects = "fixed"))
    expect_identical(.fx$term, c("tka", "tcl", "tv"))
    expect_identical(rownames(confint(fit)), c("tka", "tcl", "tv", "add.sd"))
  })
})
