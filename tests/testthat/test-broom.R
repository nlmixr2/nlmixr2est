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
    expect_equal(exp(td$conf.low[1:3]), .confLow[1:3], tolerance = tol)

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
      expect_equal(exp(td$conf.low[1:3]), .confLow[1:3], tolerance = tol)

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
      expect_equal(exp(td$conf.low[1:3]), .confLow[1:3], tolerance = tol)
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
    .sd <- c("sd__eta.ka", "sd__eta.cl", "sd__eta.v")
    # the variability terms follow the thetas, as lme4 lists them
    expect_identical(rownames(.ci), c(.all, .sd))
    expect_equal(.ci[.all, "2.5 %"], .lo)
    expect_equal(.ci[.all, "97.5 %"], .hi)
    .rpt <- .ranParTable(fit)
    expect_identical(.rpt$term, .sd)
    expect_equal(.ci[.sd, "estimate"], .rpt$estimate)
    expect_equal(.ci[.sd, "2.5 %"], .rpt$conf.low)
    expect_equal(.ci[.sd, "97.5 %"], .rpt$conf.high)
    expect_identical(rownames(confint(fit, parm = c("tv", "sd__eta.v"))), c("tv", "sd__eta.v"))
    # exponentiate = NA back-transforms the exp() thetas only; the eta.ka row
    # of $parFixedDf, with no estimate, is not among them
    .isExp <- c(TRUE, TRUE, TRUE, FALSE)
    .ciE <- confint(fit, exponentiate = NA)
    expect_identical(rownames(.ciE), c(.all, .sd))
    expect_equal(.ciE[.all, "2.5 %"], ifelse(.isExp, exp(.lo), .lo))
    expect_equal(.ciE[.all, "97.5 %"], ifelse(.isExp, exp(.hi), .hi))
    expect_equal(.ciE[.sd, "2.5 %"], .rpt$conf.low)
    .fxE <- as.data.frame(broom.mixed::tidy(fit, effects = "fixed", exponentiate = NA))
    expect_identical(.fxE$term, .th)
    expect_equal(.fxE$std.error, unname(exp(.pf[.th, "Estimate"]) * .pf[.th, "SE"]))
    .fxC <- as.data.frame(broom.mixed::tidy(fit, effects = "fixed", conf.int = TRUE))
    expect_identical(.fxC$term, .th)
    expect_equal(.fxC$conf.low, .lo[1:3])
    expect_equal(.fxC$conf.high, .hi[1:3])
    # the residual-error row keeps its interval under ran_pars
    .rpC <- as.data.frame(broom.mixed::tidy(fit, effects = "ran_pars", conf.int = TRUE))
    expect_equal(.rpC$conf.low[4], .lo[4])
    expect_equal(.rpC$conf.high[4], .hi[4])
    # broom.mixed's spelling, ran_coefs, gives the same coefficients
    expect_warning(
      .rc <- as.data.frame(broom.mixed::tidy(fit, effects = "ran_coefs")),
      "not mu-referenced"
    )
    expect_warning(
      .rc0 <- as.data.frame(broom.mixed::tidy(fit, effects = "ran_coef")),
      "not mu-referenced"
    )
    expect_identical(unique(.rc$effect), "ran_coefs")
    expect_equal(.rc$estimate, .rc0$estimate)
  })

  test_that("ran_pars SEs are the delta method over the Omega rows of $cov", {
    .om <- matrix(c(0.4, 0, 0, 0, 0.09, 0.02, 0, 0.02, 0.05), 3, 3, dimnames = list(c("a", "b", "c"), c("a", "b", "c")))
    .cn <- c("tka", "om.a", "om.b", "cov.c.b", "om.c")
    .cov <- matrix(0.001, 5, 5, dimnames = list(.cn, .cn))
    diag(.cov) <- c(0.01, 0.02, 0.003, 0.0005, 0.002)
    .omR <- .corFromCov(.om)
    .terms <- function(v) {
      .s <- sqrt(v[c(1, 2, 4)])
      c(.s, v[3] / (.s[2] * .s[3]))
    }
    .v0 <- c(0.4, 0.09, 0.02, 0.05)
    .h <- 1e-6
    .J <- vapply(
      1:4,
      function(k) {
        .d <- replace(numeric(4), k, .h)
        (.terms(.v0 + .d) - .terms(.v0 - .d)) / (2 * .h)
      },
      numeric(4)
    )
    .vc <- .cov[-1, -1]
    expect_identical(names(.getR(.omR, TRUE)), c("sd__a", "sd__b", "sd__c", "cor__c.b"))
    .pr <- cbind(row = 3L, col = 2L)
    expect_equal(.ranParSe(list(cov = .cov), .om, .pr), sqrt(diag(.J %*% .vc %*% t(.J))), tolerance = 1e-6)
    # a theta-only $cov has no Omega rows to give an SE
    expect_identical(.ranParSe(list(cov = .cov[1, 1, drop = FALSE]), .om, .pr), rep(NA_real_, 4))
  })

  test_that("Wishart ran_pars intervals: scaled chi-square SDs, Fisher-z correlations", {
    .w <- .ranParWishart(c(0.5, 0.3), c(TRUE, FALSE), 20, 0.9)
    expect_equal(.w$low[1], 0.5 * sqrt(20 / stats::qchisq(0.95, 20)))
    expect_equal(.w$high[1], 0.5 * sqrt(20 / stats::qchisq(0.05, 20)))
    expect_equal(.w$se, c(0.5 / sqrt(40), (1 - 0.09) / sqrt(20)))
    expect_equal(.w$low[2], tanh(atanh(0.3) - stats::qnorm(0.95) / sqrt(17)))
    expect_equal(.w$high[2], tanh(atanh(0.3) + stats::qnorm(0.95) / sqrt(17)))
    # three subjects give no Fisher-z interval
    expect_identical(.ranParWishart(0.3, FALSE, 3, 0.9)$low, NA_real_)
  })

  test_that("a fixed Omega element gives its terms no interval", {
    .fx <- function(ui) {
      .ranParFixed(ui$iniDf, ui$omega, .ranParPairs(ui$iniDf, dimnames(ui$omega)[[2]]))
    }
    .ui <- rxode2::rxode2(.bsvCorBlockFit)
    # sd__eta.ka, sd__eta.cl, sd__eta.v, cor__eta.v.eta.cl
    expect_identical(.fx(.ui), c(FALSE, FALSE, FALSE, FALSE))
    expect_identical(.fx(rxode2::ini(.ui, eta.ka ~ fix(0.6))), c(TRUE, FALSE, FALSE, FALSE))
    expect_identical(.fx(rxode2::ini(.ui, eta.cl + eta.v ~ fix(0.3, 0.05, 0.1))), c(FALSE, TRUE, TRUE, TRUE))
  })

  test_that("auto ran_pars intervals use $cov when it has Omega, else Wishart", {
    .ui <- rxode2::rxode2(.bsvCorBlockFit)
    .om <- .ui$omega
    .cn <- c("tka", "om.eta.ka", "om.eta.cl", "cov.eta.v.eta.cl", "om.eta.v")
    .cov <- diag(c(0.01, 0.02, 0.003, 0.0005, 0.002))
    dimnames(.cov) <- list(.cn, .cn)
    .x <- list(ui = .ui, omega = .om, omegaR = .corFromCov(.om), cov = .cov, eta = data.frame(ID = 1:12), est = "focei")
    expect_identical(.ranParTable(.x)$method, rep("cov", 4))
    expect_identical(.ranParTable(.x, method = "wishart")$method, rep("wishart", 4))
    .w <- .ranParWishart(.ranParTable(.x)$estimate, c(TRUE, TRUE, TRUE, FALSE), 12, 0.95)
    .x$cov <- .cov[1, 1, drop = FALSE]
    .t <- .ranParTable(.x)
    expect_identical(.t$method, rep("wishart", 4))
    expect_equal(.t$conf.low, .w$low)
    expect_identical(.ranParTable(.x, method = "cov")$conf.low, rep(NA_real_, 4))
    # a posthoc fit takes Omega from ini(); nothing about it was estimated
    .x$est <- "posthoc"
    expect_identical(.ranParTable(.x)$std.error, rep(NA_real_, 4))
  })

  test_that("a declared covariance estimated at zero keeps its cor__ term", {
    .ui <- rxode2::rxode2(.bsvCorBlockFit)
    .om <- .ui$omega
    .om["eta.cl", "eta.v"] <- .om["eta.v", "eta.cl"] <- 0
    .x <- list(ui = .ui, omega = .om, omegaR = .corFromCov(.om), cov = NULL, eta = data.frame(ID = 1:12), est = "focei")
    .t <- .ranParTable(.x)
    expect_identical(.t$term, c("sd__eta.ka", "sd__eta.cl", "sd__eta.v", "cor__eta.v.eta.cl"))
    expect_identical(.t$estimate[4], 0)
    # a ui rebuilt from that omega has lost the covariance row; the model as
    # estimated (iniDf0) still has it
    .rebuilt <- rxode2::rxode2(.ui)
    .rebuilt$iniDf <- .ui$iniDf[.ui$iniDf$name != "(eta.cl,eta.v)", , drop = FALSE]
    .env <- new.env(parent = emptyenv())
    .env$iniDf0 <- .ui$iniDf
    .x$ui <- .rebuilt
    expect_identical(.ranParTable(.x)$term[1:3], c("sd__eta.ka", "sd__eta.cl", "sd__eta.v"))
    expect_length(.ranParTable(.x)$term, 3L)
    .x$env <- .env
    expect_identical(.ranParTable(.x)$term[4], "cor__eta.v.eta.cl")
    # an undeclared pair (eta.ka with eta.cl) is not a term
    expect_identical(.ranParPairs(.ui$iniDf, c("eta.ka", "eta.cl", "eta.v")), cbind(row = 3L, col = 2L))
  })

  test_that("ran_pars intervals use log(sd) and atanh(cor), with an SD above 1", {
    .z <- stats::qnorm(0.975)
    # an SD above 1 must not reach atanh()
    .ci <- withr::with_options(
      list(warn = 2),
      .ranParCi(c(2, 1, 0.5), c(0.4, 0.1, 0.2), c(TRUE, TRUE, FALSE), 0.95)
    )
    expect_equal(.ci$low, c(2 * exp(-.z * 0.2), exp(-.z * 0.1), tanh(atanh(0.5) - .z * 0.2 / 0.75)))
    expect_equal(.ci$high, c(2 * exp(.z * 0.2), exp(.z * 0.1), tanh(atanh(0.5) + .z * 0.2 / 0.75)))
    expect_identical(.ranParCi(1, NA_real_, TRUE, 0.95)$low, NA_real_)
  })

  test_that("tidy() and confint() leave out an IOV fit's occasion variability row", {
    skip_if_not_installed("broom.mixed")
    fit <- .iovFit("focei", control = list(maxOuterIterations = 0, calcTables = FALSE))
    expect_identical(rownames(fit$parFixedDf), c("tka", "tcl", "tv", "add.sd", "iov.cl"))
    .fx <- as.data.frame(broom.mixed::tidy(fit, effects = "fixed"))
    expect_identical(.fx$term, c("tka", "tcl", "tv"))
    # only the subject-level variability is listed
    .rn <- rownames(confint(fit))
    expect_identical(.rn, c("tka", "tcl", "tv", "add.sd", .ranParTable(fit)$term))
    expect_false(any(grepl("iov", .rn)))
  })

  test_that("broom's theta rows come from the model's thetas, and none is an error", {
    .ui <- rxode2::rxode2(.broomEtaRowMod)
    .df <- data.frame(
      Estimate = c(0.45, 0.7, NA),
      SE = c(0.1, 0.05, NA),
      row.names = c("tka", "add.sd", "eta.ka")
    )
    .got <- .broomParFixedDf(list(parFixedDf = .df, ui = .ui))
    expect_identical(rownames(.got), c("tka", "add.sd"))
    expect_identical(names(.got), c("model.est", "std.error"))
    # a table whose rows name no theta of the model is not a fit's table
    expect_error(
      .broomParFixedDf(list(parFixedDf = .df["eta.ka", , drop = FALSE], ui = .ui)),
      "no row of $parFixedDf is a theta of the fit's model",
      fixed = TRUE
    )
    # an empty table has nothing to report
    expect_identical(nrow(.broomParFixedDf(list(parFixedDf = .df[0, , drop = FALSE], ui = .ui))), 0L)
  })
})
