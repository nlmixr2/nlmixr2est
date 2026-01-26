nmTest({

  tol <- 1e-5
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
      expect_true(all(exp.names %fin% colnames(o)))
    }
  }
  tmpDir <- tempdir()
  if (!dir.exists(tmpDir)) dir.create(tmpDir)
  options(
    nlmixr.save = TRUE,
    nlmixr.save.dir = tmpDir
  )

  # Use centralized model from helper-models.R
  fitS <- .nlmixr(one.compartment.add.err, theo_sd, est = "saem", control = saemControlFast)

  test_that("tidy works on nlmixr fit SAEM fits", {

    td <- broom.mixed::tidy(fitS, exponentiate = NA)
    check_tidy(td, 7, 7, c("effect", "group", "term", "estimate", "std.error", "statistic", "p.value"))
    expect_equal(
      td$term,
      c(
        "tka", "tcl", "tv", "sd__eta.ka", "sd__eta.cl", "sd__eta.v",
        "add.err"
      )
    )
    td <- broom.mixed::tidy(fitS, conf.level = 0.9, exponentiate = NA)
    check_tidy(
      td, 7, 9,
      c(
        "effect", "group", "term", "estimate", "std.error", "statistic", "p.value",
        "conf.low", "conf.high"
      )
    )
    expect_equal(
      td$term,
      c(
        "tka", "tcl", "tv", "sd__eta.ka", "sd__eta.cl", "sd__eta.v",
        "add.err"
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
    check_tidy(td, 4, 4, c("effect", "group", "term", "estimate"))

    td <- broom.mixed::tidy(fitS, effects = "ran_pars", exponentiate = FALSE)
    td2 <- td$estimate
    check_tidy(td, 4, 4, c("effect", "group", "term", "estimate"))

    td <- broom.mixed::tidy(fitS, effects = "ran_pars", exponentiate = TRUE)
    td3 <- td$estimate
    check_tidy(td, 4, 4, c("effect", "group", "term", "estimate"))

    expect_equal(td1, td2)
    expect_equal(td2, td3)

  })

  skip_on_cran()

  for (f in c("focei", "foce")) {
    fitF <- .nlmixr(one.compartment, theo_sd, est = f, control = foceiControlFast)

    test_that(sprintf("tidy works on nlmixr fit %s fits", f), {

      td <- broom.mixed::tidy(fitF, exponentiate = NA)
      check_tidy(td, 7, 7, c("effect", "group", "term", "estimate", "std.error", "statistic", "p.value"))
      expect_equal(
        td$term,
        c(
          "tka", "tcl", "tv", "sd__eta.ka", "sd__eta.cl", "sd__eta.v",
          "add.err"
        )
      )
      td <- broom.mixed::tidy(fitF, conf.level = 0.9, exponentiate = NA)
      check_tidy(td, 7, 9, c(
        "effect", "group", "term", "estimate", "std.error", "statistic", "p.value",
        "conf.low", "conf.high"
      ))
      expect_equal(
        td$term,
        c(
          "tka", "tcl", "tv", "sd__eta.ka", "sd__eta.cl", "sd__eta.v",
          "add.err"
        )
      )
      .est <- td$estimate
      .stdErr <- td$std.error
      .confLow <- td$conf.low

      td <- broom.mixed::tidy(fitF, conf.level = 0.9, exponentiate = FALSE)
      check_tidy(td)
      expect_equal(td$estimate[1:3], log(.est)[1:3],
                   tolerance = tol
                   )
      expect_equal(td$estimate[-(1:3)], .est[-(1:3)],
                   tolerance = tol
                   )
      expect_equal(setNames(exp(td$estimate[1:3]) * td$std.error[1:3], NULL),
                   .stdErr[1:3],
                   tolerance = tol
                   )
      expect_equal(exp(td$conf.low), .confLow,
                   tolerance = tol
                   )

      for (ef in c("ran_vals", "random")) {
        td <- broom.mixed::tidy(fitF, effects = ef, exponentiate = NA)
        td1 <- td$estimate
        check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))
        ##
        td <- broom.mixed::tidy(fitF, effects = ef, exponentiate = FALSE)
        td2 <- td$estimate
        check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))
        ##
        td <- broom.mixed::tidy(fitF, effects = ef, exponentiate = TRUE)
        td3 <- td$estimate
        check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))
        ##
        expect_equal(td1, td2)
        expect_equal(td2, td3)
      }

      ##
      td <- broom.mixed::tidy(fitF, effects = "ran_coef", exponentiate = NA)
      td1 <- td$estimate
      check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))
      ##
      td <- broom.mixed::tidy(fitF, effects = "ran_coef", exponentiate = FALSE)
      td2 <- td$estimate
      check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))
      ##
      td <- broom.mixed::tidy(fitF, effects = "ran_coef", exponentiate = TRUE)
      td3 <- td$estimate
      check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))
      ##
      expect_equal(log(td1), td2)
      expect_equal(td2, log(td3))
      ##
      td <- broom.mixed::tidy(fitF, effects = "ran_pars", exponentiate = NA)
      td1 <- td$estimate
      check_tidy(td, 4, 4, c("effect", "group", "term", "estimate"))
      ##
      td <- broom.mixed::tidy(fitF, effects = "ran_pars", exponentiate = FALSE)
      td2 <- td$estimate
      check_tidy(td, 4, 4, c("effect", "group", "term", "estimate"))
      ##
      td <- broom.mixed::tidy(fitF, effects = "ran_pars", exponentiate = TRUE)
      td3 <- td$estimate
      check_tidy(td, 4, 4, c("effect", "group", "term", "estimate"))
      ##
      expect_equal(td1, td2, tolerance = tol)
      expect_equal(td2, td3, tolerance = tol)
    })
  }

  for (f in c("foi", "fo")) {
    fitF <- suppressMessages(suppressWarnings(nlmixr(one.compartment, theo_sd, est = f, control=list(print=0))))
    test_that(sprintf("tidy works on nlmixr fit %s fits", f), {
      td <- broom.mixed::tidy(fitF, exponentiate = NA)
      check_tidy(td, 7, 7, c("effect", "group", "term", "estimate", "std.error", "statistic", "p.value"))
      expect_equal(
        td$term,
        c(
          "tka", "tcl", "tv", "sd__eta.ka", "sd__eta.cl", "sd__eta.v",
          "add.err"
        )
      )
      td <- broom.mixed::tidy(fitF, conf.level = 0.9, exponentiate = NA)
      check_tidy(td, 7, 9, c(
        "effect", "group", "term", "estimate", "std.error", "statistic", "p.value",
        "conf.low", "conf.high"
      ))
      expect_equal(
        td$term,
        c(
          "tka", "tcl", "tv", "sd__eta.ka", "sd__eta.cl", "sd__eta.v",
          "add.err"
        )
      )
      .est <- td$estimate
      .stdErr <- td$std.error
      .confLow <- td$conf.low

      td <- broom.mixed::tidy(fitF, conf.level = 0.9, exponentiate = FALSE)
      check_tidy(td)
      expect_equal(td$estimate[1:3], log(.est)[1:3],
                   tolerance = tol
                   )
      expect_equal(td$estimate[-(1:3)], .est[-(1:3)],
                   tolerance = tol
                   )
      ## exp(.df$model.est[.exp])*.df$std.error[.exp]
      expect_equal(setNames(exp(td$estimate[1:3]) * td$std.error[1:3], NULL),
                   .stdErr[1:3],
                   tolerance = tol
                   )
      expect_equal(exp(td$conf.low), .confLow,
                   tolerance = tol
                   )
      for (ef in c("ran_vals", "random")) {
        td <- broom.mixed::tidy(fitF, effects = ef, exponentiate = NA)
        td1 <- td$estimate
        check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))
        ##
        td <- broom.mixed::tidy(fitF, effects = ef, exponentiate = FALSE)
        td2 <- td$estimate
        check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))
        ##
        td <- broom.mixed::tidy(fitF, effects = ef, exponentiate = TRUE)
        td3 <- td$estimate
        check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))
        ##
        expect_equal(td1, td2, tolerance = tol)
        expect_equal(td2, td3, tolerance = tol)
      }
      ##
      td <- broom.mixed::tidy(fitF, effects = "ran_coef", exponentiate = NA)
      td1 <- td$estimate
      check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))
      ##
      td <- broom.mixed::tidy(fitF, effects = "ran_coef", exponentiate = FALSE)
      td2 <- td$estimate
      check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))
      ##
      td <- broom.mixed::tidy(fitF, effects = "ran_coef", exponentiate = TRUE)
      td3 <- td$estimate
      check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))
      ##
      expect_equal(log(td1), td2, tolerance = tol)
      expect_equal(td2, log(td3), tolerance = tol)
      ##
      td <- broom.mixed::tidy(fitF, effects = "ran_pars", exponentiate = NA)
      td1 <- td$estimate
      check_tidy(td, 4, 4, c("effect", "group", "term", "estimate"))
      ##
      td <- broom.mixed::tidy(fitF, effects = "ran_pars", exponentiate = FALSE)
      td2 <- td$estimate
      check_tidy(td, 4, 4, c("effect", "group", "term", "estimate"))
      ##
      td <- broom.mixed::tidy(fitF, effects = "ran_pars", exponentiate = TRUE)
      td3 <- td$estimate
      check_tidy(td, 4, 4, c("effect", "group", "term", "estimate"))
      ##
      expect_equal(td1, td2, tolerance = tol)
      expect_equal(td2, td3, tolerance = tol)
    })
  }

  fitP <- suppressMessages(suppressWarnings(nlmixr(one.compartment, theo_sd, est = "posthoc",
                                                   control=posthocControl(covMethod=0,
                                                                          calcTables=FALSE))))

  test_that("tidy works on posthoc fit fits", {

    td <- broom.mixed::tidy(fitP, exponentiate = NA)
    check_tidy(td, 7, 7, c("effect", "group", "term", "estimate", "std.error", "statistic", "p.value"))
    expect_equal(td$term, c(
      "tka", "tcl", "tv", "sd__eta.ka", "sd__eta.cl", "sd__eta.v",
      "add.err"
    ))
    td <- broom.mixed::tidy(fitP, conf.level = 0.9, exponentiate = NA)
    check_tidy(td, 7, 9, c(
      "effect", "group", "term", "estimate", "std.error", "statistic", "p.value",
      "conf.low", "conf.high"
    ))
    expect_equal(td$term, c(
      "tka", "tcl", "tv", "sd__eta.ka", "sd__eta.cl", "sd__eta.v",
      "add.err"
    ))
    expect_equal(td$estimate, c(
      1.56831218549017, 2.71828182845905, 31.5003923087479, 0.774596669241483,
      0.547722557505166, 0.316227766016838, 0.7
    ),
    tolerance = tol
    )
    expect_equal(td$std.error, c(
      NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_,
      NA_real_
    ))
    expect_equal(td$conf.low, c(
      NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_,
      NA_real_
    ))

    td <- broom.mixed::tidy(fitP, conf.level = 0.9, exponentiate = FALSE)
    check_tidy(td)
    expect_equal(td$estimate, c(
      0.45, 1, 3.45,
      0.774596669241483, 0.547722557505166, 0.316227766016838, 0.7
    ), tolerance = tol)

    expect_equal(setNames(td$std.error, NULL),
                 c(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_),
                 tolerance = tol
                 )

    expect_equal(td$conf.low, c(
      NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_,
      NA_real_
    ), tolerance = tol)

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

    expect_equal(td1, c(
      1.75611262206464, 1.92886734474446, 2.36872361414073, 1.18776870402381,
      1.48421770349779, 1.1603466406079, 0.728548560764397, 1.37496750947703,
      6.63504345280614, 0.728113370661147, 3.56192283422217, 0.87667797314966,
      1.62338691901516, 3.22874000491057, 2.81088025117325, 2.7068918775575,
      2.37701456492688, 4.04231470515951, 3.23451222726721, 3.26106788182658,
      2.88069603771061, 1.86803407932834, 3.73237924876953, 2.49388025078765,
      29.2305913871308, 31.8379070004274, 33.9077584685152, 31.2454399174589,
      27.0633487903421, 40.6859508443866, 33.6487106457298, 35.5118803877647,
      31.940988239621, 26.0728299589223, 37.2139995031868, 24.7515885649664
    ), tolerance = tol)

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
    expect_equal(td1, c(0.774596669241483, 0.547722557505166, 0.316227766016838, 0.7),
                 tolerance = tol
                 )

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
