# Model-generation side of M2/M3/M4 censoring for a llik-forced
# t()/cauchy()/dnorm() endpoint under the FOCEi family (#992).  The
# log-density in rx_pred_ hides the location/scale/degrees of freedom
# censEst.h needs, so .fixCensRNuLine() restores a real rx_r_/rx_nu_ and
# .rxFinalizeInner()/.rxFinalizePred() emit them (plus d(f)/d(eta)) as real
# lhs.  Cheap (no fits) -- the end-to-end numeric validation lives in
# test-focei-cens-t-fit.R.

nmTest({

  .withCensNuFix <- function(code) {
    .old <- nlmixr2global$rxCensNuFix
    .oldLlik <- nlmixr2global$rxPredLlik
    on.exit({
      nlmixr2global$rxCensNuFix <- .old
      nlmixr2global$rxPredLlik <- .oldLlik
    })
    nlmixr2global$rxCensNuFix <- TRUE
    nlmixr2global$rxPredLlik <- TRUE
    force(code)
  }

  .cauchyLines <- function() {
    list(quote(rx_yj_ ~ 122),
         quote(rx_lambda_ ~ 1),
         quote(rx_pred_f_ ~ cp),
         quote(rx_pred_ ~ rx_pred_f_),
         quote(rx_rll_ ~ sqrt((add.sd)^2)),
         quote(rx_pred_ ~ llikCauchy(DV, rx_pred_, rx_rll_)),
         quote(rx_r_ ~ 0))
  }

  .tLines <- function() {
    list(quote(rx_pred_f_ ~ cp),
         quote(rx_pred_ ~ rx_pred_f_),
         quote(rx_rll_ ~ sqrt((add.sd)^2)),
         quote(rx_pred_ ~ llikT(DV, nu, rx_pred_, rx_rll_)),
         quote(rx_r_ ~ 0))
  }

  test_that(".fixCensRNuLine only fires when the censoring flag is set", {
    .lines <- .cauchyLines()
    .old <- nlmixr2global$rxCensNuFix
    on.exit(nlmixr2global$rxCensNuFix <- .old)
    nlmixr2global$rxCensNuFix <- FALSE
    expect_identical(.fixCensRNuLine(.lines), .lines)
  })

  test_that(".fixCensRNuLine restores rx_r_ and adds rx_nu_ for cauchy()", {
    .out <- .withCensNuFix(.fixCensRNuLine(.cauchyLines()))
    # cauchy is Student-t with nu=1
    expect_true(any(vapply(.out, identical, logical(1), quote(rx_nu_ ~ 1))))
    expect_true(any(vapply(.out, identical, logical(1), quote(rx_r_ ~ rx_rll_^2))))
    # ...and the hardcoded zero is gone
    expect_false(any(vapply(.out, identical, logical(1), quote(rx_r_ ~ 0))))
  })

  test_that(".fixCensRNuLine picks up t()'s own degrees of freedom", {
    .out <- .withCensNuFix(.fixCensRNuLine(.tLines()))
    expect_true(any(vapply(.out, identical, logical(1), quote(rx_nu_ ~ nu))))
    expect_true(any(vapply(.out, identical, logical(1), quote(rx_r_ ~ rx_rll_^2))))
  })

  test_that(".fixCensRNuLine leaves an ar() endpoint alone", {
    # rx_rll_ is the MARGINAL sd for an AR(1) endpoint; the conditional scale
    # handed to llikCauchy() is rx_rll_*sqrt(1-phi^2), so squaring rx_rll_ back
    # into rx_r_ would feed the correction the wrong scale.
    .lines <- c(.cauchyLines()[1:5],
                list(quote(rx_arPhi_cp <- rx_arNf_cp * cor^rx_arDt_cp)),
                .cauchyLines()[6:7])
    expect_identical(.withCensNuFix(.fixCensRNuLine(.lines)), .lines)
  })

  test_that(".fixCensRNuLine leaves a distribution with no rx_rll_ alone", {
    .lines <- list(quote(rx_pred_f_ ~ cp), quote(rx_pred_ ~ rx_pred_f_),
                   quote(rx_r_ ~ 0))
    expect_identical(.withCensNuFix(.fixCensRNuLine(.lines)), .lines)
  })

  .cauchyUi <- function() {
    .f <- function() {
      ini({tka <- 0.45; tcl <- 1; tv <- 3.45; add.sd <- c(0, 0.7)
        eta.ka ~ 0.2; eta.cl ~ 0.1})
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd) + cauchy()
      })
    }
    rxode2::rxUiDecompress(rxode2::rxode2(.f))
  }

  .innerText <- function(censFix) {
    .old <- nlmixr2global$rxCensNuFix
    .oldLlik <- nlmixr2global$rxPredLlik
    on.exit({
      nlmixr2global$rxCensNuFix <- .old
      nlmixr2global$rxPredLlik <- .oldLlik
    })
    nlmixr2global$rxCensNuFix <- censFix
    nlmixr2global$rxPredLlik <- TRUE
    .s <- suppressMessages(rxUiGet.foceiEnv(list(.cauchyUi(), TRUE)))
    .s$..inner
  }

  test_that("the llik inner model exposes the censoring inputs (#992)", {
    .on <- suppressMessages(.innerText(TRUE))
    # the location and the degrees of freedom become real (read-back) lhs
    expect_match(.on, "rx_pred_f_=", fixed = TRUE)
    expect_match(.on, "rx_nu_=", fixed = TRUE)
    # ...and so does d(f)/d(eta), which the eta-gradient correction chains through
    expect_match(.on, "rx__sens_rx_pred_f__BY_ETA_1___=", fixed = TRUE)
    expect_match(.on, "rx__sens_rx_pred_f__BY_ETA_2___=", fixed = TRUE)
    # rx_r_ is no longer the hardcoded zero rxode2 emits for a llik endpoint
    expect_false(any(strsplit(.on, "\n")[[1]] == "rx_r_=0"))
    # they are APPENDED after the FOCEi eta block, which likInner0 reads
    # arithmetically from predOffset -- rx_r_'s last eta column must still
    # precede them
    expect_lt(regexpr("rx__sens_rx_r__BY_ETA_2___=", .on, fixed = TRUE),
              regexpr("rx_pred_f_=", .on, fixed = TRUE))
  })

  test_that("without the flag the llik inner model is unchanged (#992)", {
    .off <- suppressMessages(.innerText(FALSE))
    expect_true(any(strsplit(.off, "\n")[[1]] == "rx_r_=0"))
    expect_false(grepl("rx_nu_=", .off, fixed = TRUE))
    expect_false(grepl("rx__sens_rx_pred_f__BY_ETA_1___=", .off, fixed = TRUE))
  })

  test_that("censoring is no longer reported as ignored for focei t()/cauchy()", {
    .f <- function() {
      ini({tka <- 0.45; tcl <- 1; tv <- 3.45; add.sd <- c(0, 0.7); eta.ka ~ 0.2})
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd) + cauchy()
      })
    }
    .ui <- rxode2::rxUiDecompress(rxode2::rxode2(.f))
    .d <- nlmixr2data::theo_sd
    .d <- .d[.d$ID <= 2, ]
    .d$CENS <- ifelse(.d$EVID == 0 & .d$DV < 3, 1L, 0L)
    for (.est in c("focei", "foce", "laplace", "agq", "posthoc", "nlm")) {
      expect_no_warning(.preProcessCensDistWarn(.ui, .est, .d, NULL))
    }
    # a method on a kernel that has no t()/cauchy() correction still warns
    for (.est in c("saem", "imp", "npag")) {
      expect_warning(.preProcessCensDistWarn(.ui, .est, .d, NULL),
                     "censoring ignored")
    }
  })
})
