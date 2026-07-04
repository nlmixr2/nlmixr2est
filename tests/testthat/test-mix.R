nmTest({

  test_that("test mixture models -- focei fit", {

    one.cmt <- function() {
      ini({
        tka <- 0.45
        tcl1 <- log(c(0, 2.7, 100))
        tcl2 <- log(c(0, 0.1, 120))
        tv <- 3.45
        p1 <- 0.3
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- mix(exp(tcl1 + eta.cl), p1, exp(tcl2 + eta.cl))
        v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }

    withr::with_seed(42, {
      fit <- nlmixr2(one.cmt, nlmixr2data::theo_sd, "focei",
                     control = foceiControl(print = 0, maxOuterIterations = 5))
    })

    # ranef should have ID + ETA columns only (no MIXEST)
    expect_true(all(c("ID", "eta.ka", "eta.cl", "eta.v") %in% names(fit$ranef)))
    expect_false("MIXEST" %in% names(fit$ranef))
    expect_equal(nrow(fit$ranef), length(unique(nlmixr2data::theo_sd$ID)))

    # mixNum: one row per subject
    mn <- fit$mixNum
    expect_true(is.data.frame(mn))
    expect_true(all(c("ID", "mixnum") %in% names(mn)))
    expect_equal(nrow(mn), length(unique(nlmixr2data::theo_sd$ID)))
    expect_true(all(mn$mixnum %in% c(1L, 2L)))

    # mixList: one component per mixture
    ml <- fit$mixList
    expect_true(is.list(ml))
    expect_equal(length(ml), 2L)
    expect_true(all(c("mix1", "mix2") %in% names(ml)))
    for (k in seq_along(ml)) {
      expect_true(all(c("ID", "eta.ka", "eta.cl", "eta.v", "prob") %in% names(ml[[k]])))
      expect_equal(nrow(ml[[k]]), length(unique(nlmixr2data::theo_sd$ID)))
      expect_true(all(ml[[k]]$prob >= 0 & ml[[k]]$prob <= 1))
    }

    # Residuals/table step completed (CWRES exists)
    expect_true("CWRES" %in% names(fit))

    # "Back-Transformed" rows in parHistData must show probabilities (in (0,1)).
    # "Unscaled" rows (fit$parHist) legitimately stay as mlogit (unbounded).
    phd <- fit$parHistData
    bt_rows <- as.character(phd$type) == "Back-Transformed"
    expect_true(any(bt_rows),
                label = "parHistData has Back-Transformed rows")
    expect_true(all(phd$p1[bt_rows] > 0 & phd$p1[bt_rows] < 1),
                label = "Back-Transformed rows: p1 is a valid probability in (0,1)")

    # Posterior mixture probabilities must sum to 1 per subject
    ml <- fit$mixList
    prob_sums <- ml[["mix1"]]$prob + ml[["mix2"]]$prob
    expect_equal(prob_sums, rep(1, nrow(ml[["mix1"]])), tolerance = 1e-8,
                 label = "posterior mixture probabilities sum to 1 per subject")
  })

  if (rxode2hasLlik()) {
    test_that("test mixture models -- focei fit with llik (dnorm) error", {
      one.cmt.ll <- function() {
        ini({
          tka <- 0.45
          tcl1 <- log(c(0, 2.7, 100))
          tcl2 <- log(c(0, 0.1, 120))
          tv <- 3.45
          p1 <- 0.3
          eta.ka ~ 0.6
          eta.cl ~ 0.3
          eta.v ~ 0.1
          add.sd <- 0.7
        })
        model({
          ka <- exp(tka + eta.ka)
          cl <- mix(exp(tcl1 + eta.cl), p1, exp(tcl2 + eta.cl))
          v <- exp(tv + eta.v)
          linCmt() ~ add(add.sd) + dnorm()
        })
      }

      withr::with_seed(42, {
        fit <- nlmixr2(one.cmt.ll, nlmixr2data::theo_sd, "focei",
                       control = foceiControl(print = 0, maxOuterIterations = 5))
      })

      # Should complete without error; p1 should be a valid probability
      expect_true("CWRES" %in% names(fit))
      expect_true(fit$fixef["p1"] > 0 & fit$fixef["p1"] < 1,
                  label = "mixture probability p1 in (0,1) for llik model")

      # Back-Transformed rows in parHistData must show probabilities
      phd <- fit$parHistData
      bt_rows <- as.character(phd$type) == "Back-Transformed"
      expect_true(any(bt_rows))
      expect_true(all(phd$p1[bt_rows] > 0 & phd$p1[bt_rows] < 1),
                  label = "parHistData Back-Transformed p1 in (0,1) for llik model")
    })
  }

  test_that("test mixture models -- ui components", {

    one.cmt <- function() {
      ini({
        ## You may label each parameter with a comment
        tka <- 0.45 # Log Ka
        tcl1 <- log(c(0, 2.7, 100)) # Log Cl
        tcl2 <- log(c(0, 0.1, 120)) # Log Cl
        ## This works with interactive models
        ## You may also label the preceding line with label("label text")
        tv <- 3.45; label("log V")
        p1 <- 0.3
        ## the label("Label name") works with all models
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- mix(exp(tcl1 + eta.cl), p1, exp(tcl2 + eta.cl))
        v <- exp(tv + eta.v)
        me <- mixest
        mn <- mixnum
        mu <- mixunif
        linCmt() ~ add(add.sd)
      })
    }

    ui <- rxode2::rxode2(one.cmt())

    expect_equal(ui$thetaIniMix,
                 c(tka = 0.45, tcl1 = 0.993251773010283, tcl2 = -2.30258509299405,
                   tv = 3.45, p1 = -0.847297860387204, add.sd = 0.7))

    expect_equal(ui$thetaMixIndex, 5L)


    one.cmt <- function() {
      ini({
        ## You may label each parameter with a comment
        tka <- 0.45 # Log Ka
        tcl1 <- log(c(0, 2.7, 100)) # Log Cl
        ## This works with interactive models
        ## You may also label the preceding line with label("label text")
        tv <- 3.45; label("log V")
        ## the label("Label name") works with all models
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl1 + eta.cl)
        v <- exp(tv + eta.v)
        me <- mixest
        mn <- mixnum
        mu <- mixunif
        linCmt() ~ add(add.sd)
      })
    }

    ui <- rxode2::rxode2(one.cmt())

    expect_equal(ui$thetaIniMix,
                 c(tka = 0.45, tcl1 = log(2.7), tv = 3.45, add.sd = 0.7))

    expect_equal(ui$thetaMixIndex, integer(0))

  })

  test_that("nlme estimation of a mix() model errors instead of silently dropping the mixture probability", {
    one.compartment.mix.nlme <- function() {
      ini({
        tka <- log(1.5)
        tcl1 <- log(1.0)
        tcl2 <- log(5.0)
        tv <- log(20)
        p1 <- 0.5
        eta.cl ~ 0.01
        eta.v ~ 0.01
        eta.ka ~ 0.01
        add.sd <- 0.05
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- mix(exp(tcl1 + eta.cl), p1, exp(tcl2 + eta.cl))
        v <- exp(tv + eta.v)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }
    d <- data.frame(ID = 1, TIME = c(0, 1), AMT = c(100, 0), EVID = c(1, 0),
                     DV = c(0, 1), CMT = c(1, 2))
    expect_error(
      .nlmixr(one.compartment.mix.nlme, d, est = "nlme", control = list(verbose = FALSE)),
      "mix\\(\\) models are not supported"
    )
  })

  test_that("saemMixProb validates initial mixture probabilities", {
    # Test that valid mixture probabilities work (should fit without error)
    valid.mix <- function() {
      ini({
        tka <- 0.45
        tcl1 <- log(c(0, 2.7, 100))
        tcl2 <- log(c(0, 0.1, 120))
        tv <- 3.45
        p1 <- 0.3  # Valid: in [0, 1] and sum <= 1
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- mix(exp(tcl1 + eta.cl), p1, exp(tcl2 + eta.cl))
        v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    # This should complete without error (valid probabilities pass validation)
    withr::with_seed(42, {
      fit <- .nlmixr(valid.mix, nlmixr2data::theo_sd, est = "saem",
                     saemControl(print = 0, nBurn = 1, nEm = 1))
    })
    expect_true(!is.null(fit))
    # Verify the mixture probability is in valid range
    expect_true(fit$fixef["p1"] > 0 & fit$fixef["p1"] < 1)
  })

  test_that("rxUiGet.saemMixProb rejects invalid mixture probabilities directly", {
    # Mock UI objects use mixProbs as a character vector of theta names,
    # matching the real UI (ui$mixProbs is character, e.g. c("p1","p2"),
    # not integer positions -- see R/focei.R's %in% .theta$name usage).
    # Test case: probability > 1
    ui.obj <- list(list(mixProbs = "p1", theta = c(p1 = 1.5)))
    expect_error(
      rxUiGet.saemMixProb(ui.obj),
      "mixture probabilities must each be in \\[0, 1\\]"
    )

    # Test case: probability < 0
    ui.obj <- list(list(mixProbs = "p1", theta = c(p1 = -0.1)))
    expect_error(
      rxUiGet.saemMixProb(ui.obj),
      "mixture probabilities must each be in \\[0, 1\\]"
    )

    # Test case: sum > 1
    ui.obj <- list(list(mixProbs = c("p1", "p2"), theta = c(p1 = 0.7, p2 = 0.6)))
    expect_error(
      rxUiGet.saemMixProb(ui.obj),
      "mixture probabilities must each be in \\[0, 1\\]"
    )

    # Test case: valid probabilities should not error
    ui.obj <- list(list(mixProbs = "p1", theta = c(p1 = 0.5)))
    result <- rxUiGet.saemMixProb(ui.obj)
    expect_equal(result, c(p1 = 0.5, 0.5))

    # Test case: valid with multiple probabilities
    ui.obj <- list(list(mixProbs = c("p1", "p2"), theta = c(p1 = 0.3, p2 = 0.2)))
    result <- rxUiGet.saemMixProb(ui.obj)
    expect_equal(result, c(p1 = 0.3, p2 = 0.2, 0.5))
  })

  test_that("rxUiGet.thetaIniMix errors on invalid initial mixture probabilities", {
    # Mock UI objects use mixProbs as a character vector of theta names,
    # matching the real UI (ui$mixProbs is character, e.g. c("p1","p2")).
    # Same mocking pattern as the rxUiGet.saemMixProb tests above -- this
    # bypasses rxode2::mix()'s own UDF check (which only rejects at UI-build
    # time when the *sum* of probabilities is outside (0, 1); a real UI whose
    # probabilities are individually invalid but sum to something plausible,
    # e.g. p1 = 1.5, p2 = -0.6 summing to 0.9, builds successfully and only
    # rxUiGet.thetaIniMix can catch it before theta is left on the wrong
    # scale).

    # Single probability > 1
    ui.obj <- list(list(mixProbs = "p1", theta = c(tka = 0.45, p1 = 1.2)))
    expect_error(rxUiGet.thetaIniMix(ui.obj), "invalid")

    # Single probability < 0
    ui.obj <- list(list(mixProbs = "p1", theta = c(tka = 0.45, p1 = -0.1)))
    expect_error(rxUiGet.thetaIniMix(ui.obj), "invalid")

    # Two probabilities, each individually invalid, but summing to a
    # plausible-looking value (0.9) -- the real-world gap this guards against.
    ui.obj <- list(list(mixProbs = c("p1", "p2"),
                        theta = c(tka = 0.45, p1 = 1.5, p2 = -0.6)))
    expect_error(rxUiGet.thetaIniMix(ui.obj), "invalid")

    # Two valid-individually probabilities whose sum exceeds 1
    ui.obj <- list(list(mixProbs = c("p1", "p2"),
                        theta = c(tka = 0.45, p1 = 0.7, p2 = 0.6)))
    expect_error(rxUiGet.thetaIniMix(ui.obj), "invalid")

    # Valid probabilities should not error, and should be mlogit-transformed
    ui.obj <- list(list(mixProbs = "p1", theta = c(tka = 0.45, p1 = 0.5)))
    result <- expect_no_error(rxUiGet.thetaIniMix(ui.obj))
    expect_equal(unname(result["p1"]), rxode2::mlogit(0.5))
    expect_equal(unname(result["tka"]), 0.45)
  })

  test_that("rxUiGet.thetaIniMix errors (not just warns) on a real UI with individually invalid mixture probabilities", {
    # This mirrors a real 3-population mixture model (two explicit
    # probabilities) where rxode2::mix()'s own parse-time check only looks at
    # the sum (0.9, which passes) and so does not catch that p1 and p2 are
    # each individually out of [0, 1]. If rxUiGet.thetaIniMix silently left
    # these untransformed, the raw, wrong-scale values would flow into focei's
    # initial parameter vector and crash later with a confusing, unrelated
    # error (e.g. "infinite while evaluating initial objective function")
    # instead of a clear diagnostic -- verified empirically while implementing
    # this fix.
    threePop <- function() {
      ini({
        tka   <- log(1.5)
        tcl1  <- log(1.0)
        tcl2  <- log(3.0)
        tcl3  <- log(6.0)
        tv    <- log(20)
        p1    <- 1.5   # invalid: > 1
        p2    <- -0.6  # invalid: < 0 (sum = 0.9, looks plausible)
        eta.cl1 ~ 0.01
        eta.cl2 ~ 0.01
        eta.cl3 ~ 0.01
        eta.v  ~ 0.01
        add.sd <- 0.05
      })
      model({
        ka <- exp(tka)
        cl <- mix(exp(tcl1 + eta.cl1), p1, exp(tcl2 + eta.cl2), p2, exp(tcl3 + eta.cl3))
        v  <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    # UI construction succeeds -- rxode2's own check only validates the sum
    ui <- expect_no_error(rxode2::rxode2(threePop))
    expect_error(
      rxUiGet.thetaIniMix(list(ui)),
      "invalid"
    )
  })

  test_that(".mixFix guards against a zero row total (all components underflow for a subject)", {
    .mixFix <- nlmixr2est:::.mixFix
    env <- new.env()
    assign("mixIdx", 1L, envir = env)  # non-empty so the early-return guards pass
    assign("etaObfFull", data.frame(
      ID = c(1L, 2L, 1L, 2L),
      MIXEST = c(1L, 1L, 2L, 2L),
      `ETA[1]` = c(0.1, 0.2, 0.1, 0.2),
      OBJI = c(1e6, 1e6, 1e6, 1e6),  # underflows exp(-0.5*OBJI) to exactly 0 for BOTH components
      check.names = FALSE
    ), envir = env)
    assign("etaObf", data.frame(ID = c(1L, 2L), MIXEST = c(1L, 1L),
                                 `ETA[1]` = c(0.1, 0.2), check.names = FALSE), envir = env)
    assign("ranef", data.frame(ID = c(1L, 2L), `ETA[1]` = c(0.1, 0.2), check.names = FALSE), envir = env)
    assign("fixef", c(p1 = 0), envir = env)  # mlogit(0.5) prior when back-transformed

    expect_warning(
      nlmixr2est:::.mixFix(env, ui = NULL),
      "underflow|zero.*likelihood|posterior"
    )
    .mixList <- get("mixList", envir = env)
    expect_false(any(vapply(.mixList, function(m) any(is.nan(m$prob)), logical(1))))
  })

})
