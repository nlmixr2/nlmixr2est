# Censored NPDE is calibrated: data simulated from the model and censored
# below an LOQ, NPDE at the true parameters, should be N(0, 1)
nmTest({
  test_that("NPDE of censored data simulated from the model is N(0, 1)", {
    skip_on_cran()
    .mod <- function() {
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
    .d0 <- nlmixr2data::theo_sd
    .d <- do.call(
      rbind,
      lapply(1:10, function(r) {
        .x <- .d0
        .x$ID <- .x$ID + 100 * r
        .x
      })
    )
    .s <- rxode2::rxSolve(.mod, .d, returnType = "data.frame", seed = 42)
    .d$DV[.d$EVID == 0] <- .s$sim
    # about half of the observations
    .d$CENS <- ifelse(.d$EVID == 0 & .d$DV < 5, 1, 0)
    .d$DV[.d$CENS == 1] <- 5
    expect_gt(mean(.d$CENS[.d$EVID == 0]), 0.4)
    for (.cm in c("truncated-normal", "ipred", "pred")) {
      .fit <- suppressWarnings(suppressMessages(nlmixr2(
        .mod,
        .d,
        "posthoc",
        table = tableControl(npde = TRUE, censMethod = .cm, nsim = 300)
      )))
      .z <- .fit$NPDE
      expect_true(all(is.finite(.z)), label = .cm)
      expect_lt(abs(mean(.z)), 0.1, label = .cm)
      expect_lt(abs(sd(.z) - 1), 0.08, label = .cm)
      expect_gt(suppressWarnings(stats::ks.test(.z, "pnorm")$p.value), 0.01, label = .cm)
    }
  })
})
