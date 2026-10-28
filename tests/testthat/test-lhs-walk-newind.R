# A model reads `newind` as 1 at a subject's first record and 2 after it, as
# rxode2's output gives.  Each estimator's own walk over the records gives the
# same, so a model that reads newind fits exactly as the model that reads a
# data column flagging each subject's first record.

# a one-compartment model with an offset after the first record, read from
# newind or from the data's `first` column
.newindMod <- function(flag, eta = TRUE) {
  .after <- if (flag) "if (first == 0) e0 <- 0.5" else "if (newind > 1) e0 <- 0.5"
  eval(parse(
    text = sprintf(
      "function() {
        ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; %s add.sd <- 0.7 })
        model({
          ka <- exp(tka)
          cl <- exp(tcl%s)
          v <- exp(tv)
          e0 <- 0
          %s
          d/dt(depot) <- -ka * depot
          d/dt(center) <- ka * depot - cl / v * center
          cp <- center / v + e0
          cp ~ add(add.sd)
        })
      }",
      if (eta) "eta.cl ~ 0.1;" else "",
      if (eta) " + eta.cl" else "",
      .after
    )
  ))
}

# theo_sd with each dose 0.01 h after the time-0 sample, so a subject's first
# record is an observation, flagged by `first`
.newindData <- function() {
  .d <- nlmixr2data::theo_sd
  .d$TIME[.d$EVID != 0] <- 0.01
  .d <- .d[order(.d$ID, .d$TIME), ]
  .d$first <- as.integer(!duplicated(.d$ID))
  .d
}

.newindFit <- function(flag, est, control, eta = TRUE) {
  suppressWarnings(suppressMessages(nlmixr2(.newindMod(flag, eta), .newindData(), est, control = control)))
}

nmTest({
  test_that("a model reading newind fits as one reading the first-record flag", {
    .ests <- list(
      focei = foceiControl(print = 0, covMethod = "", maxOuterIterations = 2L),
      foce = foceControl(print = 0, covMethod = "", maxOuterIterations = 2L),
      saem = saemControl(print = 0, nBurn = 5, nEm = 5, nmc = 1, seed = 42, covMethod = "")
    )
    for (.est in names(.ests)) {
      .flag <- .newindFit(TRUE, .est, .ests[[.est]])
      .new <- .newindFit(FALSE, .est, .ests[[.est]])
      expect_equal(.new$objf, .flag$objf, tolerance = 1e-8, label = .est)
      expect_equal(.new$theta, .flag$theta, tolerance = 1e-8, label = .est)
      expect_equal(.new$IPRED, .flag$IPRED, tolerance = 1e-8, label = .est)
    }
    # nlm fits a model without etas
    .flag <- .newindFit(TRUE, "nlm", list(print = 0, covMethod = ""), eta = FALSE)
    .new <- .newindFit(FALSE, "nlm", list(print = 0, covMethod = ""), eta = FALSE)
    expect_equal(.new$objf, .flag$objf, tolerance = 1e-8)
    expect_equal(.new$theta, .flag$theta, tolerance = 1e-8)
  })
})
