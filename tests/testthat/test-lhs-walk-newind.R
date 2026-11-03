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

# theo_sd as shipped, each subject's dose its first record, or with each dose
# 0.01 h after the time-0 sample, so a subject's first record is an
# observation; `first` flags the first record
.newindData <- function(doseFirst) {
  .d <- nlmixr2data::theo_sd
  if (!doseFirst) {
    .d$TIME[.d$EVID != 0] <- 0.01
  }
  .d <- .d[order(.d$ID, .d$TIME), ]
  .d$first <- as.integer(!duplicated(.d$ID))
  .d
}

.newindFit <- function(flag, est, control, doseFirst, eta = TRUE) {
  suppressWarnings(suppressMessages(nlmixr2(.newindMod(flag, eta), .newindData(doseFirst), est, control = control)))
}

# each estimator's fit of the newind model equals the flag model's
.newindExpectSame <- function(est, control, doseFirst, label) {
  .flag <- .newindFit(TRUE, est, control, doseFirst)
  .new <- .newindFit(FALSE, est, control, doseFirst)
  expect_equal(.new$objf, .flag$objf, tolerance = 1e-8, label = label)
  expect_equal(.new$theta, .flag$theta, tolerance = 1e-8, label = label)
  expect_equal(.new$IPRED, .flag$IPRED, tolerance = 1e-8, label = label)
}

test_that(".rxode2Optional() reports setIndNewind from rxode2's pointer table", {
  .p <- rxode2::.rxode2ptrs()
  .has <- length(.p) > 104L && names(.p)[105L] == "rxode2setIndNewind"
  expect_identical(.rxode2Optional(), c(setIndNewind = .has))
  on.exit(.Call(`_nlmixr2est_iniRxode2Optional`, .p, PACKAGE = "nlmixr2est"), add = TRUE)
  # a table that ends before slot 104, or names another entry there, loads none
  .Call(`_nlmixr2est_iniRxode2Optional`, .p[1:99], PACKAGE = "nlmixr2est")
  expect_identical(.rxode2Optional(), c(setIndNewind = FALSE))
  if (.has) {
    .q <- .p
    names(.q)[105L] <- "rxode2other"
    .Call(`_nlmixr2est_iniRxode2Optional`, .q, PACKAGE = "nlmixr2est")
    expect_identical(.rxode2Optional(), c(setIndNewind = FALSE))
    .Call(`_nlmixr2est_iniRxode2Optional`, .p, PACKAGE = "nlmixr2est")
    expect_identical(.rxode2Optional(), c(setIndNewind = TRUE))
  }
})

nmTest({
  test_that("a model reading newind fits as one reading the first-record flag", {
    # rxode2's setIndNewind() (slot 104) sets the newind the walks give; the
    # skip goes with test-dependency-workarounds.R's slot-104 block
    skip_if_not(.rxode2Optional()[["setIndNewind"]], "the installed rxode2 has no setIndNewind()")
    .ests <- list(
      list("focei", foceiControl(print = 0, covMethod = "", maxOuterIterations = 2L), "focei"),
      list("foce", foceControl(print = 0, covMethod = "", maxOuterIterations = 2L), "foce"),
      # these move the thetas in two iterations, and the second walks the
      # observations alone for its inner Hessian
      list("foceif", foceiControl(print = 0, covMethod = "", maxOuterIterations = 2L), "foceif"),
      list(
        "focei",
        foceiControl(print = 0, covMethod = "", maxOuterIterations = 2L, fast = TRUE, innerHessian = "conditional"),
        "focei conditional"
      ),
      list("saem", saemControl(print = 0, nBurn = 5, nEm = 5, nmc = 1, seed = 42, covMethod = ""), "saem")
    )
    for (.doseFirst in c(TRUE, FALSE)) {
      for (.e in .ests) {
        .newindExpectSame(.e[[1]], .e[[2]], .doseFirst, paste(.e[[3]], if (.doseFirst) "dose first" else "obs first"))
      }
      # nlm fits a model without etas
      .flag <- .newindFit(TRUE, "nlm", list(print = 0, covMethod = ""), .doseFirst, eta = FALSE)
      .new <- .newindFit(FALSE, "nlm", list(print = 0, covMethod = ""), .doseFirst, eta = FALSE)
      expect_equal(.new$objf, .flag$objf, tolerance = 1e-8)
      expect_equal(.new$theta, .flag$theta, tolerance = 1e-8)
    }
  })
})
