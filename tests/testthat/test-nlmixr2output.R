test_that(".updateParFixedGetEtaRow returns correct values", {
  envPrep <- new.env()
  expect_equal(
    .updateParFixedGetEtaRow(
      .eta = "iivemax",
      .env = envPrep,
      .ome = matrix(25, nrow = 1, dimnames = list("iivemax", "iivemax")),
      .omegaFix = c(iivemax = FALSE),
      .muRefCurEval = data.frame(parameter = "iivemax", curEval = "", low = NA_real_, hi = NA_real_),
      .sigdig = 3L
    ),
    data.frame(ch = "5.00", v = 5)
  )
  expect_false(envPrep$.cvOnly)

  envPrep <- new.env()
  expect_equal(
    .updateParFixedGetEtaRow(
      .eta = "iivemax",
      .env = envPrep,
      .ome = matrix(0.4, nrow = 1, dimnames = list("iivemax", "iivemax")),
      .omegaFix = c(iivemax = FALSE),
      .muRefCurEval = data.frame(parameter = "iivemax", curEval = "exp", low = NA_real_, hi = NA_real_),
      .sigdig = 3L
    ),
    data.frame(ch = "70.1", v = sqrt(exp(0.4) - 1) * 100)
  )
  expect_false(envPrep$.sdOnly)
})

test_that("formatMinWidth", {
  # Special values
  expect_equal(
    formatMinWidth(x = c(NA, 0, Inf, -Inf, NaN)),
    c("NA", "0", "Inf", "-Inf", "NaN")
  )
  # Rounding occurs to requested significant digits
  expect_equal(
    formatMinWidth(x = -123456*10^(-10:10)),
    c("-1.23e-5", "-1.23e-4", "-0.00123", "-0.0123", "-0.123", "-1.23",
      "-12.3", "-123", "-1230", "-12300", "-123000", "-1.23e6", "-1.23e7",
      "-1.23e8", "-1.23e9", "-1.23e10", "-1.23e11", "-1.23e12", "-1.23e13",
      "-1.23e14", "-1.23e15")
  )
  # Rounding up works as expected; scientific notation values drop extraneous
  # zeros in the exponent
  expect_equal(
    formatMinWidth(x = -9999*10^(-10:10)),
    c("-1.00e-6", "-1.00e-5", "-1.00e-4", "-0.00100", "-0.0100",
      "-0.100", "-1.00", "-10.0", "-100", "-1000", "-10000", "-100000",
      "-1.00e6", "-1.00e7", "-1.00e8", "-1.00e9", "-1.00e10", "-1.00e11",
      "-1.00e12", "-1.00e13", "-1.00e14")
  )
  # Planned significant digits are shown, including when digits are added
  expect_equal(
    formatMinWidth(x = 12*10^(-10:10)),
    c("1.20e-9", "1.20e-8", "1.20e-7", "1.20e-6", "1.20e-5", "1.20e-4",
      "0.00120", "0.0120", "0.120", "1.20", "12.0", "120", "1200",
      "12000", "120000", "1.20e6", "1.20e7", "1.20e8", "1.20e9", "1.20e10",
      "1.20e11")
  )
  # Negative values
  expect_equal(
    formatMinWidth(x = -12*10^(-10:10)),
    c("-1.20e-9", "-1.20e-8", "-1.20e-7", "-1.20e-6", "-1.20e-5",
      "-1.20e-4", "-0.00120", "-0.0120", "-0.120", "-1.20", "-12.0",
      "-120", "-1200", "-12000", "-120000", "-1.20e6", "-1.20e7", "-1.20e8",
      "-1.20e9", "-1.20e10", "-1.20e11")
  )
  # input must be numeric
  expect_error(
    formatMinWidth("A"),
    regexp = "Assertion on 'x' failed: Must be of type 'numeric', not 'character'.",
    fixed = TRUE
  )
})

test_that("formatMinWidth in parFixed", {
  one.compartment <- function() {
    ini({
      tka <- log(1.57)
      tcl <- log(2.72)
      bsvCl ~ 0.1
      tv <- log(31.5)
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl + bsvCl)
      v <- exp(tv)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }

  # Simple ----
  fit <- .nlmixr(one.compartment, theo_sd, est="focei", control = foceiControlFast)
  # The table is formatted with sigdigTable (the dedicated "final output table"
  # digits), which defaults to sigdig.  Deriving the expected digits from the
  # control keeps the test robust to the foceiControl() sigdig default.
  .digs <- fit$control$sigdigTable
  expect_equal(
    fit$parFixed,
    structure(
      list(
        Est. = formatMinWidth(fit$parFixedDf$Estimate, digits = .digs),
        SE = formatMinWidth(fit$parFixedDf$SE, digits = .digs, naValue = ""),
        `%RSE` = formatMinWidth(fit$parFixedDf$`%RSE`, digits = .digs, naValue = ""),
        `Back-transformed(95%CI)` =
          c(
            sprintf(
              "%s (%s, %s)",
              formatMinWidth(fit$parFixedDf$`Back-transformed`, digits = .digs),
              formatMinWidth(fit$parFixedDf$`CI Lower`, digits = .digs),
              formatMinWidth(fit$parFixedDf$`CI Upper`, digits = .digs)
            )[!is.na(fit$parFixedDf$`CI Upper`)],
            formatMinWidth(fit$parFixedDf$`Back-transformed`[is.na(fit$parFixedDf$`CI Upper`)], digits = .digs)
          ),
        `BSV(CV%)` = formatMinWidth(fit$parFixedDf$`BSV(CV%)`, digits = .digs, naValue = ""),
        `Shrink(SD)%` = paste0(formatMinWidth(fit$parFixedDf$`Shrink(SD)%`, digits = .digs, naValue = ""), c("", ">", "", ""))
      ),
      class = c("nlmixr2ParFixed", "data.frame"),
      row.names = c("tka", "tcl", "tv", "add.sd")
    )
  )

  # Fixed parameter ----
  one.compartment.fixed <- function() {
    ini({
      tka <- fixed(log(1.57))
      tcl <- log(2.72)
      tv <- log(31.5)
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl)
      v <- exp(tv)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }

  fitFixed <- .nlmixr(one.compartment.fixed, theo_sd, est="focei", control = foceiControlFast)
  .digsF <- fitFixed$control$sigdigTable
  expect_equal(
    fitFixed$parFixed,
    structure(
      list(
        Est. = formatMinWidth(fitFixed$parFixedDf$Estimate, digits = .digsF),
        SE = c("FIXED", formatMinWidth(fitFixed$parFixedDf$SE[2:3], digits = .digsF), ""),
        `%RSE` = c("FIXED", formatMinWidth(fitFixed$parFixedDf$`%RSE`[2:3], digits = .digsF), ""),
        `Back-transformed(95%CI)` = c(
          formatMinWidth(fitFixed$parFixedDf$`Back-transformed`[1], digits = .digsF),
          sprintf(
            "%s (%s, %s)",
            formatMinWidth(fitFixed$parFixedDf$`Back-transformed`, digits = .digsF),
            formatMinWidth(fitFixed$parFixedDf$`CI Lower`, digits = .digsF),
            formatMinWidth(fitFixed$parFixedDf$`CI Upper`, digits = .digsF)
          )[2:3],
          formatMinWidth(fitFixed$parFixedDf$`Back-transformed`[4], digits = .digsF)
        ),
        `BSV(SD)` = c("", "", "", ""),
        `Shrink(SD)%` = c("", "", "", "")
      ),
      class = c("nlmixr2ParFixed", "data.frame"),
      row.names = c("tka", "tcl", "tv", "add.sd")
    )
  )

  # Fixed parameter, labeled ----
  one.compartment.labeled <- function() {
    ini({
      tka <- fixed(log(1.57)); label("ka")
      tcl <- log(2.72); label("clearance")
      tv <- log(31.5)
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl)
      v <- exp(tv)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }

  fitFixedLabel <- .nlmixr(one.compartment.labeled, theo_sd,  est="focei", control = list(print = 0))
  .digsL <- fitFixedLabel$control$sigdigTable
  expect_equal(
    fitFixedLabel$parFixed,
    structure(
      list(
        Parameter = c("ka", "clearance", "", ""),
        Est. = formatMinWidth(fitFixedLabel$parFixedDf$Estimate, digits = .digsL),
        SE = c("FIXED", formatMinWidth(fitFixedLabel$parFixedDf$SE[2:3], digits = .digsL), ""),
        `%RSE` = c("FIXED", formatMinWidth(fitFixedLabel$parFixedDf$`%RSE`[2:3], digits = .digsL), ""),
        `Back-transformed(95%CI)` = c(
          formatMinWidth(fitFixedLabel$parFixedDf$`Back-transformed`[1], digits = .digsL),
          sprintf(
            "%s (%s, %s)",
            formatMinWidth(fitFixedLabel$parFixedDf$`Back-transformed`, digits = .digsL),
            formatMinWidth(fitFixedLabel$parFixedDf$`CI Lower`, digits = .digsL),
            formatMinWidth(fitFixedLabel$parFixedDf$`CI Upper`, digits = .digsL)
          )[2:3],
          formatMinWidth(fitFixedLabel$parFixedDf$`Back-transformed`[4], digits = .digsL)
        ),
        `BSV(SD)` = c("", "", "", ""),
        `Shrink(SD)%` = c("", "", "", "")
      ),
      class = c("nlmixr2ParFixed", "data.frame"),
      row.names = c("tka", "tcl", "tv", "add.sd")
    )
  )

  # Works with .ret$control$ci and .ret$control$sigdig ----
  fitFixedLabelCI <- .nlmixr(one.compartment.labeled, theo_sd,  est="focei", control = list(print = 0, ci = 0.9, sigdig = 4))
  expect_equal(
    fitFixedLabelCI$parFixed,
    structure(
      list(
        Parameter = c("ka", "clearance", "", ""),
        Est. = formatMinWidth(fitFixedLabelCI$parFixedDf$Estimate, digits = 4),
        SE = c("FIXED", formatMinWidth(fitFixedLabelCI$parFixedDf$SE[2:3], digits = 4), ""),
        `%RSE` = c("FIXED", formatMinWidth(fitFixedLabelCI$parFixedDf$`%RSE`[2:3], digits = 4), ""),
        `Back-transformed(90%CI)` =
          c(
            formatMinWidth(fitFixedLabelCI$parFixedDf$`Back-transformed`[1], digits = 4),
            sprintf(
              "%s (%s, %s)",
              formatMinWidth(fitFixedLabelCI$parFixedDf$`Back-transformed`, digits = 4),
              formatMinWidth(fitFixedLabelCI$parFixedDf$`CI Lower`, digits = 4),
              formatMinWidth(fitFixedLabelCI$parFixedDf$`CI Upper`, digits = 4)
            )[2:3],
            formatMinWidth(fitFixedLabelCI$parFixedDf$`Back-transformed`[4], digits = 4)),
        `BSV(SD)` = c("", "", "", ""),
        `Shrink(SD)%` = c("", "", "", "")
      ),
      class = c("nlmixr2ParFixed", "data.frame"),
      row.names = c("tka", "tcl", "tv", "add.sd")
    )
  )
})

test_that(".updateParFixedGetEtaRow returns a numeric `v` when there is no BSV", {
  # When the omega is NULL (no BSV in the model) or the eta is absent from the
  # omega matrix (e.g. a fixed/zero BSV that the estimator drops), the row must
  # still be a one-row data.frame with a *numeric* `v`.  Returning a bare ""
  # coerced the `v` column to character in .updateParFixedAddBsv() (see the
  # regression test below).
  envPrep <- new.env()
  expect_equal(
    .updateParFixedGetEtaRow(
      .eta = "eta.ka",
      .env = envPrep,
      .ome = NULL,
      .omegaFix = c(eta.ka = FALSE),
      .muRefCurEval = data.frame(parameter = "eta.ka", curEval = "exp", low = NA_real_, hi = NA_real_),
      .sigdig = 3L
    ),
    data.frame(ch = "", v = NA_real_)
  )

  envPrep <- new.env()
  .row <-
    .updateParFixedGetEtaRow(
      .eta = "eta.ka",
      .env = envPrep,
      .ome = matrix(0.1, nrow = 1, dimnames = list("eta.cl", "eta.cl")),
      .omegaFix = c(eta.ka = FALSE, eta.cl = FALSE),
      .muRefCurEval = data.frame(parameter = "eta.ka", curEval = "exp", low = NA_real_, hi = NA_real_),
      .sigdig = 3L
    )
  expect_equal(.row, data.frame(ch = "", v = NA_real_))
  expect_true(is.numeric(.row$v))
})

test_that(".updateParFixedAddBsv keeps the BSV column numeric when a mu-referenced eta is absent from omega", {
  # Regression test: a mu-referenced eta that is missing from the omega matrix
  # used to make .updateParFixedGetEtaRow() return a bare "", which coerced the
  # whole BSV `v` column to character in do.call("rbind", .cvp).  That left the
  # numeric BSV values as full-precision strings (e.g. "59.1488636895083") and
  # printed the missing values as "<NA>" instead of "".
  ui <- one.compartment

  popDf <- data.frame(
    Estimate = c(0.45, 1, 3.45, 0.7),
    row.names = c("tka", "tcl", "tv", "add.sd"),
    check.names = FALSE
  )
  # omega is missing eta.ka (tka is still mu-referenced to it)
  omega <- diag(c(0.3, 0.1))
  dimnames(omega) <- list(c("eta.cl", "eta.v"), c("eta.cl", "eta.v"))

  res <-
    .updateParFixedAddBsv(
      popDf, iniDf = ui$iniDf, omega = omega, .sigdig = 3L,
      .muRefDataFrame = ui$muRefDataFrame, .muRefCurEval = ui$muRefCurEval
    )

  # The BSV column must stay numeric (not coerced to character)
  expect_true(is.numeric(res$popDf[["BSV(CV%)"]]))
  expect_equal(
    res$popDf[["BSV(CV%)"]],
    c(NA_real_, sqrt(exp(0.3) - 1) * 100, sqrt(exp(0.1) - 1) * 100, NA_real_)
  )

  # And the formatted $parFixed shows 3 significant figures with "" (not "<NA>")
  # for the parameters that have no BSV.
  fmt <-
    .updateParFixedApplySig(
      res$popDf, digits = 3L, ci = 0.95,
      fixedNames = character(), bsvFixedNames = res$bsvFixedNames
    )
  expect_equal(fmt[["BSV(CV%)"]], c("", "59.1", "32.4", ""))
})

test_that("$parFixed honors sigdigTable and ci from the control for models with FIXED parameters", {
  # Regression test: .updateParFixed() read sigdig/ci via rxGetControl(.ui, ...),
  # but for models with fixed parameters .ui is the unfixed model, which carries
  # no control -- so the user's sigdigTable/ci were silently dropped and the
  # table used the defaults (3 significant digits, 95% CI).  They must come from
  # the fit's control object instead.
  mod <- function() {
    ini({
      tka <- fixed(log(1.57))
      tcl <- log(2.72)
      tv <- log(31.5)
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka); cl <- exp(tcl); v <- exp(tv)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }
  fit <- .nlmixr(mod, theo_sd, est = "focei", control = list(print = 0, ci = 0.9, sigdig = 4))
  # ci propagates to the back-transformed column name
  expect_true(any(grepl("90%CI", names(fit$parFixed), fixed = TRUE)))
  # sigdigTable (= sigdig = 4 here) propagates to the formatted estimate
  expect_equal(fit$parFixed["tcl", "Est."], formatMinWidth(fit$parFixedDf["tcl", "Estimate"], digits = 4))
})
