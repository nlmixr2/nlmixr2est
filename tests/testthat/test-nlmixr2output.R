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
    data.frame(ch = "5.00", v = 25)
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
  expect_equal(
    fit$parFixed,
    structure(
      list(
        Est. = formatMinWidth(fit$parFixedDf$Estimate),
        SE = formatMinWidth(fit$parFixedDf$SE, naValue = ""),
        `%RSE` = formatMinWidth(fit$parFixedDf$`%RSE`, naValue = ""),
        `Back-transformed(95%CI)` =
          c(
            sprintf(
              "%s (%s, %s)",
              formatMinWidth(fit$parFixedDf$`Back-transformed`),
              formatMinWidth(fit$parFixedDf$`CI Lower`),
              formatMinWidth(fit$parFixedDf$`CI Upper`)
            )[!is.na(fit$parFixedDf$`CI Upper`)],
            formatMinWidth(fit$parFixedDf$`Back-transformed`[is.na(fit$parFixedDf$`CI Upper`)])
          ),
        `BSV(CV%)` = formatMinWidth(fit$parFixedDf$`BSV(CV%)`, naValue = ""),
        `Shrink(SD)%` = paste0(formatMinWidth(fit$parFixedDf$`Shrink(SD)%`, naValue = ""), c("", ">", "", ""))
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
  expect_equal(
    fitFixed$parFixed,
    structure(
      list(
        Est. = formatMinWidth(fitFixed$parFixedDf$Estimate),
        SE = c("FIXED", formatMinWidth(fitFixed$parFixedDf$SE[2:3]), ""),
        `%RSE` = c("FIXED", formatMinWidth(fitFixed$parFixedDf$`%RSE`[2:3]), ""),
        `Back-transformed(95%CI)` = c(
          formatMinWidth(fitFixed$parFixedDf$`Back-transformed`[1]),
          sprintf(
            "%s (%s, %s)",
            formatMinWidth(fitFixed$parFixedDf$`Back-transformed`),
            formatMinWidth(fitFixed$parFixedDf$`CI Lower`),
            formatMinWidth(fitFixed$parFixedDf$`CI Upper`)
          )[2:3],
          formatMinWidth(fitFixed$parFixedDf$`Back-transformed`[4])
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
  expect_equal(
    fitFixedLabel$parFixed,
    structure(
      list(
        Parameter = c("ka", "clearance", "", ""),
        Est. = formatMinWidth(fitFixedLabel$parFixedDf$Estimate),
        SE = c("FIXED", formatMinWidth(fitFixedLabel$parFixedDf$SE[2:3]), ""),
        `%RSE` = c("FIXED", formatMinWidth(fitFixedLabel$parFixedDf$`%RSE`[2:3]), ""),
        `Back-transformed(95%CI)` = c(
          formatMinWidth(fitFixedLabel$parFixedDf$`Back-transformed`[1]),
          sprintf(
            "%s (%s, %s)",
            formatMinWidth(fitFixedLabel$parFixedDf$`Back-transformed`),
            formatMinWidth(fitFixedLabel$parFixedDf$`CI Lower`),
            formatMinWidth(fitFixedLabel$parFixedDf$`CI Upper`)
          )[2:3],
          formatMinWidth(fitFixedLabel$parFixedDf$`Back-transformed`[4])
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
