nmTest({
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

  test_that(".updateParFixedGetEtaRow returns a numeric `v` when there is no BSV", {
    # The rows are rbound and their `v` column taken as the BSV values, so a
    # row of a different type there turns the whole column character.
    expect_equal(
      .updateParFixedGetEtaRow(
        .eta = "iivemax",
        .env = new.env(),
        .ome = NULL,
        .omegaFix = c(iivemax = FALSE),
        .muRefCurEval = data.frame(parameter = "iivemax", curEval = "", low = NA_real_, hi = NA_real_),
        .sigdig = 3L
      ),
      data.frame(ch = "", v = NA_real_)
    )
    # an eta that is not in the omega matrix, e.g. a fixed BSV parameter
    .row <-
      .updateParFixedGetEtaRow(
        .eta = "iivemax",
        .env = new.env(),
        .ome = matrix(0.4, nrow = 1, dimnames = list("iivcl", "iivcl")),
        .omegaFix = c(iivemax = FALSE),
        .muRefCurEval = data.frame(parameter = "iivemax", curEval = "", low = NA_real_, hi = NA_real_),
        .sigdig = 3L
      )
    expect_equal(.row, data.frame(ch = "", v = NA_real_))
    expect_true(is.numeric(.row$v))
  })

  test_that(".updateParFixedAddBsv keeps the BSV column numeric when a mu-referenced eta is absent from omega", {
    ui <- one.compartment

    popDf <- data.frame(
      Estimate = c(0.45, 1, 3.45, 0.7),
      row.names = c("tka", "tcl", "tv", "add.sd"),
      check.names = FALSE
    )
    # omega is missing eta.ka, and tka is still mu-referenced to it
    omega <- diag(c(0.3, 0.1))
    dimnames(omega) <- list(c("eta.cl", "eta.v"), c("eta.cl", "eta.v"))

    res <-
      .updateParFixedAddBsv(
        popDf,
        iniDf = ui$iniDf,
        omega = omega,
        .sigdig = 3L,
        .muRefDataFrame = ui$muRefDataFrame,
        .muRefCurEval = ui$muRefCurEval
      )

    expect_true(is.numeric(res$popDf[["BSV(CV%)"]]))
    expect_equal(
      res$popDf[["BSV(CV%)"]],
      c(NA_real_, sqrt(exp(0.3) - 1) * 100, sqrt(exp(0.1) - 1) * 100, NA_real_)
    )

    # and the formatted table shows the requested significant figures, with ""
    # for the parameters that have no BSV
    fmt <-
      .updateParFixedApplySig(
        res$popDf,
        digits = 3L,
        ci = 0.95,
        fixedNames = character(),
        bsvFixedNames = res$bsvFixedNames
      )
    expect_equal(fmt[["BSV(CV%)"]], c("", "59.1", "32.4", ""))
  })

  test_that("VarCorr() reports the omega correlations of a non-nlme fit (issue 1140)", {
    one.compartment.block <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        eta.ka ~ 0.6
        eta.cl + eta.v ~ c(0.3, 0.05, 0.1)
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    fit <- .nlmixr(
      one.compartment.block,
      theo_sd,
      "focei",
      control = foceiControl(print = 0, maxOuterIterations = 0, covMethod = "", calcTables = FALSE)
    )
    .vc <- VarCorr(fit)
    .eta <- c("eta.ka", "eta.cl", "eta.v")
    expect_equal(rownames(.vc), .eta)
    expect_equal(names(.vc), c("Variance", "StdDev", "eta.ka", "eta.cl"))
    expect_equal(.vc$Variance, unname(diag(fit$omega)))
    expect_equal(.vc$StdDev, sqrt(unname(diag(fit$omega))))
    .cor <- stats::cov2cor(fit$omega)
    expect_equal(.vc$eta.cl, c(NA, NA, .cor[["eta.v", "eta.cl"]]))
    expect_equal(.vc$eta.ka, c(NA, 0, 0))
  })

  test_that("VarCorr()'s omega table spans every level of an IOV omega (issue 1140)", {
    .nm <- c("eta.cl", "eta.v")
    .om <- list(
      id = matrix(c(0.1, 0.02, 0.02, 0.2), 2, 2, dimnames = list(.nm, .nm)),
      occ = matrix(0.05, 1, 1, dimnames = list("iov.cl", "iov.cl"))
    )
    .vc <- .varCorrOmega(.om)
    expect_equal(rownames(.vc), c(.nm, "iov.cl"))
    expect_equal(.vc$Variance, c(0.1, 0.2, 0.05))
    expect_equal(.vc$StdDev, sqrt(c(0.1, 0.2, 0.05)))
    expect_equal(.vc$eta.cl, c(NA, 0.02 / sqrt(0.02), 0))
    expect_equal(.vc$eta.v, c(NA, NA, 0))
    # a diagonal omega keeps the variance-only table
    .d <- diag(c(0.1, 0.2))
    dimnames(.d) <- list(.nm, .nm)
    expect_equal(.varCorrOmega(.d), data.frame(Variance = c(0.1, 0.2), StdDev = sqrt(c(0.1, 0.2)), row.names = .nm))
  })
})

test_that("formatMinWidth", {
  # Special values
  expect_equal(
    formatMinWidth(x = c(NA, 0, Inf, -Inf, NaN)),
    c("NA", "0", "Inf", "-Inf", "NaN")
  )
  # Rounding occurs to requested significant digits
  expect_equal(
    formatMinWidth(x = -123456 * 10^(-10:10)),
    c(
      "-1.23e-5",
      "-1.23e-4",
      "-0.00123",
      "-0.0123",
      "-0.123",
      "-1.23",
      "-12.3",
      "-123",
      "-1230",
      "-12300",
      "-123000",
      "-1.23e6",
      "-1.23e7",
      "-1.23e8",
      "-1.23e9",
      "-1.23e10",
      "-1.23e11",
      "-1.23e12",
      "-1.23e13",
      "-1.23e14",
      "-1.23e15"
    )
  )
  # Rounding up works as expected; scientific notation values drop extraneous
  # zeros in the exponent
  expect_equal(
    formatMinWidth(x = -9999 * 10^(-10:10)),
    c(
      "-1.00e-6",
      "-1.00e-5",
      "-1.00e-4",
      "-0.00100",
      "-0.0100",
      "-0.100",
      "-1.00",
      "-10.0",
      "-100",
      "-1000",
      "-10000",
      "-100000",
      "-1.00e6",
      "-1.00e7",
      "-1.00e8",
      "-1.00e9",
      "-1.00e10",
      "-1.00e11",
      "-1.00e12",
      "-1.00e13",
      "-1.00e14"
    )
  )
  # Planned significant digits are shown, including when digits are added
  expect_equal(
    formatMinWidth(x = 12 * 10^(-10:10)),
    c(
      "1.20e-9",
      "1.20e-8",
      "1.20e-7",
      "1.20e-6",
      "1.20e-5",
      "1.20e-4",
      "0.00120",
      "0.0120",
      "0.120",
      "1.20",
      "12.0",
      "120",
      "1200",
      "12000",
      "120000",
      "1.20e6",
      "1.20e7",
      "1.20e8",
      "1.20e9",
      "1.20e10",
      "1.20e11"
    )
  )
  # Negative values
  expect_equal(
    formatMinWidth(x = -12 * 10^(-10:10)),
    c(
      "-1.20e-9",
      "-1.20e-8",
      "-1.20e-7",
      "-1.20e-6",
      "-1.20e-5",
      "-1.20e-4",
      "-0.00120",
      "-0.0120",
      "-0.120",
      "-1.20",
      "-12.0",
      "-120",
      "-1200",
      "-12000",
      "-120000",
      "-1.20e6",
      "-1.20e7",
      "-1.20e8",
      "-1.20e9",
      "-1.20e10",
      "-1.20e11"
    )
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
  fit <- .nlmixr(one.compartment, theo_sd, est = "focei", control = foceiControlFast)
  ## $parFixed formats with control$sigdigTable (see .updateParFixedApplySig in
  ## R/nlmixr2output.R), which defaults to foceiControl(sigdig=).  That default
  ## moved 4 -> 3 in 7d3c7b62d and this expectation kept its hardcoded 4, so the
  ## file failed on every run.  Read it off the fit: the point of the test is
  ## that $parFixed agrees with formatMinWidth at the fit's OWN digits, not that
  ## the digits take any particular value.
  .digits <- fit$control$sigdigTable
  expect_equal(
    fit$parFixed,
    structure(
      list(
        Est. = formatMinWidth(fit$parFixedDf$Estimate, digits = .digits),
        SE = formatMinWidth(fit$parFixedDf$SE, digits = .digits, naValue = ""),
        `%RSE` = formatMinWidth(fit$parFixedDf$`%RSE`, digits = .digits, naValue = ""),
        `Back-transformed(95%CI)` = c(
          sprintf(
            "%s (%s, %s)",
            formatMinWidth(fit$parFixedDf$`Back-transformed`, digits = .digits),
            formatMinWidth(fit$parFixedDf$`CI Lower`, digits = .digits),
            formatMinWidth(fit$parFixedDf$`CI Upper`, digits = .digits)
          )[!is.na(fit$parFixedDf$`CI Upper`)],
          formatMinWidth(fit$parFixedDf$`Back-transformed`[is.na(fit$parFixedDf$`CI Upper`)], digits = .digits)
        ),
        `BSV(CV%)` = formatMinWidth(fit$parFixedDf$`BSV(CV%)`, digits = .digits, naValue = ""),
        `Shrink(SD)%` = paste0(
          formatMinWidth(fit$parFixedDf$`Shrink(SD)%`, digits = .digits, naValue = ""),
          c("", ">", "", "")
        )
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

  fitFixed <- .nlmixr(one.compartment.fixed, theo_sd, est = "focei", control = foceiControlFast)
  expect_equal(
    fitFixed$parFixed,
    structure(
      list(
        Est. = formatMinWidth(fitFixed$parFixedDf$Estimate, digits = .digits),
        SE = c("FIXED", formatMinWidth(fitFixed$parFixedDf$SE[2:4], digits = .digits)),
        `%RSE` = c("FIXED", formatMinWidth(fitFixed$parFixedDf$`%RSE`[2:4], digits = .digits)),
        `Back-transformed(95%CI)` = c(
          formatMinWidth(fitFixed$parFixedDf$`Back-transformed`[1], digits = .digits),
          sprintf(
            "%s (%s, %s)",
            formatMinWidth(fitFixed$parFixedDf$`Back-transformed`, digits = .digits),
            formatMinWidth(fitFixed$parFixedDf$`CI Lower`, digits = .digits),
            formatMinWidth(fitFixed$parFixedDf$`CI Upper`, digits = .digits)
          )[2:4]
        )
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

  fitFixedLabel <- .nlmixr(one.compartment.labeled, theo_sd, est = "focei", control = list(print = 0))
  expect_equal(
    fitFixedLabel$parFixed,
    structure(
      list(
        Parameter = c("ka", "clearance", "", ""),
        Est. = formatMinWidth(fitFixedLabel$parFixedDf$Estimate, digits = .digits),
        SE = c("FIXED", formatMinWidth(fitFixedLabel$parFixedDf$SE[2:4], digits = .digits)),
        `%RSE` = c("FIXED", formatMinWidth(fitFixedLabel$parFixedDf$`%RSE`[2:4], digits = .digits)),
        `Back-transformed(95%CI)` = c(
          formatMinWidth(fitFixedLabel$parFixedDf$`Back-transformed`[1], digits = .digits),
          sprintf(
            "%s (%s, %s)",
            formatMinWidth(fitFixedLabel$parFixedDf$`Back-transformed`, digits = .digits),
            formatMinWidth(fitFixedLabel$parFixedDf$`CI Lower`, digits = .digits),
            formatMinWidth(fitFixedLabel$parFixedDf$`CI Upper`, digits = .digits)
          )[2:4]
        )
      ),
      class = c("nlmixr2ParFixed", "data.frame"),
      row.names = c("tka", "tcl", "tv", "add.sd")
    )
  )

  # Works with .ret$control$ci and .ret$control$sigdig ----
  fitFixedLabelCI <- .nlmixr(
    one.compartment.labeled,
    theo_sd,
    est = "focei",
    control = list(print = 0, ci = 0.9, sigdig = 4)
  )
  # digits from THIS fit, not the one at the top of the file: fitFixedLabelCI is
  # built with its own sigdig, so formatting the expectation with the other fit's
  # sigdigTable compares a table rendered at one precision against one checked at
  # another.  Same principle as the .digits comment above -- read it off the fit.
  .digitsCI <- fitFixedLabelCI$control$sigdigTable
  expect_equal(
    fitFixedLabelCI$parFixed,
    structure(
      list(
        Parameter = c("ka", "clearance", "", ""),
        Est. = formatMinWidth(fitFixedLabelCI$parFixedDf$Estimate, digits = .digitsCI),
        SE = c("FIXED", formatMinWidth(fitFixedLabelCI$parFixedDf$SE[2:4], digits = .digitsCI)),
        `%RSE` = c("FIXED", formatMinWidth(fitFixedLabelCI$parFixedDf$`%RSE`[2:4], digits = .digitsCI)),
        `Back-transformed(90%CI)` = c(
          formatMinWidth(fitFixedLabelCI$parFixedDf$`Back-transformed`[1], digits = .digitsCI),
          sprintf(
            "%s (%s, %s)",
            formatMinWidth(fitFixedLabelCI$parFixedDf$`Back-transformed`, digits = .digitsCI),
            formatMinWidth(fitFixedLabelCI$parFixedDf$`CI Lower`, digits = .digitsCI),
            formatMinWidth(fitFixedLabelCI$parFixedDf$`CI Upper`, digits = .digitsCI)
          )[2:4]
        )
      ),
      class = c("nlmixr2ParFixed", "data.frame"),
      row.names = c("tka", "tcl", "tv", "add.sd")
    )
  )

  # The FIXED marker must appear even when the fixed parameter is not
  # literally fixed out of the model (literalFix = FALSE), where the fixed
  # theta remains in $popDf with an NA standard error.
  fitNoLiteralFix <-
    .nlmixr(
      one.compartment.fixed,
      theo_sd,
      est = "focei",
      control = foceiControl(
        print = 0,
        maxInnerIterations = 1,
        maxOuterIterations = 1,
        eval.max = 1,
        literalFix = FALSE
      )
    )
  expect_equal(unname(fitNoLiteralFix$parFixed["tka", "SE"]), "FIXED")
  expect_equal(unname(fitNoLiteralFix$parFixed["tka", "%RSE"]), "FIXED")
  expect_equal(
    unname(fitNoLiteralFix$parFixed["tcl", "SE"]),
    formatMinWidth(fitNoLiteralFix$parFixedDf["tcl", "SE"], digits = .digits)
  )

  # nlmixr2est#355: without etas the BSV/shrinkage columns are always blank
  expect_false(any(startsWith(names(fitFixed$parFixed), "BSV(")))
  expect_false("Shrink(SD)%" %in% names(fitFixed$parFixed))
  expect_false(any(startsWith(names(fitFixed$parFixedDf), "BSV(")))
  expect_false("Shrink(SD)%" %in% names(fitFixed$parFixedDf))
  # ...but a model with etas keeps them
  expect_true(any(startsWith(names(fit$parFixed), "BSV(")))
  expect_true("Shrink(SD)%" %in% names(fit$parFixed))
})

test_that("the output builder writes a method's estimates into the ui and keeps the original model", {
  .ui <- rxode2::rxode2(function() {
    ini({
      tka <- 0.45
      tcl <- 1
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      y <- ka + cl
      y ~ add(add.sd)
    })
  })
  .om <- matrix(c(0.2, 0.05, 0.05, 0.1), 2, 2, dimnames = list(c("eta.ka", "eta.cl"), c("eta.ka", "eta.cl")))
  .old <- nlmixr2global$nlmixr2EstEnv$iniDf0
  .depth <- nlmixr2global$nlmixr2EstEnv$estDepth
  withr::defer({
    nlmixr2global$nlmixr2EstEnv$iniDf0 <- .old
    nlmixr2global$nlmixr2EstEnv$estDepth <- .depth
  })
  # called outside an nlmixr2() run: the ui passed in is the original model
  nlmixr2global$nlmixr2EstEnv$estDepth <- 0L
  nlmixr2global$nlmixr2EstEnv$iniDf0 <- data.frame(name = "stale")
  .env <- new.env(parent = emptyenv())
  .env$fullTheta <- c(tka = 0.5, tcl = 1.2, add.sd = 0.4)
  .env$omega <- .om
  .u <- .nlmixr2OutputContract(rxode2::rxUiDecompress(.ui), .env)
  expect_equal(.u$omega, .om)
  expect_equal(.u$theta[c("tka", "tcl", "add.sd")], .env$fullTheta)
  expect_equal(.env$iniDf0$est, .ui$iniDf$est)
  # a method that already updated its ui is left alone
  .env2 <- new.env(parent = emptyenv())
  .env2$iniDf0 <- data.frame(.ui$iniDf)
  .env2$omega <- .om
  .u2 <- .nlmixr2OutputContract(rxode2::rxUiDecompress(.ui), .env2)
  expect_equal(.u2$omega, .ui$omega)
  expect_equal(nlmixr2global$nlmixr2EstEnv$iniDf0, .env2$iniDf0)
})

# A 2x2 block repeated (eta.c, eta.d ~ same()), its names in either order
.sameBlockAB <- function() {
  ini({
    ta <- 1
    tb <- 2
    tc <- 3
    td <- 4
    add.sd <- 0.7
    eta.a + eta.b ~ c(0.3, 0.01, 0.1)
    eta.c + eta.d ~ same()
  })
  model({
    a <- exp(ta + eta.a)
    b <- exp(tb + eta.b)
    cc <- exp(tc + eta.c)
    d <- exp(td + eta.d)
    y <- a + b + cc + d
    y ~ add(add.sd)
  })
}

.sameBlockBA <- function() {
  ini({
    ta <- 1
    tb <- 2
    tc <- 3
    td <- 4
    add.sd <- 0.7
    eta.b + eta.a ~ c(0.1, 0.01, 0.3)
    eta.d + eta.c ~ same()
  })
  model({
    a <- exp(ta + eta.a)
    b <- exp(tb + eta.b)
    cc <- exp(tc + eta.c)
    d <- exp(td + eta.d)
    y <- a + b + cc + d
    y ~ add(add.sd)
  })
}

test_that("the output builder keeps every row of a repeated (same()) block", {
  .old <- nlmixr2global$nlmixr2EstEnv$iniDf0
  .depth <- nlmixr2global$nlmixr2EstEnv$estDepth
  withr::defer({
    nlmixr2global$nlmixr2EstEnv$iniDf0 <- .old
    nlmixr2global$nlmixr2EstEnv$estDepth <- .depth
  })
  for (.mod in list(.sameBlockAB, .sameBlockBA)) {
    .ui <- rxode2::rxode2(.mod)
    # an estimate: the block's values scaled, the copy still equal to it
    .om <- unclass(.ui$omega) * 1.5
    attributes(.om) <- list(dim = dim(.om), dimnames = dimnames(.ui$omega))
    nlmixr2global$nlmixr2EstEnv$estDepth <- 0L
    nlmixr2global$nlmixr2EstEnv$iniDf0 <- NULL
    .env <- new.env(parent = emptyenv())
    .env$fullTheta <- c(ta = 1.1, tb = 2.1, tc = 3.1, td = 4.1, add.sd = 0.5)
    .env$omega <- .om
    .u <- .nlmixr2OutputContract(rxode2::rxUiDecompress(.ui), .env)
    .label <- paste(rownames(.om), collapse = ",")
    # every eta row, the covariance rows included, keeps the model's condition
    .want <- .ui$iniDf[!is.na(.ui$iniDf$neta1), c("name", "condition")]
    .got <- .u$iniDf[!is.na(.u$iniDf$neta1), c("name", "condition")]
    expect_setequal(.got$name, .want$name)
    expect_identical(.got$condition[match(.want$name, .got$name)], .want$condition, label = .label)
    expect_identical(.u$omegaSameMap, .ui$omegaSameMap, label = .label)
    expect_equal(unclass(.u$omega)[rownames(.om), colnames(.om)], .om, ignore_attr = TRUE, label = .label)
  }
})

# Models for the BSV correlation rows; the ones not ending in "Fit" only build a ui
.bsvCor2 <- function() {
  ini({
    tcl <- 1
    tv <- 3.45
    eta.cl + eta.v ~ c(0.3, 0.05, 0.1)
    add.sd <- 0.7
  })
  model({
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    y <- cl + v
    y ~ add(add.sd)
  })
}

.bsvCor3 <- function() {
  ini({
    tka <- 0.45
    tcl <- 1
    tv <- 3.45
    eta.ka + eta.cl + eta.v ~ c(0.6, 0.01, 0.3, 0.02, 0.03, 0.1)
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    y <- ka + cl + v
    y ~ add(add.sd)
  })
}

.bsvCor4 <- function() {
  ini({
    ta <- 1
    tb <- 2
    tc <- 3
    td <- 4
    eta.a + eta.b + eta.c + eta.d ~ c(0.4, 0.01, 0.3, 0.02, 0.03, 0.2, 0.04, 0.05, 0.06, 0.1)
    add.sd <- 0.7
  })
  model({
    a <- exp(ta + eta.a)
    b <- exp(tb + eta.b)
    cc <- exp(tc + eta.c)
    d <- exp(td + eta.d)
    y <- a + b + cc + d
    y ~ add(add.sd)
  })
}

.bsvCorTwoBlocks <- function() {
  ini({
    tka <- 0.45
    tcl <- 1
    tv <- 3.45
    tq <- 1.2
    eta.ka + eta.cl ~ c(0.6, 0.01, 0.3)
    eta.v + eta.q ~ c(0.1, 0.02, 0.2)
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    q <- exp(tq + eta.q)
    y <- ka + cl + v + q
    y ~ add(add.sd)
  })
}

# eta.ka is not mu-referenced; tka carries a label
.bsvCorMixed <- function() {
  ini({
    tka <- 0.45
    label("Ka")
    tcl <- 1
    tv <- 3.45
    eta.ka + eta.cl ~ c(0.6, 0.01, 0.3)
    eta.v ~ 0.1
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka) * exp(eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    y <- ka + cl + v
    y ~ add(add.sd)
  })
}

.bsvCorFix <- function() {
  ini({
    tcl <- 1
    tv <- 3.45
    eta.cl + eta.v ~ fix(0.3, 0.05, 0.1)
    add.sd <- 0.7
  })
  model({
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    y <- cl + v
    y ~ add(add.sd)
  })
}

.bsvCorIov <- function() {
  ini({
    tcl <- 1
    tv <- 3.45
    add.sd <- 0.7
    eta.cl + eta.v ~ c(0.3, 0.05, 0.1)
    iov.cl + iov.v ~ c(0.1, 0.01, 0.1) | occ
  })
  model({
    cl <- exp(tcl + eta.cl + iov.cl)
    v <- exp(tv + eta.v + iov.v)
    y <- cl + v
    y ~ add(add.sd)
  })
}

.bsvCorBlockFit <- function() {
  ini({
    tka <- 0.45
    tcl <- 1
    tv <- 3.45
    eta.ka ~ 0.6
    eta.cl + eta.v ~ c(0.3, 0.05, 0.1)
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    linCmt() ~ add(add.sd)
  })
}

.bsvCorMixedFit <- function() {
  ini({
    tka <- 0.45
    label("Ka")
    tcl <- 1
    tv <- 3.45
    eta.ka + eta.cl ~ c(0.6, 0.01, 0.3)
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

.bsvCorIovFit <- function() {
  ini({
    tka <- 0.45
    tcl <- 1
    tv <- 3.45
    add.sd <- 0.7
    eta.cl + eta.v ~ c(0.3, 0.05, 0.1)
    iov.cl + iov.v ~ c(0.1, 0.01, 0.1) | occ
  })
  model({
    ka <- exp(tka)
    cl <- exp(tcl + eta.cl + iov.cl)
    v <- exp(tv + eta.v + iov.v)
    linCmt() ~ add(add.sd)
  })
}

.bsvCorOmega <- function(m) {
  .n <- c("eta.cl", "eta.v")
  matrix(m, 2, 2, dimnames = list(.n, .n))
}

nmTest({
  test_that(".updateParFixedBsvCor reports a 2x2 block's correlation from the fit's omega", {
    .ui <- rxode2::rxode2(.bsvCor2)
    # not the ini() values, so the correlation is read off the omega given
    .om <- .bsvCorOmega(c(0.2, -0.03, -0.03, 0.4))
    .r <- .updateParFixedBsvCor(.ui$iniDf, .om, .ui$muRefDataFrame)
    expect_identical(.r$name, "bsvcorr(tcl,tv)")
    expect_identical(.r$label, "BSV correlation between tcl and tv")
    expect_equal(.r$cor, stats::cov2cor(.om)[["eta.cl", "eta.v"]], tolerance = 1e-12)
    expect_equal(.r$cor, -0.03 / sqrt(0.2 * 0.4), tolerance = 1e-12)
    expect_identical(.r$fix, FALSE)
    # the row is found by its name
    expect_identical(
      .parFixedIsBsvCor(c("tcl", "eta.cl", "bsvcorr(tcl,tv)", "BSV correlation between tcl and tv")),
      c(FALSE, FALSE, TRUE, FALSE)
    )
  })

  test_that(".updateParFixedBsvCor orders the pairs of a block by the first eta, then the second", {
    .ui <- rxode2::rxode2(.bsvCor3)
    .om <- .ui$omega
    .cor <- stats::cov2cor(.om)
    .r <- .updateParFixedBsvCor(.ui$iniDf, .om, .ui$muRefDataFrame)
    expect_identical(
      .r$name,
      c(
        "bsvcorr(tka,tcl)",
        "bsvcorr(tka,tv)",
        "bsvcorr(tcl,tv)"
      )
    )
    expect_equal(
      .r$cor,
      c(.cor[["eta.ka", "eta.cl"]], .cor[["eta.ka", "eta.v"]], .cor[["eta.cl", "eta.v"]]),
      tolerance = 1e-12
    )
    expect_identical(.r$fix, rep(FALSE, 3L))
    # ini() lists a block's covariances row by row (b-a, c-a, c-b, d-a, ...);
    # the table groups them by the first eta instead
    .ui4 <- rxode2::rxode2(.bsvCor4)
    .r4 <- .updateParFixedBsvCor(.ui4$iniDf, .ui4$omega, .ui4$muRefDataFrame)
    expect_identical(
      .r4$name,
      c("bsvcorr(ta,tb)", "bsvcorr(ta,tc)", "bsvcorr(ta,td)", "bsvcorr(tb,tc)", "bsvcorr(tb,td)", "bsvcorr(tc,td)")
    )
    .cor4 <- stats::cov2cor(.ui4$omega)
    expect_equal(
      .r4$cor,
      c(
        .cor4[["eta.a", "eta.b"]],
        .cor4[["eta.a", "eta.c"]],
        .cor4[["eta.a", "eta.d"]],
        .cor4[["eta.b", "eta.c"]],
        .cor4[["eta.b", "eta.d"]],
        .cor4[["eta.c", "eta.d"]]
      ),
      tolerance = 1e-12
    )
  })

  test_that(".updateParFixedBsvCor gives no row to the zeros between separate blocks", {
    .ui <- rxode2::rxode2(.bsvCorTwoBlocks)
    .om <- .ui$omega
    .cor <- stats::cov2cor(.om)
    .r <- .updateParFixedBsvCor(.ui$iniDf, .om, .ui$muRefDataFrame)
    expect_identical(.r$name, c("bsvcorr(tka,tcl)", "bsvcorr(tv,tq)"))
    expect_equal(.r$cor, c(.cor[["eta.ka", "eta.cl"]], .cor[["eta.v", "eta.q"]]), tolerance = 1e-12)
    # the etas of one block left out, as the IOV rewrite's occasion etas are
    .r <- .updateParFixedBsvCor(.ui$iniDf, .om, .ui$muRefDataFrame, exclude = c("eta.v", "eta.q"))
    expect_identical(.r$name, "bsvcorr(tka,tcl)")
    expect_equal(.r$cor, .cor[["eta.ka", "eta.cl"]], tolerance = 1e-12)
  })

  test_that(".updateParFixedBsvCor has no rows without a covariance", {
    .none <- data.frame(name = character(0), label = character(0), cor = numeric(0), fix = logical(0))
    # diagonal omega
    .ui <- rxode2::rxode2(one.compartment)
    expect_identical(.updateParFixedBsvCor(.ui$iniDf, .ui$omega, .ui$muRefDataFrame), .none)
    # no omega
    .ui <- rxode2::rxode2(.bsvCor2)
    expect_identical(.updateParFixedBsvCor(.ui$iniDf, NULL, .ui$muRefDataFrame), .none)
  })

  test_that(".updateParFixedBsvCor names a non-mu-referenced eta by itself", {
    .ui <- rxode2::rxode2(.bsvCorMixed)
    expect_false("eta.ka" %in% .ui$muRefDataFrame$eta)
    .om <- .ui$omega
    .r <- .updateParFixedBsvCor(.ui$iniDf, .om, .ui$muRefDataFrame)
    expect_identical(.r$name, "bsvcorr(eta.ka,tcl)")
    expect_equal(.r$cor, stats::cov2cor(.om)[["eta.ka", "eta.cl"]], tolerance = 1e-12)
  })

  test_that(".updateParFixedBsvCor flags a fixed covariance and keeps a zero estimate", {
    .ui <- rxode2::rxode2(.bsvCorFix)
    .r <- .updateParFixedBsvCor(.ui$iniDf, .ui$omega, .ui$muRefDataFrame)
    expect_identical(.r$name, "bsvcorr(tcl,tv)")
    expect_identical(.r$fix, TRUE)
    expect_equal(.r$cor, 0.05 / sqrt(0.3 * 0.1), tolerance = 1e-12)
    # a covariance estimated at exactly zero inside its block keeps its row
    .ui <- rxode2::rxode2(.bsvCor2)
    .r <- .updateParFixedBsvCor(.ui$iniDf, .bsvCorOmega(c(0.2, 0, 0, 0.4)), .ui$muRefDataFrame)
    expect_identical(.r$name, "bsvcorr(tcl,tv)")
    expect_identical(.r$cor, 0)
    .popDf <- data.frame(Estimate = c(1, 3.45, 0.7), row.names = c("tcl", "tv", "add.sd"))
    .res <- .updateParFixedAddBsv(
      .popDf,
      .ui$iniDf,
      .bsvCorOmega(c(0.2, 0, 0, 0.4)),
      3L,
      .ui$muRefDataFrame,
      .ui$muRefCurEval
    )
    .fmt <- .updateParFixedApplySig(.res$popDf, 3L, 0.95, character(), .res$bsvFixedNames)
    expect_identical(.fmt["bsvcorr(tcl,tv)", "BSV(CV%)"], "0")
  })

  test_that(".updateParFixedBsvCor leaves the correlation of an eta absent from the omega empty", {
    .ui <- rxode2::rxode2(.bsvCorFix)
    .om <- matrix(0.3, 1, 1, dimnames = list("eta.cl", "eta.cl"))
    .r <- .updateParFixedBsvCor(.ui$iniDf, .om, .ui$muRefDataFrame)
    expect_identical(.r$name, "bsvcorr(tcl,tv)")
    expect_identical(.r$cor, NA_real_)
    expect_identical(.r$fix, TRUE)
    # an empty cell is not wrapped as fix()
    .popDf <- data.frame(Estimate = c(1, 3.45, 0.7), row.names = c("tcl", "tv", "add.sd"))
    .res <- .updateParFixedAddBsv(.popDf, .ui$iniDf, .om, 3L, .ui$muRefDataFrame, .ui$muRefCurEval)
    expect_false("bsvcorr(tcl,tv)" %in% .res$bsvFixedNames)
    .fmt <- .updateParFixedApplySig(.res$popDf, 3L, 0.95, character(), .res$bsvFixedNames)
    expect_identical(.fmt["bsvcorr(tcl,tv)", "BSV(CV%)"], "")
  })

  test_that(".updateParFixedBsvCor reports only the subject level of an IOV omega", {
    .ui <- rxode2::rxode2(.bsvCorIov)
    .om <- .ui$omega
    expect_true(is.list(.om))
    .r <- .updateParFixedBsvCor(.ui$iniDf, .om, .ui$muRefDataFrame)
    expect_identical(.r$name, "bsvcorr(tcl,tv)")
    expect_equal(.r$cor, stats::cov2cor(.om$id)[["eta.cl", "eta.v"]], tolerance = 1e-12)
  })

  test_that("the BSV correlation rows end the parameter table with their description and BSV cell", {
    .ui <- rxode2::rxode2(.bsvCorMixed)
    .popDf <- data.frame(Estimate = c(0.45, 1, 3.45, 0.7), row.names = c("tka", "tcl", "tv", "add.sd"))
    .popDf <- .updateParFixedAddParameterLabel(.popDf, .ui$iniDf)
    .n <- c("eta.ka", "eta.cl", "eta.v")
    .om <- matrix(c(0.5, 0.02, 0, 0.02, 0.2, 0, 0, 0, 0.1), 3, 3, dimnames = list(.n, .n))
    .res <- .updateParFixedAddBsv(.popDf, .ui$iniDf, .om, 3L, .ui$muRefDataFrame, .ui$muRefCurEval)
    .cor <- "bsvcorr(eta.ka,tcl)"
    # after the theta rows and the non-mu-referenced eta's row
    expect_identical(rownames(.res$popDf), c("tka", "tcl", "tv", "add.sd", "eta.ka", .cor))
    # the header is the one the variances alone choose
    .resNoCor <- .updateParFixedAddBsv(
      .popDf,
      .ui$iniDf,
      .om,
      3L,
      .ui$muRefDataFrame,
      .ui$muRefCurEval,
      exclude = c("eta.ka", "eta.cl")
    )
    expect_identical(rownames(.resNoCor$popDf), c("tka", "tcl", "tv", "add.sd", "eta.ka"))
    expect_identical(names(.res$popDf), names(.resNoCor$popDf))
    expect_identical(names(.res$popDf), c("Parameter", "Estimate", "BSV(CV%)"))
    # a model with labels describes the row in its label column
    expect_identical(.res$popDf[.cor, "Parameter"], "BSV correlation between eta.ka and tcl")
    expect_identical(.res$popDf[.cor, "Estimate"], NA_real_)
    expect_equal(.res$popDf[.cor, "BSV(CV%)"], stats::cov2cor(.om)[["eta.ka", "eta.cl"]], tolerance = 1e-12)
    expect_equal(.res$popDf[.cor, "BSV(CV%)"], 0.02 / sqrt(0.5 * 0.2), tolerance = 1e-12)
    .sh <- .updateParFixedAddShrinkage(
      .res$popDf,
      shrink = data.frame(eta.ka = 1:7, eta.cl = 11:17, eta.v = 21:27, IWRES = 31:37),
      ui = .ui
    )
    expect_identical(.sh[["Shrink(SD)%"]], c(NA, 17, 27, 37, 7, NA))
    .fmt <- .updateParFixedApplySig(.sh, 3L, 0.95, character(), .res$bsvFixedNames)
    expect_identical(
      unlist(.fmt[.cor, ]),
      c(Parameter = "BSV correlation between eta.ka and tcl", Est. = "", `BSV(CV%)` = "0.0632", `Shrink(SD)%` = "")
    )
    # the label column of the theta rows is untouched
    expect_identical(.fmt[["Parameter"]][1:4], c("Ka", "", "", ""))
  })

  test_that("a fixed covariance's correlation is printed as fix() and survives an SE refresh", {
    .ui <- rxode2::rxode2(.bsvCorFix)
    .popDf <- data.frame(
      Estimate = c(1, 3.45, 0.7),
      SE = c(0.1, 0.2, 0.05),
      `%RSE` = c(10, 0.2 / 3.45 * 100, 0.05 / 0.7 * 100),
      `Back-transformed` = c(exp(1), exp(3.45), 0.7),
      `CI Lower` = NA_real_,
      `CI Upper` = NA_real_,
      row.names = c("tcl", "tv", "add.sd"),
      check.names = FALSE
    )
    .res <- .updateParFixedAddBsv(.popDf, .ui$iniDf, .ui$omega, 3L, .ui$muRefDataFrame, .ui$muRefCurEval)
    .cor <- "bsvcorr(tcl,tv)"
    expect_identical(.res$bsvFixedNames, c("tcl", "tv", .cor))
    .fmt <- .updateParFixedApplySig(.res$popDf, 3L, 0.95, character(), .res$bsvFixedNames)
    .want <- c(
      Est. = "",
      SE = "",
      `%RSE` = "",
      `Back-transformed(95%CI)` = "",
      `BSV(CV%)` = "fix(0.289)"
    )
    expect_identical(unlist(.fmt[.cor, ]), .want)
    .env <- new.env(parent = emptyenv())
    .env$ui <- .ui
    .env$control <- list(ci = 0.95, sigdigTable = 3L)
    .env$parFixedDf <- .res$popDf
    class(.fmt) <- c("nlmixr2ParFixed", "data.frame")
    .env$parFixed <- .fmt
    .cov <- diag(c(0.01, 0.04))
    dimnames(.cov) <- list(c("tcl", "tv"), c("tcl", "tv"))
    .updateParFixedRefreshSeFromCov(.env, .cov)
    expect_identical(rownames(.env$parFixedDf), c("tcl", "tv", "add.sd", .cor))
    expect_equal(.env$parFixedDf[.cor, "BSV(CV%)"], 0.05 / sqrt(0.3 * 0.1), tolerance = 1e-12)
    expect_identical(unlist(.env$parFixed[.cor, ]), .want)
    expect_identical(.env$parFixed[["SE"]][1:2], c("0.100", "0.200"))
  })

  test_that("broom leaves the BSV correlation rows out of the fixed-effect table", {
    .df <- data.frame(
      Estimate = c(1, 3.45, NA),
      SE = c(0.1, 0.2, NA),
      `Back-transformed` = c(exp(1), exp(3.45), NA),
      `BSV(CV%)` = c(59.1, 32.4, 0.289),
      check.names = FALSE,
      row.names = c("tcl", "tv", "bsvcorr(tcl,tv)")
    )
    .got <- .broomParFixedDf(list(parFixedDf = .df))
    expect_identical(rownames(.got), c("tcl", "tv"))
    expect_identical(names(.got), c("model.est", "std.error", "estimate", "bsv"))
    expect_identical(.got$model.est, c(1, 3.45))
  })

  test_that("a fit reports its estimated BSV correlations in $parFixed and $parFixedDf", {
    skip_if_not_installed("broom.mixed")
    .ctl <- foceiControl(print = 0, maxOuterIterations = 0, calcTables = FALSE)
    fit <- .nlmixr(.bsvCorBlockFit, theo_sd, "focei", control = .ctl)
    .cor <- "bsvcorr(tcl,tv)"
    .rows <- c("tka", "tcl", "tv", "add.sd", .cor)
    .cols <- c("Estimate", "SE", "%RSE", "Back-transformed", "CI Lower", "CI Upper", "BSV(CV%)", "Shrink(SD)%")
    expect_identical(rownames(fit$parFixedDf), .rows)
    # no labels, so no label column: the row name alone names the row
    expect_identical(names(fit$parFixedDf), .cols)
    expect_false("Parameter" %in% names(fit$parFixed))
    expect_equal(
      fit$parFixedDf[.cor, "BSV(CV%)"],
      stats::cov2cor(fit$omega)[["eta.cl", "eta.v"]],
      tolerance = 1e-12
    )
    for (.c in setdiff(.cols, "BSV(CV%)")) {
      expect_identical(unname(fit$parFixedDf[.cor, .c]), NA_real_, label = .c)
      expect_false(anyNA(names(fit$parFixedDf[[.c]])), label = .c)
    }
    expect_identical(rownames(fit$parFixed), .rows)
    expect_identical(
      names(fit$parFixed),
      c("Est.", "SE", "%RSE", "Back-transformed(95%CI)", "BSV(CV%)", "Shrink(SD)%")
    )
    # maxOuterIterations = 0 keeps the omega at its ini() values
    expect_identical(
      unlist(fit$parFixed[.cor, ]),
      c(Est. = "", SE = "", `%RSE` = "", `Back-transformed(95%CI)` = "", `BSV(CV%)` = "0.289", `Shrink(SD)%` = "")
    )
    expect_identical(
      fit$parFixed[.cor, "BSV(CV%)"],
      formatMinWidth(stats::cov2cor(fit$omega)[["eta.cl", "eta.v"]], digits = fit$control$sigdigTable)
    )
    # broom keeps the correlation in the random-effect terms only
    .td <- as.data.frame(broom.mixed::tidy(fit))
    expect_identical(.td$term[.td$effect == "fixed"], c("tka", "tcl", "tv"))
    expect_identical(
      .td$term[.td$effect == "ran_pars"],
      c("sd__eta.ka", "sd__eta.cl", "sd__eta.v", "cor__eta.v.eta.cl", "add.sd")
    )
    expect_identical(rownames(confint(fit)), c("tka", "tcl", "tv", "add.sd"))
  })

  test_that("a fit puts its non-mu-referenced eta's row before the BSV correlation rows", {
    .ctl <- foceiControl(print = 0, maxOuterIterations = 0, calcTables = FALSE, covMethod = "")
    fit <- .nlmixr(.bsvCorMixedFit, theo_sd, "focei", control = .ctl)
    .cor <- "bsvcorr(eta.ka,tcl)"
    expect_identical(rownames(fit$parFixedDf), c("tka", "tcl", "tv", "add.sd", "eta.ka", .cor))
    # the model has labels, so the label column describes the row
    .desc <- "BSV correlation between eta.ka and tcl"
    expect_identical(unname(fit$parFixedDf[.cor, "Parameter"]), .desc)
    expect_equal(
      fit$parFixedDf[.cor, "BSV(CV%)"],
      stats::cov2cor(fit$omega)[["eta.ka", "eta.cl"]],
      tolerance = 1e-12
    )
    expect_identical(fit$parFixed[.cor, "Parameter"], .desc)
    expect_identical(fit$parFixed[["Parameter"]][1:4], c("Ka", "", "", ""))
    expect_identical(fit$parFixed[.cor, "BSV(CV%)"], "0.0236")
  })

  test_that("an IOV fit reports only its subject-level BSV correlations", {
    .d <- theo_sd
    .d$occ <- ifelse(.d$TIME < 6, 1L, 2L)
    .ctl <- foceiControl(print = 0, maxOuterIterations = 0, calcTables = FALSE, covMethod = "", iovMethod = "omega")
    fit <- .nlmixr(.bsvCorIovFit, .d, "focei", control = .ctl)
    .cor <- "bsvcorr(tcl,tv)"
    expect_identical(rownames(fit$parFixedDf), c("tka", "tcl", "tv", "add.sd", "iov.cl", "iov.v", .cor))
    expect_identical(rownames(fit$parFixed), rownames(fit$parFixedDf))
    expect_equal(
      fit$parFixedDf[.cor, "BSV(CV% or SD)"],
      stats::cov2cor(fit$omega$id)[["eta.cl", "eta.v"]],
      tolerance = 1e-12
    )
    expect_identical(fit$parFixed[.cor, "BSV(CV% or SD)"], "0.289")
  })
})

# tcl is mu-referenced by two etas, so neither eta is named by it
.bsvCorSharedFit <- function() {
  ini({
    tcl <- 1
    tv <- 3.45
    eta.cl + eta.cl2 + eta.v ~ c(0.3, 0.02, 0.2, 0.03, 0.04, 0.1)
    add.sd <- 0.7
  })
  model({
    cl <- exp(tcl + eta.cl)
    ka <- exp(tcl + eta.cl2)
    v <- exp(tv + eta.v)
    linCmt() ~ add(add.sd)
  })
}

.bsvCorSameFit <- function() {
  ini({
    tka <- 0.45
    tcl <- 1
    tv <- 3.45
    tq <- 1
    tvp <- 3
    add.sd <- 0.7
    eta.cl + eta.v ~ c(0.3, 0.05, 0.1)
    eta.q + eta.vp ~ same()
  })
  model({
    ka <- exp(tka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    q <- exp(tq + eta.q)
    vp <- exp(tvp + eta.vp)
    linCmt() ~ add(add.sd)
  })
}

nmTest({
  test_that(".updateParFixedBsvCorLabel names an eta by its theta only when each has one mu-reference", {
    .mu <- data.frame(
      theta = c("tcl", "tcl", "tv", "tka", "tq"),
      eta = c("eta.cl", "eta.cl2", "eta.v", "eta.ka", "eta.ka"),
      level = "id"
    )
    expect_identical(
      .updateParFixedBsvCorLabel(c("eta.cl", "eta.cl2", "eta.v", "eta.ka", "eta.x"), .mu),
      c("eta.cl", "eta.cl2", "tv", "eta.ka", "eta.x")
    )
    expect_identical(.updateParFixedBsvCorLabel(c("eta.cl", "eta.v"), NULL), c("eta.cl", "eta.v"))
  })

  test_that("the etas of a theta mu-referenced twice keep their own names in the correlation rows", {
    .ui <- rxode2::rxode2(.bsvCorSharedFit)
    expect_identical(as.character(.ui$muRefDataFrame$theta), c("tcl", "tcl", "tv"))
    .om <- .ui$omega
    .cor <- stats::cov2cor(.om)
    .r <- .updateParFixedBsvCor(.ui$iniDf, .om, .ui$muRefDataFrame)
    expect_identical(.r$name, c("bsvcorr(eta.cl,eta.cl2)", "bsvcorr(eta.cl,tv)", "bsvcorr(eta.cl2,tv)"))
    expect_identical(
      .r$label,
      paste0("BSV correlation between ", c("eta.cl and eta.cl2", "eta.cl and tv", "eta.cl2 and tv"))
    )
    expect_equal(
      .r$cor,
      c(.cor[["eta.cl", "eta.cl2"]], .cor[["eta.cl", "eta.v"]], .cor[["eta.cl2", "eta.v"]]),
      tolerance = 1e-12
    )
    .popDf <- data.frame(Estimate = c(1, 3.45, 0.7), row.names = c("tcl", "tv", "add.sd"))
    .res <- .updateParFixedAddBsv(.popDf, .ui$iniDf, .om, 3L, .ui$muRefDataFrame, .ui$muRefCurEval)
    expect_identical(rownames(.res$popDf), c("tcl", "tv", "add.sd", "eta.cl", "eta.cl2", .r$name))
  })

  test_that("an eta without variance leaves its correlation empty", {
    .om <- .bsvCorOmega(c(0.3, 0, 0, 0))
    .ui <- rxode2::rxode2(.bsvCor2)
    expect_no_warning(.r <- .updateParFixedBsvCor(.ui$iniDf, .om, .ui$muRefDataFrame))
    expect_identical(.r$cor, NA_real_)
    # a fixed covariance (its ini() value over a zero variance) is not wrapped as fix()
    .ui <- rxode2::rxode2(.bsvCorFix)
    expect_no_warning(.r <- .updateParFixedBsvCor(.ui$iniDf, .om, .ui$muRefDataFrame))
    expect_identical(.r$cor, NA_real_)
    expect_identical(.r$fix, TRUE)
    .popDf <- data.frame(Estimate = c(1, 3.45, 0.7), row.names = c("tcl", "tv", "add.sd"))
    .res <- .updateParFixedAddBsv(.popDf, .ui$iniDf, .om, 3L, .ui$muRefDataFrame, .ui$muRefCurEval)
    .fmt <- .updateParFixedApplySig(.res$popDf, 3L, 0.95, character(), .res$bsvFixedNames)
    expect_identical(.fmt["bsvcorr(tcl,tv)", "BSV(CV%)"], "")
  })

  test_that("a fixed covariance is its ini() value over the fit's variances", {
    .ui <- rxode2::rxode2(.bsvCorFix)
    .cor <- "bsvcorr(tcl,tv)"
    # the omega a fit returns for a covariance fixed at zero
    .om <- .bsvCorOmega(c(0.3, 1.7e-11, 1.7e-11, 0.1))
    .ini <- .ui$iniDf
    .ini$est[.ini$name == "(eta.cl,eta.v)"] <- 0
    .r <- .updateParFixedBsvCor(.ini, .om, .ui$muRefDataFrame)
    expect_identical(.r$cor, 0)
    .popDf <- data.frame(Estimate = c(1, 3.45, 0.7), row.names = c("tcl", "tv", "add.sd"))
    .res <- .updateParFixedAddBsv(
      .popDf,
      .ui$iniDf,
      .om,
      3L,
      .ui$muRefDataFrame,
      .ui$muRefCurEval,
      corIniDf = .ini
    )
    .fmt <- .updateParFixedApplySig(.res$popDf, 3L, 0.95, character(), .res$bsvFixedNames)
    expect_identical(.fmt[.cor, "BSV(CV%)"], "fix(0)")
    # nonzero: the fixed 0.05, over variances that moved
    .om <- .bsvCorOmega(c(0.2, 0.04, 0.04, 0.4))
    .r <- .updateParFixedBsvCor(.ui$iniDf, .om, .ui$muRefDataFrame)
    expect_equal(.r$cor, 0.05 / sqrt(0.2 * 0.4), tolerance = 1e-12)
    # an estimated covariance is the omega's
    .ui2 <- rxode2::rxode2(.bsvCor2)
    .r <- .updateParFixedBsvCor(.ui2$iniDf, .om, .ui2$muRefDataFrame)
    expect_equal(.r$cor, stats::cov2cor(.om)[["eta.cl", "eta.v"]], tolerance = 1e-12)
  })

  test_that("the covariances reported come from the model as estimated when the fit carries it", {
    .ui <- rxode2::rxode2(.bsvCor2)
    .env <- new.env(parent = emptyenv())
    expect_identical(.updateParFixedBsvCorIni(.env, .ui), .ui$iniDf)
    .ini0 <- .ui$iniDf
    .ini0$est[.ini0$name == "tcl"] <- 2
    .env$iniDf0 <- .ini0
    expect_identical(.updateParFixedBsvCorIni(.env, .ui), .ini0)
    # iniDf0 can be a compressed ui, which is not a frame to read
    .env$iniDf0 <- rxode2::rxUiCompress(.ui)
    expect_identical(.updateParFixedBsvCorIni(.env, .ui), .ui$iniDf)
    # a ui rebuilt from an omega estimated at zero has lost the covariance row;
    # the model as estimated still has it
    .rebuilt <- .ui$iniDf[.ui$iniDf$name != "(eta.cl,eta.v)", , drop = FALSE]
    .om <- .bsvCorOmega(c(0.2, 0, 0, 0.4))
    .popDf <- data.frame(Estimate = c(1, 3.45, 0.7), row.names = c("tcl", "tv", "add.sd"))
    .res <- .updateParFixedAddBsv(
      .popDf,
      .rebuilt,
      .om,
      3L,
      .ui$muRefDataFrame,
      .ui$muRefCurEval,
      corIniDf = .ui$iniDf
    )
    expect_identical(rownames(.res$popDf), c("tcl", "tv", "add.sd", "bsvcorr(tcl,tv)"))
    expect_identical(.res$popDf["bsvcorr(tcl,tv)", "BSV(CV%)"], 0)
  })

  test_that("a fit with a theta mu-referenced by two correlated etas reports every pair", {
    .ctl <- foceiControl(print = 0, maxOuterIterations = 0, calcTables = FALSE, covMethod = "")
    fit <- .nlmixr(.bsvCorSharedFit, theo_sd, "focei", control = .ctl)
    .rows <- c("bsvcorr(eta.cl,eta.cl2)", "bsvcorr(eta.cl,tv)", "bsvcorr(eta.cl2,tv)")
    expect_identical(rownames(fit$parFixedDf), c("tcl", "tv", "add.sd", "eta.cl", "eta.cl2", .rows))
    expect_identical(rownames(fit$parFixed), rownames(fit$parFixedDf))
    .cor <- stats::cov2cor(fit$omega)
    expect_equal(
      unname(fit$parFixedDf[.rows, "BSV(CV%)"]),
      c(.cor[["eta.cl", "eta.cl2"]], .cor[["eta.cl", "eta.v"]], .cor[["eta.cl2", "eta.v"]]),
      tolerance = 1e-12
    )
    expect_identical(unname(fit$parFixed[.rows, "BSV(CV%)"]), c("0.0816", "0.173", "0.283"))
  })

  test_that("a fit reports the copy of a repeated (same()) block as its own correlation", {
    .ctl <- foceiControl(print = 0, maxOuterIterations = 0, calcTables = FALSE, covMethod = "")
    fit <- .nlmixr(.bsvCorSameFit, theo_sd, "focei", control = .ctl)
    .rows <- c("bsvcorr(tcl,tv)", "bsvcorr(tq,tvp)")
    expect_identical(rownames(fit$parFixedDf), c("tka", "tcl", "tv", "tq", "tvp", "add.sd", .rows))
    .cor <- stats::cov2cor(fit$omega)
    expect_equal(
      unname(fit$parFixedDf[.rows, "BSV(CV%)"]),
      c(.cor[["eta.cl", "eta.v"]], .cor[["eta.q", "eta.vp"]]),
      tolerance = 1e-12
    )
    expect_identical(unname(fit$parFixed[.rows, "BSV(CV%)"]), c("0.289", "0.289"))
  })

  test_that("a saem fit reports its BSV correlation, even one estimated at zero", {
    .ctl <- saemControl(print = 0, nBurn = 1, nEm = 1, nmc = 1, nu = c(1, 1, 1), calcTables = FALSE)
    fit <- .nlmixr(.bsvCorBlockFit, theo_sd, "saem", control = .ctl)
    .row <- "bsvcorr(tcl,tv)"
    # two iterations leave the covariance at exactly zero; saem's ui, rebuilt
    # from that omega, has no row for it, so the row comes from iniDf0
    expect_identical(fit$omega[["eta.cl", "eta.v"]], 0)
    expect_false("(eta.cl,eta.v)" %in% fit$ui$iniDf$name)
    expect_identical(rownames(fit$parFixedDf), c("tka", "tcl", "tv", "add.sd", .row))
    expect_identical(rownames(fit$parFixed), rownames(fit$parFixedDf))
    .value <- stats::cov2cor(fit$omega)[["eta.cl", "eta.v"]]
    expect_equal(fit$parFixedDf[.row, "BSV(CV%)"], .value, tolerance = 1e-12)
    expect_identical(fit$parFixed[.row, "BSV(CV%)"], formatMinWidth(.value, digits = fit$control$sigdigTable))
    expect_identical(fit$parFixed[.row, "Est."], "")
  })

  test_that("the BSV correlation is printed to the fit's sigdigTable", {
    .ctl <- foceiControl(print = 0, maxOuterIterations = 0, calcTables = FALSE, covMethod = "", sigdigTable = 5)
    fit <- .nlmixr(.bsvCorBlockFit, theo_sd, "focei", control = .ctl)
    .row <- "bsvcorr(tcl,tv)"
    expect_equal(fit$control$sigdigTable, 5)
    expect_identical(fit$parFixed[.row, "BSV(CV%)"], "0.28868")
    expect_identical(fit$parFixed["tcl", "BSV(CV%)"], "59.149")
  })
})
