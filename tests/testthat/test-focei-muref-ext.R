nmTest({
  # Extended mu-referencing (.foceiAnalyticMuRefExtDirMap): a theta that acts only through an
  # eta-carrying individual parameter (dE/dtheta = c*dE/deta for every model expression,
  # c = dmu/dtheta free of etas, states and time) reuses the eta's sensitivities scaled by c,
  # whatever the (nonlinear) covariate model.

  .extMod <- function(clLine, extraIni = "", tcl = 1, err = "add(add.sd)", vLine = "v <- exp(tv + eta.v)") {
    .errIni <- if (grepl("^prop", err)) "prop.sd <- 0.1" else "add.sd <- 0.7"
    eval(parse(
      text = paste0(
        "function() { ini({ tka <- 0.45; tcl <- ",
        tcl,
        "; tv <- 3.45; ",
        extraIni,
        if (nzchar(extraIni)) "; " else "",
        "eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1; ",
        .errIni,
        " });",
        " model({ ka <- exp(tka + eta.ka); ",
        clLine,
        "; ",
        vLine,
        ";",
        " d/dt(depot) <- -ka * depot; d/dt(center) <- ka * depot - cl / v * center;",
        " cp <- center / v; cp ~ ",
        err,
        " }) }"
      )
    ))
  }
  .extUi <- function(f, constCovs = c("PMA", "WT")) {
    .ui <- rxode2::rxUiDecompress(f())
    rxode2::rxAssignControlValue(.ui, "foceiConstCovs", constCovs)
    .ui
  }
  # every eta plus every estimated structural theta that is not a plain mu-referenced typical value
  .extDirs <- function(ui) {
    .ini <- ui$iniDf
    .th <- .ini[!is.na(.ini$ntheta) & !.ini$fix & is.na(.ini$err), ]
    .eta <- .ini[!is.na(.ini$neta1) & .ini$neta1 == .ini$neta2, ]
    c(
      paste0("ETA_", .eta$neta1, "_"),
      paste0("THETA_", .th$ntheta[!(.th$name %in% ui$muRefDataFrame$theta)], "_")
    )
  }
  .extMap <- function(ui) {
    .s <- ui$loadPruneSens
    .foceiAnalyticMuRefExtDirMap(
      ui,
      .s,
      rxode2::rxStateOde(.s),
      .extDirs(ui),
      list(),
      rxode2::rxGetControl(ui, "foceiConstCovs", NULL)
    )
  }
  .mat <- "cl <- exp(tcl + eta.cl) * PMA^exp(lhill) / (PMA^exp(lhill) + exp(ltm50)^exp(lhill))"
  .matIni <- "ltm50 <- 3.7; lhill <- 1.1"

  test_that("nonlinear covariate parameters and natural-scale typical values reuse the eta", {
    .m <- .extMap(.extUi(.extMod(.mat, .matIni)))
    expect_setequal(names(.m), c("THETA_4_", "THETA_5_"))
    expect_true(all(vapply(.m, function(x) x$etaDir, "") == "ETA_2_"))
    .m <- .extMap(.extUi(.extMod("cl <- exp(tcl + eta.cl) * (1 + emax * WT / (ec50 + WT))", "emax <- 0.5; ec50 <- 60")))
    expect_setequal(names(.m), c("THETA_4_", "THETA_5_"))
    .m <- .extMap(.extUi(.extMod("cl <- exp(tcl + eta.cl) * (WT / 70)^wt.cl", "wt.cl <- 0.75")))
    expect_named(.m, "THETA_4_")
    .m <- .extMap(.extUi(.extMod("cl <- exp(tcl + eta.cl) * (1 + crcl.cl * (WT - 70))", "crcl.cl <- 0.01")))
    expect_named(.m, "THETA_4_")
    .m <- .extMap(.extUi(.extMod("cl <- tcl * exp(eta.cl)", tcl = 2.7)))
    expect_named(.m, "THETA_2_")
    expect_equal(.m$THETA_2_$etaDir, "ETA_2_")
    # foceiControl(muRefExt = FALSE) turns the detection off
    .ui <- .extUi(.extMod(.mat, .matIni))
    rxode2::rxAssignControlValue(.ui, "muRefExt", FALSE)
    expect_length(.extMap(.ui), 0L)
    # a fixed theta may appear in the scale
    .m <- .extMap(.extUi(.extMod(.mat, "ltm50 <- 3.7; lhill <- fix(1.1)")))
    expect_named(.m, "THETA_4_")
  })

  test_that("reuse is refused when the identity cannot hold", {
    # time-varying covariate in the scale
    expect_length(.extMap(.extUi(.extMod(.mat, .matIni), constCovs = "WT")), 0L)
    # d(cl)/d(tcl) = 1 is not a constant multiple of d(cl)/d(eta.cl) = exp(eta.cl)
    expect_length(.extMap(.extUi(.extMod("cl <- exp(eta.cl) + tcl"))), 0L)
    # eta.cl shared by cl and v: d/deta.cl moves v too, d/dlhill does not
    .m <- .extMap(.extUi(.extMod(.mat, .matIni, vLine = "v <- exp(tv + eta.v + 0.5 * eta.cl)")))
    expect_false(any(c("THETA_4_", "THETA_5_") %in% names(.m)))
    # one exponent on two parameters: d/dwt moves cl and v, d/deta.cl only cl
    .m <- .extMap(.extUi(.extMod(
      "cl <- exp(tcl + eta.cl) * (WT / 70)^wt",
      "wt <- 0.75",
      vLine = "v <- exp(tv + eta.v) * (WT / 70)^wt"
    )))
    expect_length(.m, 0L)
  })

  # solve the augmented model with and without the reuse at the same parameters
  .extCompare <- function(f, d, nStates = NULL, constCovs = c("PMA", "WT")) {
    .ui <- .extUi(f, constCovs)
    .dirs <- .extDirs(.ui)
    .ini <- .ui$iniDf
    .th <- .ini[!is.na(.ini$ntheta), ]
    .eta <- .ini[!is.na(.ini$neta1) & .ini$neta1 == .ini$neta2, ]
    .par <- c(
      stats::setNames(.th$est * 1.05 + 0.02, paste0("THETA_", .th$ntheta, "_")),
      stats::setNames(0.1 * (-1)^seq_len(nrow(.eta)) * seq_len(nrow(.eta)), paste0("ETA_", .eta$neta1, "_"))
    )
    .solve <- function(flag) {
      rxode2::rxAssignControlValue(.ui, "muRefExt", flag)
      .am <- .foceiAnalyticAugModelDirs(.ui, .dirs)
      list(
        n = length(rxode2::rxState(.am$augMod)),
        s = rxode2::rxSolve(.am$augMod, .par, d[d$ID == 3, ], returnType = "data.frame", atol = 1e-12, rtol = 1e-12)
      )
    }
    .on <- .solve(TRUE)
    .off <- .solve(FALSE)
    if (!is.null(nStates)) {
      expect_equal(c(.on$n, .off$n), nStates)
    }
    expect_lt(.on$n, .off$n)
    # the reuse changes how columns are computed, never which columns exist
    expect_setequal(
      grep("^rx_(f1|f2|rvar1|rvar2|rsig1)_", names(.on$s), value = TRUE),
      grep("^rx_(f1|f2|rvar1|rvar2|rsig1)_", names(.off$s), value = TRUE)
    )
    .cols <- grep("^rx_(f1|f2|rvar1|rvar2|rsig1)_.*THETA", names(.on$s), value = TRUE)
    # eta-only columns too: the reuse must leave them untouched
    for (.c in grep("^rx_(f1|f2|rvar1|rvar2|rsig1)_", names(.on$s), value = TRUE)) {
      expect_equal(.on$s[[.c]], .off$s[[.c]], tolerance = 1e-8, info = .c)
    }
    .cols
  }

  test_that("reused sensitivity columns equal the full symbolic build", {
    skip_if_not_installed("nlmixr2data")
    .d <- nlmixr2data::theo_sd
    .d$PMA <- 30 + 3 * .d$ID
    .d$WT <- 70 + 10 * sin(.d$ID)
    # maturation pair (cross curvature) with proportional error (variance columns)
    .cols <- .extCompare(.extMod(.mat, .matIni, err = "prop(prop.sd)"), .d, nStates = c(26L, 62L))
    expect_true("rx_f2_THETA_4__THETA_5_" %in% .cols)
    expect_true(any(grepl("^rx_rvar2_.*THETA_4_", .cols)))
    # natural-scale typical value: c = 1/tcl, self curvature
    .cols <- .extCompare(.extMod("cl <- tcl * exp(eta.cl)", tcl = 2.7, err = "prop(prop.sd)"), .d)
    expect_true("rx_f2_THETA_2__THETA_2_" %in% .cols)
    # a mu-referenced covariate coefficient (theta-free scale) next to extended directions on
    # the same eta: mixed scales in the pair columns
    .cols <- .extCompare(
      .extMod(
        "cl <- exp(tcl + eta.cl + wt.cl * log(WT / 70)) * PMA^exp(lhill) / (PMA^exp(lhill) + exp(ltm50)^exp(lhill))",
        "wt.cl <- 0.75; ltm50 <- 3.7; lhill <- 1.1"
      ),
      .d
    )
    expect_true(all(c("rx_f2_THETA_4__THETA_5_", "rx_f2_THETA_5__THETA_6_") %in% .cols))
  })

  test_that("reuse through a dosing modifier keeps the jump sensitivities", {
    skip_if_not_installed("nlmixr2data")
    .d <- nlmixr2data::theo_sd
    .d$WT <- 70 + 10 * sin(.d$ID)
    .f <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45; tlag <- -1; wt.lag <- 0.5
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1; eta.lag ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        lagD <- exp(tlag + eta.lag) * (1 + wt.lag * (WT - 70) / 70)
        d/dt(depot) <- -ka * depot
        alag(depot) <- lagD
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }
    .m <- .extMap(.extUi(.f, constCovs = "WT"))
    expect_named(.m, "THETA_5_")
    expect_equal(.m$THETA_5_$etaDir, "ETA_4_")
    .extCompare(.f, .d, constCovs = "WT")
  })

  test_that("reuse is refused for delay(), unmaterialized ODEs and error-model thetas", {
    # delay() with the default history: no rx__pastRhs_ variable, only the delay() call
    .ui <- .extUi(eval(parse(
      text = sub(
        "d/dt(center) <- ka * depot",
        "d/dt(center) <- ka * delay(depot, 0.5)",
        paste(deparse(.extMod(.mat, .matIni)), collapse = "\n"),
        fixed = TRUE
      )
    )))
    expect_length(.extMap(.ui), 0L)
    # an ODE right-hand side missing from the pruned env (matExp()/indLin() before .rxJacobian)
    .ui <- .extUi(.extMod(.mat, .matIni))
    .s <- .ui$loadPruneSens
    rm(list = "rx__d_dt_center__", envir = .s)
    expect_length(
      .foceiAnalyticMuRefExtDirMap(.ui, .s, rxode2::rxStateOde(.s), .extDirs(.ui), list(), c("PMA", "WT")),
      0L
    )
    # mtime() right-hand sides are outside the proof
    .ui <- .extUi(.extMod(.mat, .matIni))
    .s <- .ui$loadPruneSens
    expect_true(
      length(.foceiAnalyticMuRefExtDirMap(.ui, .s, rxode2::rxStateOde(.s), .extDirs(.ui), list(), c("PMA", "WT"))) > 0L
    )
    assign("..mtime", "mtime(rx_t1_) <- exp(ltm50)", envir = .s)
    expect_length(
      .foceiAnalyticMuRefExtDirMap(.ui, .s, rxode2::rxStateOde(.s), .extDirs(.ui), list(), c("PMA", "WT")),
      0L
    )
    # an error-model theta is never reused, even when offered as a direction
    .ui <- .extUi(.extMod("cl <- tcl * exp(eta.cl)", tcl = 2.7, err = "prop(prop.sd)"))
    .s <- .ui$loadPruneSens
    .sg <- paste0("THETA_", .ui$iniDf$ntheta[.ui$iniDf$name == "prop.sd"], "_")
    .m <- .foceiAnalyticMuRefExtDirMap(.ui, .s, rxode2::rxStateOde(.s), c(.extDirs(.ui), .sg), list(), character(0))
    expect_false(.sg %in% names(.m))
    expect_true("THETA_2_" %in% names(.m))
  })

  test_that("the analytic outer gradient with extended directions matches central differences", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    .d <- nlmixr2data::theo_sd
    .d$PMA <- 30 + 3 * .d$ID
    .ctl <- foceiControl(
      fast = TRUE,
      maxOuterIterations = 0L,
      maxInnerIterations = 1000L,
      epsilon = 1e-10,
      print = 0,
      covMethod = "",
      calcTables = FALSE,
      compress = FALSE,
      rxControl = rxode2::rxControl(atol = 1e-11, rtol = 1e-11)
    )
    .f <- .extMod(.mat, .matIni)
    .fit <- .nlmixr(.f, .d, "focei", control = .ctl)
    expect_gt(.fit$env$nAnalyticGradDirect, 0)
    # the fit saw PMA as subject-constant, so the maturation thetas were reused
    .cc <- .fit$control$foceiConstCovs
    expect_true("PMA" %in% .cc)
    expect_gt(length(.extMap(.extUi(.f, constCovs = .cc))), 0L)
    .g <- .foceiGradDirect(.fit)
    # the full build (no reuse) gives the same analytic gradient
    .ctlOff <- .ctl
    .ctlOff$muRefExt <- FALSE
    .fitOff <- .nlmixr(.f, .d, "focei", control = .ctlOff)
    expect_false(.fitOff$control$muRefExt)
    expect_equal(.foceiGradDirect(.fitOff), .g, tolerance = 1e-8)
    .ui <- .fit$finalUi
    .objAt <- function(ui) .nlmixr(ui, .d, "focei", control = .ctl)$objf
    .h <- 1e-4
    expect_equal(
      .g[["ltm50"]],
      (.objAt(rxode2::ini(.ui, ltm50 = 3.7 + .h)) - .objAt(rxode2::ini(.ui, ltm50 = 3.7 - .h))) / (2 * .h),
      tolerance = 1e-3
    )
    expect_equal(
      .g[["lhill"]],
      (.objAt(rxode2::ini(.ui, lhill = 1.1 + .h)) - .objAt(rxode2::ini(.ui, lhill = 1.1 - .h))) / (2 * .h),
      tolerance = 1e-3
    )
  })
})
