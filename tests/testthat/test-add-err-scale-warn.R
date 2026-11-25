nmTest({
  # #1383: warn when a combined error model's additive SD starts negligibly
  # small relative to the observed DV scale
  .mod <- function(addSd = 0.01, err = "comb", fixAdd = FALSE) {
    .f <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        add.sd <- 0.01
        prop.sd <- 0.2
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl)
        v <- exp(tv)
        cp <- linCmt()
        cp ~ add(add.sd) + prop(prop.sd)
      })
    }
    .ui <- rxode2::rxode2(.f)
    .ui <- rxode2::ini(.ui, add.sd = addSd)
    if (fixAdd) {
      .ui <- rxode2::ini(.ui, add.sd = fix(addSd))
    }
    if (err == "add") {
      .ui <- rxode2::model(.ui, cp ~ add(add.sd))
    } else if (err == "prop") {
      .ui <- rxode2::model(.ui, cp ~ prop(prop.sd))
    }
    .ui
  }

  .dat <- nlmixr2data::theo_sd
  .dat <- .dat[.dat$ID <= 4, ]
  .big <- .dat
  .big$DV <- .big$DV * 50 # median |DV| in the hundreds

  .warns <- function(ui, data, est = "focei") {
    .w <- character(0)
    withCallingHandlers(
      .preProcessAddErrScaleWarn(ui, est, data, list()),
      warning = function(w) {
        .w <<- c(.w, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    .w
  }

  test_that("tiny additive start on a combined model warns", {
    .w <- .warns(.mod(0.01), .big)
    expect_length(.w, 1L)
    expect_match(.w, "'cp' add sd 'add.sd'=0.01 < 1% median |DV|", fixed = TRUE)
    expect_lt(nchar(.w), 75)
  })

  test_that("a data-scaled start does not warn", {
    .med <- stats::median(abs(.big$DV[.big$EVID == 0]))
    expect_length(.warns(.mod(0.1 * .med), .big), 0L)
  })

  test_that("additive-only, proportional-only and fixed additive do not warn", {
    expect_length(.warns(.mod(0.01, "add"), .big), 0L)
    expect_length(.warns(.mod(0.01, "prop"), .big), 0L)
    expect_length(.warns(.mod(0.01, fixAdd = TRUE), .big), 0L)
  })

  test_that("no positive DV and simulation do not warn", {
    .neg <- .big
    .neg$DV <- -abs(.neg$DV)
    expect_length(.warns(.mod(0.01), .neg), 0L)
    expect_length(.warns(.mod(0.01), .big, "rxSolve"), 0L)
  })

  test_that("each endpoint is judged on its own DV scale", {
    .f <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        add.sd <- 0.5
        prop.sd <- 0.2
        add.sd2 <- 0.5
        prop.sd2 <- 0.2
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl)
        v <- exp(tv)
        d / dt(depot) <- -ka * depot
        d / dt(central) <- ka * depot - cl / v * central
        cp <- central / v
        cp2 <- cp
        cp ~ add(add.sd) + prop(prop.sd)
        cp2 ~ add(add.sd2) + prop(prop.sd2)
      })
    }
    .ui <- rxode2::rxode2(.f)
    .dose <- .dat[.dat$EVID != 0, ]
    .d1 <- .dat[.dat$EVID == 0, ]
    .d2 <- .d1
    .d2$DV <- .d2$DV * 1000
    .med1 <- stats::median(abs(.d1$DV))
    .med2 <- stats::median(abs(.d2$DV))
    .stack <- function(code1, code2, col) {
      .dose[[col]] <- if (col == "CMT") "depot" else NA
      .d1[[col]] <- code1
      .d2[[col]] <- code2
      .d <- rbind(.dose, .d1, .d2)
      .d[order(.d$ID, .d$TIME, -.d$EVID), ]
    }
    for (.d in list(.stack(1L, 2L, "DVID"), .stack("cp", "cp2", "CMT"))) {
      # 0.5 is fine against cp's median but under 1% of cp2's
      .w <- .warns(.ui, .d)
      expect_length(.w, 1L)
      expect_match(.w, paste0("'cp2' add sd 'add.sd2'=0.5 < 1% median |DV| ", signif(.med2, 3)), fixed = TRUE)
      # the small endpoint warns on its own median, not the pooled one
      .w <- .warns(rxode2::ini(.ui, add.sd = 0.001, add.sd2 = 0.5 * .med2), .d)
      expect_length(.w, 1L)
      expect_match(.w, paste0("'cp' add sd 'add.sd'=0.001 < 1% median |DV| ", signif(.med1, 3)), fixed = TRUE)
    }
  })

  test_that("warfarin cp and pca are each checked against their own median", {
    .f <- function() {
      ini({
        tktr <- log(1)
        tka <- log(1)
        tcl <- log(0.1)
        tv <- log(10)
        prop.err <- 0.1
        pkadd.err <- 0.001
        temax <- logit(0.8)
        tec50 <- log(0.5)
        tkout <- log(0.05)
        te0 <- log(100)
        pdadd.err <- 0.001
        pdprop.err <- 0.1
      })
      model({
        ktr <- exp(tktr)
        ka <- exp(tka)
        cl <- exp(tcl)
        v <- exp(tv)
        emax <- expit(temax)
        ec50 <- exp(tec50)
        kout <- exp(tkout)
        e0 <- exp(te0)
        DCP <- center / v
        PD <- 1 - emax * DCP / (ec50 + DCP)
        effect(0) <- e0
        kin <- e0 * kout
        d / dt(depot) <- -ktr * depot
        d / dt(gut) <- ktr * depot - ka * gut
        d / dt(center) <- ka * gut - cl / v * center
        d / dt(effect) <- kin * PD - kout * effect
        cp <- center / v
        cp ~ prop(prop.err) + add(pkadd.err)
        effect ~ add(pdadd.err) + prop(pdprop.err) | pca
      })
    }
    .w <- nlmixr2data::warfarin
    .med <- tapply(.w$dv[.w$evid == 0], .w$dvid[.w$evid == 0], function(x) median(abs(x)))
    expect_false(isTRUE(all.equal(.med[["cp"]], .med[["pca"]])))
    .ws <- .warns(rxode2::rxode2(.f), .w)
    expect_length(.ws, 2L)
    expect_true(any(grepl(
      paste0("'cp' add sd 'pkadd.err'=0.001 < 1% median |DV| ", signif(.med[["cp"]], 3), ";"),
      .ws,
      fixed = TRUE
    )))
    expect_true(any(grepl(
      paste0("'pca' add sd 'pdadd.err'=0.001 < 1% median |DV| ", signif(.med[["pca"]], 3), ";"),
      .ws,
      fixed = TRUE
    )))
  })

  test_that("the warning reaches $runInfo and leaves the initial estimates alone", {
    skip_on_cran()
    .ui <- .mod(0.01)
    fit <- suppressMessages(suppressWarnings(
      nlmixr2(.ui, .big, est = "focei", control = foceiControl(print = 0, maxOuterIterations = 0, covMethod = ""))
    ))
    expect_true(any(grepl("'cp' add sd 'add.sd'=0.01 <", fit$runInfo, fixed = TRUE)))
    # no outer iterations: the estimate is the (untouched) initial value
    expect_equal(unname(fit$theta["add.sd"]), 0.01)
  })

  test_that("a cached simulationModel does not break rebuilding the ui", {
    .ui <- rxode2::rxUiDecompress(.mod(0.01))
    invisible(.ui$simulationModel)
    .f <- .getUiFunFromIniAndModel(.ui, quote(ini({ a <- 1 })), quote(model({ y <- a })))
    expect_true(is.function(.f))
  })
})
