# eta.v repeats eta.cl's subject-level block ("id:same:eta.cl")
.npdeSameMod <- function() {
  ini({
    tka <- 0.45
    tcl <- 1
    tv <- 3.45
    add.sd <- 0.7
    eta.cl ~ 0.3
    eta.v ~ same()
  })
  model({
    ka <- exp(tka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    linCmt() ~ add(add.sd)
  })
}

# A one-compartment model, and theo_sd with subject 3's observations all NA,
# so the fit drops it and its table reinserts it
.npdeNoObsMod <- function() {
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

.npdeNoObsData <- function() {
  .d <- nlmixr2data::theo_sd
  .d$DV[.d$ID == 3 & .d$EVID == 0] <- NA
  .d
}

# The source to load in a child process when the tests run on a loaded (not
# installed) nlmixr2est, or "" for the installed one
.npdeChildPkg <- function() {
  # a namespace devtools::load_all() made carries .__DEVTOOLS__
  .ns <- asNamespace("nlmixr2est")
  if (exists(".__DEVTOOLS__", envir = .ns, inherits = FALSE)) {
    getNamespaceInfo(.ns, "path")
  } else {
    ""
  }
}

# In a child process: npdeCalc() given a simulation with 3 observations per
# sample for 2 observations, which used to index past them
.npdeCalcMismatchChild <- function(pkg) {
  if (nzchar(pkg)) {
    devtools::load_all(pkg, quiet = TRUE)
  } else {
    library(nlmixr2est)
  }
  .sim <- data.frame(
    sim.id = rep(1:2, each = 3),
    id = rep(c(1L, 1L, 2L), 2),
    sim = c(1, 2, 3, 1, 2, 3),
    rxLambda = 1,
    rxYj = 0,
    rxLow = 0,
    rxHi = 1
  )
  tryCatch(
    .Call("_nlmixr2est_npdeCalc", .sim, c(1, 2), c(0L, 0L), NULL, NULL, tableControl(), PACKAGE = "nlmixr2est"),
    error = conditionMessage
  )
}

# In a child process, so an abort fails the test: the NPDE of a fit with an
# observation-less subject, from the fit's table or from addNpde()
.npdeNoObsChild <- function(pkg, est, mode, mod, dat) {
  if (nzchar(pkg)) {
    devtools::load_all(pkg, quiet = TRUE)
  } else {
    library(nlmixr2est)
  }
  .ctl <- list(print = 0, covMethod = "")
  .fit <- if (mode == "table") {
    suppressWarnings(suppressMessages(nlmixr2(mod, dat, est, control = .ctl, table = tableControl(npde = TRUE))))
  } else {
    suppressWarnings(suppressMessages(addNpde(nlmixr2(mod, dat, est, control = .ctl), updateObject = FALSE)))
  }
  list(id = as.character(.fit$ID), npde = .fit$NPDE)
}

nmTest({
  test_that("npde", {
    one.cmt <- function() {
      ini({
        ## You may label each parameter with a comment
        tka <- 0.45 # Ka
        tcl <- log(c(0, 2.7, 100)) # Log Cl
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
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }

    # Don't use saemControlFast because numeric results are tested below
    fit <- .nlmixr(one.cmt, theo_sd, est = "saem")

    expect_false(all(c("EPRED", "ERES", "NPDE", "NPD", "PDE", "PD") %in% names(fit)))
    # npde's uniforms are seeded per observation, so the thread count cannot move them
    .cols <- c("EPRED", "ERES", "NPDE", "NPD", "PDE", "PD")
    .old <- rxode2::getRxThreads(verbose = FALSE)
    rxode2::setRxThreads(1L)
    .t1 <- as.data.frame(suppressMessages(addNpde(fit, updateObject = FALSE)))[, .cols]
    rxode2::setRxThreads(2L)
    .t2 <- as.data.frame(suppressMessages(addNpde(fit, updateObject = FALSE)))[, .cols]
    rxode2::setRxThreads(.old)
    expect_identical(.t1, .t2)
    suppressMessages(expect_error(addNpde(fit), NA))
    expect_true(all(c("EPRED", "ERES", "NPDE", "NPD", "PDE", "PD") %in% names(fit)))
    .range1 <- range(fit$PDE)
    .range2 <- range(fit$PD)
    expect_true(.range1[1] > 0)
    expect_true(.range2[1] > 0)
    expect_true(.range1[2] < 1)
    expect_true(.range2[2] < 1)

    .range1 <- range(fit$NPDE)
    .range2 <- range(fit$NPD)
    .range3 <- range(fit$ERES)

    expect_true(.range1[1] < -1.7)
    expect_true(.range2[1] < -1.7)
    expect_true(.range3[1] < -1.7)
    expect_true(.range1[2] > 1.7)
    expect_true(.range2[2] > 1.7)
    expect_true(.range3[2] > 1.7)

    .range4 <- range(fit$EPRED)
    expect_true(.range4[1] > -0.1)
    expect_true(.range4[2] > 7)

    fit <- .nlmixr(one.cmt, theo_sd, est = "saem")

    expect_false(all(c("EPRED", "ERES", "NPDE", "NPD", "PDE", "PD") %in% names(fit)))
    fit2 <- suppressMessages(addNpde(fit, updateObject = FALSE))
    expect_false(all(c("EPRED", "ERES", "NPDE", "NPD", "PDE", "PD") %in% names(fit)))
    expect_true(all(c("EPRED", "ERES", "NPDE", "NPD", "PDE", "PD") %in% names(fit2)))

    fit <- .nlmixr(one.cmt, theo_sd, est = "saem", control = saemControlFast, table = tableControl(npde = TRUE))

    expect_true(all(c("EPRED", "ERES", "NPDE", "NPD", "PDE", "PD") %in% names(fit)))
  })

  test_that("pheno", {
    pheno <- function() {
      ini({
        tcl <- log(0.008) # typical value of clearance
        tv <-  log(0.6)   # typical value of volume
        ## var(eta.cl)
        eta.cl + eta.v ~ c(1,
                           0.01, 1) ## cov(eta.cl, eta.v), var(eta.v)
        # interindividual variability on clearance and volume
        add.err <- 0.1    # residual variability
      })
      model({
        cl <- exp(tcl + eta.cl) # individual value of clearance
        v <- exp(tv + eta.v)    # individual value of volume
        ke <- cl / v            # elimination rate constant
        d/dt(A1) = - ke * A1    # model differential equation
        cp = A1 / v             # concentration in plasma
        cp ~ add(add.err)       # define error model
      })
    }

    fit <- .nlmixr(pheno, pheno_sd, est = "saem", control = saemControlFast, table = list(npde = TRUE))

    # Since there is a correlation here the npde and npd

    expect_false(isTRUE(all.equal(fit$NPDE, fit$NPD)))
    expect_false(isTRUE(all.equal(fit$PDE, fit$PD)))
  })

  test_that(".uiHasIov() finds occasion etas, not repeated subject blocks", {
    expect_true(.uiHasIov(rxode2::rxode2(.iovFitMod)))
    expect_false(.uiHasIov(rxode2::rxode2(.npdeSameMod)))
    # the rows: eta.cl and eta.v ("id:same:eta.cl") are subject-level, iov.cl
    # is the occasion eta, and the thetas and add.sd have no condition or are
    # residual
    .iniDf <- rxode2::rxode2(.iovFitMuMod)$iniDf
    expect_identical(.iniDf$name[.iniDfIovRows(.iniDf)], "iov.cl")
    .iniDf <- rxode2::rxode2(.npdeSameMod)$iniDf
    expect_identical(.iniDf$condition[.iniDf$name == "eta.v"], "id:same:eta.cl")
    expect_false(any(.iniDfIovRows(.iniDf)))
  })

  test_that("npdeCalc() stops when the simulation does not match the observations", {
    skip_if_not_installed("callr")
    expect_identical(
      callr::r(.npdeCalcMismatchChild, args = list(pkg = .npdeChildPkg())),
      "npdeCalc: the simulation has 3 observations per sample (2 samples, 6 rows), but the data 2"
    )
  })

  test_that("NPDE of a fit with a subject without observations", {
    # the NPDE simulation solves the fit's data, the dropped subject included;
    # its rows misaligned every later subject's and aborted R
    skip_if_not_installed("callr")
    .pkg <- .npdeChildPkg()
    for (.est in c("focei", "saem")) {
      for (.mode in c("table", "add")) {
        .lab <- paste(.est, .mode)
        .r <- callr::r(
          .npdeNoObsChild,
          args = list(pkg = .pkg, est = .est, mode = .mode, mod = .npdeNoObsMod, dat = .npdeNoObsData())
        )
        expect_identical(length(.r$npde), 132L, label = .lab)
        expect_identical(sum(.r$id == "3"), 11L, label = .lab)
        expect_true(all(is.na(.r$npde[.r$id == "3"])), label = .lab)
        expect_true(all(is.finite(.r$npde[.r$id != "3"])), label = .lab)
      }
    }
  })
})
