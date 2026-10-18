# The table step's C routines (resCalc, cwresCalc, iresCalc) must survive a
# garbage collection at any allocation.  Each builds its table from three
# getDfSubsetVars() results (the state, lhs and covariate columns), and
# getDfSubsetVars() returns an unprotected SEXP: one left unprotected while the
# next allocation runs is collected, and the table loses those columns or R's
# heap is corrupted.  A normal run needs a GC to land in that window of a few
# allocations, so the defect shows rarely; gctorture() puts a GC at every
# allocation, so it shows within a few calls.  ASAN does not report it: the
# dangling SEXP is only touched through libR, which is not instrumented.

# Call a table routine on a deep copy of `args` (the routines write into some
# of their inputs), with a GC at every allocation when `torture` is TRUE.
.residGcCall <- function(routine, args, torture) {
  args <- unserialize(serialize(args, NULL))
  if (torture) {
    gctorture(TRUE)
    on.exit(gctorture(FALSE))
  }
  do.call(.Call, c(list(routine), args))
}

# The routine must give the table the fit got, GC or not.  `captured` holds
# the routine's arguments and the result the table step got from them.
.residGcExpect <- function(routine, captured) {
  # the arguments are the ones the table step used: same result without a GC
  expect_identical(.residGcCall(routine, captured$args, FALSE), captured$ret)
  # a torture pass loses columns only when the freed node is reused before it
  # is read, so try three times
  for (.i in 1:3) {
    expect_identical(.residGcCall(routine, captured$args, TRUE), captured$ret)
  }
}

nmTest({
  test_that("the table routines keep every column when a GC runs at every allocation", {
    skip_on_cran()
    .ns <- asNamespace("nlmixr2est")
    .cap <- new.env(parent = emptyenv())
    # record each routine's arguments and result as the table step calls it.
    # Copied at once: the table step goes on to change both in place (it drops
    # the result's dim attributes, for one).
    trace(
      ".calcCwres0",
      where = .ns,
      print = FALSE,
      exit = bquote({
        if (!npde) {
          assign(
            if (predOnly) "res" else "cwres",
            unserialize(serialize(
              list(
                args = list(
                  .prdLst,
                  fit$omega,
                  fit$eta,
                  .prdLst$ipred$dv,
                  .prdLst$ipred$evid,
                  .prdLst$ipred$cens,
                  .prdLst$ipred$limit,
                  .lhs,
                  .state,
                  .params,
                  fit$IDlabel,
                  table
                ),
                ret = returnValue()
              ),
              NULL
            )),
            envir = .(.cap)
          )
        }
      })
    )
    withr::defer(untrace(".calcCwres0", where = .ns))
    trace(
      ".calcIres",
      where = .ns,
      print = FALSE,
      exit = bquote(assign(
        "ires",
        unserialize(serialize(
          list(
            args = list(.ipred, dv, .ipred$evid, .ipred$cens, .ipred$limit, .lhs, .state, .params, fit$IDlabel, table),
            ret = .ret
          ),
          NULL
        )),
        envir = .(.cap)
      ))
    )
    withr::defer(untrace(".calcIres", where = .ns))

    # a covariate and the lhs make all three column groups non-empty
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
        cl <- exp(tcl + eta.cl) * (WT / 70)^0.75
        v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    .popMod <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl) * (WT / 70)^0.75
        v <- exp(tv)
        linCmt() ~ add(add.sd)
      })
    }
    .d <- nlmixr2data::theo_sd
    # saem has no inner model: resCalc
    .nlmixr(.mod, .d, "saem", saemControl(print = 0, nBurn = 1, nEm = 1, covMethod = ""))
    # focei: cwresCalc
    .nlmixr(.mod, .d, "focei", foceiControl(print = 0, maxOuterIterations = 0, covMethod = ""))
    # no etas: iresCalc
    .nlmixr(.popMod, .d, "focei", foceiControl(print = 0, maxOuterIterations = 0, covMethod = ""))

    # premise: each routine ran and built its column groups (the population
    # fit's solve does not carry the covariates, so iresCalc has no WT)
    expect_true(all(c("res", "cwres", "ires") %in% ls(.cap)))
    expect_true(all(c("depot", "central", "ka", "cl", "v", "WT") %in% names(.cap$res$ret[[3]])))
    expect_true(all(c("depot", "central", "ka", "cl", "v", "WT") %in% names(.cap$cwres$ret[[3]])))
    expect_true(all(c("depot", "central", "ka", "cl", "v") %in% names(.cap$ires$ret)))

    .residGcExpect(`_nlmixr2est_resCalc`, .cap$res)
    .residGcExpect(`_nlmixr2est_cwresCalc`, .cap$cwres)
    .residGcExpect(`_nlmixr2est_iresCalc`, .cap$ires)
  })
})

test_that("popResFinal keeps every column when a GC runs at every allocation", {
  skip_on_cran()
  # the two shapes the table step hands it: residuals only, and residuals with
  # the ipred/npde list beside them
  .res <- list(c(1, 2, 3), data.frame(a = 1:3, b = 4:6), data.frame(c = 7:9), list(shrink = 0.1))
  .ipred <- list(c(1.5, 2.5, 3.5), data.frame(d = 10:12))
  for (.in in list(list(.res, NULL), list(.res, .ipred), list(.ipred, .res))) {
    .want <- .residGcCall("_nlmixr2est_popResFinal", list(.in), FALSE)
    expect_named(.want, c("resid", "shrink"))
    for (.i in 1:5) {
      expect_identical(.residGcCall("_nlmixr2est_popResFinal", list(.in), TRUE), .want)
    }
  }
  expect_identical(
    names(.residGcCall("_nlmixr2est_popResFinal", list(list(.res, NULL)), FALSE)$resid),
    c("a", "b", "DV", "c")
  )
})

# The arguments .calcCwres0() hands resCalc() (cwres = FALSE) or cwresCalc()
# (cwres = TRUE), with `eta` as the eta table
.residCall <- function(fit, eta, cwres) {
  .prd <- .foceiPredIpredList(
    fit,
    data = fit$dataSav,
    keep = "nlmixrRowNums",
    thetaEtaParameters = fit$foceiThetaEtaParameters,
    predOnly = !cwres,
    addDosing = FALSE,
    subsetNonmem = TRUE
  )
  .state <- c(fit$predOnlyModel$state, fit$predOnlyModel$stateExtra)
  .lhs <- setdiff(unique(.getRelevantLhs(fit, "nlmixrRowNums", .prd$ipred)), .state)
  .params <- setdiff(intersect(names(fit$dataSav), fit$predOnlyModel$params), c("CMT", "cmt", "Cmt", .state, .lhs))
  .tab <- tableControl()
  .tab$doSim <- TRUE
  .Call(
    if (cwres) `_nlmixr2est_cwresCalc` else `_nlmixr2est_resCalc`,
    .prd,
    fit$omega,
    eta,
    .prd$ipred$dv,
    .prd$ipred$evid,
    .prd$ipred$cens,
    .prd$ipred$limit,
    .lhs,
    .state,
    .params,
    fit$IDlabel,
    .tab,
    PACKAGE = "nlmixr2est"
  )
}

# A model with a normal endpoint (cp) and a count (cnt), on theo_sd with
# counts for every subject; `cpless` subjects have counts only
.residTwoEndpointMod <- function() {
  ini({
    tka <- 0.45
    tcl <- 1
    tv <- 3.45
    tl <- log(3)
    eta.ka ~ 0.6
    eta.cl ~ 0.3
    eta.v ~ 0.1
    eta.l ~ 0.1
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    lam <- exp(tl + eta.l)
    d / dt(depot) <- -ka * depot
    d / dt(center) <- ka * depot - cl / v * center
    cp <- center / v
    cp ~ add(add.sd)
    cnt ~ dpois(lam)
  })
}

.residTwoEndpointData <- function(cpless = integer(0), drop = integer(0)) {
  .d <- nlmixr2data::theo_sd
  .d <- .d[!(.d$ID %in% drop), ]
  .obs <- .d[.d$EVID == 0, ]
  .cnt <- .obs
  .cnt$DV <- (seq_len(nrow(.cnt)) %% 5) + 1
  .cnt$CMT <- "cnt"
  .obs$CMT <- "cp"
  .dose <- .d[.d$EVID != 0, ]
  .dose$CMT <- "depot"
  .d <- rbind(.dose, .obs[!(.obs$ID %in% cpless), ], .cnt)
  .d[order(.d$ID, .d$TIME, -.d$EVID), ]
}

nmTest({
  test_that("residuals take each row's etas from its own subject", {
    # CWRES uses only the normally distributed rows, so a subject with counts
    # only has none there, and its absence shifted the etas by position.  At
    # maxOuterIterations = 0 the other subjects' etas and residuals do not
    # depend on that subject, so they equal those of the same data without
    # it, to the inner problem's convergence (the shifted ones were 0.3-0.9
    # off).
    # With subject 5 count-only, every subject before it took the next
    # subject's etas and the first one's were left unset; with subject 1,
    # subject 2's were.
    .ctl <- foceiControl(print = 0, covMethod = "", maxOuterIterations = 0L)
    for (.k in c(1L, 5L)) {
      .with <- suppressWarnings(suppressMessages(
        nlmixr2(.residTwoEndpointMod, .residTwoEndpointData(cpless = .k), "focei", control = .ctl)
      ))
      .without <- suppressWarnings(suppressMessages(
        nlmixr2(.residTwoEndpointMod, .residTwoEndpointData(drop = .k), "focei", control = .ctl)
      ))
      .a <- as.data.frame(.with)
      .b <- as.data.frame(.without)
      .a <- .a[.a$CMT == "cp", ]
      .b <- .b[.b$CMT == "cp", ]
      expect_identical(nrow(.a), 121L, label = .k)
      expect_identical(as.character(.a$ID), as.character(.b$ID), label = .k)
      for (.c in c("CWRES", "CRES", "CPRED", "eta.cl")) {
        expect_equal(.a[[.c]], .b[[.c]], tolerance = 1e-4, label = paste(.k, .c))
      }
    }
  })

  test_that("residuals stop when a subject of the rows has no eta row", {
    .fit <- suppressWarnings(suppressMessages(nlmixr2(
      .residTwoEndpointMod,
      .residTwoEndpointData(),
      "focei",
      control = foceiControl(print = 0, covMethod = "", maxOuterIterations = 0L)
    )))
    .eta <- .fit$eta
    expect_true(is.list(.residCall(.fit, .eta, FALSE)))
    expect_true(is.list(.residCall(.fit, .eta, TRUE)))
    .bad <- .eta[.eta$ID != "2", ]
    for (.cw in c(FALSE, TRUE)) {
      .e <- expect_error(.residCall(.fit, .bad, .cw))
      expect_identical(
        conditionMessage(.e),
        "subject code 2 of the solved rows has no row in the eta table",
        label = .cw
      )
    }
  })
})
