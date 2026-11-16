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

nmTest({
  test_that("addTable() takes the table's keep and drop as its arguments", {
    # as a fit's own table step does.  It solves the fit's processed data
    # (`dataSav`), which has only the model's columns, so rxode2 warns that a
    # data column can't be kept; the IOV route refits from the data, so it can
    .fit <- one.compartment.fit.focei
    expect_true("IWRES" %in% names(.fit))
    .tab <- tableControl(keep = "WT", drop = "IWRES")
    .w <- capture_warnings(.t <- suppressMessages(addTable(.fit, updateObject = FALSE, table = .tab)))
    expect_identical(unique(.w), "Cannot keep missing columns: WT")
    expect_false("IWRES" %in% names(.t))
    .w <- capture_warnings(.t <- suppressMessages(addTable(.fit, updateObject = FALSE, keep = "WT", drop = "IWRES")))
    expect_identical(unique(.w), "Cannot keep missing columns: WT")
    expect_false("IWRES" %in% names(.t))
  })
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

# A cwresCalc() input with known sensitivities, without a fit.  `occ` > 1 makes
# each eta active on only one block of each subject's rows, like IOV etas.
.cwresSynth <- function(nid, nper, neta, zeroFrac = 0.5, occ = 1L) {
  n <- nid * nper
  id <- rep(seq_len(nid), each = nper)
  .sens <- function() {
    m <- matrix(rnorm(n * neta), n, neta)
    m[runif(n * neta) < zeroFrac] <- 0
    if (occ > 1L) {
      .blk <- rep(rep(seq_len(occ), length.out = nper), nid)
      m[cbind(seq_len(n), 1L)] <- m[cbind(seq_len(n), 1L)] * (.blk == 1L)
      for (j in seq_len(neta)[-1]) {
        m[.blk != ((j - 1L) %% occ) + 1L, j] <- 0
      }
    }
    m
  }
  fp <- .sens()
  fi <- .sens()
  om <- crossprod(matrix(rnorm(neta * neta), neta)) / neta + diag(neta)
  etas <- matrix(rnorm(nid * neta, sd = 0.3), nid, neta)
  .mk <- function(f, pr, r) {
    d <- data.frame(ID = id, time = rep(seq_len(nper), nid), rx_pred_ = pr)
    for (j in seq_len(neta)) {
      d[[paste0("rx__sens_rx_pred__BY_ETA_", j, "___")]] <- f[, j]
    }
    d$rx_r_ <- r
    d$rxLambda <- 1
    d$rxYj <- 0
    d$rxLow <- 0
    d$rxHi <- 1
    d
  }
  pr <- rnorm(n, 10)
  ipr <- pr + rnorm(n, sd = 0.2)
  rp <- runif(n, 0.5, 1)
  ri <- runif(n, 0.5, 1)
  ipred <- .mk(fi, ipr, ri)
  etasDf <- data.frame(ID = seq_len(nid), etas)
  names(etasDf) <- c("ID", paste0("eta", seq_len(neta)))
  list(
    ipred = ipred,
    pred = .mk(fp, pr, rp),
    ebe = ipred[, c("ID", "time")],
    om = om,
    etasDf = etasDf,
    etaLst = lapply(seq_len(neta + 1), function(j) double(10)),
    dv = ipr + rnorm(n, sd = 0.5),
    fp = fp,
    fi = fi,
    etas = etas,
    id = id,
    pr = pr,
    ipr = ipr,
    rp = rp,
    ri = ri
  )
}

.cwresSynthCall <- function(s) {
  .Call(
    `_nlmixr2est_cwresCalc`,
    list(s$ipred, s$pred, s$etaLst, s$ebe),
    s$om,
    s$etasDf,
    s$dv,
    NULL,
    NULL,
    NULL,
    character(0),
    character(0),
    character(0),
    NULL,
    tableControl()
  )[[3]]
}

test_that("cwresCalc() WRES/CWRES equal the dense F Omega F' formula", {
  set.seed(1356)
  for (.occ in c(1L, 3L)) {
    s <- .cwresSynth(30, 9, 3, occ = .occ)
    r <- .cwresSynthCall(s)
    # rxLambda = 1, rxYj = 0 is a Box-Cox transform x - 1; predictions are on that scale
    .vp <- diag(s$fp %*% s$om %*% t(s$fp))
    .vi <- diag(s$fi %*% s$om %*% t(s$fi))
    expect_equal(r$WRES, (s$dv - 1 - s$pr) / sqrt(abs(.vp + s$rp)))
    .cp <- s$ipr - rowSums(s$fi * s$etas[s$id, ])
    expect_equal(r$CWRES, (s$dv - 1 - .cp) / sqrt(.vi + s$ri))
    expect_equal(r$CPRED, .cp + 1)
  }
})

test_that("cwresCalc() never forms a dense nobs x nobs matrix", {
  skip_on_cran()
  set.seed(1356)
  # 60000 observations: one dense V_fo would need ~29 GB
  s <- .cwresSynth(6000, 10, 3, zeroFrac = 0.3)
  r <- .cwresSynthCall(s)
  expect_equal(nrow(r), 60000L)
  expect_true(all(is.finite(r$CWRES)))
})

nmTest({
  test_that("pooled table solves equal the rxSolve() table solves (#1356)", {
    one.cmt <- function() {
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
        v <- exp(tv + eta.v + 0.01 * (WT - 70))
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }
    .d <- nlmixr2data::theo_sd
    .d$WT <- .d$WT + .d$TIME / 10
    fit <- suppressMessages(suppressWarnings(nlmixr2(
      one.cmt,
      .d,
      "focei",
      foceiControl(print = 0, maxOuterIterations = 0, calcTables = FALSE)
    )))
    fit$innerModelForce
    .te <- fit$foceiThetaEtaParameters
    .refI <- .foceiSolvePars(fit, fit$innerModel, .te$ipred, returnType = "data.frame.TBS", what = "ipred")
    .refP <- .foceiSolvePars(fit, fit$innerModel, .te$pred, what = "pred")
    .refE <- .foceiSolvePars(fit, fit$predOnlyModel, .te$ipred, what = "ebe", addCov = TRUE)

    .ok <- .tablePoolSetup(
      fit,
      fit$innerModel,
      fit$predOnlyModel,
      .te$ipred,
      .residOdeMethodName(fit$methodOde),
      fit$atol[1],
      fit$rtol[1]
    )
    on.exit(tablePoolFree_(), add = TRUE)
    expect_equal(.ok, c(inner = TRUE, predOnly = TRUE))
    .n0 <- .odeSwapInfo()$tablePoolSolveN

    .mvI <- rxode2::rxModelVars(fit$innerModel)
    .sens <- grep("^rx__sens_rx_pred__BY_ETA_", .mvI$lhs, value = TRUE)
    .i <- .tablePoolSolve(
      fit$innerModel,
      "inner",
      .te$ipred,
      lhs = c("rx_pred_", "rx_r_"),
      state = .mvI$state,
      sens = .sens,
      tbs = TRUE
    )
    expect_true(.i$rowsOk)
    expect_equal(.i$time, .refI$time)
    for (.n in names(.i$vals)) {
      expect_equal(.i$vals[[.n]], .refI[[.n]], tolerance = 0, label = .n)
    }
    # the eta sensitivities come back sparse, and only their nonzero entries
    .dense <- matrix(0, nrow(.refI), length(.sens))
    .dense[cbind(.i$sens$i, .i$sens$k)] <- .i$sens$x
    expect_equal(.dense, unname(as.matrix(.refI[, .sens])), tolerance = 0)
    expect_true(all(.i$sens$x != 0))

    .p <- .tablePoolSolve(fit$innerModel, "inner", .te$pred, lhs = c("rx_pred_", "rx_r_"))
    expect_equal(.p$vals$rx_pred_, .refP$rx_pred_, tolerance = 0)
    expect_equal(.p$vals$rx_r_, .refP$rx_r_, tolerance = 0)

    .mvE <- rxode2::rxModelVars(fit$predOnlyModel)
    .e <- .tablePoolSolve(
      fit$predOnlyModel,
      "predOnly",
      .te$ipred,
      lhs = .mvE$lhs,
      state = .mvE$state,
      cov = "WT"
    )
    for (.n in names(.e$vals)) {
      expect_equal(.e$vals[[.n]], .refE[[.n]], tolerance = 0, label = .n)
    }
    # the solve's records map back onto the fit's events
    .rows <- unlist(lapply(split(seq_len(nrow(fit$dataSav)), fit$dataSav$ID), function(.r) .r[1] - 1L))
    expect_equal(fit$dataSav$nlmixrRowNums[.rows[.e$id] + .e$rec + 1L], .refE$nlmixrRowNums)
    expect_equal(.odeSwapInfo()$tablePoolSolveN, .n0 + 3)

    tablePoolFree_()
    expect_equal(.odeSwapInfo()$nLoaded, 0L)
    # a fit after the table pool is unaffected by it
    .fit2 <- suppressMessages(suppressWarnings(nlmixr2(
      one.cmt,
      .d,
      "focei",
      foceiControl(print = 0, maxOuterIterations = 0, calcTables = FALSE)
    )))
    expect_equal(.fit2$objf, fit$objf)
  })
})

test_that(".residSolveWithFallback() loosens tolerances, then switches method", {
  .fit <- list(
    methodOde = "liblsoda",
    atol = 1e-6,
    rtol = 1e-6,
    foceiControl = list(rxControl = list(maxAtolRtolFactor = 0.1), odeRecalcFactor = 10, stickyRecalcN = 2)
  )
  .calls <- list()
  .solve <- function(okWhen) {
    function(method, atol, rtol) {
      .calls[[length(.calls) + 1L]] <<- c(method, atol)
      list(rx_pred_ = if (okWhen(method, atol)) 1 else NA_real_)
    }
  }
  # recovers after one loosening
  expect_equal(.residSolveWithFallback(.fit, "ipred", .solve(function(m, a) a > 1e-6))$rx_pred_, 1)
  expect_equal(.calls, list(c("liblsoda", "1e-06"), c("liblsoda", "1e-05")))
  # the next method, with a warning naming the failed one
  .calls <- list()
  expect_warning(
    .r <- .residSolveWithFallback(.fit, "pred", .solve(function(m, a) m == "dop853")),
    "Problems solving pred with liblsoda, returning results from dop853"
  )
  expect_equal(.r$rx_pred_, 1)
  expect_equal(vapply(.calls, `[`, character(1), 1), c("liblsoda", "liblsoda", "dop853"))
  # every method fails: the first method's result, with a warning
  expect_warning(
    .r <- .residSolveWithFallback(.fit, "ebe", .solve(function(m, a) FALSE)),
    "returning results from the first method"
  )
  expect_true(is.na(.r$rx_pred_))
})

nmTest({
  .tablePoolModel <- function() {
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
      cl <- exp(tcl + eta.cl) * (WT / 70)
      v <- exp(tv + eta.v)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }
  .tablePoolData <- function() {
    .d <- nlmixr2data::theo_sd
    .d$WT <- .d$WT + .d$TIME / 5
    .d
  }
  .tablePoolFit <- function(model = .tablePoolModel, data = .tablePoolData()) {
    suppressMessages(suppressWarnings(nlmixr2(
      model,
      data,
      "focei",
      foceiControl(print = 0, maxOuterIterations = 0, calcTables = FALSE)
    )))
  }
  # the table of the rxSolve() route
  .tableOld <- function(fit, ...) {
    withr::with_options(list(nlmixr2.tablePool = FALSE), suppressMessages(addTable(fit, ...)))
  }

  test_that("addTable() solves in the pool and equals the rxSolve() table (#1356)", {
    fit <- .tablePoolFit()
    .ref <- .tableOld(fit)
    .i0 <- .tablePoolInfo()
    .tab <- suppressMessages(addTable(fit))
    .i1 <- .tablePoolInfo()
    expect_equal(.i1$pooledN, .i0$pooledN + 1L)
    expect_equal(.i1$declineN, .i0$declineN)
    expect_equal(as.data.frame(.tab), as.data.frame(.ref), tolerance = 0)
    # the slots are cleared after the table step
    expect_equal(.odeSwapInfo()$nLoaded, 0L)
    # cwres = FALSE solves only the predOnly model
    .ref2 <- .tableOld(fit, table = tableControl(cwres = FALSE))
    .tab2 <- suppressMessages(addTable(fit, table = tableControl(cwres = FALSE)))
    expect_equal(as.data.frame(.tab2), as.data.frame(.ref2), tolerance = 0)
  })

  test_that("the pooled table switches ODE method like the rxSolve() table", {
    fit <- .tablePoolFit()
    .both <- function(...) {
      .w1 <- character(0)
      .old <- withCallingHandlers(.tableOld(fit), warning = function(w) {
        .w1 <<- c(.w1, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
      .w2 <- character(0)
      .i0 <- .tablePoolInfo()
      .new <- withCallingHandlers(suppressMessages(addTable(fit)), warning = function(w) {
        .w2 <<- c(.w2, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
      .i1 <- .tablePoolInfo()
      expect_equal(.i1$pooledN, .i0$pooledN + 1L)
      expect_equal(as.data.frame(.new), as.data.frame(.old), tolerance = 0)
      expect_equal(.w2, .w1)
      .i1$rebuildN - .i0$rebuildN
    }
    # liblsoda "fails" for every solve, so each one is redone with dop853
    withr::with_options(list(nlmixr2.tableFailMethods = "liblsoda"), {
      expect_gt(.both(), 0L)
    })
    # every method fails: the first method's result, with the warnings
    withr::with_options(list(nlmixr2.tableFailMethods = c("liblsoda", "dop853")), {
      expect_gt(.both(), 0L)
    })
  })

  test_that("addTable() gets the right table once the fit's pool is gone", {
    fit <- .tablePoolFit()
    .ref <- .tableOld(fit)
    # another model's fit and a plain rxSolve() replace the registry and the solve
    .other <- function() {
      ini({
        tcl <- 1
        tv <- 3.45
        eta.cl ~ 0.3
        prop.sd <- 0.2
      })
      model({
        cl <- exp(tcl + eta.cl)
        v <- exp(tv)
        cp <- linCmt()
        cp ~ prop(prop.sd)
      })
    }
    .fit2 <- .tablePoolFit(.other, nlmixr2data::theo_sd)
    invisible(rxode2::rxSolve(.fit2, nlmixr2data::theo_sd))
    .tab <- suppressMessages(addTable(fit))
    expect_equal(as.data.frame(.tab), as.data.frame(.ref), tolerance = 0)
    # twice in a row
    .tab2 <- suppressMessages(addTable(fit))
    expect_equal(as.data.frame(.tab2), as.data.frame(.ref), tolerance = 0)
    expect_equal(.odeSwapInfo()$nLoaded, 0L)
    # a saved and reloaded fit
    .f <- tempfile(fileext = ".rds")
    on.exit(unlink(.f), add = TRUE)
    saveRDS(fit, .f)
    .fit3 <- readRDS(.f)
    .i0 <- .tablePoolInfo()
    .tab3 <- suppressMessages(addTable(.fit3))
    expect_equal(.tablePoolInfo()$pooledN, .i0$pooledN + 1L)
    expect_equal(as.data.frame(.tab3), as.data.frame(.ref), tolerance = 0)
  })
})
