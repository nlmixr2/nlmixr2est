nmTest({
  # mtime() in a model used to break every fit: etTrans() materializes the
  # modeled times as EVID 10-99 records (TIME=0, AMT=NA) and $dataSav kept them,
  # so the estimation solves saw them as input doses with a missing amt.  On top
  # of that rxode2::rxS() drops the mtime() assignment, so the generated models
  # lost the modeled times and left the mtime variable undefined (#919).
  .mkMtime <- function(useMtime) {
    .bdy <- c("ka <- exp(tka + eta.ka)", "cl <- exp(tcl)", "v <- exp(tv)",
              if (useMtime) "mtime(t5) <- 5",
              if (useMtime) "kmult <- ifelse(t < t5, 1.0, 2.0)"
              else "kmult <- ifelse(t < 5, 1.0, 2.0)",
              "d/dt(depot) <- -ka * depot",
              "d/dt(center) <- ka * depot - kmult * cl / v * center",
              "cp <- center / v", "cp ~ add(add.sd)")
    eval(parse(text=paste0(
      "function() {\n ini({tka <- 0.45; tcl <- -3.2; tv <- -1; eta.ka ~ 0.1; add.sd <- 0.7})\n",
      " model({\n", paste(.bdy, collapse="\n"), "\n })\n}")))
  }

  .theo <- nlmixr2data::theo_sd

  test_that(".rxMtimeRhs() pulls the mtime declarations out of model text", {
    expect_equal(.rxMtimeRhs("d/dt(depot)=-ka*depot;\n"), character(0))
    expect_equal(.rxMtimeRhs("mtime(t5)=5;\nd/dt(depot)=-ka*depot;\n"),
                 c(t5="5"))
    expect_equal(.rxMtimeRhs("mtime(t5)~5;\nmtime(tx)=2*THETA[1];\n"),
                 c(t5="5", tx="2*THETA[1]"))
  })

  test_that("$dataSav keeps no mtime (EVID 10-99) records", {
    .ui <- rxode2::rxode2(.mkMtime(TRUE))
    # the mtime records DO come out of etTrans(), which is what used to be saved
    .et <- as.data.frame(rxode2::etTrans(.theo, .ui, addCmt=TRUE, dropUnits=TRUE,
                                         allTimeVar=TRUE, keepDosingOnly=FALSE))
    expect_true(any(.et$EVID >= 10 & .et$EVID <= 99))

    .env <- new.env(parent=emptyenv())
    .env$table <- tableControl()
    .foceiPreProcessData(.theo, .env, .ui, rxode2::rxControl())
    expect_false(any(.env$dataSav$EVID >= 10 & .env$dataSav$EVID <= 99))
    # EVID 9 (system init) is below the mtime range and must survive
    expect_equal(sum(.env$dataSav$EVID == 0), sum(.theo$EVID == 0))
  })

  test_that("every generated model keeps the mtime() declaration", {
    .ui <- rxode2::rxUiDecompress(.mkMtime(TRUE)())
    .f <- .ui$focei
    for (.m in c("inner", "predOnly", "predNoLhs")) {
      expect_equal(rxode2::rxModelVars(.f[[.m]])$nMtime, 1L,
                   info=paste0("focei ", .m, " keeps mtime()"))
    }
    expect_equal(rxode2::rxModelVars(.ui$saemModel)$nMtime, 1L)
    expect_equal(rxode2::rxModelVars(.ui$saemModelPred$predOnly)$nMtime, 1L)
    # it is emitted suppressed, so it adds no output column that would shift
    # the positional lhs layout inner.cpp reads
    expect_false("t5" %in% rxode2::rxModelVars(.f$inner)$lhs)
  })

  test_that("focei fits a model with mtime() and matches an independent solve", {
    .ctl <- foceiControl(print=0, maxOuterIterations=0, maxInnerIterations=0,
                         covMethod="", calcTables=TRUE, sigdig=8,
                         rxControl=rxode2::rxControl(atol=1e-10, rtol=1e-10))
    .fit <- nlmixr2(.mkMtime(TRUE), .theo, est="focei", control=.ctl)
    expect_true(inherits(.fit, "nlmixr2FitData"))
    # the mtime records are model output, not data: no extra DV=NA rows
    .df <- as.data.frame(.fit)
    expect_equal(nrow(.df), sum(.theo$EVID == 0))
    expect_false(anyNA(.df$DV))

    # PRED at eta=0 must reproduce an independent rxode2 solve.  The reference is
    # plain rxode2 model text writing the same switch with a literal 5 instead of
    # the mtime variable, so it needs none of the machinery under test, emits no
    # mtime records, and goes nowhere near the UI/estimation pipeline.
    .mod <- rxode2::rxode2("
      ka=exp(tka);cl=exp(tcl);v=exp(tv);
      kmult=ifelse(t<5,1,2);
      d/dt(depot)=-ka*depot;
      d/dt(center)=ka*depot-kmult*cl/v*center;
      cp=center/v;")
    .s <- suppressWarnings(
      rxode2::rxSolve(.mod, c(tka=0.45, tcl=-3.2, tv=-1), .theo,
                      atol=1e-10, rtol=1e-10, returnType="data.frame"))
    expect_equal(nrow(.s), nrow(.df))
    expect_equal(.df$PRED, .s$cp, tolerance=1e-5)
  })

  test_that("an mtime right hand side that is not a constant is carried", {
    .mk <- function(lines, useTsw=FALSE) {
      eval(parse(text=paste0(
        "function() {\n ini({tka <- 0.45; tcl <- -3.2; tv <- -1;",
        if (useTsw) " tsw <- 1.6;" else "", " eta.ka ~ 0.1; add.sd <- 0.7})\n",
        " model({\n", paste(lines, collapse="\n"), "\n })\n}")))
    }
    .base <- c("ka <- exp(tka + eta.ka)", "cl <- exp(tcl)", "v <- exp(tv)")
    .tail <- c("d/dt(depot) <- -ka * depot",
               "d/dt(center) <- ka * depot - kmult * cl / v * center",
               "cp <- center / v", "cp ~ add(add.sd)")
    .mtimeLines <- function(mod) {
      grep("^mtime", strsplit(rxode2::rxModelVars(mod)$model["normModel"], "\n")[[1]],
           value=TRUE)
    }
    # two mtimes: both survive, in declaration order
    .ui <- rxode2::rxUiDecompress(.mk(c(.base, "mtime(t5) <- 5", "mtime(t8) <- 8",
      "kmult <- ifelse(t < t5, 1.0, ifelse(t < t8, 2.0, 3.0))", .tail))())
    expect_equal(.mtimeLines(.ui$focei$inner), c("mtime(t5)~5;", "mtime(t8)~8;"))

    # one mtime referencing an earlier one: rxS() leaves the earlier variable
    # unbound, so the expansion has to bind it to itself
    .ui <- rxode2::rxUiDecompress(.mk(c(.base, "mtime(ta) <- 5", "mtime(tb) <- ta + 3",
      "kmult <- ifelse(t < ta, 1.0, ifelse(t < tb, 2.0, 3.0))", .tail))())
    expect_equal(.mtimeLines(.ui$focei$inner), c("mtime(ta)~5;", "mtime(tb)~3+ta;"))

    # a theta comes back in the generated model's own namespace: THETA[#] for
    # focei, the natural name for saem
    .ui <- rxode2::rxUiDecompress(.mk(c(.base, "mtime(tsw5) <- exp(tsw)",
      "kmult <- ifelse(t < tsw5, 1.0, 2.0)", .tail), useTsw=TRUE)())
    expect_equal(.mtimeLines(.ui$focei$inner), "mtime(tsw5)~exp(THETA[4]);")
    expect_equal(grep("^mtime", strsplit(.ui$saemModel, "\n")[[1]], value=TRUE),
                 "mtime(tsw5)~exp(tsw)")

    # a covariate stays a covariate
    .ui <- rxode2::rxUiDecompress(.mk(c(.base, "mtime(tw) <- WT / 10",
      "kmult <- ifelse(t < tw, 1.0, 2.0)", .tail))())
    expect_equal(.mtimeLines(.ui$focei$inner), "mtime(tw)~0.1*WT;")
  })

  test_that("mtime() does not cost the table its ADDL doses", {
    # the extra mtime rows are dropped from a table solve; an ADDL-expanded dose
    # also has no source data row and must NOT be dropped with them
    .mkAddl <- function(useMtime) {
      .bdy <- c("ka <- exp(tka + eta.ka)", "cl <- exp(tcl)", "v <- exp(tv)",
                if (useMtime) "mtime(t5) <- 5",
                if (useMtime) "kmult <- ifelse(t < t5, 1.0, 2.0)"
                else "kmult <- ifelse(t < 5, 1.0, 2.0)",
                "d/dt(depot) <- -ka * depot",
                "d/dt(center) <- ka * depot - kmult * cl / v * center",
                "cp <- center / v", "cp ~ add(add.sd)")
      eval(parse(text=paste0(
        "function() {\n ini({tka <- 0.45; tcl <- -3.2; tv <- -1; eta.ka ~ 0.1; add.sd <- 0.7})\n",
        " model({\n", paste(.bdy, collapse="\n"), "\n })\n}")))
    }
    .d <- do.call(rbind, lapply(1:3, function(.i) {
      rbind(data.frame(ID=.i, TIME=0, DV=NA_real_, AMT=100, EVID=1, CMT=1, II=12, ADDL=2),
            data.frame(ID=.i, TIME=seq(1, 36, by=4), DV=1, AMT=0, EVID=0, CMT=2,
                       II=0, ADDL=0))
    }))
    .ctl <- foceiControl(print=0, maxOuterIterations=0, maxInnerIterations=0,
                         covMethod="", calcTables=TRUE)
    .n <- vapply(c(TRUE, FALSE), function(.u) {
      vapply(c(FALSE, TRUE), function(.a) {
        nrow(as.data.frame(nlmixr2(.mkAddl(.u), .d, est="focei", control=.ctl,
                                   table=tableControl(addDosing=.a))))
      }, numeric(1))
    }, numeric(2))
    # same shape with and without the mtime: 27 observations, plus 9 doses
    # (3 subjects x 1 dose + 2 ADDL) when addDosing is on
    expect_equal(.n[, 1], .n[, 2])
    expect_equal(.n[, 1], c(27, 36))
  })

  test_that("saem fits a model with mtime()", {
    .fit <- nlmixr2(.mkMtime(TRUE), .theo, est="saem",
                    control=saemControl(print=0, nBurn=5, nEm=5, seed=42))
    expect_true(inherits(.fit, "nlmixr2FitData"))
    expect_true(is.finite(.fit$objf))
    .df <- as.data.frame(.fit)
    expect_equal(nrow(.df), sum(.theo$EVID == 0))
    expect_false(anyNA(.df$IPRED))
  })
})
