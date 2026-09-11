nmTest({
  # mtime() in a model used to break every fit: etTrans() materializes the
  # modeled times as EVID 10-99 records (TIME=0, AMT=NA) and $dataSav kept them,
  # so the estimation solves saw them as input doses with a missing amt.  On top
  # of that rxode2::rxS() drops the mtime() assignment, so the generated models
  # lost the modeled times and left the mtime variable undefined (#919).
  # The fits that exercise the other estimation routines live in
  # test-mtime-fit.R (a weekly batch).
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

  test_that(".addMtimeLines() splices after the leading declarations", {
    .s <- new.env(parent=emptyenv())
    .s$..mtime <- "mtime(t5)~5"
    expect_equal(.addMtimeLines("param(a)\ncmt(x)\nd/dt(x)=-a*x\ncmt(cp)", .s),
                 "param(a)\ncmt(x)\nmtime(t5)~5\nd/dt(x)=-a*x\ncmt(cp)")
    # no leading declaration at all
    expect_equal(.addMtimeLines("d/dt(x)=-a*x", .s), "mtime(t5)~5\nd/dt(x)=-a*x")
    # nothing to add
    .s$..mtime <- character(0)
    expect_equal(.addMtimeLines("d/dt(x)=-a*x", .s), "d/dt(x)=-a*x")
  })

  test_that(".rxMtimeToAssign() loads the declaration as a suppressed assignment", {
    # no mtime: the text is handed to rxS() untouched
    expect_equal(.rxMtimeToAssign("d/dt(depot)=-ka*depot;\n"),
                 "d/dt(depot)=-ka*depot;\n")
    expect_equal(.rxMtimeToAssign("mtime(t5)=5;\nd/dt(x)=-a*x;"),
                 "t5~5;\nd/dt(x)=-a*x;")
    # every declaration form, and `~` so no generated model gains an output column
    expect_equal(.rxMtimeToAssign("mtime(t5)~5;\nmtime(tx)=2*THETA[1];"),
                 "t5~5;\ntx~2*THETA[1];")
  })

  test_that("a modeled time that moves with a parameter is differentiated", {
    # A boundary that moves with an estimated parameter contributes to the
    # sensitivities.  rxS() drops the mtime() right hand side, so the switch time
    # used to reach symengine as a free symbol and every derivative through it was
    # zero -- silently, and only for the mtime() spelling.  The reference is the
    # same switch written out in place, which has always been differentiated.
    .mkSw <- function(useMtime) {
      .bdy <- c("ka <- exp(tka)", "cl <- exp(tcl)", "v <- exp(tv)",
                if (useMtime) c("mtime(tsw5) <- exp(tsw + eta.sw)",
                                "kmult <- ifelse(t < tsw5, 1.0, 2.0)")
                else "kmult <- ifelse(t < exp(tsw + eta.sw), 1.0, 2.0)",
                "d/dt(depot) <- -ka * depot",
                "d/dt(center) <- ka * depot - kmult * cl / v * center",
                "cp <- center / v", "cp ~ add(add.sd)")
      eval(parse(text=paste0(
        "function() {\n ini({tka <- 0.45; tcl <- -3.2; tv <- -1; tsw <- 1.386;",
        " eta.sw ~ 0.1; add.sd <- 0.7})\n model({\n",
        paste(.bdy, collapse="\n"), "\n })\n}")))
    }
    .ddt <- function(useMtime) {
      .s <- rxode2::rxUiDecompress(.mkSw(useMtime)())$loadPruneSens
      get("rx__d_dt_center__", .s)
    }
    .mt <- .ddt(TRUE)
    # the loaded equation is the in-place one: the switch time is its expansion,
    # not an opaque name
    expect_equal(paste(.mt), paste(.ddt(FALSE)))
    expect_match(paste(.mt), "rxLt(t, exp(ETA_1_ + THETA_4_))", fixed=TRUE)
    # ...so the eta reaches the branch and the derivative is not identically zero
    expect_false(paste(symengine::D(.mt, symengine::S("ETA_1_"))) == "0")
  })

  test_that(".rxMtimeDeps() walks back through the preceding assignments", {
    .lines <- c("a=1;", "b=a+2;", "mtime(tx)=b;", "c=3;")
    .lhs <- .rxLineLhs(.lines)
    expect_equal(.lhs, c("a", "b", NA, "c"))
    # tx depends on b, and b on a -- but not on c, which comes after it
    expect_setequal(.rxMtimeDeps(.lines, .lhs, 3L, "b"), c("a", "b"))
    # a name with no preceding assignment is a leaf (a parameter or covariate)
    expect_setequal(.rxMtimeDeps(.lines, .lhs, 3L, "WT/10"), "WT")
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
    expect_equal(rxode2::rxModelVars(.ui$nlmRxModel$predOnly)$nMtime, 1L)
    # it is emitted suppressed, so it adds no output column that would shift
    # the positional lhs layout inner.cpp reads
    expect_false("t5" %in% rxode2::rxModelVars(.f$inner)$lhs)
    # nlme off a FRESH ui: building the focei bundle caches focei's own predDf
    # on the ui, which the nlme builder does not accept
    .ui2 <- rxode2::rxUiDecompress(.mkMtime(TRUE)())
    expect_equal(rxode2::rxModelVars(.ui2$nlmeRxModel)$nMtime, 1L)
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
})
