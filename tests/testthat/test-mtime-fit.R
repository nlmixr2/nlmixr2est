nmTest({
  # The rest of the #919 coverage: one fit per estimation routine that builds
  # its own rxode2 model, plus the table-shape cases.  test-mtime.R holds the
  # cheap half that runs on every push.
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
  .nObs <- sum(.theo$EVID == 0)

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
    .ui <- rxode2::rxUiDecompress(.mk(c(.base, "mtime(tsw5) <- exp(tsw)",
      "kmult <- ifelse(t < tsw5, 1.0, 2.0)", .tail), useTsw=TRUE)())
    expect_equal(grep("^mtime", strsplit(.ui$saemModel, "\n")[[1]], value=TRUE),
                 "mtime(tsw5)~exp(tsw)")
    # saem's predOnly (residuals/tables) body is in the NATURAL names, which
    # only the mu-reference replacement block defines, so the declaration has
    # to come after it -- not with the other prologue lines
    .ui <- rxode2::rxUiDecompress(.mk(c(.base, "mtime(tsw5) <- exp(tsw)",
      "kmult <- ifelse(t < tsw5, 1.0, 2.0)", .tail), useTsw=TRUE)())
    .txt <- strsplit(rxode2::rxModelVars(.ui$saemModelPred$predOnly)$model["normModel"],
                     "\n")[[1]]
    expect_true(grep("^mtime\\(tsw5\\)", .txt) > grep("^tsw=", .txt))

    # a covariate stays a covariate
    .ui <- rxode2::rxUiDecompress(.mk(c(.base, "mtime(tw) <- WT / 10",
      "kmult <- ifelse(t < tw, 1.0, 2.0)", .tail))())
    expect_equal(.mtimeLines(.ui$focei$inner), "mtime(tw)~0.1*WT;")

    # a branch in the right hand side: .foceiPrune() rewrites ifelse() to
    # arithmetic before symengine sees it, so it expands like anything else
    # and agrees with what plain rxode2 computes in place
    .ui <- rxode2::rxUiDecompress(.mk(c(.base, "kbase <- 5",
      "mtime(tb) <- ifelse(kbase > 1, 3, 4)",
      "kmult <- ifelse(t < tb, 1.0, 2.0)", .tail))())
    expect_equal(.mtimeLines(.ui$focei$inner), "mtime(tb)~3;")
    .ref <- rxode2::rxode2("kbase=5;mtime(tb)=ifelse(kbase>1,3,4);d/dt(x)=-x;")
    expect_equal(unique(as.data.frame(
      rxode2::rxSolve(.ref, rxode2::et(amt=1) |> rxode2::et(0:5)))$tb), 3)

    # rxS() keeps only a variable's FINAL value, so a right hand side whose
    # dependency is assigned again after the declaration would expand to the
    # later value -- rxode2 evaluates the declaration in place, so refuse
    .ui <- rxode2::rxUiDecompress(.mk(c(.base, "t1 <- 10", "mtime(tx) <- t1",
      "t1 <- 20", "kmult <- ifelse(t < tx, 1.0, 2.0) * t1 / 20", .tail))())
    expect_error(.ui$focei, "assigns again after it")
  })

  test_that("an emitted mtime() line resolves in every generated model", {
    # Structural invariant: every name in an emitted mtime() line must be a
    # declared parameter of that model, or assigned on an EARLIER line of it.
    # The declaration is re-emitted near the top, so a right hand side naming
    # something the model only defines further down would read it unset.
    .resolves <- function(txt, label) {
      .lines <- strsplit(paste(txt, collapse="\n"), "\n")[[1]]
      .w <- grep("^mtime", .lines)
      expect_true(length(.w) > 0L, info=paste(label, "emits an mtime"))
      .pw <- grep("^params?\\(", .lines)
      .pars <- if (length(.pw) > 0L) {
        all.vars(str2lang(sub(";$", "", .lines[.pw[1]])))
      } else character(0)
      .lhs <- .rxLineLhs(.lines)
      for (.i in .w) {
        .rhs <- sub(";$", "", sub("^mtime\\([^)]*\\)[=~]", "", .lines[.i]))
        for (.v in setdiff(all.vars(str2lang(.rhs)), c("t", "time"))) {
          .ok <- .v %in% .pars ||
            any(!is.na(.lhs) & .lhs == .v & seq_along(.lhs) < .i) ||
            any(grepl(paste0("^mtime\\(", .v, "\\)"), .lines[seq_len(.i - 1L)]))
          expect_true(.ok, info=paste0(label, ": '", .v, "' in ", .lines[.i]))
        }
      }
    }
    .mkTheta <- function() {
      eval(parse(text=paste0(
        "function() {\n ini({tka <- 0.45; tcl <- -3.2; tv <- -1; tsw <- 1.6;",
        " eta.ka ~ 0.1; add.sd <- 0.7})\n model({\n",
        paste(c("ka <- exp(tka + eta.ka)", "cl <- exp(tcl)", "v <- exp(tv)",
                "mtime(tsw5) <- exp(tsw)",
                "kmult <- ifelse(t < tsw5, 1.0, 2.0)",
                "d/dt(depot) <- -ka * depot",
                "d/dt(center) <- ka * depot - kmult * cl / v * center",
                "cp <- center / v", "cp ~ add(add.sd)"), collapse="\n"),
        "\n })\n}")))
    }
    .norm <- function(m) rxode2::rxModelVars(m)$model["normModel"]
    # a FRESH ui per builder: building one bundle caches its own predDf on the ui
    .f <- rxode2::rxUiDecompress(.mkTheta()())$focei
    for (.m in c("inner", "predOnly", "predNoLhs")) {
      .resolves(.norm(.f[[.m]]), paste0("focei$", .m))
    }
    .resolves(rxode2::rxUiDecompress(.mkTheta()())$saemModel, "saemModel")
    .resolves(.norm(rxode2::rxUiDecompress(.mkTheta()())$saemModelPred$predOnly),
              "saemModelPred")
    .resolves(.norm(rxode2::rxUiDecompress(.mkTheta()())$nlmRxModel$predOnly),
              "nlmRxModel")
    .resolves(.norm(rxode2::rxUiDecompress(.mkTheta()())$nlsRxModel$predOnly),
              "nlsRxModel")
    .resolves(.norm(rxode2::rxUiDecompress(.mkTheta()())$nlmeRxModel), "nlmeRxModel")
    # the impmap theta-sensitivity model builds off a LIGHTWEIGHT list, not the
    # symengine environment (that is deliberate, it used to leak), so the mtime
    # lines have to be carried in that list
    .uiI <- rxode2::rxUiDecompress(.mkTheta()())
    rxode2::rxAssignControlValue(.uiI, "combSens", TRUE)
    .ts <- .impmapThetaSensModel(.uiI)
    expect_false(is.null(.ts))
    .resolves(.norm(.ts), "impmapThetaSens")
  })

  test_that("mtime() does not cost the table its ADDL doses", {
    # the extra mtime rows are dropped from a table solve; an ADDL-expanded dose
    # also has no source data row and must NOT be dropped with them
    .d <- do.call(rbind, lapply(1:3, function(.i) {
      rbind(data.frame(ID=.i, TIME=0, DV=NA_real_, AMT=100, EVID=1, CMT=1, II=12, ADDL=2),
            data.frame(ID=.i, TIME=seq(1, 36, by=4), DV=1, AMT=0, EVID=0, CMT=2,
                       II=0, ADDL=0))
    }))
    .ctl <- foceiControl(print=0, maxOuterIterations=0, maxInnerIterations=0,
                         covMethod="", calcTables=TRUE)
    .n <- vapply(c(TRUE, FALSE), function(.u) {
      vapply(c(FALSE, TRUE), function(.a) {
        nrow(as.data.frame(nlmixr2(.mkMtime(.u), .d, est="focei", control=.ctl,
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
    expect_equal(nrow(.df), .nObs)
    expect_false(anyNA(.df$IPRED))
  })

  test_that("nlm fits a model with mtime()", {
    # .rxFinalizeNlm() assembles and OPTIMIZES its own model text; rxOptExpr()
    # cannot parse mtime(), so the lines have to be spliced in afterwards
    .pop <- function() {
      ini({tka <- 0.45; tcl <- -3.2; tv <- -1; add.sd <- 0.7})
      model({
        ka <- exp(tka)
        cl <- exp(tcl)
        v <- exp(tv)
        mtime(t5) <- 5
        kmult <- ifelse(t < t5, 1.0, 2.0)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - kmult * cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }
    .fit <- nlmixr2(.pop, .theo, est="nlm", control=nlmControl(print=0))
    expect_true(inherits(.fit, "nlmixr2FitData"))
    expect_equal(nrow(as.data.frame(.fit)), .nObs)
  })

  test_that("nlme fits a model with mtime()", {
    # the nlme objective returns one prediction per row of the solve, so the
    # extra mtime rows made the returned vector longer than nlme's response
    .mk <- function(useMtime) {
      .bdy <- c("ka <- exp(tka + eta.ka)", "cl <- exp(tcl + eta.cl)",
                "v <- exp(tv + eta.v)",
                if (useMtime) "mtime(t5) <- 5",
                "d/dt(depot) = -ka * depot",
                "d/dt(center) = ka * depot - cl / v * center",
                "cp = center / v", "cp ~ add(add.sd)")
      eval(parse(text=paste0(
        "function() {\n ini({tka <- 0.45; tcl <- 1; tv <- 3.45; eta.ka ~ 0.6;",
        " eta.cl ~ 0.3; eta.v ~ 0.1; add.sd <- 0.7})\n",
        " model({\n", paste(.bdy, collapse="\n"), "\n })\n}")))
    }
    .ctl <- nlmeControl(maxIter=5, verbose=FALSE, returnObject=TRUE)
    .fit <- nlmixr2(.mk(TRUE), .theo, est="nlme", control=.ctl)
    expect_true(inherits(.fit, "nlmixr2FitData"))
    expect_equal(nrow(as.data.frame(.fit)), .nObs)
    # and the modeled time has not moved the answer much off the plain model
    .ref <- nlmixr2(.mk(FALSE), .theo, est="nlme", control=.ctl)
    expect_equal(.fit$objf, .ref$objf, tolerance=1e-2)
  })

  test_that("a moving boundary: mtime() sensitivities match the in-place switch", {
    skip_on_cran()
    # An mtime() whose time depends on an estimated parameter is a moving
    # discontinuity in the right hand side.  rxS() drops the declaration, so the
    # switch time reached symengine as a free symbol and the boundary term was
    # silently lost -- the eta sensitivity came back identically zero and every
    # EBE for that eta stayed pinned at its initial value.  The reference is the
    # same switch written in place, which has always carried the term.
    .mkSw <- function(useMtime) {
      .bdy <- c("ka <- exp(tka)", "cl <- exp(tcl)", "v <- exp(tv)",
                if (useMtime) c("mtime(tsw5) <- exp(tsw + eta.sw)",
                                "kmult <- ifelse(t < tsw5, 1.0, 2.0)")
                else "kmult <- ifelse(t < exp(tsw + eta.sw), 1.0, 2.0)",
                "d/dt(depot) <- -ka * depot",
                "d/dt(center) <- ka * depot - kmult * cl / v * center",
                "cp <- center / v", "cp ~ add(add.sd)")
      eval(parse(text=paste0(
        "function() {\n ini({tka <- 0.45; tcl <- -3.2; tv <- -1; tsw <- ", log(4),
        "; eta.sw ~ 0.2; add.sd <- 0.7})\n model({\n",
        paste(.bdy, collapse="\n"), "\n })\n}")))
    }
    .grid <- seq(0, 12, by=0.25)
    .ev <- rxode2::et(rxode2::et(amt=4.02, cmt="depot"), .grid)
    .th <- c(0.45, -3.2, -1, log(4), 0.7)
    # d(pred)/d(eta) off the generated inner model, against central differences
    # of the same model's prediction
    .sens <- function(useMtime, h=1e-4) {
      .inner <- rxode2::rxUiDecompress(.mkSw(useMtime)())$focei$inner
      # loading the declaration as an assignment must not turn the modeled time
      # into an output column: that would shift the positional lhs layout
      # inner.cpp reads
      if (useMtime) expect_false("tsw5" %in% rxode2::rxModelVars(.inner)$lhs)
      .at <- function(.e) {
        .p <- stats::setNames(c(.th, .e),
                              c(paste0("THETA[", 1:5, "]"), "ETA[1]"))
        .s <- suppressWarnings(
          rxode2::rxSolve(.inner, .p, .ev, returnType="data.frame",
                          atol=1e-10, rtol=1e-10))
        # the mtime record is an extra row whose time MOVES with the eta, so
        # keep the requested grid to compare like with like
        .s <- .s[.s$time %in% .grid, ]
        .s[!duplicated(.s$time, fromLast=TRUE), ]
      }
      .s0 <- .at(0); .sp <- .at(h); .sm <- .at(-h)
      data.frame(time=.s0$time,
                 analytic=.s0[["rx__sens_rx_pred__BY_ETA_1___"]],
                 fd=(.sp$rx_pred_ - .sm$rx_pred_) / (2 * h))
    }
    .mt <- .sens(TRUE)
    .ref <- .sens(FALSE)
    # the boundary term is there at all
    expect_gt(max(abs(.mt$analytic)), 1)
    # ...and is the same term the in-place switch gets
    expect_equal(.mt$analytic, .ref$analytic, tolerance=1e-8)
    # ...and it is the right term: agrees with central differences away from the
    # smoothing window symengine puts around the branch
    .w <- .mt$time >= 4 + 0.25          # the switch is at exp(tsw) = 4
    expect_equal(.mt$analytic[.w], .mt$fd[.w], tolerance=0.01)
  })

  test_that("a moving boundary: the EBEs match the in-place switch", {
    skip_on_cran()
    .mkSw <- function(useMtime) {
      .bdy <- c("ka <- exp(tka + eta.ka)", "cl <- exp(tcl)", "v <- exp(tv)",
                if (useMtime) c("mtime(tsw5) <- exp(tsw + eta.sw)",
                                "kmult <- ifelse(t < tsw5, 1.0, 2.0)")
                else "kmult <- ifelse(t < exp(tsw + eta.sw), 1.0, 2.0)",
                "d/dt(depot) <- -ka * depot",
                "d/dt(center) <- ka * depot - kmult * cl / v * center",
                "cp <- center / v", "cp ~ add(add.sd)")
      eval(parse(text=paste0(
        "function() {\n ini({tka <- 0.45; tcl <- -3.2; tv <- -1; tsw <- ", log(4),
        "; eta.ka ~ 0.1; eta.sw ~ 0.2; add.sd <- 0.7})\n model({\n",
        paste(.bdy, collapse="\n"), "\n })\n}")))
    }
    .ctl <- foceiControl(print=0L, covMethod="", sigdig=4, calcTables=FALSE,
                         maxOuterIterations=0L, maxInnerIterations=300L)
    .fit <- nlmixr2(.mkSw(TRUE), .theo, est="focei", control=.ctl)
    .ref <- nlmixr2(.mkSw(FALSE), .theo, est="focei", control=.ctl)
    # the inner problem can move the eta the boundary depends on: a zero
    # sensitivity leaves every one of these at 0
    expect_gt(max(abs(.fit$eta$eta.sw)), 0.05)
    expect_equal(.fit$eta$eta.sw, .ref$eta$eta.sw, tolerance=0.01)
    expect_equal(.fit$objf, .ref$objf, tolerance=1e-5)
  })
})
