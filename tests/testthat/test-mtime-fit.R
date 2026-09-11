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
    expect_equal(grep("^mtime", strsplit(.ui$saemModel, "\n")[[1]], value=TRUE),
                 "mtime(tsw5)~exp(tsw)")

    # a covariate stays a covariate
    .ui <- rxode2::rxUiDecompress(.mk(c(.base, "mtime(tw) <- WT / 10",
      "kmult <- ifelse(t < tw, 1.0, 2.0)", .tail))())
    expect_equal(.mtimeLines(.ui$focei$inner), "mtime(tw)~0.1*WT;")

    # rxS() keeps only a variable's FINAL value, so a right hand side whose
    # dependency is assigned again after the declaration would expand to the
    # later value -- rxode2 evaluates the declaration in place, so refuse
    .ui <- rxode2::rxUiDecompress(.mk(c(.base, "t1 <- 10", "mtime(tx) <- t1",
      "t1 <- 20", "kmult <- ifelse(t < tx, 1.0, 2.0) * t1 / 20", .tail))())
    expect_error(.ui$focei, "assigns again after it")
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
})
