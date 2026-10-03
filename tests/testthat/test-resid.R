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
