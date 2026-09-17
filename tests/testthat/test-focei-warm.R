nmTest({
  test_that("foceiControl(warm=) option mapping", {
    expect_equal(foceiControl()$warm, 1L)
    expect_equal(foceiControl(warm="calc")$warm, 1L)
    expect_equal(foceiControl(warm="save")$warm, 0L)
    expect_equal(foceiControl(warm="none")$warm, 2L)
    expect_equal(foceiControl(warm=0L)$warm, 0L)
    expect_equal(foceiControl(warm=1L)$warm, 1L)
    expect_equal(foceiControl(warm=2L)$warm, 2L)
    expect_error(foceiControl(warm="bogus"))
    expect_error(foceiControl(warm=3L))
    for (.w in c("calc", "save", "none")) {
      .ctl <- foceiControl(warm=.w)
      expect_equal(do.call(foceiControl, .ctl)$warm, .ctl$warm, info=.w)
    }
  })

  one.cmt <- function() {
    ini({
      tka <- 0.45
      tcl <- log(c(0, 2.7, 100))
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

  test_that("focei warm='calc' matches warm='save'", {
    # innerOpt is pinned to n1qn1 on purpose: warm= seeds the n1qn1 inner Hessian
    # and the default innerOpt="auto" resolves to "trust", which builds its own
    # exact Hessian at every trial point and ignores warm= entirely -- comparing
    # the arms under "auto" compares a fit against itself.
    .fit <- function(warm, maxOuterIterations) {
      suppressWarnings(suppressMessages(
        nlmixr(one.cmt, nlmixr2data::theo_sd, "focei",
               foceiControl(maxOuterIterations=maxOuterIterations,
                            covMethod="", calcTables=FALSE, print=0,
                            innerOpt="n1qn1", warm=warm))))
    }

    # posthoc: same inner problems converged to the same etas/objective
    .p1 <- .fit("calc", 0L)
    .p2 <- .fit("save", 0L)
    .p3 <- .fit("none", 0L)
    expect_equal(.p1$objf, .p2$objf, tolerance=1e-4)
    expect_equal(.p1$objf, .p3$objf, tolerance=1e-4)
    expect_equal(as.data.frame(.p1$eta), as.data.frame(.p2$eta),
                 tolerance=1e-4)
    expect_equal(as.data.frame(.p1$eta), as.data.frame(.p3$eta),
                 tolerance=1e-4)

    # short optimization run finishes and agrees
    .f1 <- .fit("calc", 5L)
    .f2 <- .fit("save", 5L)
    .f3 <- .fit("none", 5L)
    expect_true(inherits(.f1, "nlmixr2FitCore"))
    expect_true(inherits(.f2, "nlmixr2FitCore"))
    expect_true(inherits(.f3, "nlmixr2FitCore"))
    expect_equal(.f1$objf, .f2$objf, tolerance=1e-2)
    expect_equal(.f1$objf, .f3$objf, tolerance=1e-2)
  })

  test_that("warm='save' actually restarts from the saved curvature (#1043)", {
    # updateZm() zeroed zm before reading it back, so mode=2 always got an
    # all-zero factorization and n1qn1 self-initialized on every inner solve:
    # warm="save" reused nothing for the whole life of the option.  Assert the
    # MECHANISM, not just that the numbers are reasonable -- the old code
    # produced perfectly reasonable numbers.
    .fit <- function(warm) {
      suppressWarnings(suppressMessages(
        nlmixr(one.cmt, nlmixr2data::theo_sd, "focei",
               foceiControl(maxOuterIterations=5L, covMethod="", calcTables=FALSE,
                            print=0, innerOpt="n1qn1", warm=warm))))
    }
    .save <- .fit("save")
    .ws <- .save$env$nWarmSave
    expect_false(is.null(.ws))
    # every inner solve after the first per subject reuses the previous one
    expect_gt(.ws[["reused"]], 0L)
    # the reconstructed curvature is a usable Hessian, so no solve falls back
    expect_equal(.ws[["selfInit"]], 0L)
    # ... and the counter only exists for warm="save"
    expect_null(.fit("none")$env$nWarmSave)
    expect_null(.fit("calc")$env$nWarmSave)
  })

  test_that("warm='save' and warm='none' are different runs (#1043)", {
    # Before the fix these two were bit-identical -- "save" fell through to
    # n1qn1's self-init.  They must now be genuinely different inner solves
    # (while still converging to the same neighbourhood).
    .fit <- function(warm) {
      suppressWarnings(suppressMessages(
        nlmixr(one.cmt, nlmixr2data::theo_sd, "focei",
               foceiControl(maxOuterIterations=5L, covMethod="", calcTables=FALSE,
                            print=0, innerOpt="n1qn1", warm=warm))))
    }
    .s <- .fit("save")
    .n <- .fit("none")
    expect_false(isTRUE(all.equal(.s$objf, .n$objf, tolerance=1e-10)))
    expect_equal(.s$objf, .n$objf, tolerance=1e-2)
  })

  test_that("the mceta eta=0 floor pass keeps the warm='save' seed (#1043)", {
    # mceta>=1 runs a second, eta=0 "floor" pass whose job is to be exactly the
    # run mceta=0 would have made.  n1qn1 overwrites zm in place, so that pass
    # has to be handed the seed the first pass got back; it used to self-init
    # while the mceta=0 run it must reproduce got the saved curvature.
    # $nWarmSave["floorReseed"] counts the floor passes that really got a
    # mode=2 seed.
    .f <- suppressWarnings(suppressMessages(
      nlmixr(one.cmt, nlmixr2data::theo_sd, "focei",
             foceiControl(maxOuterIterations=5L, covMethod="", calcTables=FALSE,
                          print=0, innerOpt="n1qn1", warm="save", mceta=5L))))
    expect_true(is.finite(.f$objf))
    # a sampled eta beat eta=0 somewhere, so the floor pass ran at all
    expect_gt(.f$env$nMcetaStart[["sample"]], 0L)
    expect_gt(.f$env$nWarmSave[["floorReseed"]], 0L)
  })

  test_that("warm='save' seeds a single-eta model too (#1043)", {
    # The reconstruction also had `if (n == 1) H = D` with D still zeroed, so a
    # one-eta model could not have been seeded even with the fill removed.
    .m <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        add.sd <- 0.7
        eta.cl ~ 0.3
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv)
        linCmt() ~ add(add.sd)
      })
    }
    .f <- suppressWarnings(suppressMessages(
      nlmixr(.m, nlmixr2data::theo_sd, "focei",
             foceiControl(maxOuterIterations=5L, covMethod="", calcTables=FALSE,
                          print=0, innerOpt="n1qn1", warm="save"))))
    expect_true(is.finite(.f$objf))
    expect_gt(.f$env$nWarmSave[["reused"]], 0L)
    expect_equal(.f$env$nWarmSave[["selfInit"]], 0L)
  })
})
