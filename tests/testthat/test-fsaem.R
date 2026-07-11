nmTest({
  test_that("saemControl fast* options validate and round-trip", {
    .c0 <- saemControl()
    expect_false(.c0$fast)
    expect_equal(.c0$fastKernel, "firstN")
    expect_equal(.c0$fastCov, "auto")
    expect_equal(.c0$fastIter, 20L)
    expect_equal(.c0$fastLik, "focei")

    .c1 <- saemControl(fast=TRUE, fastKernel="throughout", fastCov="hessian",
                       fastIter=5, fastLik="foce")
    expect_true(.c1$fast)
    expect_equal(.c1$fastKernel, "throughout")
    expect_equal(.c1$fastCov, "hessian")
    expect_equal(.c1$fastIter, 5L)
    expect_equal(.c1$fastLik, "foce")

    # feeding a control back through saemControl() (as getValidNlmixrCtl does)
    # preserves the fast settings
    .c1b <- do.call(saemControl, .c1)
    expect_true(.c1b$fast)
    expect_equal(.c1b$fastKernel, "throughout")
    expect_equal(.c1b$fastCov, "hessian")
    expect_equal(.c1b$fastIter, 5L)
    expect_equal(.c1b$fastLik, "foce")

    expect_error(saemControl(fastIter=0), "fastIter")
  })

  test_that("fsaemControl forces fast=TRUE and validates as an fsaem control", {
    .fc <- fsaemControl(fast=FALSE, nBurn=7)
    expect_s3_class(.fc, "fsaemControl")
    expect_true(.fc$fast)
    expect_equal(.fc$mcmc$niter[1], 7)

    .gv <- getValidNlmixrControl(NULL, "fsaem")
    expect_s3_class(.gv, "fsaemControl")
    expect_true(.gv$fast)
  })

  test_that("est='fsaem' runs the live IMH kernel and converges to the SAEM MLE", {
    one.cmt <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v  <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    .ctl <- saemControl(nBurn=200, nEm=100, nmc=3, seed=42, print=0L, calcTables=FALSE)
    .fs <- suppressMessages(nlmixr2(one.cmt, nlmixr2data::theo_sd, est="fsaem", control=.ctl))
    .ss <- suppressMessages(nlmixr2(one.cmt, nlmixr2data::theo_sd, est="saem", control=.ctl))
    # fast flag flows through to the stored control (mechanism is wired up)
    expect_true(.fs$saemControl$fast)
    expect_false(.ss$saemControl$fast)
    # the fast kernel changes the simulation trajectory (it fired) -- so fsaem is
    # NOT bit-identical to saem, but converges to the same MLE
    expect_false(isTRUE(all.equal(unname(fixef(.fs)), unname(fixef(.ss)))))
    expect_lt(max(abs(fixef(.fs) - fixef(.ss))), 0.05)
  })
})
