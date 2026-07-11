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

  test_that("est='fsaem' converges faster than saem from poor starts (the point)", {
    # deliberately poor initial estimates; with a short burn-in the f-SAEM IMH
    # kernel should land closer to the MLE than the standard random-walk SAEM.
    poor <- function() {
      ini({
        tka <- 0.1; tcl <- 0.5; tv <- 3.0
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    .d <- nlmixr2data::theo_sd
    .ref <- fixef(suppressMessages(nlmixr2(poor, .d, est = "saem",
      control = saemControl(nBurn = 400, nEm = 150, nmc = 3, seed = 7, print = 0L, calcTables = FALSE))))
    .short <- function(est, nb) {
      .f <- suppressMessages(nlmixr2(poor, .d, est = est,
        control = saemControl(nBurn = nb, nEm = 20, nmc = 3, seed = 7, print = 0L,
                              calcTables = FALSE, fastIter = nb)))
      sqrt(sum((fixef(.f) - .ref)^2))
    }
    expect_lt(.short("fsaem", 20L), .short("saem", 20L))
  })

  test_that("est='fsaem' degrades to standard SAEM for unsupported models", {
    # proportional error is outside the current fast-kernel envelope, so fsaem
    # must fall back to standard SAEM bit-for-bit (with a message).
    prop <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
        prop.sd <- 0.2
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        linCmt() ~ prop(prop.sd)
      })
    }
    .ctl <- saemControl(nBurn = 40, nEm = 30, nmc = 2, seed = 3, print = 0L, calcTables = FALSE)
    expect_message(
      .f <- nlmixr2(prop, nlmixr2data::theo_sd, est = "fsaem", control = .ctl),
      "not yet supported")
    .s <- suppressMessages(nlmixr2(prop, nlmixr2data::theo_sd, est = "saem", control = .ctl))
    expect_equal(unname(fixef(.f)), unname(fixef(.s)))
  })
})
