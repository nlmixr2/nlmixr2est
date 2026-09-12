nmTest({

  ## Cheap, fit-free checks of the vae mixture plumbing.  Deliberately NOT in
  ## .slowBatches: these are the ones that must run on every push, because each
  ## of them guards a defect that produced a silently wrong number rather than
  ## an error.

  .mixMod <- function() {
    ini({
      tka <- log(1.5)
      tcl1 <- log(1.0)
      tcl2 <- log(5.0)
      p1 <- 0.3
      tv <- log(20)
      eta.cl ~ 0.01
      eta.v ~ 0.01
      add.sd <- 0.05
    })
    model({
      ka <- exp(tka)
      cl <- mix(exp(tcl1 + eta.cl), p1, exp(tcl2 + eta.cl))
      v <- exp(tv + eta.v)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }

  test_that("vae prep carries the mixture proportion on the estimation scale", {
    ## The inner problem reads that theta slot through .getMixFromLog (mexpit),
    ## so prep$th has to hold the MLOGIT value.  Handing it the raw ini()
    ## probability meant p1 = 0.3 was read as mexpit(0.3) = 0.574 -- a wrong
    ## number with no error anywhere.
    .ui <- rxode2::assertRxUi(.mixMod())
    .prep <- .vaeDataPrep(.ui, nlmixr2data::theo_sd)
    expect_equal(unname(.prep$th[.ui$thetaMixIndex]), rxode2::mlogit(0.3))
    ## and it round trips back to the simplex the combine step needs
    expect_equal(.getMixFromLog(.prep$th, .ui$thetaMixIndex), c(0.3, 0.7))
  })

  test_that("a mixture proportion is not a vae regress theta", {
    ## p1 appears in mix(a, p1, b) so it survived the "is it in a model
    ## expression" filter and became a bobyqa regress theta, moved against the
    ## (-Inf, Inf) it carries in iniDf -- producing e.g. p1 <- 10.6 in ini().
    ## npag excludes them for the same reason (.npMuExpand).
    .ui <- rxode2::assertRxUi(.mixMod())
    expect_false("p1" %in% .vaeNonMuThetas(.ui))
    .prep <- .vaeDataPrep(.ui, nlmixr2data::theo_sd)
    expect_false("p1" %in% .prep$regressNames)
    ## the ordinary non-mu structural thetas are still candidates
    expect_true(all(c("tka", "tcl1", "tcl2") %in% .vaeNonMuThetas(.ui)))
  })

  test_that("vae writes the fitted mixture proportion back as a probability", {
    .ui <- rxode2::assertRxUi(.mixMod())
    .setIni <- function(u, expr) do.call(rxode2::ini, list(u, str2lang(expr)))
    .ui2 <- .vaeSetIniMixProb(.ui, .ui, list(mixProb = c(0.42, 0.58)), .setIni)
    expect_equal(unname(.ui2$iniDf$est[.ui2$iniDf$name == "p1"]), 0.42)
    ## and the result is a model the ui itself accepts
    expect_error(rxode2::assertRxUi(.ui2$fun), NA)

    ## a proportion that reached the boundary is clamped rather than written as
    ## a 0/1 that ini() validation rejects
    .ui3 <- .vaeSetIniMixProb(.ui, .ui, list(mixProb = c(1, 0)), .setIni)
    .p3 <- unname(.ui3$iniDf$est[.ui3$iniDf$name == "p1"])
    expect_true(.p3 > 0 && .p3 < 1)
    expect_error(rxode2::assertRxUi(.ui3$fun), NA)

    ## a non-mixture fit, or a malformed mixProb, leaves the ui alone
    expect_equal(.vaeSetIniMixProb(.ui, .ui, list(mixProb = NULL), .setIni)$iniDf$est,
                 .ui$iniDf$est)
    expect_equal(.vaeSetIniMixProb(.ui, .ui, list(mixProb = c(0.5, 0.3, 0.2)), .setIni)$iniDf$est,
                 .ui$iniDf$est)
  })

})
