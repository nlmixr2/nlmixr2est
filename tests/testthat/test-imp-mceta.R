# `mceta` sets how the FOCEi inner MAP picks its starting eta, and imp reaches
# that same inner problem -- so `impmapControl(mceta=)` is a real lever on an
# imp fit.  Two things made it look otherwise, and both are asserted here.
#
# 1. `fit$env$nMcetaStart` -- the only evidence the omega draws were explored --
#    was stamped inside the non-EM branch of foceiFinalizeTables().  est="imp"
#    takes the isImpmap branch, so the counter was absent from every imp fit
#    whatever mceta was set to.  Absence was not evidence of inertness, which
#    is exactly how a genuinely inert setting would also look.
# 2. `.impmapControlToFoceiControl()` stripped mceta as though it were an
#    IS-only control name.  Estimation reads mceta off the impmapControl
#    directly so it was unaffected, but the consumers of the down-converted
#    control saw foceiControl()'s default instead of the user's value.
#
# The estimates are no evidence either way: the EM is stochastic, so mceta=0 and
# mceta=10 differ a little whether or not the draws were ever explored.
nmTest({
  test_that("impmapControl(mceta=) survives the down-conversion to foceiControl", {
    # A pure control-level assertion, so it fails on the mapping rather than on
    # anything a fit happens to do.
    expect_false("mceta" %in% nlmixr2est:::.impmapIsControlNames)
    .e <- new.env(parent = emptyenv())
    .e$impmapControl <- impmapControl(mceta = 7L)
    .fc <- nlmixr2est:::.impmapControlToFoceiControl(.e, assign = FALSE)
    expect_true("mceta" %in% names(.fc))
    expect_equal(.fc$mceta, 7L)
    # and the down-converted list is still something foceiControl() accepts
    expect_silent(do.call(foceiControl, .fc))
  })

  test_that("no foceiControl() user knob is stripped by the down-conversion", {
    # The general form of the mceta defect.  .impmapIsControlNames exists to drop
    # names foceiControl() does not accept; listing one it DOES accept silently
    # replaces the user's value with foceiControl()'s default in every consumer
    # of the down-converted control.  That is invisible whenever the two
    # defaults agree and wrong whenever they do not -- etaDistCorSuff is TRUE on
    # foceiControl() and FALSE on impmapControl(), so stripping it handed the
    # output path the opposite of what an imp user asked for.
    .strip <- nlmixr2est:::.impmapIsControlNames
    .fc <- names(formals(nlmixr2est::foceiControl))
    # flatEtaIdx is the one deliberate exception: it is a foceiControl argument,
    # but the value on an impmap control is imp's runtime index map rather than
    # a user setting, so it is dropped on purpose.
    .bad <- setdiff(intersect(.strip, .fc), "flatEtaIdx")
    expect_equal(.bad, character(0))
  })

  test_that("est=\"imp\" actually reaches the inner MAP with mceta >= 1", {
    .mod <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv)
        linCmt() ~ add(add.sd)
      })
    }
    .f <- suppressMessages(
      nlmixr2(.mod, nlmixr2data::theo_sd, est = "imp",
              control = impmapControl(nIter = 2L, isample = 50L, print = 0L,
                                      mceta = 10L, covMethod = "")))
    .ns <- .f$env$nMcetaStart
    expect_false(is.null(.ns))
    expect_true(sum(.ns) > 0)
    # the draws must be EXPLORED, not merely counted: eta=0 winning every
    # subject is what the old inert path looked like from the outside
    expect_true(.ns[["sample"]] > 0)

    # control: the counter is gated on mceta >= 1, so mceta=0 must not stamp it
    .f0 <- suppressMessages(
      nlmixr2(.mod, nlmixr2data::theo_sd, est = "imp",
              control = impmapControl(nIter = 2L, isample = 50L, print = 0L,
                                      mceta = 0L, covMethod = "")))
    expect_null(.f0$env$nMcetaStart)
  })
})
