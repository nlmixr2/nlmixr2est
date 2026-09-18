# `mceta` sets how the FOCEi inner MAP picks its starting eta, and imp reaches
# that same inner problem -- so `impmapControl(mceta=)` is a real lever on an
# imp fit.  Nothing reported it.
#
# `fit$env$nMcetaStart` -- the only evidence the omega draws were explored --
# was stamped deep inside the non-EM branch of `foceiFinalizeTables()`, the one
# that builds the `extra` details string.  est="imp"/"impmap" takes the
# isImpmap branch instead, so the counter was absent from every imp fit
# whatever mceta was set to.  Absence was not evidence of inertness, which is
# exactly how a genuinely inert setting would also look -- measured on
# origin/main with three etas, mceta=0 and mceta=10 gave objf 116.8385 and
# 116.8386 (the lever IS working) while both fits reported no counter at all.
#
# The estimates are no evidence either way: the EM is stochastic, so mceta=0
# and mceta=10 differ a little whether or not the draws were ever explored.
# That is why this asserts the counter rather than the numbers.
nmTest({
  test_that("impmapControl(mceta=) survives the down-conversion to foceiControl", {
    # A regression guard rather than a fix: mceta is a genuine foceiControl()
    # argument, and stripping it on the way down would leave the consumers of
    # the down-converted control (.setOfvFo(), the general-likelihood tables)
    # on the package default rather than the user's value.  Estimation reads
    # mceta off the impmapControl directly, so a fit would not notice.
    expect_false("mceta" %in% .impmapIsControlNames)
    .e <- new.env(parent = emptyenv())
    .e$impmapControl <- impmapControl(mceta = 7L)
    .fc <- .impmapControlToFoceiControl(.e, assign = FALSE)
    expect_true("mceta" %in% names(.fc))
    expect_equal(.fc$mceta, 7L)
    # and the down-converted list is still something foceiControl() accepts
    expect_silent(do.call(foceiControl, .fc))
  })

  test_that("est=\"imp\" reports which mceta candidate its inner solves started from", {
    skip_on_cran()
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
    # subject is what an inert setting would also look like from the outside
    expect_true(.ns[["sample"]] > 0)

    # control: the counter is gated on mceta >= 1, so mceta=0 must not stamp it
    .f0 <- suppressMessages(
      nlmixr2(.mod, nlmixr2data::theo_sd, est = "imp",
              control = impmapControl(nIter = 2L, isample = 50L, print = 0L,
                                      mceta = 0L, covMethod = "")))
    expect_null(.f0$env$nMcetaStart)
  })
})
