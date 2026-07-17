nmTest({

  # Regression tests for issue #517: the FO/FOI estimation path returned the
  # intermediate (est="none") linearization fit env with an empty control, so
  # .updateParFixed() saw a NULL $control and silently fell back to defaults.

  test_that("nmObjGetControl.default surfaces a raw control binding (issue #517)", {
    # An intermediate fit whose est has no method-specific nmObjGetControl
    # should still return its stored control instead of NULL.
    .e <- new.env(parent = emptyenv())
    assign("est", "none", envir = .e)
    assign("control", foceiControl(), envir = .e)
    .obj <- list(.e)
    class(.obj) <- "none"
    expect_true(inherits(nlmixr2est:::nmObjGetControl(.obj), "foceiControl"))

    # With no control binding at all it must still return NULL.
    .e2 <- new.env(parent = emptyenv())
    assign("est", "none", envir = .e2)
    .obj2 <- list(.e2)
    class(.obj2) <- "none"
    expect_null(nlmixr2est:::nmObjGetControl(.obj2))
  })

  test_that("FO/FOI fits never see a NULL control in .updateParFixed (issue #517)", {
    .seen <- new.env(parent = emptyenv())
    .seen$nullControl <- logical(0)
    trace(
      nlmixr2est:::.updateParFixed,
      tracer = bquote({
        .env517 <- .(.seen)
        .env517$nullControl <- c(.env517$nullControl, is.null(.ret$control))
      }),
      print = FALSE
    )
    on.exit(suppressMessages(untrace(nlmixr2est:::.updateParFixed)), add = TRUE)

    fitFo <- .nlmixr(one.compartment, theo_sd, est = "fo",
                     control = foControl(print = 0))
    fitFoi <- .nlmixr(one.compartment, theo_sd, est = "foi",
                      control = foiControl(print = 0))

    # .updateParFixed ran on both the intermediate est="none" objects and the
    # final fits; none should have carried a NULL control.
    expect_gt(length(.seen$nullControl), 0L)
    expect_false(any(.seen$nullControl))

    # The final fits resolve their method control and produce a parFixed table.
    expect_true(inherits(fitFo$control, "foControl"))
    expect_true(inherits(fitFoi$control, "foiControl"))
    expect_s3_class(fitFo$parFixed, "data.frame")
    expect_s3_class(fitFoi$parFixed, "data.frame")
  })

  test_that("a built control with internal residThetaIdx round-trips for fo/foi", {
    # .foceiFamilyReturn stashes ui$foceiResidTheta on the control as
    # $residThetaIdx.  The fo/foi post-fit table step re-validates the stored
    # control by round-tripping it through foceiControl() (nmObjGetControl.fo /
    # .foi), so residThetaIdx must be an accepted internal dot-argument.  When it
    # was not, the round-trip errored with "unused argument: 'residThetaIdx'" and
    # the fit died with "cannot find fo/foi related control object".
    expect_true("residThetaIdx" %in% nlmixr2est:::.foceiControlInternal)

    .ctl <- foceiControl()
    .ctl$residThetaIdx <- c(1L, 2L, 3L)
    # foceiControl() accepts the round-tripped internal field
    expect_silent(do.call(foceiControl, .ctl))
    # and the fo/foi validators return the method control rather than erroring
    expect_true(inherits(nlmixr2est:::getValidNlmixrCtl.fo(list(.ctl)), "foControl"))
    expect_true(inherits(nlmixr2est:::getValidNlmixrCtl.foi(list(.ctl)), "foiControl"))
  })

})
