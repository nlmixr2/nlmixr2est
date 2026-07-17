nmTest({
  # Issue #286: a model mixing linCmt() with ODEs.  rxode2 requires the
  # linCmt() compartments to be the last states of the solve, so the FOCEi
  # inner model's eta-sensitivity states push depot/central past the
  # compartment numbers the data was translated against.  The linear part is
  # solved as ODEs instead for the FOCEi family.

  .pure <- function() {
    ini({
      tka <- 0.5; tcl <- 1; tv <- 3.5
      eta.ka ~ 0.2
      p <- 0.1
    })
    model({
      ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv)
      cp <- linCmt()
      cp ~ prop(p)
    })
  }
  .mixed <- function() {
    ini({
      tka <- 0.5; tcl <- 1; tv <- 3.5; tke0 <- 0
      eta.ka ~ 0.2
      p <- 0.1
    })
    model({
      ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv); ke0 <- exp(tke0)
      C2 <- linCmt()
      d/dt(ce) <- ke0 * (C2 - ce)
      ce ~ add(p)
    })
  }
  .ode <- function() {
    ini({
      tka <- 0.5; tcl <- 1; tv <- 3.5
      eta.ka ~ 0.2
      p <- 0.1
    })
    model({
      ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv)
      d/dt(depot) <- -ka * depot
      d/dt(central) <- ka * depot - (cl / v) * central
      cp <- central / v
      cp ~ prop(p)
    })
  }

  test_that("mixed linCmt()/ODE models are detected", {
    expect_true(.uiIsMixedLinCmtOde(.mixed()))
    # a linCmt() model with no other ODE keeps the analytic solution
    expect_false(.uiIsMixedLinCmtOde(.pure()))
    # depot/central written as ODEs is not a linCmt() model
    expect_false(.uiIsMixedLinCmtOde(.ode()))
  })

  test_that("the linCmt()/ODE hook only applies to the FOCEi family", {
    expect_null(.preProcessLinCmtOde(.mixed(), "saem", NULL, NULL))
    expect_null(.preProcessLinCmtOde(.mixed(), "nlme", NULL, NULL))
    expect_null(.preProcessLinCmtOde(.pure(), "focei", NULL, NULL))
    expect_true(is.list(.preProcessLinCmtOde(.mixed(), "focei", NULL, NULL)))
  })

  test_that("the mixed model is translated to ODEs without renumbering compartments", {
    .ui <- .mixed()
    .r <- .preProcessLinCmtOde(.ui, "focei", NULL, NULL)$ui
    # linCmt() is gone -- the compartments are real ODE states now
    expect_true(is.null(.r$mvL) || !.uiIsMixedLinCmtOde(.r))
    expect_false(any(vapply(.r$lstExpr, function(e) {
      any(all.vars(e) == "linCmt") ||
        (is.call(e) && length(e) > 2L && is.call(e[[3]]) &&
           identical(e[[3]][[1]], quote(linCmt)))
    }, logical(1))))
    # the data's numeric cmt must keep meaning the same compartment; linToOde()
    # on its own would return depot,central,ce
    expect_equal(.r$state, .ui$state)
    expect_equal(.r$state, c("ce", "depot", "central"))
  })

  test_that("the translated model keeps the linCmt() output defined before it is used", {
    .r <- .preProcessLinCmtOde(.mixed(), "focei", NULL, NULL)$ui
    .lines <- .r$lstExpr
    .isDdtCe <- vapply(.lines, function(e) {
      is.call(e) && is.call(e[[2]]) && identical(e[[2]][[2]], quote(d)) &&
        as.character(e[[2]][[3]][[2]]) == "ce"
    }, logical(1))
    .isC2 <- vapply(.lines, function(e) {
      is.call(e) && is.name(e[[2]]) && identical(e[[2]], quote(C2))
    }, logical(1))
    expect_true(any(.isC2))
    expect_true(any(.isDdtCe))
    # C2 <- central/v must precede d/dt(ce) <- ke0*(C2 - ce)
    expect_lt(which(.isC2)[1], which(.isDdtCe)[1])
  })
})
