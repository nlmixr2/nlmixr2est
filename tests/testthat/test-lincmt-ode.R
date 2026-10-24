nmTest({
  # Issue #286: a model mixing linCmt() with ODEs.  rxode2 requires the
  # linCmt() compartments to be the last states of the solve, so the
  # sensitivity states these methods add (one per eta for FOCEi, one per theta
  # for nlm) push depot/central past the compartment numbers the data was
  # translated against.  The linear part is solved as ODEs instead, which is
  # warned about since the model then no longer mixes a solved system with
  # ODEs at all.

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

  # dose each numeric cmt in both models: the same compartment must receive it
  .sameCmtNumbers <- function(orig, translated) {
    .sig <- function(m, cmt) {
      .s <- as.data.frame(rxode2::rxSolve(m, rxode2::et(amt = 100, cmt = cmt) |> rxode2::et(c(0.5, 2)), omega = NA))
      unlist(.s[, c("C2", .linCmtOdeDdtStates(orig$lstExpr))])
    }
    # every compartment of the translation, peripherals included (orig$state
    # does not list those)
    for (.k in seq_along(translated$state)) {
      expect_equal(.sig(translated, .k), .sig(orig, .k), tolerance = 1e-5, info = paste("cmt", .k))
    }
  }

  test_that("mixed linCmt()/ODE models are detected", {
    expect_true(.uiIsMixedLinCmtOde(.mixed()))
    # a linCmt() model with no other ODE keeps the analytic solution
    expect_false(.uiIsMixedLinCmtOde(.pure()))
    # depot/central written as ODEs is not a linCmt() model
    expect_false(.uiIsMixedLinCmtOde(.ode()))
  })

  test_that("the linCmt()/ODE hook applies to the methods that add sensitivity states", {
    # saem/nlme build no sensitivity compartments, so they keep the analytic
    # linCmt() solution -- and must not warn
    expect_silent(expect_null(.preProcessLinCmtOde(.mixed(), "saem", NULL, NULL)))
    expect_silent(expect_null(.preProcessLinCmtOde(.mixed(), "nlme", NULL, NULL)))
    # a linCmt() model with no other ODE is never shifted
    expect_silent(expect_null(.preProcessLinCmtOde(.pure(), "focei", NULL, NULL)))
    expect_silent(expect_null(.preProcessLinCmtOde(.pure(), "nlm", NULL, NULL)))
    # FOCEi family (eta sensitivities) and nlm family (theta sensitivities)
    for (.e in c("focei", "foce", "fo", "laplace", "agq", "nlm", "nlminb", "bobyqa", "optim", "n1qn1")) {
      expect_true(is.list(suppressWarnings(.preProcessLinCmtOde(.mixed(), .e, NULL, NULL))), info = .e)
    }
  })

  test_that("solving the linear compartments as ODEs is warned about, not silent", {
    # the model no longer mixes a solved system with ODEs, so the user has to
    # be told the analytic linCmt() is not being used
    expect_warning(.preProcessLinCmtOde(.mixed(), "focei", NULL, NULL), "cannot use the analytic 'linCmt\\(\\)'")
    expect_warning(.preProcessLinCmtOde(.mixed(), "nlm", NULL, NULL), "solved as ODEs")
    # the warning names the routine that could not use it
    expect_warning(.preProcessLinCmtOde(.mixed(), "focei", NULL, NULL), "focei")
    expect_warning(.preProcessLinCmtOde(.mixed(), "nlm", NULL, NULL), "nlm")
  })

  test_that("the mixed model is translated to ODEs without renumbering compartments", {
    .ui <- .mixed()
    .r <- suppressWarnings(.preProcessLinCmtOde(.ui, "focei", NULL, NULL))$ui
    # linCmt() is gone -- the compartments are real ODE states now
    expect_true(is.null(.r$mvL) || !.uiIsMixedLinCmtOde(.r))
    expect_false(any(vapply(
      .r$lstExpr,
      function(e) {
        any(all.vars(e) == "linCmt") ||
          (is.call(e) && length(e) > 2L && is.call(e[[3]]) && identical(e[[3]][[1]], quote(linCmt)))
      },
      logical(1)
    )))
    # the data's numeric cmt must keep meaning the same compartment: linCmt()'s
    # depot/central are numbered first, then the ODE states
    expect_equal(.r$state, c("depot", "central", "ce"))
    .sameCmtNumbers(.ui, .r)
  })

  test_that("the compartment numbers survive for other linCmt() shapes", {
    .translate <- function(ui) {
      .w <- NULL
      .r <- withCallingHandlers(
        .preProcessLinCmtOde(ui, "focei", NULL, NULL)$ui,
        warning = function(w) {
          .w <<- c(.w, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )
      # the numbering was kept, so nothing may claim it was renumbered
      expect_false(any(grepl("renumbered from", .w)))
      .r
    }
    # an ODE declared before linCmt()
    .before <- function() {
      ini({
        tka <- 0.5; tcl <- 1; tv <- 3.5; tke0 <- 0
        eta.ka ~ 0.2
        p <- 0.1
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv); ke0 <- exp(tke0)
        d/dt(eff) <- -ke0 * eff
        C2 <- linCmt()
        d/dt(ce) <- ke0 * (C2 - ce) + eff
        ce ~ add(p)
      })
    }
    .ui <- rxode2::rxode2(.before)
    .r <- .translate(.ui)
    expect_equal(.r$state, c("depot", "central", "eff", "ce"))
    .sameCmtNumbers(.ui, .r)
    # a one-compartment IV linCmt()
    .iv1 <- function() {
      ini({
        tcl <- 1; tv <- 3.5; tke0 <- 0
        eta.cl ~ 0.2
        p <- 0.1
      })
      model({
        cl <- exp(tcl + eta.cl); v <- exp(tv); ke0 <- exp(tke0)
        d/dt(eff) <- -ke0 * eff
        C2 <- linCmt()
        d/dt(ce) <- ke0 * (C2 - ce) + eff
        ce ~ add(p)
      })
    }
    .ui <- rxode2::rxode2(.iv1)
    .r <- .translate(.ui)
    expect_equal(.r$state, c("central", "eff", "ce"))
    .sameCmtNumbers(.ui, .r)
    # a two-compartment IV linCmt(): its peripheral compartment goes last
    .iv2 <- function() {
      ini({
        tcl <- 1; tv <- 3.5; tq <- 0; tvp <- 4; tke0 <- 0
        eta.cl ~ 0.2
        p <- 0.1
      })
      model({
        cl <- exp(tcl + eta.cl); v <- exp(tv); q <- exp(tq); vp <- exp(tvp); ke0 <- exp(tke0)
        C2 <- linCmt()
        d/dt(eff) <- ke0 * (C2 - eff)
        eff ~ add(p)
      })
    }
    .ui <- rxode2::rxode2(.iv2)
    .r <- .translate(.ui)
    expect_equal(.r$state, c("central", "eff", "peripheral1"))
    .sameCmtNumbers(.ui, .r)
    # a two-compartment oral linCmt() with an ODE declared before it
    .oral2 <- function() {
      ini({
        tka <- 0.5; tcl <- 1; tv <- 3.5; tq <- 0; tvp <- 4; tke0 <- 0
        eta.ka ~ 0.2
        p <- 0.1
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv); q <- exp(tq); vp <- exp(tvp); ke0 <- exp(tke0)
        d/dt(eff) <- -ke0 * eff
        C2 <- linCmt()
        d/dt(ce) <- ke0 * (C2 - ce) + eff
        ce ~ add(p)
      })
    }
    .ui <- rxode2::rxode2(.oral2)
    .r <- .translate(.ui)
    expect_equal(.r$state, c("depot", "central", "eff", "ce", "peripheral1"))
    .sameCmtNumbers(.ui, .r)
    # a three-compartment oral linCmt(): both peripherals go last
    .oral3 <- function() {
      ini({
        tka <- 0.5; tcl <- 1; tv <- 3.5; tq <- 0; tvp <- 4; tq2 <- -1; tvp2 <- 5; tke0 <- 0
        eta.ka ~ 0.2
        p <- 0.1
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv); q <- exp(tq); vp <- exp(tvp)
        q2 <- exp(tq2); vp2 <- exp(tvp2); ke0 <- exp(tke0)
        d/dt(eff) <- -ke0 * eff
        C2 <- linCmt()
        d/dt(ce) <- ke0 * (C2 - ce) + eff
        ce ~ add(p)
      })
    }
    .ui <- rxode2::rxode2(.oral3)
    .r <- .translate(.ui)
    expect_equal(.r$state, c("depot", "central", "eff", "ce", "peripheral1", "peripheral2"))
    .sameCmtNumbers(.ui, .r)
    # the model's only ODE inside an if () block is still an ODE of the model
    .ifOde <- function() {
      ini({
        tka <- 0.5; tcl <- 1; tv <- 3.5; tke0 <- 0
        eta.ka ~ 0.2
        p <- 0.1
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv); ke0 <- exp(tke0)
        C2 <- linCmt()
        if (t > -1) {
          d/dt(ce) <- ke0 * (C2 - ce)
        } else {
          d/dt(ce) <- 0
        }
        ce ~ add(p)
      })
    }
    .ui <- rxode2::rxode2(.ifOde)
    expect_true(.uiIsMixedLinCmtOde(.ui))
    .r <- .translate(.ui)
    expect_equal(.r$state, c("depot", "central", "ce"))
    .sameCmtNumbers(.ui, .r)
    # the model's own cmt() declaration does not move linCmt()'s numbers
    .ownCmt <- function() {
      ini({
        tka <- 0.5; tcl <- 1; tv <- 3.5; tke0 <- 0
        eta.ka ~ 0.2
        p <- 0.1
      })
      model({
        cmt(eff)
        ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv); ke0 <- exp(tke0)
        d/dt(eff) <- -ke0 * eff
        C2 <- linCmt()
        d/dt(ce) <- ke0 * (C2 - ce) + eff
        ce ~ add(p)
      })
    }
    .ui <- rxode2::rxode2(.ownCmt)
    .r <- .translate(.ui)
    expect_equal(.r$state, c("depot", "central", "eff", "ce"))
    .sameCmtNumbers(.ui, .r)
    # an IV linCmt() next to the model's own ODE named depot: that depot is not
    # linCmt()'s, so it keeps its place after central
    .ivDepot <- function() {
      ini({
        tcl <- 1; tv <- 3.5; tke0 <- 0
        eta.cl ~ 0.2
        p <- 0.1
      })
      model({
        cl <- exp(tcl + eta.cl); v <- exp(tv); ke0 <- exp(tke0)
        d/dt(depot) <- -ke0 * depot
        C2 <- linCmt()
        d/dt(ce) <- ke0 * (C2 - ce) + depot
        ce ~ add(p)
      })
    }
    .ui <- rxode2::rxode2(.ivDepot)
    .r <- .translate(.ui)
    expect_equal(.r$state, c("central", "depot", "ce"))
    .sameCmtNumbers(.ui, .r)
  })

  test_that("renumbering a translation declares cmt() and moves no model line", {
    # what an older linToOde() produced for an ODE declared before linCmt():
    # eff numbered first.  A variable reassigned between two d/dt() lines has
    # to keep being read where it was, so no line may move.
    .old <- function() {
      ini({
        tka <- 0.5; tcl <- 1; tv <- 3.5; tke0 <- 0
        eta.ka ~ 0.2
        p <- 0.1
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv); ke0 <- exp(tke0)
        r <- 1
        d/dt(eff) <- -ke0 * eff * r
        d/dt(depot) <- -ka * depot
        d/dt(central) <- ka * depot - cl / v * central
        C2 <- central / v
        r <- 2
        d/dt(ce) <- ke0 * (C2 - ce) * r
        ce ~ add(p)
      })
    }
    .ui <- rxode2::rxode2(.old)
    expect_equal(.ui$state, c("eff", "depot", "central", "ce"))
    .target <- c("depot", "central", "eff", "ce")
    .r <- .linCmtOdeRestoreStateOrder(.ui, .target)
    expect_equal(.r$state, .target)
    # the model lines are untouched, only cmt() declarations lead
    .n <- length(.target)
    expect_equal(.r$lstExpr[-seq_len(.n)], .ui$lstExpr)
    expect_equal(vapply(.r$lstExpr[seq_len(.n)], deparse1, ""), paste0("cmt(", .target, ")"))
    # and a dose to each compartment by name solves exactly as before
    for (.c in .target) {
      .ev <- rxode2::et(amt = 100, cmt = .c) |> rxode2::et(c(0.5, 2))
      .a <- as.data.frame(rxode2::rxSolve(.ui, .ev, omega = NA))
      .b <- as.data.frame(rxode2::rxSolve(.r, .ev, omega = NA))
      expect_equal(.b[, c("C2", "eff", "ce")], .a[, c("C2", "eff", "ce")], info = .c)
    }
    # already in order: returned as is
    expect_identical(.linCmtOdeRestoreStateOrder(.r, .target), .r)
  })

  test_that("the translated model keeps the linCmt() output defined before it is used", {
    .r <- suppressWarnings(.preProcessLinCmtOde(.mixed(), "focei", NULL, NULL))$ui
    .lines <- .r$lstExpr
    .isDdtCe <- vapply(
      .lines,
      function(e) {
        is.call(e) && is.call(e[[2]]) && identical(e[[2]][[2]], quote(d)) && as.character(e[[2]][[3]][[2]]) == "ce"
      },
      logical(1)
    )
    .isC2 <- vapply(
      .lines,
      function(e) {
        is.call(e) && is.name(e[[2]]) && identical(e[[2]], quote(C2))
      },
      logical(1)
    )
    expect_true(any(.isC2))
    expect_true(any(.isDdtCe))
    # C2 <- central/v must precede d/dt(ce) <- ke0*(C2 - ce)
    expect_lt(which(.isC2)[1], which(.isDdtCe)[1])
  })
})
