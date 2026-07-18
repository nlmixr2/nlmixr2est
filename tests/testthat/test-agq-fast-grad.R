# Analytic AGQ outer gradient (agqControl(fast=TRUE)): AGQ is FOCEI with the objective's
# data term replaced by log(sum_k a_k) over the quadrature grid, so at nAGQ=1 the kernel
# must reproduce the FOCEI gradient exactly, and for nAGQ>1 it must agree with central
# differences of the AGQ objective.  Out-of-scope models fall back to finite differences.

nmTest({
  .agq_one_cmt <- function() {
    ini({
      tka <- log(1.5); tcl <- log(2.7); tv <- log(31.5)
      eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
      d/dt(depot)  <- -ka * depot
      d/dt(center) <-  ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }
  .ctl <- function(...) {
    agqControl(maxOuterIterations = 0L, maxInnerIterations = 500L, covMethod = "",
               calcTables = FALSE, print = 0L, ...)
  }

  test_that("agqf/magqf/iagqf are registered and default to fast=TRUE", {
    for (.m in c("agqf", "magqf", "iagqf")) {
      expect_true(!is.null(getS3method("nlmixr2Est", .m, optional = TRUE)))
      expect_true(!is.null(getS3method("getValidNlmixrCtl", .m, optional = TRUE)))
    }
    expect_true(getValidNlmixrCtl.agqf(list(agqControl()))$fast)
    expect_true(getValidNlmixrCtl.magqf(list(magqControl()))$fast)
    expect_true(getValidNlmixrCtl.iagqf(list(iagqControl()))$fast)
    # a defaulted outer optimizer re-defaults under fast; an explicit one is kept
    expect_equal(getValidNlmixrCtl.agqf(list(agqControl()))$outerOptTxt, "lbfgsb3c")
    expect_equal(getValidNlmixrCtl.agqf(list(agqControl(outerOpt = "nlminb")))$outerOptTxt, "nlminb")
    # agqControl forwards ... to foceiControl, so fast is accepted directly too
    expect_true(agqControl(fast = TRUE)$fast)
    expect_false(agqControl()$fast)
  })

  test_that("nAGQ=1: the AGQ dispatch reproduces the Laplace analytic gradient", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    # At one node (x=0, pi=1) the quadrature collapses onto the mode, so the AGQ objective
    # equals the Laplace objective and the outer gradients must agree (the envelope term
    # Phi_eta(etahat)'etaP is zero at a converged EBE).  Both fits fix theta
    # (maxOuterIterations=0), so the two analytic gradients are evaluated at the same point;
    # this is the identity that guards the duplicated eta-hat block.
    .c0 <- list(maxOuterIterations = 0L, maxInnerIterations = 500L,
                covMethod = "", calcTables = FALSE, print = 0L)
    .fl <- suppressMessages(nlmixr2(.agq_one_cmt, nlmixr2data::theo_sd, "laplace",
      do.call(laplaceControl, c(list(fast = TRUE), .c0))))
    .fa <- suppressMessages(nlmixr2(.agq_one_cmt, nlmixr2data::theo_sd, "agq",
      do.call(agqControl, c(list(fast = TRUE, nAGQ = 1L), .c0))))
    .gl <- .foceiGradAnalyticCalc(.fl)
    .ga <- .foceiGradAnalyticCalc(.fa)
    expect_false(is.null(.gl))
    expect_false(is.null(.ga))
    expect_true(all(is.finite(.ga)))
    expect_equal(unname(.ga), unname(.gl), tolerance = 1e-4)
  })

  test_that("analytic AGQ gradient matches central differences (nAGQ=2 and 3)", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    for (.n in c(2L, 3L)) {
      .f <- suppressMessages(nlmixr2(.agq_one_cmt, nlmixr2data::theo_sd, "agq",
                                     .ctl(nAGQ = .n, fast = TRUE, sigdig = 7)))
      .g <- .foceiGradAnalyticCalc(.f)
      expect_false(is.null(.g))
      .base <- fixef(.f)
      .ofvAt <- function(nm, val) {
        .ui <- do.call(rxode2::ini, c(list(.f$finalUi), setNames(list(val), nm)))
        suppressMessages(suppressWarnings(nlmixr2(.ui, nlmixr2data::theo_sd, "agq",
                                                  .ctl(nAGQ = .n, sigdig = 7))))$objf
      }
      # NB h: the AGQ objective's central-difference error bottoms out around 3e-3..1e-2;
      # 1e-4 sits on the noisy side of the V and reads ~1e-3 relative even for an exact
      # gradient, so do not tighten this.
      .fd <- vapply(names(.base), function(nm) {
        h <- 3e-3 * max(abs(.base[[nm]]), 1)
        (.ofvAt(nm, .base[nm] + h) - .ofvAt(nm, .base[nm] - h)) / (2 * h)
      }, numeric(1))
      expect_equal(unname(.g[names(.base)]), unname(.fd), tolerance = 0.02)
    }
  })

  test_that("out-of-scope AGQ models fall back to the finite-difference gradient", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    # fast=FALSE: no analytic gradient at all
    .f0 <- suppressMessages(nlmixr2(.agq_one_cmt, nlmixr2data::theo_sd, "agq", .ctl(nAGQ = 2L)))
    expect_null(.foceiGradAnalyticCalc(.f0))
    # an active agqLow/agqHi clamp kinks the objective (both default to +/-Inf)
    .fc <- suppressMessages(nlmixr2(.agq_one_cmt, nlmixr2data::theo_sd, "agq",
                                    .ctl(nAGQ = 2L, fast = TRUE, agqLow = -1e6)))
    expect_null(.foceiGradAnalyticCalc(.fc))
    # cholSEOpt uses a different Cholesky factor than chol(), and the factor places the
    # quadrature nodes -- differentiating chol() would be the wrong function
    .fs <- suppressMessages(nlmixr2(.agq_one_cmt, nlmixr2data::theo_sd, "agq",
                                    .ctl(nAGQ = 2L, fast = TRUE, cholSEOpt = TRUE)))
    expect_null(.foceiGradAnalyticCalc(.fs))
  })

  test_that("fast=TRUE AGQ fit matches the finite-difference fit", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    # outerOpt forced on BOTH arms: fast=TRUE otherwise re-defaults nlminb -> lbfgsb3c,
    # which would make this compare optimizers rather than gradients.
    .fit <- function(fast) {
      suppressMessages(suppressWarnings(nlmixr2(.agq_one_cmt, nlmixr2data::theo_sd, "agq",
        agqControl(nAGQ = 2L, fast = fast, outerOpt = "lbfgsb3c", covMethod = "",
                   calcTables = FALSE, print = 0L))))
    }
    .fd <- .fit(FALSE); .an <- .fit(TRUE)
    expect_equal(.an$objf, .fd$objf, tolerance = 1e-3)
    expect_equal(unname(fixef(.an)[names(fixef(.fd))]), unname(fixef(.fd)), tolerance = 0.02)
    expect_equal(unname(diag(.an$omega)), unname(diag(.fd$omega)), tolerance = 0.05)
    expect_equal(.an$objDf[["OBJF"]], .fd$objDf[["OBJF"]], tolerance = 1e-3)
  })

  test_that("est='agqf' equals est='agq' with fast=TRUE", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    .a <- suppressMessages(suppressWarnings(nlmixr2(.agq_one_cmt, nlmixr2data::theo_sd, "agqf",
      agqControl(nAGQ = 2L, outerOpt = "lbfgsb3c", covMethod = "", calcTables = FALSE, print = 0L))))
    .b <- suppressMessages(suppressWarnings(nlmixr2(.agq_one_cmt, nlmixr2data::theo_sd, "agq",
      agqControl(nAGQ = 2L, fast = TRUE, outerOpt = "lbfgsb3c", covMethod = "",
                 calcTables = FALSE, print = 0L))))
    expect_equal(.a$objf, .b$objf, tolerance = 1e-8)
    expect_equal(unname(fixef(.a)), unname(fixef(.b)), tolerance = 1e-6)
  })

  test_that("analytic AGQ gradient matches central differences for a covariate model", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    .cov <- function() {
      ini({ tka <- log(1.5); tcl <- log(2.7); tv <- log(31.5); wt.cl <- 0.1
            eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1; add.sd <- 0.7 })
      model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl + wt.cl * (WT - 70))
              v <- exp(tv + eta.v)
              d/dt(depot) <- -ka * depot; d/dt(center) <- ka * depot - cl / v * center
              cp <- center / v; cp ~ add(add.sd) })
    }
    .f <- suppressMessages(nlmixr2(.cov, nlmixr2data::theo_sd, "agq",
                                   .ctl(nAGQ = 2L, fast = TRUE, sigdig = 7)))
    .g <- .foceiGradAnalyticCalc(.f)
    expect_false(is.null(.g))
    .base <- fixef(.f)
    .ofvAt <- function(nm, val) {
      .ui <- do.call(rxode2::ini, c(list(.f$finalUi), setNames(list(val), nm)))
      suppressMessages(suppressWarnings(nlmixr2(.ui, nlmixr2data::theo_sd, "agq",
                                                .ctl(nAGQ = 2L, sigdig = 7))))$objf
    }
    .fd <- vapply(names(.base), function(nm) {
      h <- 3e-3 * max(abs(.base[[nm]]), 1)
      (.ofvAt(nm, .base[nm] + h) - .ofvAt(nm, .base[nm] - h)) / (2 * h)
    }, numeric(1))
    expect_equal(unname(.g[names(.base)]), unname(.fd), tolerance = 0.02)
  })
})
