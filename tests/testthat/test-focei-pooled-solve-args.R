nmTest({
  test_that(".foceiPoolCores resolves the rxControl default to rxode2's threads", {
    # rxControl(cores = 0) is the default and MEANS "use rxode2's thread setting"; the
    # pooled augmented solve read it as a literal 1 and so never went parallel over
    # subjects (outerSolveFill only forks when cores > 1).  Asserted here rather than
    # through a fit because a host with getRxThreads() == 1 cannot tell the two apart.
    .n <- as.integer(rxode2::getRxThreads())
    expect_equal(.foceiPoolCores(0L), .n)
    expect_equal(.foceiPoolCores(0), .n)
    expect_equal(.foceiPoolCores(NULL), .n)
    expect_equal(.foceiPoolCores(NA_integer_), .n)
    # an explicit request is taken literally, including a deliberate cores = 1
    expect_equal(.foceiPoolCores(1L), 1L)
    expect_equal(.foceiPoolCores(3L), 3L)
    # nothing usable -> rxode2's threads, silently: this runs inside a fit, and a
    # coercion warning here would be collected into the fit's $runInfo
    expect_silent(expect_equal(.foceiPoolCores("nope"), .n))
    expect_equal(.foceiPoolCores(c(2L, 4L)), 2L)
  })

  test_that("covSolveTol beats the sigdig default for the augmented solve tolerance", {
    # The tolerance the pooled route was ignoring.  Both sources reach it through the same
    # `solveTol` argument, so there is no branch that could honour one and drop the other
    # -- what needs pinning is which of the two wins, and that it tracks sigdig rather
    # than a frozen literal.  No fit: this is control plumbing.
    .m <- function() {
      ini({ tka <- 0.45; add.sd <- 0.7; eta.ka ~ 0.6 })
      model({ ka <- exp(tka + eta.ka); cp <- ka; cp ~ add(add.sd) })
    }
    .ui <- suppressMessages(rxode2::rxUiDecompress(nlmixr2(.m)))
    .tol <- function(...) { .u <- .ui; .u$control <- foceiControl(...)
      .foceiAnalyticSolveTol(.u) }
    expect_equal(.tol(sigdig = 3), 1e-9)
    expect_equal(.tol(sigdig = 4), 1e-10)
    expect_equal(.tol(sigdig = 6), 1e-12)
    # an explicit covSolveTol wins, and does not move with sigdig
    expect_equal(.tol(covSolveTol = 1e-7), 1e-7)
    expect_equal(.tol(covSolveTol = 1e-7, sigdig = 6), 1e-7)
  })
})
