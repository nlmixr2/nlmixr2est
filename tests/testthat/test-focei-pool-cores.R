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
})
