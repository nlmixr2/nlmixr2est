nmTest({
  test_that("saem E-step combined-error SD honors addProp (#912)", {
    # saemFormG() is the exact function every _scratch_g E-step/simulation
    # site in src/saem.cpp calls; pin it directly against the M-step
    # objective's own combined1/combined2 formulas (obj()/objH()/objI()'s
    # _saemAddProp branch) rather than comparing fitted estimates.
    a <- c(0.7, 0.7, 0.7)
    b <- c(0.2, 0.2, 0.2)
    f <- c(2, -3, 0)

    # combined1: g = a + b*|f|
    g1 <- as.vector(nlmixr2est:::saemFormGTest(a, b, abs(f), c(1L, 1L, 1L)))
    expect_equal(g1, a + b * abs(f))

    # combined2: g = sqrt(a^2 + b^2*f^2)
    g2 <- as.vector(nlmixr2est:::saemFormGTest(a, b, abs(f), c(2L, 2L, 2L)))
    expect_equal(g2, sqrt(a^2 + b^2 * f^2))

    # a multi-endpoint fit mixes both per observation in one call
    gm <- as.vector(nlmixr2est:::saemFormGTest(a, b, abs(f), c(1L, 2L, 1L)))
    expect_equal(gm, c(g1[1], g2[2], g1[3]))

    # the two formulas must actually disagree when both terms are present,
    # otherwise this test would still pass under the old hardcoded combined1
    expect_false(isTRUE(all.equal(g1, g2)))
  })
})
