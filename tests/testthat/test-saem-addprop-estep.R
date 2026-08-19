nmTest({
  test_that("saem E-step combined-error SD honors addProp (#912)", {
    # saemFormG() is the exact function every _scratch_g E-step/simulation
    # site in src/saem.cpp calls; pin it directly against the M-step
    # objective's own combined1/combined2 formulas (obj()/objH()/objI()'s
    # _saemAddProp branch) rather than comparing fitted estimates.  saemFormG()
    # takes the signed ft (it applies fabs() itself, matching every call site).
    a <- c(0.7, 0.7, 0.7)
    b <- c(0.2, 0.2, 0.2)
    f <- c(2, -3, 0)

    # combined1: g = a + b*|f|
    g1 <- as.vector(saemFormGTest(a, b, f, c(1L, 1L, 1L)))
    expect_equal(g1, a + b * abs(f))

    # combined2: g = sqrt(a^2 + b^2*f^2)
    g2 <- as.vector(saemFormGTest(a, b, f, c(2L, 2L, 2L)))
    expect_equal(g2, sqrt(a^2 + b^2 * f^2))

    # a multi-endpoint fit mixes both per observation in one call
    gm <- as.vector(saemFormGTest(a, b, f, c(1L, 2L, 1L)))
    expect_equal(gm, c(g1[1], g2[2], g1[3]))

    # the two formulas must actually disagree when both terms are present,
    # otherwise this test would still pass under the old hardcoded combined1
    expect_false(isTRUE(all.equal(g1, g2)))
  })

  test_that("saem E-step actually runs with addProp plumbed through (#912)", {
    # the unit test above only pins the standalone saemFormG() math; run a
    # minimal real fit so a deleted/mis-sized vecaddProp (uninitialized uvec,
    # wrong length vs. ntotal) fails here instead of only at runtime for a
    # user's model.
    mod <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
        add.sd <- 0.3; prop.sd <- 0.2
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd) + prop(prop.sd)
      })
    }
    f1 <- .nlmixr(mod, theo_sd, est = "saem",
                  control = saemControl(nBurn = 10, nEm = 10, print = 0, nmc = 2,
                                        addProp = "combined1"))
    f2 <- .nlmixr(mod, theo_sd, est = "saem",
                  control = saemControl(nBurn = 10, nEm = 10, print = 0, nmc = 2,
                                        addProp = "combined2"))
    expect_true(all(is.finite(f1$parFixedDf$Estimate)))
    expect_true(all(is.finite(f2$parFixedDf$Estimate)))
  })

  test_that("saem M-step add+pow combined2 objective takes sqrt() (objC)", {
    # objC() (case rmAddPow: add()+pow(), no boxCox/yeoJohnson) formed the
    # combined2 branch as a^2+b^2*f^(2*pw) and used it directly as the SD --
    # missing the sqrt() every sibling combined2 objective (obj()/objH()/
    # objI()) applies -- found reviewing the addProp branches for #912.
    # Regression: a finite, sane fit is impossible under the squared (wrong)
    # objective for typical theo_sd-scale add/pow starting values.
    mod <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
        add.sd <- 0.3; prop.sd <- 0.2; pw <- 1
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd) + pow(prop.sd, pw)
      })
    }
    fit <- .nlmixr(mod, theo_sd, est = "saem",
                   control = saemControl(nBurn = 20, nEm = 20, print = 0, nmc = 2))
    expect_true(all(is.finite(fit$parFixedDf$Estimate)))
    expect_true(fit$parFixedDf["add.sd", "Estimate"] > 0)
    expect_true(fit$parFixedDf["prop.sd", "Estimate"] > 0)
  })
})
