# A model with an occasion eta, and the same model with a subject-level eta in
# its place
.npdeIovMod <- function() {
  ini({
    tka <- 0.45
    tcl <- 1
    tv <- 3.45
    add.sd <- 0.7
    eta.v ~ 0.1
    iov.cl ~ 0.1 | occ
  })
  model({
    ka <- exp(tka)
    cl <- exp(tcl + iov.cl)
    v <- exp(tv + eta.v)
    linCmt() ~ add(add.sd)
  })
}

# eta.v repeats eta.cl's subject-level block ("id:same:eta.cl")
.npdeSameMod <- function() {
  ini({
    tka <- 0.45
    tcl <- 1
    tv <- 3.45
    add.sd <- 0.7
    eta.cl ~ 0.3
    eta.v ~ same()
  })
  model({
    ka <- exp(tka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    linCmt() ~ add(add.sd)
  })
}

.npdeIovData <- function() {
  .d <- nlmixr2data::theo_sd
  .d$occ <- ifelse(.d$TIME < 6, 1L, 2L)
  .d
}

.npdeIovFit <- function(est, table = tableControl()) {
  suppressWarnings(suppressMessages(
    nlmixr2(.npdeIovMod, .npdeIovData(), est, control = list(print = 0, covMethod = ""), table = table)
  ))
}

nmTest({
  test_that("npde", {
    one.cmt <- function() {
      ini({
        ## You may label each parameter with a comment
        tka <- 0.45 # Ka
        tcl <- log(c(0, 2.7, 100)) # Log Cl
        ## This works with interactive models
        ## You may also label the preceding line with label("label text")
        tv <- 3.45; label("log V")
        ## the label("Label name") works with all models
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }

    # Don't use saemControlFast because numeric results are tested below
    fit <- .nlmixr(one.cmt, theo_sd, est = "saem")

    expect_false(all(c("EPRED", "ERES", "NPDE", "NPD", "PDE", "PD") %in% names(fit)))
    # npde's uniforms are seeded per observation, so the thread count cannot move them
    .cols <- c("EPRED", "ERES", "NPDE", "NPD", "PDE", "PD")
    .old <- rxode2::getRxThreads(verbose = FALSE)
    rxode2::setRxThreads(1L)
    .t1 <- as.data.frame(suppressMessages(addNpde(fit, updateObject = FALSE)))[, .cols]
    rxode2::setRxThreads(2L)
    .t2 <- as.data.frame(suppressMessages(addNpde(fit, updateObject = FALSE)))[, .cols]
    rxode2::setRxThreads(.old)
    expect_identical(.t1, .t2)
    suppressMessages(expect_error(addNpde(fit), NA))
    expect_true(all(c("EPRED", "ERES", "NPDE", "NPD", "PDE", "PD") %in% names(fit)))
    .range1 <- range(fit$PDE)
    .range2 <- range(fit$PD)
    expect_true(.range1[1] > 0)
    expect_true(.range2[1] > 0)
    expect_true(.range1[2] < 1)
    expect_true(.range2[2] < 1)

    .range1 <- range(fit$NPDE)
    .range2 <- range(fit$NPD)
    .range3 <- range(fit$ERES)

    expect_true(.range1[1] < -1.7)
    expect_true(.range2[1] < -1.7)
    expect_true(.range3[1] < -1.7)
    expect_true(.range1[2] > 1.7)
    expect_true(.range2[2] > 1.7)
    expect_true(.range3[2] > 1.7)

    .range4 <- range(fit$EPRED)
    expect_true(.range4[1] > -0.1)
    expect_true(.range4[2] > 7)

    fit <- .nlmixr(one.cmt, theo_sd, est = "saem")

    expect_false(all(c("EPRED", "ERES", "NPDE", "NPD", "PDE", "PD") %in% names(fit)))
    fit2 <- suppressMessages(addNpde(fit, updateObject = FALSE))
    expect_false(all(c("EPRED", "ERES", "NPDE", "NPD", "PDE", "PD") %in% names(fit)))
    expect_true(all(c("EPRED", "ERES", "NPDE", "NPD", "PDE", "PD") %in% names(fit2)))

    fit <- .nlmixr(one.cmt, theo_sd, est = "saem", control = saemControlFast, table = tableControl(npde = TRUE))

    expect_true(all(c("EPRED", "ERES", "NPDE", "NPD", "PDE", "PD") %in% names(fit)))
  })

  test_that("pheno", {
    pheno <- function() {
      ini({
        tcl <- log(0.008) # typical value of clearance
        tv <-  log(0.6)   # typical value of volume
        ## var(eta.cl)
        eta.cl + eta.v ~ c(1,
                           0.01, 1) ## cov(eta.cl, eta.v), var(eta.v)
        # interindividual variability on clearance and volume
        add.err <- 0.1    # residual variability
      })
      model({
        cl <- exp(tcl + eta.cl) # individual value of clearance
        v <- exp(tv + eta.v)    # individual value of volume
        ke <- cl / v            # elimination rate constant
        d/dt(A1) = - ke * A1    # model differential equation
        cp = A1 / v             # concentration in plasma
        cp ~ add(add.err)       # define error model
      })
    }

    fit <- .nlmixr(pheno, pheno_sd, est = "saem", control = saemControlFast, table = list(npde = TRUE))

    # Since there is a correlation here the npde and npd

    expect_false(isTRUE(all.equal(fit$NPDE, fit$NPD)))
    expect_false(isTRUE(all.equal(fit$PDE, fit$PD)))
  })

  test_that("addNpde() on a fit with IOV gives the NPDE the fit computes itself", {
    # the finished fit's thetas, etas and data are the user's model's, so
    # addNpde() refits; its columns match those of a fit asked for NPDE
    .cols <- c("EPRED", "ERES", "NPDE", "NPD", "PDE", "PD")
    for (.est in c("focei", "saem")) {
      .with <- .npdeIovFit(.est, tableControl(npde = TRUE))
      .fit <- .npdeIovFit(.est)
      .new <- suppressMessages(addNpde(.fit, updateObject = FALSE))
      expect_identical(setdiff(names(.new), names(.fit)), .cols, label = .est)
      expect_identical(nrow(.new), nrow(.fit), label = .est)
      for (.c in .cols) {
        expect_equal(.new[[.c]], .with[[.c]], tolerance = 1e-10, label = paste(.est, .c))
      }
      expect_true(all(is.finite(.new$NPDE)), label = .est)
    }
  })

  test_that(".uiHasIov() finds occasion etas, not repeated subject blocks", {
    expect_true(.uiHasIov(rxode2::rxode2(.npdeIovMod)))
    expect_false(.uiHasIov(rxode2::rxode2(.npdeSameMod)))
  })
})
