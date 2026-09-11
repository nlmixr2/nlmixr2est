nmTest({
  # #1047: saem parameterizes a random effect by the population parameter it is
  # added to, so an eta paired with none owns no Gamma2_phi1 column.  It was
  # then silently dropped from the kernel's `model$omega` -- `m[NA, NA] <- 1` is
  # a no-op in R -- and never sampled, so the fit ran on a model without that
  # random effect and only died at the very end, assembling the reported omega,
  # with "subscript out of bounds".

  # A non-mu eta on a non-"id" condition: rxode2's mu-ref downgrade only
  # records "id" etas into `nonMuEtas`, so this one is left with no phi.
  .noPhiMod <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      eta.ka ~ 0.6
      eta.occ ~ 0.1 | occ
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + 0.3 * logit(pnorm(eta.occ)))
      v <- exp(tv)
      linCmt() ~ add(add.sd)
    })
  }

  .muMod <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
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

  test_that("an eta with no phi parameter is named, not given a slot", {
    .ui <- rxode2::rxUiDecompress(rxode2::rxode2(.noPhiMod))
    # the eta maps to no saem population parameter at all
    expect_equal(.ui$saemEtaTrans, c(1L, NA_integer_))
    # ... and therefore owns no Gamma2_phi1 column.  Handing it the dense rank 2
    # is what indexed past a 1x1 matrix.
    expect_equal(.ui$saemOmegaTrans, c(1L, NA_integer_))
    expect_equal(.saemEtaNoPhi(.ui), "eta.occ")

    .ui2 <- rxode2::rxUiDecompress(rxode2::rxode2(.muMod))
    expect_equal(.ui2$saemOmegaTrans, c(1L, 2L))
    expect_equal(.saemEtaNoPhi(.ui2), character(0))
  })

  test_that("saem refuses such a model up front, naming the random effect", {
    .ui <- rxode2::rxUiDecompress(rxode2::rxode2(.noPhiMod))
    expect_error(.saemAssertEtaPhi(.ui), "eta.occ")
    expect_error(.saemAssertEtaPhi(.ui), "cannot sample")
    .ok <- rxode2::rxUiDecompress(rxode2::rxode2(.muMod))
    expect_silent(.saemAssertEtaPhi(.ok))

    # the check is wired into the estimation method, so it fires before any
    # iteration rather than after the whole run
    .env <- new.env(parent=emptyenv())
    .env$ui <- .ui
    expect_error(nlmixr2Est.saem(.env), "eta.occ")
  })

  # Two random effects mu-referenced to the SAME population parameter: saem
  # gives a phi one Gamma2_phi1 column, so only the first is sampled.
  .sharedPhiMod <- function() {
    ini({
      tka <- 0.45
      tv <- 3.45
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tka + eta.cl)
      v <- exp(tv)
      linCmt() ~ add(add.sd)
    })
  }

  test_that("two etas on one population parameter are refused, not silently merged", {
    .ui <- rxode2::rxUiDecompress(rxode2::rxode2(.sharedPhiMod))
    # both map to tka, so the model gets ONE phi1 column for the two of them
    expect_equal(.ui$saemEtaTrans, c(1L, 1L))
    expect_equal(sum(diag(.ui$saemModelOmega)), 1)
    # the second is the one with no column of its own
    expect_equal(.saemEtaNoPhi(.ui), "eta.cl")
    expect_error(.saemAssertEtaPhi(.ui), "eta.cl")
  })

  test_that("the gate does not refuse a model whose eta simply is not mu-referenced", {
    # Every one of these fits today: rxode2 records an eta that is not
    # `theta + eta` in `nonMuEtas`, which `saemParamsToEstimate()` appends to
    # the phi list, so the eta DOES own a phi1 column.  Refusing any of them
    # would be a regression, not a fix.
    .shared <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        eta.x ~ 0.6
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.x)
        cl <- exp(tcl + eta.x)
        v <- exp(tv)
        linCmt() ~ add(add.sd)
      })
    }
    .mixShared <- function() {
      ini({
        tka1 <- 0.45
        tka2 <- 0.8
        tcl <- 1
        tv <- 3.45
        p1 <- 0.5
        eta.ka ~ 0.6
        add.sd <- 0.7
      })
      model({
        ka <- mix(exp(tka1 + eta.ka), p1, exp(tka2 + eta.ka))
        cl <- exp(tcl)
        v <- exp(tv)
        linCmt() ~ add(add.sd)
      })
    }
    .muCov <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        cl.wt <- 0.1
        eta.cl ~ 0.3
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl + cl.wt * WT + eta.cl)
        v <- exp(tv)
        linCmt() ~ add(add.sd)
      })
    }
    for (.f in list(.shared, .mixShared, .muCov)) {
      .ui <- rxode2::rxUiDecompress(rxode2::rxode2(.f))
      expect_equal(.saemEtaNoPhi(.ui), character(0))
      expect_false(anyNA(.ui$saemOmegaTrans))
    }
  })

  # Two mixture components' random effects on ONE population parameter.  This
  # is the shape that reproduced #1047 end to end: saem ran every iteration and
  # then died with "subscript out of bounds" assembling the reported omega.
  .mixOnePhiMod <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      p1 <- 0.5
      eta.cl1 ~ 0.3
      eta.cl2 ~ 0.3
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka)
      cl <- mix(exp(tcl + eta.cl1), p1, exp(tcl + eta.cl2))
      v <- exp(tv)
      linCmt() ~ add(add.sd)
    })
  }

  # The spelling saem DOES support: one population parameter per component, so
  # each random effect owns a phi1 column and the pooled reporting omega
  # (Gamma2_phi1Report) has something to pool.
  .mixSplitMod <- function() {
    ini({
      tka <- 0.45
      tcl1 <- 1
      tcl2 <- 1.6
      tv <- 3.45
      p1 <- 0.5
      eta.cl1 ~ 0.3
      eta.cl2 ~ 0.3
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka)
      cl <- mix(exp(tcl1 + eta.cl1), p1, exp(tcl2 + eta.cl2))
      v <- exp(tv)
      linCmt() ~ add(add.sd)
    })
  }

  test_that("the #1047 model is refused before it fits, and its working twin is not", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    .ctl <- saemControl(nBurn = 5, nEm = 5, print = 0, calcTables = FALSE,
                        covMethod = "")
    # saemEtaNames() collapses the two onto one slot, so the kernel only ever
    # knew about the second -- the first was never sampled
    .ui <- rxode2::rxUiDecompress(rxode2::rxode2(.mixOnePhiMod))
    expect_equal(.ui$saemEtaTrans, c(2L, 2L))
    expect_equal(.ui$saemEtaNames, "eta.cl2")
    expect_error(
      suppressMessages(nlmixr2(.mixOnePhiMod, nlmixr2data::theo_sd, "saem", .ctl)),
      "eta.cl2")

    .ui2 <- rxode2::rxUiDecompress(rxode2::rxode2(.mixSplitMod))
    expect_equal(.ui2$saemEtaTrans, c(2L, 3L))
    expect_equal(.saemEtaNoPhi(.ui2), character(0))
    .fit <- suppressMessages(nlmixr2(.mixSplitMod, nlmixr2data::theo_sd, "saem", .ctl))
    expect_s3_class(.fit, "nlmixr2FitCore")
  })

  test_that(".getSaemOmega() reports the mismatch instead of running off the end", {
    # Gamma2_phi1 one column short of the UI's etas: the backstop for a
    # disagreement the up-front check did not catch.
    .env <- new.env(parent=emptyenv())
    .env$ui <- rxode2::rxUiDecompress(rxode2::rxode2(.muMod))
    .env$saem <- list(Gamma2_phi1=matrix(0.6, 1, 1))
    expect_error(.getSaemOmega(.env), "eta.cl")
    expect_error(.getSaemOmega(.env), "no variance")

    .env$ui <- rxode2::rxUiDecompress(rxode2::rxode2(.noPhiMod))
    expect_error(.getSaemOmega(.env), "eta.occ")

    # A missing Gamma2_phi1 has to count as zero columns.  `x > nrow(NULL)` is
    # logical(0), so comparing the ranks against it makes every eta look in
    # range and the message is lost again.
    .env$ui <- rxode2::rxUiDecompress(rxode2::rxode2(.muMod))
    .env$saem <- list()
    expect_error(.getSaemOmega(.env), "eta.ka, eta.cl")
  })
})
