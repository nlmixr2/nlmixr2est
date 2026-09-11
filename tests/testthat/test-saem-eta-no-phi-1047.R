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
  })
})
