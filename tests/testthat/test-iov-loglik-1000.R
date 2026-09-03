# issue #1000: SAEM with IOV and a general-likelihood (ll()) endpoint returned
# an astronomically wrong objf (1.28e13).  The cause was NOT the objective
# calculation: refinePhi0Lik() froze the ODE for every distribution==4 fit on
# the premise that a general-likelihood phi0 parameter is always a likelihood
# SD the solve never sees.  The IOV magnitude theta is a phi0 parameter that
# DOES drive the structural model, so with the solve frozen the objective was
# exactly constant in it and the bounded optimizer walked to its upper bound --
# the magnitude diverged geometrically (1e17 at the issue's iteration counts)
# and the objf followed.
test_that("saem: an ll() endpoint with IOV gives a sane objf (#1000)", {
  skip_on_cran()

  .theoIov <- nlmixr2data::theo_md
  .theoIov$occ <- 1
  .theoIov$occ[.theoIov$TIME >= 144] <- 2

  .ll <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      lsd <- 0.7
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      iov.cl ~ 0.1 | occ
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl + iov.cl)
      v <- exp(tv + eta.v)
      sd <- exp(lsd)
      cp <- linCmt()
      ll(err) ~ -lsd - 0.5 * log(2 * pi) - 0.5 * ((DV - cp) / sd)^2
    })
  }

  .add <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      add.sd <- exp(0.7)
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      iov.cl ~ 0.1 | occ
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl + iov.cl)
      v <- exp(tv + eta.v)
      cp <- linCmt()
      cp ~ add(add.sd)
    })
  }

  .ctl <- saemControl(nBurn = 60, nEm = 80, seed = 42L, print = 0L,
                      covMethod = "", calcTables = FALSE)
  .fL <- suppressWarnings(nlmixr2(.ll, .theoIov, est = "saem", control = .ctl))
  .fA <- suppressWarnings(nlmixr2(.add, .theoIov, est = "saem", control = .ctl))

  # the defect: an objf eleven orders of magnitude too large
  expect_true(is.finite(.fL$objf))
  expect_true(.fL$objf > 0)
  expect_true(.fL$objf < 10 * .fA$objf)

  # the two forms are the same model written two ways, so the structural
  # estimates have to land in the same place
  expect_equal(.fL$fixef[["tka"]], .fA$fixef[["tka"]], tolerance = 0.2)
  expect_equal(.fL$fixef[["tcl"]], .fA$fixef[["tcl"]], tolerance = 0.05)
  expect_equal(.fL$fixef[["tv"]], .fA$fixef[["tv"]], tolerance = 0.02)
  expect_equal(exp(.fL$fixef[["lsd"]]), .fA$fixef[["add.sd"]], tolerance = 0.2)

  # the mechanism, not just the symptom: the IOV magnitude is a phi0 parameter
  # refined by refinePhi0Lik().  Assert its whole iteration history stayed
  # bounded -- unfixed it grows geometrically past 1e50 within 60 iterations.
  # (the magnitude's own column cannot be named from the finished fit -- by then
  # .uiFinalizeIov() has turned it back into an eta -- so bound the whole
  # history; every other entry is a theta or a variance of order 1)
  .ph <- .fL$env$saem$par_hist
  expect_true(all(is.finite(.ph)))
  expect_true(max(abs(.ph)) < 1e3)

  # and the reported inter-occasion variance is a plausible variance
  expect_true(is.finite(.fL$omega$occ[1, 1]))
  expect_true(.fL$omega$occ[1, 1] < 1)
})
