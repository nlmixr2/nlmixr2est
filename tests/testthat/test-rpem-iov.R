# IOV support for est="rpem", enabled by the "iov" attribute on nlmixr2Est.rpem:
# the shared .uiApplyIov hook materializes each occasion-level random effect into
# unit-variance (fix=TRUE) per-occasion etas plus a structural magnitude theta.
# RPEM then fits it as a plain diagonal model where the occasion etas are
# non-mu-referenced (typical value fixed at 0) with their omega held at 1.  The goal
# here is that IOV models RUN end-to-end and are structured correctly, not that the
# occasion variance is recovered (that needs shared-omega occasion blocks -- future).

.rpemIovMod <- function() {
  ini({ lka <- log(1.8); lke <- log(0.086); lV <- log(32)
    eta.ka ~ 0.3; eta.ke ~ 0.03; eta.V ~ 0.03
    iov.ka ~ 0.1 | occ
    add.err <- 0.7 })
  model({ ka <- exp(lka + eta.ka + iov.ka); ke <- exp(lke + eta.ke); V <- exp(lV + eta.V)
    d/dt(depot) = -ka * depot; d/dt(central) = ka * depot - ke * central
    cp <- central / V; cp ~ add(add.err) })
}

test_that("est=rpem declares IOV support and classifies occasion etas as centered", {
  expect_true(isTRUE(attr(nlmixr2est:::nlmixr2Est.rpem, "iov")))
  dat <- nlmixr2data::theo_md; dat$occ <- 1L; dat$occ[dat$TIME >= 144] <- 2L
  res <- nlmixr2est:::.uiApplyIov(rxode2::assertRxUi(.rpemIovMod()), "rpem", dat, rpemControl())
  cl <- nlmixr2est:::.rpemClassify(res$ui)
  # id-level etas are mu-referenced; the materialized occasion etas are centered
  expect_equal(cl$muRef, c(TRUE, TRUE, TRUE, FALSE, FALSE))
  expect_true(all(is.na(cl$muNames[!cl$muRef])))                 # no typical-value theta
  expect_true(all(cl$etaFix[!cl$muRef]))                         # unit variance held
  # the IOV magnitude is a non-mu-ref structural theta
  expect_true("iov.ka" %in% cl$thetaNames[cl$structIdx + 1L])
})

test_that("est=rpem runs an IOV model end-to-end", {
  skip_on_cran()
  skip_on_ci()  # heavy: ODE solve + multi-iteration RPEM + FOCEI eval

  dat <- nlmixr2data::theo_md; dat$occ <- 1L; dat$occ[dat$TIME >= 144] <- 2L
  fit <- suppressMessages(suppressWarnings(nlmixr2(.rpemIovMod(), dat, est = "rpem",
    control = rpemControl(nGauss = 200L, nMH = 40000L, mhBurn = 4000L,
                          niter = 12L, collect = 6L, seed = 1L, cores = 4L))))
  # a full fit object comes back
  expect_s3_class(fit, "nlmixr2FitData")
  .pf <- fit$parFixedDf
  expect_true(all(is.finite(.pf[c("lka", "lke", "lV", "add.err"), "Estimate"])))
  # the occasion-eta variance is held at 1 (magnitude carried by iov.ka), and the
  # id-level omegas are estimated (finite, positive)
  .om <- diag(fit$omega)
  .iov <- grepl("^rx\\.iov", names(.om))
  expect_true(all(.om[.iov] == 1))
  expect_true(all(is.finite(.om[!.iov]) & .om[!.iov] > 0))
})
