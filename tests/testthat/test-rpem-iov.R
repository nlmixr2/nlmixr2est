# IOV support for est="rpem", enabled by the "iov" attribute on nlmixr2Est.rpem: the
# shared .uiApplyIov hook materializes each occasion effect into unit-variance
# per-occasion etas (rx.iov.<v>.<occ>) scaled by a magnitude theta (iov.<v>).  RPEM
# recovers the IOV variance the SAEM way -- a single shared omega across a parameter's
# occasion etas (a "shared-omega occasion block") -- by fixing the magnitude theta at 1
# and estimating one pooled omega for the block, rather than a magnitude fixed effect.

.rpemIovMod <- function() {
  ini({ lka <- log(1.8); lke <- log(0.086); lV <- log(32)
    eta.ka ~ 0.3; eta.ke ~ 0.03; eta.V ~ 0.03
    iov.ka ~ 0.1 | occ
    add.err <- 0.7 })
  model({ ka <- exp(lka + eta.ka + iov.ka); ke <- exp(lke + eta.ke); V <- exp(lV + eta.V)
    d/dt(depot) = -ka * depot; d/dt(central) = ka * depot - ke * central
    cp <- central / V; cp ~ add(add.err) })
}

test_that("est=rpem declares IOV support and builds a shared-omega occasion block", {
  expect_true(isTRUE(attr(nlmixr2est:::nlmixr2Est.rpem, "iov")))
  dat <- nlmixr2data::theo_md; dat$occ <- 1L; dat$occ[dat$TIME >= 144] <- 2L
  res <- nlmixr2est:::.uiApplyIov(rxode2::assertRxUi(.rpemIovMod()), "rpem", dat, rpemControl())
  cl <- nlmixr2est:::.rpemClassify(res$ui)
  # id-level etas are mu-referenced; the materialized occasion etas are centered
  expect_equal(cl$muRef, c(TRUE, TRUE, TRUE, FALSE, FALSE))
  expect_true(all(is.na(cl$muNames[!cl$muRef])))
  # the two occasion etas share one omega group (shared-omega occasion block)
  .iov <- !cl$muRef
  expect_equal(length(unique(cl$omGroup[.iov])), 1L)
  expect_false(any(cl$etaFix[.iov]))                       # the shared omega is estimated
  # the IOV magnitude theta (iov.ka) is fixed at 1 (variance moves into the omega)
  expect_true("iov.ka" %in% cl$thetaNames[cl$iovMagIdx + 1L])
  expect_equal(cl$base[cl$iovMagIdx + 1L], 1)
})

test_that("est=rpem recovers a shared IOV omega (SAEM parity)", {
  skip_on_cran()

  # simulate IOV: id-level eta.ka ~ N(0, 0.09), per-occasion deviation ~ N(0, 0.2)
  simu <- rxode2::rxode2({ ka <- exp(lka + eka + iovk); cl <- exp(lcl); v <- exp(lv); cp <- linCmt() })
  set.seed(11); nsub <- 60L; obsT1 <- seq(0.5, 20, by = 2)
  dat <- do.call(rbind, lapply(seq_len(nsub), function(i) {
    eka <- stats::rnorm(1, 0, sqrt(0.09))
    do.call(rbind, lapply(1:2, function(oc) {
      iovk <- stats::rnorm(1, 0, sqrt(0.2))
      ev <- rxode2::et(amt = 100, cmt = "depot"); ev <- rxode2::et(ev, obsT1)
      s <- rxode2::rxSolve(simu, params = c(lka = 0.5, lcl = 1, lv = 3.45, eka = eka, iovk = iovk),
                           events = ev, returnType = "data.frame", addDosing = FALSE)
      d <- as.data.frame(ev); d$id <- i; d$occ <- oc; d$time <- d$time + (oc - 1) * 24
      o <- d$evid == 0; d$DV <- 0; d$DV[o] <- s$cp + stats::rnorm(nrow(s), 0, 0.1)
      d
    }))
  }))
  mod <- function() {
    ini({ lka <- 0.5; lcl <- fix(1.0); lv <- fix(3.45); eta.ka ~ 0.09; iov.ka ~ 0.15 | occ; add.sd <- 0.2 })
    model({ ka <- exp(lka + eta.ka + iov.ka); cl <- exp(lcl); v <- exp(lv); cp <- linCmt(); cp ~ add(add.sd) })
  }
  fit <- suppressMessages(suppressWarnings(nlmixr2(mod, dat, est = "rpem",
    control = rpemControl(nGauss = 300L, nMH = 60000L, mhBurn = 6000L,
                          niter = 30L, collect = 12L, seed = 1L, cores = 4L))))
  expect_s3_class(fit, "nlmixr2FitData")
  .om <- diag(fit$omega)
  .iov <- grepl("^rx\\.iov", names(.om))
  # the occasion etas report ONE shared omega (equal across occasions), finite/positive
  expect_gt(length(unique(round(.om[.iov], 8))), 0L)
  expect_equal(length(unique(round(.om[.iov], 8))), 1L)   # shared
  expect_true(all(is.finite(.om) & .om > 0))
  # id-level omega no longer absorbs the occasion variance (SAEM ~0.083, truth 0.09;
  # the magnitude-theta approach blew this up to ~0.23)
  expect_lt(.om[!.iov], 0.2)
  # the IOV magnitude theta is held at 1 (variance now in the shared omega)
  expect_equal(unname(fit$parFixedDf["iov.ka", "Estimate"]), 1)
})
