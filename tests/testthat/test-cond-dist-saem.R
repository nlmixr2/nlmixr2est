# addCondDist(method = "saem"): a SAEM fit's own chains, or the SAEM MCMC held
# at any fit's estimates (#1378).

test_that(".condDistSaemDraws subtracts the prior mean, masks and maps etas", {
  # N=2, nmc=2, niter=3, nphi=3 with phi1 columns 1 and 3
  .phiM <- array(seq_len(2 * 2 * 3 * 3), c(2, 2, 3, 3))
  .saem <- list(
    mprior_phi = matrix(c(1, 2, 0, 0, 10, 20), 2),
    etaUe = matrix(c(1, 1, 1, 0), 2)
  )
  attr(.saem, "saem.cfg") <- list(i1 = c(0L, 2L))
  .d <- .condDistSaemDraws(.saem, .phiM, 2:3, c(2L, 1L), c("eta.b", "eta.a"))
  expect_equal(dim(.d), c(2L, 2L, 4L))
  # subject 1, kernel eta 1 (phi column 1): chain 1 iter 2, chain 2 iter 2, ...
  expect_equal(.d[1, 2, ], c(.phiM[1, 1, 2, 1], .phiM[1, 2, 2, 1], .phiM[1, 1, 3, 1], .phiM[1, 2, 3, 1]) - 1)
  expect_equal(.d[1, 1, ], c(.phiM[1, 1, 2, 3], .phiM[1, 2, 2, 3], .phiM[1, 1, 3, 3], .phiM[1, 2, 3, 3]) - 10)
  # subject 2, kernel eta 2 (phi column 3) is not sampled
  expect_equal(.d[2, 1, ], rep(NA_real_, 4))
  # a shape mismatch gives NULL
  expect_null(.condDistSaemDraws(.saem, .phiM, 3:4, c(2L, 1L), c("eta.b", "eta.a")))
  expect_null(.condDistSaemDraws(.saem, .phiM, 2:3, 3L, "eta.a"))
})

test_that("addCondDist(method = 'saem') refuses chains with temporary etas", {
  .env <- new.env()
  .env$saemPseudoEtas <- "rx.eta.sd"
  .fit <- list(
    env = .env,
    eta = data.frame(ID = 1, eta.a = 0),
    est = "saem",
    omega = matrix(0.1, 1, 1, dimnames = list("eta.a", "eta.a"))
  )
  expect_error(addCondDist(.fit, method = "saem"), "temporary etas")
})

test_that("addCondDist(method = 'saem') agrees with importance sampling", {
  skip_on_cran()

  # nolint start: object_usage_linter.
  .mod <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      add.sd <- 0.7
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }
  # nolint end
  .dat <- nlmixr2data::theo_sd
  .n <- length(unique(.dat$ID))
  .agree <- function(a, b, z = 0.35) {
    expect_lt(mean(abs(a$MEAN - b$MEAN) / b$SD), z)
    expect_equal(a$SD, b$SD, tolerance = 0.3)
  }

  fit <- .nlmixr(.mod, .dat, est = "saem", control = saemControl(print = 0, nBurn = 50, nEm = 200))
  .theta <- fit$theta
  # the fit's own chains: 100 iterations of 3 chains, nothing simulated
  .chain <- addCondDist(fit, keepSamples = TRUE)$condDist
  expect_equal(attr(.chain, "method"), "saem")
  expect_equal(dim(attr(.chain, "samples")), c(.n, 3L, 300L))
  expect_true(all(.chain$SD > 0))
  # more draws than the fit's second phase holds: a held run
  .held <- addCondDist(fit, nsamp = 900)$condDist
  expect_equal(attr(.held, "method"), "mcmc")
  expect_equal(addCondDist(fit, nsamp = 900)$condDist, .held)
  expect_equal(fit$theta, .theta)
  .imp <- addCondDist(fit, method = "imp")$condDist
  .agree(.chain, .held)
  .agree(.held, .imp)

  # a FOCEi fit: the held SAEM run starts cold
  fit2 <- .nlmixr(.mod, .dat, est = "focei", control = foceiControl(print = 0))
  .mcmc <- addCondDist(fit2, method = "saem", nsamp = 600)$condDist
  expect_equal(attr(.mcmc, "method"), "mcmc")
  .agree(.mcmc, addCondDist(fit2)$condDist)
})
