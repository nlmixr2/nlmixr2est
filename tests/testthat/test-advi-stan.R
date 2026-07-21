## Optional external cross-check of est="advi" against Stan's own ADVI
## (rstan::vb()).  Uses a simple linear random-intercept model
##   y_ij = theta + eta_i + err,  eta_i ~ N(0, omega^2), err ~ N(0, sigma^2)
## which has the SAME variational objective in both tools, so the population
## posterior means should agree.  This needs rstan + a C++ toolchain and is
## therefore SKIP-BY-DEFAULT: set NLMIXR2_ADVI_STAN=true (and have rstan
## installed) to run it.  It is not part of any CI batch.

test_that("est='advi' population posterior agrees with rstan::vb()", {
  skip_on_cran()
  skip_if_not(identical(Sys.getenv("NLMIXR2_ADVI_STAN"), "true"),
              "set NLMIXR2_ADVI_STAN=true to run the Stan cross-check")
  skip_if_not_installed("rstan")

  ## simulate a linear random-intercept dataset
  .testSeed(42)
  nsub <- 60L; nobs <- 6L
  thetaTrue <- 5; omegaTrue <- 1; sigmaTrue <- 0.7
  eta <- stats::rnorm(nsub, 0, omegaTrue)
  dat <- do.call(rbind, lapply(seq_len(nsub), function(i) {
    data.frame(ID = i, TIME = seq_len(nobs),
               DV = thetaTrue + eta[i] + stats::rnorm(nobs, 0, sigmaTrue),
               EVID = 0, AMT = 0)
  }))

  ## est="advi" (full-Bayes) on the same model
  linmod <- function() {
    ini({ theta <- 4; eta ~ 1; add.sd <- 1 })
    model({ pred <- theta + eta; pred ~ add(add.sd) })
  }
  fA <- suppressMessages(suppressWarnings(
    nlmixr2(linmod, dat, est = "advi",
            control = adviControl(iters = 800L, print = 0L, returnAdvi = TRUE,
                                  pointEstimate = FALSE))))

  ## Stan ADVI (mean-field vb) on the identical generative model
  stanCode <- "
  data { int<lower=1> N; int<lower=1> J; int<lower=1,upper=J> id[N]; vector[N] y; }
  parameters { real theta; real<lower=0> omega; real<lower=0> sigma; vector[J] eta; }
  model {
    eta ~ normal(0, omega);
    y ~ normal(theta + eta[id], sigma);
  }"
  sm <- rstan::stan_model(model_code = stanCode)
  sdat <- list(N = nrow(dat), J = nsub, id = dat$ID, y = dat$DV)
  ## (vb() may emit a Pareto-k importance-resampling diagnostic warning; benign)
  vb <- suppressWarnings(
    rstan::vb(sm, data = sdat, algorithm = "meanfield", seed = 1, output_samples = 2000))
  post <- rstan::extract(vb)

  ## population posterior means agree (theta, residual sd, between-subject sd)
  expect_equal(unname(fA$theta[1]), mean(post$theta), tolerance = 0.15)
  expect_equal(unname(fA$theta[2]), mean(post$sigma), tolerance = 0.2)   # add.sd
  expect_equal(unname(sqrt(fA$popOmega[1])), mean(post$omega), tolerance = 0.25)
})
