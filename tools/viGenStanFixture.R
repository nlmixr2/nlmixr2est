## Regenerates tests/testthat/baselines/vi-stan-golden.rds -- the Stan ADVI (rstan::vb())
## cross-check reference for est="fbvi" in test-vi-stan.R.
##
## Needs rstan and a C++ toolchain; the test itself does NOT (that is the point).  Stan's
## posterior for a fixed model, fixed data and fixed seed does not depend on anything in
## this package, so it is frozen here rather than recompiled and re-fit on every run.
##
##     Rscript tools/viGenStanFixture.R
##
## Regenerate only when the Stan model or the simulated data changes.  The saved DV vector
## is checked by the test: if the data no longer matches, the test skips with a message
## rather than asserting against a reference built from different data.

if (!requireNamespace("rstan", quietly = TRUE)) stop("rstan is required to regenerate this fixture")

## The data must be simulated EXACTLY as test-vi-stan.R does, or the reference is
## meaningless.  Keep the two in sync; the DV vector saved below is what enforces it.
set.seed(42)
nsub <- 60L
nobs <- 6L
thetaTrue <- 5
omegaTrue <- 1
sigmaTrue <- 0.7
eta <- stats::rnorm(nsub, 0, omegaTrue)
dat <- do.call(rbind, lapply(seq_len(nsub), function(i) {
  data.frame(ID = i, TIME = seq_len(nobs),
             DV = thetaTrue + eta[i] + stats::rnorm(nobs, 0, sigmaTrue),
             EVID = 0, AMT = 0)
}))

stanCode <- "
  data { int<lower=1> N; int<lower=1> J; int<lower=1,upper=J> id[N]; vector[N] y; }
  parameters { real theta; real<lower=0> omega; real<lower=0> sigma; vector[J] eta; }
  model {
    eta ~ normal(0, omega);
    y ~ normal(theta + eta[id], sigma);
  }"

sm <- rstan::stan_model(model_code = stanCode)
sdat <- list(N = nrow(dat), J = nsub, id = dat$ID, y = dat$DV)
## vb() may emit a Pareto-k importance-resampling diagnostic warning; benign
vb <- suppressWarnings(
  rstan::vb(sm, data = sdat, algorithm = "meanfield", seed = 1, output_samples = 2000))
post <- rstan::extract(vb)

golden <- list(theta = mean(post$theta),
               sigma = mean(post$sigma),
               omega = mean(post$omega),
               dv = dat$DV,                       # data fingerprint, checked by the test
               stanVersion = as.character(utils::packageVersion("rstan")),
               seed = 1L)

.out <- file.path("tests", "testthat", "baselines", "vi-stan-golden.rds")
saveRDS(golden, .out)
cat("wrote", .out, "\n")
cat(sprintf("  theta = %.6f  sigma = %.6f  omega = %.6f  (rstan %s)\n",
            golden$theta, golden$sigma, golden$omega, golden$stanVersion))
