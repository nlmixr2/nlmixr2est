# tools/vaeGenEncoderFixtures.R
# Generate golden encoder fixtures from the torch oracle for the VAE-NLME C++
# encoder (Phase 2). Implements plans/vae-encoder-spec.md exactly: single-layer
# unidirectional LSTM + FC head -> (mu, log_sigma, L), reparam z = mu + L eps.
# Dumps forward outputs and autograd gradients of a concrete scalar loss w.r.t.
# every encoder parameter so the native C++ forward/backward can be validated.
#
# Run:  Rscript tools/vaeGenEncoderFixtures.R
# Out:  tests/testthat/fixtures/vae/encoder_golden.rds

library(torch)
torch_manual_seed(20260709L)
set.seed(20260709L)

## ---- fixed tiny problem -------------------------------------------------
N     <- 2L    # subjects
Tmax  <- 3L    # max timesteps
xDim  <- 2L    # (time, DV)
hDim  <- 4L    # LSTM hidden
zDim  <- 3L    # individual params
nCov  <- 1L    # FC-head covariates
nOff  <- as.integer(zDim * (zDim - 1L) / 2L)   # strictly-lower entries
outDim <- 2L * zDim + nOff
lengths <- c(3L, 2L)                            # variable lengths (packing matters)

dataIn <- torch_randn(N, Tmax, xDim)            # standardized (time, DV)
covIn  <- torch_randn(N, nCov)                  # FC-head covariates
eps    <- torch_randn(N, zDim)                  # FIXED reparam noise (not resampled)

## ELBO surrogate targets (fixed): data term couples through z; DKL uses z_pop/omega
target <- torch_randn(N, zDim)
zPop   <- torch_tensor(c(0.2, -0.1, 0.3))
omega  <- torch_tensor(c(0.5, 0.8, 0.3))        # variances (>0)

## ---- modules (per spec) -------------------------------------------------
lstm   <- nn_lstm(input_size = xDim, hidden_size = hDim,
                  num_layers = 1L, batch_first = TRUE)
linear <- nn_linear(hDim + nCov, outDim)

## bias override: [h_inverse(mu0)=log(mu0) | sigma0(log) | zeros(nOff)]
mu0        <- torch_tensor(c(1.0, 0.5, 15.0))
sigma0log  <- torch_log(torch_tensor(c(1e-2, 5e-3, 1e-1)))   # already-logged, per reference
with_no_grad({
  linear$bias$copy_(torch_cat(list(torch_log(mu0), sigma0log, torch_zeros(nOff))))
})

## constant strictly-lower selection matrix (row-major, offset -1): [nOff, zDim*zDim]
trilFlat <- torch_zeros(nOff, zDim * zDim)
{
  idx <- 1L
  for (r in 2L:zDim) for (cc in 1L:(r - 1L)) {
    trilFlat[idx, (r - 1L) * zDim + cc] <- 1.0   # 0-based flat pos (r-1,cc-1) -> +1 for R
    idx <- idx + 1L
  }
}

## ---- forward ------------------------------------------------------------
encoderForward <- function() {
  packed <- nn_utils_rnn_pack_padded_sequence(dataIn, lengths,
                                              batch_first = TRUE, enforce_sorted = FALSE)
  res    <- lstm(packed)
  hn     <- res[[2]][[1]]            # h_n: [num_layers, N, hDim]
  hLast  <- hn[1, , ]               # [N, hDim] final hidden (last valid step via packing)
  combined <- torch_cat(list(hLast, covIn), dim = 2L)
  out    <- linear(combined)         # [N, outDim]

  mu       <- out[, 1:zDim]
  logSigma <- out[, (zDim + 1L):(2L * zDim)]
  Lmask    <- out[, (2L * zDim + 1L):outDim]                 # [N, nOff]
  Lflat    <- torch_matmul(Lmask, trilFlat)                  # [N, zDim^2]
  L        <- Lflat$view(c(N, zDim, zDim)) + torch_diag_embed(torch_exp(logSigma))
  z        <- mu + torch_matmul(L, eps$unsqueeze(-1L))$squeeze(-1L)
  list(mu = mu, logSigma = logSigma, L = L, z = z)
}

## ---- concrete scalar loss (mirrors ELBO structure) ----------------------
## loss = 0.5*sum((z - target)^2)          [data surrogate p_x_z through z]
##      + p_z(z, zPop, omega) - q_z(eps, diagL)   [DKL]
ln2pi <- log(2 * pi)
lossFun <- function(fw) {
  dataTerm <- 0.5 * torch_sum((fw$z - target)^2)
  pz <- 0.5 * torch_sum((fw$z - zPop)^2 / omega + torch_log(omega) + ln2pi)
  diagL <- torch_diagonal(fw$L, dim1 = 2L, dim2 = 3L)        # [N, zDim]
  qz <- 0.5 * torch_sum(eps^2 + ln2pi + 2 * torch_log(diagL))
  dataTerm + pz - qz
}

fw   <- encoderForward()
loss <- lossFun(fw)
loss$backward()

grad <- function(p) as.array(p$grad$to(dtype = torch_float64())$cpu())
arr  <- function(t) as.array(t$detach()$to(dtype = torch_float64())$cpu())

golden <- list(
  meta = list(N = N, Tmax = Tmax, xDim = xDim, hDim = hDim, zDim = zDim,
              nCov = nCov, nOff = nOff, outDim = outDim, lengths = lengths,
              gateOrder = "PyTorch i,f,g,o; weight_ih=[Wii;Wif;Wig;Wio]",
              seed = 20260709L, note = "diag L = exp(log_sigma); L offdiag unconstrained"),
  inputs = list(dataIn = arr(dataIn), covIn = arr(covIn), eps = arr(eps),
                target = arr(target), zPop = arr(zPop), omega = arr(omega),
                weight_ih = arr(lstm$weight_ih_l1), weight_hh = arr(lstm$weight_hh_l1),
                bias_ih   = arr(lstm$bias_ih_l1),   bias_hh   = arr(lstm$bias_hh_l1),
                fc_weight = arr(linear$weight),     fc_bias   = arr(linear$bias)),
  forward = list(mu = arr(fw$mu), logSigma = arr(fw$logSigma),
                 L = arr(fw$L), z = arr(fw$z), loss = as.numeric(arr(loss))),
  grads = list(weight_ih = grad(lstm$weight_ih_l1), weight_hh = grad(lstm$weight_hh_l1),
               bias_ih   = grad(lstm$bias_ih_l1),   bias_hh   = grad(lstm$bias_hh_l1),
               fc_weight = grad(linear$weight),     fc_bias   = grad(linear$bias))
)

outDir <- "tests/testthat/fixtures/vae"
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
saveRDS(golden, file.path(outDir, "encoder_golden.rds"))
cat("wrote", file.path(outDir, "encoder_golden.rds"), "\n")
cat("loss =", golden$forward$loss, "\n")
str(golden$forward)
