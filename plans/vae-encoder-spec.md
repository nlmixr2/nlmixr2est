# VAE-NLME encoder specification (reverse-engineered from encoder.c)

Ground truth = the compiled Cython reference `~/src/vae_nlme/VAE/encoder.c`
(embedded `.pyx` source) plus `~/src/vae_nlme/functions_theo.py` (data prep) and
`~/src/vae_nlme/functions.py` (ELBO). This spec is what the native C++ encoder
(Phase 2) must implement; the torch oracle (`tools/vaeGenEncoderFixtures.R`)
generates golden forward outputs and gradients from the SAME spec for testing.

## 1. Modules

Inner `LSTM(in_dim, h_dim, out_dim, n_cov)`:
- `nn.LSTM(input_size=in_dim, hidden_size=h_dim, num_layers=1, batch_first=True)`
  -- SINGLE LAYER, UNIDIRECTIONAL. (prototype README's "bidirectional" is WRONG.)
- `nn.Linear(h_dim + n_cov, out_dim)`.
- forward(x, covariates, lengths):
  - `x = pack_padded_sequence(x, lengths, batch_first=True, enforce_sorted=False)`
  - `out, (hidden, cell) = self.lstm(x)`
  - `combined = cat((hidden[-1], covariates), dim=1)`  # hidden[-1] = final layer h_n
  - `return self.linear(combined)`                     # [N, out_dim]
  Because of packing, `hidden[-1]` is the hidden state at each subject's LAST
  VALID timestep -- equivalent to gathering the output row at index `length`.

Outer `LSTM_Encoder(x_dim, h_dim, z_dim, nbatch, n_cov, mu0, sigma0, h_inverse)`:
- `self.lstm = LSTM(x_dim, h_dim, z_dim + z_dim*(z_dim+1)/2, n_cov)`
  so `out_dim = z_dim + z_dim*(z_dim+1)/2 = 2*z_dim + n_off`, `n_off = z_dim*(z_dim-1)/2`.
- `self.tril_indices = tril_indices(row=z_dim, col=z_dim, offset=-1)`  # strictly lower
- bias override on the FC head:
  `linear.bias.data = cat( h_inverse(mu0)[z_dim], sigma0[z_dim], zeros(n_off) )`.
  Note `sigma0` is passed already log-transformed (theophylline.py:
  `sigma0 = tensor([1e-2, 5e-3, 1e-1]).log()`), so the diagonal init
  `exp(log_sigma)` starts at the raw small SDs.
- forward(data, covariates, lengths):
  - `out = self.lstm(data, covariates, lengths)`                      # [N, out_dim]
  - `mu        = out[:, 0:z_dim]`
  - `log_sigma = out[:, z_dim:2*z_dim]`
  - `L_mask    = out[:, 2*z_dim:]`                                    # [N, n_off]
  - `L = zeros(N, z_dim, z_dim)`
  - `L[:, tril_indices[0], tril_indices[1]] = L_mask`                # strictly-lower, unconstrained
  - `L = L + diag_embed(exp(log_sigma))`                             # diagonal = exp(log_sigma) > 0
  - `eps = randn_like(mu)`                                           # [N, z_dim]
  - `z = mu + (L @ eps.unsqueeze(-1)).squeeze(-1)`                   # reparameterization
  - returns (z, mu, L, log_sigma, eps)

`tril_indices(..., offset=-1)` lists lower entries in ROW-MAJOR order; e.g.
z_dim=3 -> (1,0),(2,0),(2,1), so L_mask[0]->L[1,0], L_mask[1]->L[2,0],
L_mask[2]->L[2,1] (0-indexed).

## 2. Weight initialization

- LSTM and Linear WEIGHTS: PyTorch DEFAULT init (no override in encoder.c).
  nn.Linear/LSTM default = uniform(-sqrt(k), sqrt(k)), k = 1/hidden_size (LSTM)
  or 1/in_features (Linear).
- Only the FC BIAS is overridden (section 1).
- The R prototype's `normal_(0, 1e-3)` weight shrink is a DEVIATION from the
  reference; do not replicate it. Sharp initial posterior comes from the bias
  (small sigma0), not weight shrinkage.

C++ note: we train our own encoder from scratch (never load torch weights), so
matching torch's exact internal gate/weight LAYOUT is not required for
correctness -- only that the LSTM math and its analytic gradient are correct.
We still adopt PyTorch's gate convention (below) so the torch oracle is a direct
cross-check.

## 3. LSTM cell (PyTorch convention adopted)

For input x_t, previous h_{t-1}, c_{t-1} (h_0 = c_0 = 0):
```
i = sigmoid(W_ii x + b_ii + W_hi h + b_hi)
f = sigmoid(W_if x + b_if + W_hf h + b_hf)
g = tanh   (W_ig x + b_ig + W_hg h + b_hg)
o = sigmoid(W_io x + b_io + W_ho h + b_ho)
c_t = f * c_{t-1} + i * g
h_t = o * tanh(c_t)
```
weight_ih_l0 = [W_ii; W_if; W_ig; W_io]  (4*h_dim x x_dim)
weight_hh_l0 = [W_hi; W_hf; W_hg; W_ho]  (4*h_dim x h_dim)
bias_ih_l0, bias_hh_l0 each length 4*h_dim (PyTorch keeps both bias vectors).

## 4. Input standardization (encoder inputs only)

From functions_theo.py load_data (data columns: time, DV, dose, covariates...):
- `data_in` = (time, DV) per obs, standardized:
  - time channel: `time / max(time)`         (global max over all obs)
  - DV channel:   `(DV - mean(DV)) / sd(DV)`  (global mean/sd over all obs)
- `covariates_in` (FC-head covariates), continuous cols standardized
  `(x - mean)/sd`; categorical (sex) left raw 0/1.
- `lengths` = obs count per subject (theophylline is balanced; general case
  variable). Encoder must honor per-subject lengths.

These encoder-input scalings are SEPARATE from the covariate encoding used to
build z_pop (section 6).

## 5. ELBO / loss (functions.py) -- minimized

N_tot = sum(lengths). Per subject sigma = a + b*pred (b=0 => additive).
```
p_x_z = 0.5*sum(err^2) + 0.5*N_tot*ln(2pi) + sum(log sigma),   err=(DV-pred)/sigma
p_z   = 0.5*sum( (z - z_pop)^2 / omega + log(omega) + ln(2pi) )   # z_pop = C @ beta
q_z   = 0.5*sum( eps^2 + ln(2pi) + 2*log(diag L) )
DKL   = p_z - q_z
loss  = p_x_z + alpha_KL * DKL      # alpha_KL ramps 0.01 -> 1 during annealing
```
This "loss" is the NEGATIVE ELBO (p_x_z is a negative log-likelihood). Encoder
params are updated by Adam on `loss`; z_pop/omega/a by the closed-form M-step.

## 6. Covariate encoding for z_pop (functions_theo.py initalize_C)

Per subject i, `z_pop_i = C[i] @ beta`, `C[i] = [ I(z_dim) | covariate blocks ]`:
- continuous covariate (e.g. weight, j != sex): entry = `log(cov_ij / cov_pop_j)`
  (`cov_pop` = mean over subjects, e.g. weight_pop).
- categorical (sex, j == 1): entry = raw `cov_ij` (0/1).
- `beta` = [z_pop(z_dim) ; vec(covariate effects)]; covariate columns pruned to
  zero = deselected (covariate selection, Phase 5).

## 7. Golden fixture contract (tools/vaeGenEncoderFixtures.R)

Fixed tiny example (N=2, T=3, x_dim=2, h_dim=4, z_dim=3, n_cov=1), fixed torch
seed, FIXED eps (not random), explicit weight tensors. Dump to
tests/testthat/fixtures/vae/encoder_golden.rds:
- inputs: data_in, covariates_in, lengths, all weight/bias tensors, eps.
- forward outputs: mu, log_sigma, L, z.
- gradients: d(scalar loss)/d(every encoder parameter) from torch autograd.
Phase 2 C++ encoder must match forward to ~1e-6 and analytic backward to torch
autograd (and finite differences) to ~1e-6.
