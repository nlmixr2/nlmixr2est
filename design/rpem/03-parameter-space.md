# RPEM -- Parameter Space

## Parameter classification

Every population parameter is classified at UI-translation time into exactly one
class. The class determines how the M-step updates it.

1. **Random (eta) params** -- mu-referenced, carry between-subject variability.
   These define the mixture Gaussian: typical value -> `mu^(k)`, and their
   covariance -> `Sigma^(k)`. Updated by the conjugate formulas Eq 15-16, 18-21.
2. **Mixture weight params** -- `w^(k)`; updated by Eq 12/27. Only present when
   K>1 (deferred; see `07-mixtures.md`).
3. **Fixed effects (no BSV)** -- residual-error params (add/prop/combined/power,
   autocorrelation coefficients) and any structural or covariate coefficient
   that is NOT mu-referenced. These are the paper's `beta`. Updated numerically
   in the M-step (`05-mstep.md`), or by a conjugate form when one exists (e.g.
   the paper's `sigma^2` closed form for the pure `sigma^2*H_i` case).

The classifier is the crux of the "also support non-mu-referenced" decision
(D4). A model may mix all three classes freely.

## Transforms and the truncated-Gaussian problem (D5)

- Etas are sampled on nlmixr2's **unconstrained transformed scale** (log for
  positive params, logit for bounded, etc.), exactly as the model UI already
  declares. `theta_i = mu + eta` is formed on that scale, then back-transformed
  before the ODE solve.
- This sidesteps the paper's "Dealing with the truncated Gaussian" section
  entirely: no negative structural parameters are ever generated, so the
  normalization factor `N^(k)` is 1 and no draws are dropped.
- The paper's explicit truncation (sample untruncated, drop negatives,
  renormalize) is implemented **only if** a future spec adds raw-normal
  structural parameters without a transform. Not in scope for M1.

## Consequences for the M-step conjugate updates

The Eq 15-21 conjugate updates for `mu^(k)`/`Sigma^(k)` operate on the
transformed-scale `theta` (i.e. on the mu-referenced params). Because the
transform is deterministic, sampling and moment-matching on the transformed
scale is internally consistent -- the Gaussian mixture lives on the transformed
scale, which is also where SAEM's omega lives. This keeps RPEM's population
model aligned with SAEM/FOCEI reporting conventions.

## Covariates (mu2)

Covariate effects enter through the mu-referencing expressions (as in SAEM), read
from `ui$muRefCovariateDataFrame` (theta -> covariate -> covariateParameter). The
compiled rpem model already contains the covariate term (the covariate is a data
column in `rpemParams`, its coefficient a THETA), so the E-step solve handles the
covariate-adjusted per-subject mean with no special handling.

Estimation (D22): for a covariate on a mu-referenced (eta) param the coefficient
is estimated by generalizing the conjugate `mu` update (Eq 15) from a MEAN to a
weighted linear REGRESSION of the accepted `theta` samples on the per-subject
design matrix `X_i = [1, cov1_i, cov2_i, ...]`:
  `beta^(r+1) = (sum_acc X_i' X_i)^{-1} (sum_acc X_i' theta_ij)`,
  `omega^(r+1) = mean_acc (theta_ij - X_i beta)^2`.
The intercept-only case (no covariates) reduces to `mu = mean(theta)`, matching
the current `rpemMstepK1`. Covariate coefficients on params WITHOUT an eta are
pure fixed effects (class 3, numeric -- deferred with the other non-mu-ref
structural params, D20). The classifier builds `X_i` from
`muRefCovariateDataFrame` + the data; `.rpemBuildFit` reports the coefficients
(estimated, with SEs) rather than marking them held.

## Open items

- OI-1: Confirm the exact rxode2 transform metadata used to map eta-scale <->
  natural-scale so the classifier and back-transform match FOCEI/SAEM exactly.
- OI-2: Decide reporting scale for `Sigma^(k)` (report as omega on transformed
  scale, consistent with SAEM) -- assumed yes.
