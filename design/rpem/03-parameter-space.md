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

## Covariates

Covariate effects enter through the mu-referencing expressions (as in SAEM):
covariate coefficients that are mu-referenced ride along in class 1; those that
are not are class 3 (numeric). Covariate handling reuses
`preProcessCovariatesPresent`.

## Open items

- OI-1: Confirm the exact rxode2 transform metadata used to map eta-scale <->
  natural-scale so the classifier and back-transform match FOCEI/SAEM exactly.
- OI-2: Decide reporting scale for `Sigma^(k)` (report as omega on transformed
  scale, consistent with SAEM) -- assumed yes.
