# Draft: contacting the VAE-NLME authors about an objective deviation

Prepared for Jan Rohleff, Freya Bachmann, Uri Nahum, Dominic Braem, Britta
Steffens, Marc Pfister, Gilbert Koch, Johannes Schropp -- authors of
*Redefining Parameter Estimation and Covariate Selection via Variational
Autoencoders: One run is all you need*, CPT:PSP 2025,
doi:10.1002/psp4.70129, reference code <https://github.com/janrohleff/vae_nlme>.

This is a DRAFT for a human to review, edit and send.  It is not correspondence
that should go out unreviewed -- it makes a technical claim about someone else's
published method and asks them to weigh in.

---

## What we did

`nlmixr2` (R, CRAN) now ships `est = "vae"`, a native reimplementation of the
method, with matched default hyperparameters (`itersBurnIn = 100`,
`klWarmup = 50`, `gammaIter = 250`, `iters = 300`, `nGradStep = 5`,
`hiddenDim = 25`, Adam 8e-3 burn-in / 5e-3, KL anneal 0.01 -> 1, covariate
penalty ramp 2 -> 1).  Differences are mostly engineering: an Armadillo LSTM with
analytic BPTT instead of Torch autograd, an `rxode2`-compiled ODE decoder instead
of a hand-written per-case-study decoder, and exact branch-and-bound instead of a
GUROBI MIQP for the L0 selection (same selected subset, no commercial solver).

## The deviation we want to flag

One difference is **statistical, not engineering**, and we would value the
authors' view.

The reference trains and runs its M-step on the plain variational bound

    ELBO = p(x|z) + [ p(z) - q(z|x) ]

(`Main/theophylline.py:129-136`, `Main/neonates.py`).  There is no
Laplace/Hessian term in training, in the population M-step, or in the covariate
selection regression -- the encoder entropy `q(z|x)` fills that role, which is
exactly what a variational method is for.  `torch.linalg.slogdet` occurs once in
the whole codebase (`functions.py:106`), inside `LogLikelihood_linearization`,
which is evaluated only at the end to report the OFV/AIC/BIC.

`nlmixr2` instead scores the population M-step against the **full FOCEi outer
objective**: the Laplace determinant, `0.5*log|Omega^-1|`, and the DV-transform
Jacobian.  Two reasons:

1. **Ecosystem consistency.** The objective an `nlmixr2` fit reports is the FOCEi
   objective.  Having the M-step optimise something else means the reported OFV
   is not the quantity being optimised.
2. **Internal consistency.** We added an option that steps the population
   parameters with the exact analytic outer gradient (Almquist 2015 sensitivity
   equations, shared with our FOCEi `fast=TRUE` path).  That gradient
   differentiates the marginal (Laplace) likelihood; if the M-step optimised the
   plain ELBO the two would be optimising different functionals.

## Why we are writing

Covariate selection is scored by a BICc-ELBO criterion, so **changing the
objective can change the selected covariate set.**  On Case Study 2 (neonatal
weight, 189 neonates) we reproduce the canonical result -- gestational age on
birth weight is selected -- but our fit also retains GA on `kin` and a small sex
effect on birth weight, where the paper reports the single effect.  Indicative
coefficients from our run (short demonstration schedule):

| coefficient | estimate |
|---|---|
| `beta_lW0_GA`   | 1.858 |
| `beta_lkin_GA`  | 1.880 |
| `beta_lW0_SEX`  | 0.090 |

We do not claim this is a defect in the published method.  A variational
objective legitimately has no Laplace term, and the paper's parsimony may well be
the better behaviour.  But it does mean a reader comparing the two
implementations on the same data can see different covariate sets, and we would
rather that be documented and understood than discovered.

## Questions for the authors

1. Was scoring the L0 selection against the plain ELBO (rather than a
   Laplace-corrected marginal) a deliberate choice?  If so, was the parsimony
   difference something you observed?
2. Do you have a view on which objective *should* score the BICc penalty term?
   The penalty is `alpha * ln(N) * |S|` in both, but the goodness-of-fit term it
   is traded against differs between our implementations.
3. Would you be interested in a side-by-side on Case Study 2 -- same data, same
   schedule, the two objectives -- to see how much of the difference is the
   objective and how much is the stochastic trajectory?

## Status in our documentation

The deviation is documented for users in the `est = "vae"` article, under
*Objective -- a deliberate deviation from the reference*, including the
consequence for covariate selection, so nobody reproducing the paper is surprised
by it.

## Before sending

* Have a maintainer read and sign it -- it should go from a person, not a tool.
* Re-run Case Study 2 at the paper's full schedule (not the shortened
  demonstration settings) and quote those coefficients instead of the indicative
  ones above; the sex effect in particular is small enough that it may not
  survive a longer run, and it would be careless to raise it if it does not.
* Confirm the reference-code line numbers still match the version being cited.
