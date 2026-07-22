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
GUROBI MIQP for the L0 selection (same optimum, no commercial solver).

## The question: parameters the encoder cannot reach

Our one statistical deviation is confined to a case the case studies do not
exercise: a structural population parameter with **no random effect**.

Such a parameter does not occupy the latent space, so the encoder never sees it
and the ELBO M-step has no route to it.  In a pharmacometrics package this comes
up constantly -- a fixed allometric exponent, a shared absorption parameter, an
`Emax` estimated without IIV -- so we had to give those parameters an estimation
path of their own, sitting alongside the variational machinery rather than inside
it.

We estimate them with a regression M-step scored against the **full FOCEi outer
objective**: the frozen-eta joint likelihood *plus* the Laplace determinant,
`0.5*log|Omega^-1|` and the DV-transform Jacobian.  Everything else -- the
encoder, the ELBO training step, the covariate branch-and-bound criterion --
follows your formulation exactly.

**The Laplace term is what makes a gradient available for these parameters.**
That is the real reason we added it, and it is worth being explicit about.  We
step them with an exact analytic outer gradient (Almquist 2015 sensitivity
equations, shared with our FOCEi implementation), which differentiates the
*marginal* (Laplace) likelihood.  Under the plain variational bound there is no
Laplace term, so there is nothing for that gradient to differentiate: dropping
the term does not merely change the target, it removes the gradient option and
leaves only a derivative-free search.  With the term in place the gradient goes
to the same Adam machinery that moves the encoder weights, so an unmatched theta
is learned on the same schedule as the rest of the model -- which feels closer to
the spirit of the method than pausing each M-step to run a separate optimizer to
convergence.

It also keeps the optimized quantity equal to the objective the fit reports,
which matters for us because `nlmixr2` reports the FOCEi objective across all of
its estimation methods.

We have made this switchable: `vaeControl(mStepObjective = "elbo")` reproduces
your plain-bound M-step (and correspondingly disables the gradient option).

## A comparison we withdrew

An earlier draft of this note raised a covariate-selection difference on Case
Study 2: our fit selects several effects where the paper reports fewer.  **We
have withdrawn that entirely, because the comparison was invalid.**

The repository ships *simulated* neonatal data (189 subjects, 1120 observations)
because, as the README states, the real data cannot be shared.  The paper's
Case Study 2 results are from the real cohort of N = 2425.  We were comparing a
fit on the 189-subject simulated set against a figure produced from 2425 real
neonates -- different data, so there was never anything to reconcile.  We are
recording this here only so that the omission is deliberate rather than
accidental.

For what it is worth, the exercise was still useful to us: chasing the
non-difference turned up several genuine mismatches on our side that we have
since fixed (an off-by-one in our smoothing gain, an omega update we smoothed
twice, encoder-input standardization computed over the observed values rather
than the padded matrix, and a residual model that estimated a proportional term
where yours fixes it at zero).

## One thing we would like your read on: the padded encoder inputs

While aligning our implementation to yours we found that the encoder-input
standardization is computed across the **whole padded observation matrix**:

```python
data_mean = data[:,:,1].mean()
data_std  = data[:,:,1].std()
data_in[:,:,1] = (data_in[:,:,1] - data_mean)/data_std
```

Because `data` is padded to the longest subject, every zero-padded cell of a
short subject enters both statistics.  On the neonatal data that is 1512 cells
against 1120 real observations, so the padding is the majority of the tensor and
moves the scale substantially:

| statistic | across padded matrix | across observed values only |
|---|---|---|
| mean | 2570.6 | 3470.3 |
| SD   | 1582.4 |  505.9 |

The SD differs by roughly a factor of three.  We had originally used the
observed-only statistics; we have since made the padded form the default so we
match you, and kept the other available as an option.

Our question is simply whether that is intended.  It is defensible -- it is a
fixed, deterministic rescaling of the encoder's inputs, the LSTM is masked by
`lengths` so the padded steps do not reach the likelihood, and the encoder can
absorb an input scaling into its first-layer weights.  But it does make the input
scale depend on how ragged the dataset happens to be, which is a property of the
sampling design rather than of the data, and it would change if the same study
were padded differently.  If it is deliberate we would like to understand the
reasoning; if it is incidental, it may be worth a note for reimplementers.

## Questions for you

1. Did parameters without a random effect come up in your work?  If so, how did
   you handle them -- fix them, give them a nominal IIV, or something else?  We
   would rather match your intent than invent a convention.
2. Do you see any objection to scoring *only* those parameters against the
   Laplace-corrected marginal while the latent parameters keep the variational
   bound?  Our concern is the obvious one: it is a mixed objective, and we would
   value a second opinion on whether that is principled or merely convenient.
3. Is the padded-matrix encoder standardization above intended?
4. Both `alpha_KL` and `alpha_pen` are built as `linspace(..., kl_iter)` and then
   indexed as `alpha_KL[iter]` with `iter` running from 1.  That skips element 0
   and reaches the final value one iteration early -- the first ramp value is
   never used.  Is that intended, or an off-by-one?
5. Would you be willing to share the *selected covariate set* from the simulated
   neonatal data in the repository (not the real cohort)?  That is the one
   artifact that would let anyone reimplementing check themselves against you
   end to end, and it involves no confidential data.

## Status in our documentation

The deviation, its narrow scope, and the bit-for-bit neonatal check are all
documented for users in the `est = "vae"` article, under *Objective -- a
deliberate deviation from the reference*, so nobody reproducing the paper is
surprised by it or misreads the covariate difference as caused by it.

## Before sending

* Have a maintainer read and sign it -- it should go from a person, not a tool.
* Question 3 is the cheapest one to resolve and may dissolve the whole covariate
  section -- consider leading with it, or asking it alone first.  There is no
  point presenting an elimination trail if the premise (that they selected
  exactly one covariate) is not what the paper claims.
* Confirm the reference-code line numbers and hyperparameters still match the
  version being cited.
