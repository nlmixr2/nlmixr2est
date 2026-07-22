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

## What we checked before writing -- and a correction to ourselves

We initially suspected this deviation was behind a difference in Case Study 2
(neonatal weight, 189 neonates): we reproduce the canonical gestational-age
effect on birth weight, but our run also retains GA on `kin` and a small sex
effect on birth weight, where the paper reports the single effect.

**That suspicion was wrong, and we would rather say so plainly than send you a
false lead.**  Every structural parameter in the neonatal model is
mu-referenced, so the model contains no unmatched theta and the objective option
has nothing to act on.  Re-running the case study under
`mStepObjective = "elbo"` reproduces our default run *bit for bit* -- identical
population parameters, identical omegas, identical residual error, identical
selected covariate set.  We have pinned that identity as a regression test.

We also ruled out the obvious second explanation, a too-short schedule.  At your
full schedule (`itersBurnIn = 100`, `klWarmup = 50`, `gammaIter = 250`,
`iters = 300`) our run selects *more* terms, not fewer -- five: GA on birth
weight (1.88), GA on `kin` (2.59), sex on birth weight (0.09), mother's age on
`kin` (-0.27) and parity on `kin` (0.12).

Reading your `pop_parameter` source we found a third candidate and tested that
too: **what the L0 selection regresses.**  The penalized criterion itself is
identical on both sides -- you minimize
`sum_squares(y - X beta) + alpha_pen * ln(nbatch) * sum(z)` with `y` and `X` both
scaled by `1/sqrt(omega_pop[k])`, which is exactly our
`RSS_S/omega + alpha*ln(N)*|S|`.  But your `y` is the SAEM sufficient statistic
`s1`, an exponential moving average of the posterior means
(`s1 <- s1 + gamma*(mu - s1)`), where ours was the current posterior mean.  We
have since switched to the smoothed statistic to match you (it is now our
default, independently of this question).

It makes essentially no difference: the same five terms at the full schedule, the
same three at the short one, coefficients agreeing to about 0.05.  In hindsight
that is what the gain schedule implies -- `gamma` is exactly 1 until
`gamma_iter` (250 of 300 iterations) in both implementations, so `s1` *equals*
the posterior mean for most of a run and is averaged only over the closing tail.

**So we cannot account for the difference.**  Three candidates tested, three
eliminated.  What we have not ruled out: how the candidate covariate columns are
constructed (centering/transformation), the omega values in force at selection
time, and the possibility that the published figure reports the headline effect
rather than the complete selected set -- which would mean there is no discrepancy
at all.

Where the two objectives *do* differ is on a model that carries an unmatched
theta.  On `theo_sd` with `tv` written without a random effect (FOCEi MLE
`tv` = 3.4293), a short 80-iteration schedule gives:

| M-step | `tv` |
|---|---|
| full outer objective + analytic gradient | 3.4286 |
| full outer objective + derivative-free   | 3.4360 |
| plain variational bound (your M-step)    | 3.4214 |

## Questions for you

1. Did parameters without a random effect come up in your work?  If so, how did
   you handle them -- fix them, give them a nominal IIV, or something else?  We
   would rather match your intent than invent a convention.
2. Do you see any objection to scoring *only* those parameters against the
   Laplace-corrected marginal while the latent parameters keep the variational
   bound?  Our concern is the obvious one: it is a mixed objective, and we would
   value a second opinion on whether that is principled or merely convenient.
3. For Case Study 2, does your reported result reflect the *complete* selected
   set, or the effect you chose to highlight?  If additional terms were selected
   and not shown, there may be no discrepancy here at all, and we would like to
   stop looking for one.
4. If it is the complete set: can you see what we are doing differently?  We are
   happy to share our run.  Our remaining suspects are the construction of the
   candidate covariate columns and the omega values in force at selection time,
   but we are guessing at this point.
5. Would a shared benchmark be of interest -- your Case Study 2 at the full
   schedule, same data -- so that the difference can be pinned to something
   specific rather than guessed at?

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
