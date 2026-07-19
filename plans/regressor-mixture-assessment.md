# Assessment: regressor (fixed-membership) mixture approach for impmap / vae / advi

Context: SAEM gained `mixProbMethod="regress"` -- hard-classify each subject to
one mixture component once, hold it fixed, feed it in as a per-subject regressor
(`setIndMixest`), and solve each subject ONCE under its component instead of the
soft/marginalized per-component treatment (`nMix`-fold fewer ODE solves, more
stable, avoids the soft-EM collapse).  This assesses whether the same idea makes
sense for `est="imp"`/`"impmap"`, `est="vae"`, and `est="advi"`.

## est="imp" / "impmap" -- RECOMMENDED (best candidate)

Current mixture handling (`src/imp.cpp`): soft / marginalized.  It expands to
`nExp = nsub * Nmix` pseudo-subjects (`imp.cpp:590`), runs the importance-sampling
E/M-step per component, and forms posterior mixture responsibilities
`aMat` (nsub x Nmix, `imp.cpp:588`).  Mixtures are FORCED SERIAL
(`imp.cpp:606`, "the expanded pseudo-subjects' per-component solves race and the
proportions/Omega drift and collapse").

Fixed-membership regressor fit: strong.  Hard-classify each subject, sample/solve
only under its assigned component -> `nsub` work (not `nsub * Nmix`), drop the
responsibility E-step, and -- crucially -- REMOVE the per-component pseudo-subject
expansion that forces serial execution, so the mixture path could run PARALLEL
again.  Directly mirrors the SAEM `regress` pattern.

Payoff: ~`Nmix`-fold fewer inner solves AND restored parallelism (a second,
independent speedup the SAEM case did not have).  Same stability argument as SAEM
(avoids the documented proportion/Omega collapse).  Effort: moderate -- adapt the
imp EM loop to a fixed per-subject `mixest` and skip the responsibility update.

Verdict: WORTH DOING; clearest win of the three.

## est="vae" -- PLAUSIBLE, secondary

Current mixture handling (`R/vaeInner.R`, `R/vaeFit.R`, C++ `vaeInnerLik`):
marginalized.  It evaluates `nSub * nMix` ids and combines them with a `-2
logsumexp` over `mixProb` (`vaeInner.R:103-104`); consumes `nMix`/`mixProb` from
the ui (`vaeFit.R:191`); returns a per-subject `mixnum`.

Fixed-membership regressor fit: plausible but more entangled.  Hard-classify ->
evaluate `nSub` (assigned component) instead of `nSub * nMix`, feeding membership
as a decoder regressor; ~`nMix`-fold fewer inner evaluations.  BUT the VAE's
amortized encoder produces a per-subject posterior eta distribution and the
`logsumexp` marginalization is architecturally central -- collapsing to hard
membership changes the inference semantics (soft -> hard, like EM -> CEM) and
partly defeats the "amortized posterior" design.  It already returns `mixnum`, so
a hard-select in the driver is mechanically feasible.

Payoff: solve/eval reduction only (no parallelism story; VAE is already parallel).
Effort: moderate-high (touches the encoder/decoder combine).

Verdict: DO ONLY IF the `nSub * nMix` inner cost is a measured bottleneck;
otherwise the soft marginalization is a feature worth keeping.

## est="advi" -- NOT APPLICABLE as an optimization (no mixture support today)

Current mixture handling: NONE.  ADVI (`src/inner.cpp:10659+`) is mean-field
variational inference on the theta-sensitivity model; it does not reference
`nMix`/`mixProb`/`mixProbs` anywhere (R or C++), and `mix()` structure is not
consumed.  So there is nothing to "convert" -- ADVI simply does not do mixtures.

Fixed-membership regressor fit: this would be the way to ADD basic mixture
support, not to optimize an existing one.  A hard-classify-then-condition
approach (membership as a fixed regressor, VI conditioned on the assigned
component) is the SIMPLEST path to mixtures in ADVI -- it avoids needing a full
mixture variational family or a discrete-latent reparameterization.

Verdict: FEATURE WORK, not an optimization.  Pursue only if "mixtures under ADVI"
is a wanted capability; if so, the fixed-membership regressor is the pragmatic
first implementation.

## Summary

| estimator | today | regressor approach | recommend |
|-----------|-------|--------------------|-----------|
| imp/impmap | soft, nsub*Nmix, serial | fewer solves + restores parallelism | YES (best) |
| vae | marginalized logsumexp | fewer evals, but changes amortized semantics | MAYBE (secondary) |
| advi | no mixture support | would ADD mixtures (feature) | only if mixtures+ADVI wanted |
