# Control options for the impmap (importance-sampling EM) estimation method

A NONMEM-style Monte Carlo importance-sampling EM built on the
mu-referenced FOCEI MAP. The proposal density for each subject is
centered at the MAP mode (\`muModel="lin"\`); mu-referenced population
parameters are updated by the EM gradient, while non-mu parameters
(structural and residual error) are updated by a symbolic-sensitivity
Newton step – the importance-sampling-weighted score and Gauss-Newton
information built from the analytic \`d(f)/d(theta)\` and
\`d(V)/d(theta)\` (exact censored partials for BLQ/M2/M3/M4 points).

## Usage

``` r
impmapControl(
  sigdig = 3,
  ...,
  isample = 300L,
  nIter = 100L,
  mapIter = 1L,
  nBurn = 0L,
  burnFreezeOmega = FALSE,
  gamma = 1,
  gammaMethod = c("auto", "global", "individual"),
  gammaRule = c("target", "floor"),
  df = 0,
  proposal = c("auto", "normal", "t", "laplace", "mixture"),
  propMixScale = c(1, 9),
  propMixWeight = c(0.9, 0.1),
  auto = TRUE,
  autoNonmemSparse = FALSE,
  autoDfPatience = 2L,
  iscaleMin = 0.1,
  iscaleMax = 10,
  iaccept = 0.4,
  ctol = NULL,
  nConvWindow = 10L,
  impSeed = 42L,
  covMethod = c("imp", "analytic", "r,s", "r", "s", "sa", ""),
  qr = FALSE,
  qrShift = TRUE,
  qrRefresh = TRUE,
  qrScramble = c("none", "owen", "lms"),
  sir = FALSE,
  sirSample = NULL,
  muModel = c("lin", "none"),
  combSens = TRUE
)
```

## Arguments

- sigdig:

  Optimization significant digits. One value drives, with a single
  consistent formula, the inner/outer optimizer convergence tolerance
  (`10^-sigdig`), the boundary check tolerance (`5*10^(-sigdig+1)`), and
  the ODE solver tolerances: the `rtol` exponent IS `sigdig` and `atol`
  sits three orders below, so `rtol = 10^-sigdig`,
  `atol = 10^(-sigdig-3)` for every solver (stiff, non-stiff or
  auto-switching). The sensitivity (`atolSens`/`rtolSens`) tolerances
  match the main solve (the outer gradient and covariance are built from
  them); the steady-state (`ssAtol`/`ssRtol`) tolerances run one order
  looser. Keying the optimizer to the same `10^-sigdig` means it
  converges to exactly the precision the solve supports. At the default
  `sigdig = 3` this is `atol = 1e-6`, `rtol = 1e-3`.

- ...:

  Parameters used in the default \`foceiControl()\`

- isample:

  Number of importance samples drawn per subject per iteration (NONMEM
  ISAMPLE). Either a single count used for every subject, or a vector of
  length \`nsub\` giving a count \*\*per subject\*\*.

  Per-subject counts are the NM7 Technical Guide's own remedy for poor
  coverage (its derivation is Gaussian throughout and never mentions a t
  proposal): a subject whose weights are badly behaved can be given more
  samples without charging every other subject for them. Note this
  treats the symptom rather than the cause – more draws from a proposal
  whose tails are too light still gives weights with infinite variance,
  which \`fit\$env\$impPsisK\` will show. See \`df\` for the shape-based
  remedy.

- nIter:

  Maximum number of importance-sampling EM iterations.

- mapIter:

  MAP-assist period, in EM iterations. \`1\` (default) re-centers the
  proposal at each subject's MAP mode every iteration; \`k \> 1\`
  re-centers every \`k\`th iteration; \`0\` never re-centers after the
  startup MAP pass.

  The MAP search is the mu-referenced FOCEI inner problem and is the
  dominant per-iteration cost of \`est="impmap"\`, so raising
  \`mapIter\` trades proposal accuracy for speed. It is worth raising
  once the population parameters are moving slowly enough that the mode
  barely shifts between iterations, and not before – a proposal centered
  away from the mode costs effective sample size, which
  \`fit\$env\$impNeffFrac\` will show.

  A skipped iteration does not freeze the proposal center: the M-step
  reseeds every subject's eta with its conditional mean, so the proposal
  still moves each iteration – it is the MAP re-optimization, not the
  center, that is skipped. The covariance is still built from the inner
  Hessian there, which is what distinguishes this from \`est="imp"\`;
  \`est="imp"\` uses the running conditional variance instead and is not
  the \`mapIter = 0\` case.

  \`mapIter\` cannot change an \`est="imp"\` fit, which never re-centers
  at all.

  \*\*Measured trade-off.\*\* \`mapIter = 3\` against the default, same
  harness: theta RMSE 0.00102 vs 0.00157 (1 ETA), 0.00130 vs 0.00144 (3
  ETAs), 0.00173 vs 0.00173 (general \`ll()\`), with \`Omega\` RMSE,
  effective sample size and k-hat essentially unmoved. On those fixtures
  skipping two MAP searches in three is close to free, which is the
  point of the knob.

  At 8 ETAs it becomes a genuine trade rather than a free saving, and in
  a direction worth knowing: theta RMSE 0.0205 against 0.0146 (worse),
  but Pareto k-hat +0.88 against +0.97 and the number of subjects above
  0.7 down from 2.1 to \*\*0.25\*\* – by a wide margin the healthiest
  sampler of any setting measured on that fixture. Re-optimizing the
  mode every iteration is not automatically the safer choice at high eta
  dimension.

  \`1\` stays the default: it is the previous behaviour, and it never
  runs a stale proposal. Raise it when the per-iteration MAP cost
  matters (it is the inner FOCEI problem, so it scales with how hard
  that is), and watch \`fit\$env\$impNeffFrac\` – a proposal centred
  away from the mode costs effective sample size on a model whose mode
  moves quickly.

- nBurn:

  Number of burn-in EM iterations run before the \`nIter\` budget, to
  let the proposal-scale (\`gamma\`) and \`auto\` controllers settle
  before the estimates they influence are judged. \`0\` (default) is no
  burn-in.

  These are EXTRA iterations, not carved out of \`nIter\`, so raising
  \`nBurn\` never silently shortens the fit; a fit runs at most
  \`nBurn + nIter\` iterations. Everything runs normally during them –
  the E-step, both controllers, and the theta M-step – except that
  convergence cannot fire, and \`Omega\` is held when \`burnFreezeOmega
  = TRUE\`.

  Burn-in iterations are the first \`nBurn\` rows of \`\$parHist\` and
  of \`fit\$env\$impObjTrace\`; \`fit\$env\$impNburn\` reports how many
  there were.

  \*\*Measured trade-off.\*\* \`nBurn = 5, burnFreezeOmega = TRUE\`
  against the default, same harness as \`proposal\`: theta RMSE 0.00140
  vs 0.00157 (1 ETA), 0.00185 vs 0.00144 (3 ETAs), 0.00104 vs 0.00173
  (general \`ll()\`); \`Omega\` RMSE worse on all three (0.00404 /
  0.00261 / 0.00249 against 0.00344 / 0.00184 / 0.00170), for 2-5 extra
  iterations.

  \*\*It can be much worse than that, so do not switch it on
  speculatively.\*\* On an 8-ETA warfarin PK/PD fit the same setting
  gave theta RMSE 0.149 against 0.0146 for the default – ten times worse
  – with k-hat and the failing-subject count both degraded (+1.30 / 3.1
  against +0.97 / 2.1). Holding \`Omega\` while the thetas move is a
  coherent thing to want, but on a high-dimensional model the thetas can
  travel a long way from an \`Omega\` that is not allowed to follow, and
  the fit does not recover the difference in the iterations that remain.

  Off by default, and worth trying only when a fit's \`Omega\` visibly
  overshoots in the first rows of \`\$parHist\` and then has to come
  back – which is the specific failure it addresses.

- burnFreezeOmega:

  When \`TRUE\`, hold \`Omega\` at its starting value for the \`nBurn\`
  burn-in iterations while the structural and residual-error thetas
  update normally. Ignored when \`nBurn = 0\`.

  The case for it: the first EM iterations run under a proposal whose
  scale the controller has not yet adapted, and \`Omega\`'s update is
  the one that absorbs that – it is built from the conditional variances
  (\`Omega = mean(eta eta' + condVar)\`), which an over- or
  under-dispersed proposal estimates badly. Freezing it lets the thetas
  move the model toward the data while the proposal settles, instead of
  chasing an \`Omega\` excursion the fit then has to undo.

  Convergence is not tested until the whole trailing \`nConvWindow\`
  lies past the burn-in, so a frozen \`Omega\` cannot be mistaken for a
  settled one.

- gamma:

  Initial proposal-variance inflation factor (NONMEM ISCALE); the
  proposal covariance is \`gamma\` times the inverse of the inner
  information matrix at the mode.

- gammaMethod:

  How the proposal scale \`gamma\` is adapted during the EM.

  \`"auto"\` (default) picks per model: \`"individual"\` when the model
  is not transformably normal – a general log-likelihood (\`ll()\`)
  endpoint, or a count/categorical/time-to-event distribution – and
  \`"global"\` otherwise. That is the split the hypothesis actually
  rests on: \`gamma = 1\` is already the efficient proposal when the
  individual posterior is close to Gaussian, so a normal model gains
  nothing from per-subject adaptation and would only pay for it in
  effective sample size, while a general-likelihood model is exactly
  where the posteriors go non-Gaussian and per-subject coverage starts
  to vary. The test is \`all(ui\$predDf\$distribution == "norm")\`, the
  same line \[rxode2::assertRxUiTransformNormal()\] draws. The resolved
  value is reported in the fit's \`\$runInfo\`.

  \`"global"\` keeps one scale shared by every subject, inflated (never
  relaxed) only when the \*mean\* Kish effective-sample fraction falls
  below \`iaccept\`. It leaves \`gamma\` at its efficient starting value
  while coverage is healthy, which is the right behaviour when the
  individual posteriors are close to Gaussian.

  \`"individual"\` gives every subject its own \`gamma_i\` and adapts it
  two-sided toward a target on that subject's \`xi_i\`, clamped to
  \`\[iscaleMin, iscaleMax\]\`. This follows NONMEM's \*objective\* –
  the NM7 Technical Guide states that \`gamma\` is per-subject (eq.
  1.90) and is "continually adjusted so that xi_i approximates IACCEPT"
  (note after eq. 1.76), bounded by \`ISCALE_MIN\`/\`ISCALE_MAX\` – but
  the guide publishes no update formula, so the functional form used
  here is nlmixr2's own. Prefer this when the individual posteriors are
  heavy-tailed or the design is heterogeneous: a global scale is driven
  by the mean, so a handful of badly-covered subjects never trip it and
  their likelihood and \`Omega\` contributions end up carried by a few
  samples.

  The two modes report \*different\* efficiency statistics –
  \`"individual"\` targets \`xi\` (the mean normalized importance
  weight, NONMEM's \`IACCEPT\` quantity) while \`"global"\` targets the
  Kish effective-sample fraction. These are not comparable; the fit's
  \`\$runInfo\` says which is in force.

- gammaRule:

  How the SHARED (\`gammaMethod="global"\`) proposal scale is adapted.
  Ignored for \`gammaMethod="individual"\`, which always follows
  NONMEM's two-sided per-subject rule.

  \`"floor"\` treats \`iaccept\` as a one-sided FLOOR on the mean Kish
  effective-sample fraction: \`gamma\` stays at its efficient starting
  value while coverage is healthy and is inflated only when coverage
  drops below the floor. It never comes back down.

  \`"target"\` (default) follows the NM7 Technical Guide, which says
  \`gamma\` is "continually adjusted so that xi_i approximates IACCEPT"
  – adapted BOTH ways, on \`xi\` rather than the Kish fraction, using
  the same analytic inversion as the per-subject controller (\`gamma \*
  (xi/iaccept)^(2/neta)\`, capped 1.25x each way and clamped to
  \`\[iscaleMin, iscaleMax\]\`).

  The two rules settle at different operating points, so they are not
  interchangeable, and constants tuned against one do not carry over to
  the other.

  \*\*Measured trade-off\*\* (theophylline, 6 seeds, \`isample = 300\`,
  RMSE against an \`isample = 6000\` reference):

  - 3 ETAs – \`"floor"\`: theta RMSE 0.00185, \`Omega\` RMSE 0.00298,
    max Pareto k-hat 0.604, 0.33 subjects above 0.7, converged 100
    \`"target"\`: theta RMSE 0.00257, \`Omega\` RMSE 0.00378, max k-hat
    -0.491, 0 subjects above 0.7, converged 0

  - 8 ETAs – \`"floor"\` never adapts at all (\`gamma\` pinned at 1.0
    with \`xi\` 1.35, i.e. a proposal far too narrow), because it
    watches the mean Kish fraction, which stays healthy while \`xi\`
    says the proposal is wrong. \`"target"\` moves \`gamma\` to 1.25 and
    puts \`xi\` on 0.41.

  So \`"target"\` fixes the tail and is the only rule that reacts at
  high ETA dimension, but it pays for it: putting \`xi\` on \`iaccept\`
  deliberately widens the proposal (\`gamma\` 1.79 on the 3-ETA
  fixture), which costs effective sample size (0.96 to 0.70) and adds
  Monte-Carlo noise to the objective – enough that the fit often does
  not meet \`ctol\` within \`nIter\`.

  \`"target"\` is the default: tuned against tuned, it wins every column
  at 3 ETAs and is the only rule that adapts at all at 8 ETAs. Its cost
  is on the 1-ETA fixture, where it roughly doubles theta RMSE (0.00113
  to 0.00224) and takes about twice as many iterations. That trade
  follows the same weighting \`auto\` uses: weights with infinite
  variance are a correctness problem whose error is unbounded, while the
  extra Monte-Carlo noise is bounded and measurable.

  Choose \`"floor"\` for the previous behaviour – a proposal left at its
  efficient starting value while coverage is healthy. It is also what
  the tail-machinery tests pin, because \`"target"\` repairs the tail
  itself and leaves the \`df\` ladder nothing to fix.

  \*\*The tuned constants travel with the rule.\*\* Selecting a rule
  also selects its tuned defaults (\`nConvWindow\` 20 for \`"target"\`,
  10 for \`"floor"\`), so switching to the NONMEM method does not
  silently run NONMEM's law on the other rule's tuning. An explicitly
  supplied value always wins.

- df:

  Degrees of freedom of the importance-sampling proposal (NONMEM
  \`DF\`). \`0\` (default) uses a multivariate \*\*normal\*\* proposal;
  any value \`\> 0\` uses a multivariate \*\*t\*\* with that many
  degrees of freedom.

  This changes the proposal's SHAPE rather than its width, and that is
  the distinction that matters. \`gamma\` can only make a Gaussian
  proposal wider; it cannot give it heavier tails. When the target
  posterior has heavier tails than the proposal, the importance weights
  have infinite variance – and neither \`xi\` nor the Kish effective
  sample size can detect it, because the offending mass lies where the
  proposal rarely lands. The Pareto k-hat diagnostic
  (\`fit\$env\$impPsisK\`) does detect it: \`k \> 0.7\` means that
  subject's weights are unreliable. A t proposal has polynomial tails
  that dominate a Gaussian target's, which bounds the weights.

  NONMEM's guidance (Bauer, \*NONMEM Tutorial Part II\*) is to set a
  nonzero \`DF\` when there are fewer data points than etas, or for
  categorical data. Small values (3-8) are heavy; large values approach
  the Gaussian.

- proposal:

  Importance-sampling proposal family. \`"auto"\` (default) resolves to
  \`"t"\` when \`df \> 0\` and \`"normal"\` otherwise, so the historical
  \`df\` behaviour is unchanged.

  \* \`"normal"\` / \`"t"\` – the multivariate normal and t proposals
  \`df\` already selected. See \`df\` for the shape-versus-width
  discussion these all turn on. \* \`"laplace"\` – a spherical
  (Kotz-type) multivariate Laplace, \`f(x) ~ exp(-sqrt(x' S^-1 x))\`.
  Its exponential tail dominates the joint target's, which for a bounded
  likelihood is at worst Gaussian, so the importance weights are BOUNDED
  by construction – without having to pick a \`df\`, and without the
  low-\`df\` t's cost in effective sample size. Its scale is
  covariance-matched (\`S = Sigma/(p+1)\`), so \`gamma\`, \`iscaleMin\`
  and \`iscaleMax\` mean what they mean for \`"normal"\`. \*
  \`"mixture"\` – a defensive scale mixture about the same mode, \`sum_k
  w_k N(mode, c_k gamma Sigma)\`, with \`propMixScale\` and
  \`propMixWeight\`. A broad component covers what a narrow one misses,
  so the mixture is deliberately over-dispersed relative to \`gamma \*
  Sigma\`. Note this does NOT bound the weights: its widest component is
  still Gaussian-tailed, so use \`"laplace"\` when that is what you
  need.

  \`auto\` adapts \`df\` only on the normal/t axis – its ladder encodes
  "Gaussian" as \`df \<= 0\` and its constants were tuned against a
  Gaussian baseline, so a subject on \`"laplace"\` or \`"mixture"\`
  keeps its family. The \`isample\` budget reallocation still applies to
  every family.

  \`fit\$env\$impProposal\` reports the resolved family and
  \`fit\$env\$impPropInd\` the per-subject families actually used, so an
  \`auto\` escalation from normal to t is visible.

  \*\*Measured trade-off.\*\* 8 seeds, \`isample = 300\`, RMSE against
  an \`isample = 8000\` reference
  (\`design/qrpem/qrpem-options-bench.R\`), reported as theta RMSE /
  \`Omega\` RMSE / max Pareto k-hat / mean effective-sample fraction:

  - 1 ETA – \`"normal"\` 0.00157 / 0.00344 / -1.18 / 0.624;
    \`"laplace"\` 0.00128 / 0.00367 / -0.86 / 0.753; \`"mixture"\`
    0.00173 / 0.00510 / -1.33 / 0.592

  - 3 ETAs – \`"normal"\` 0.00144 / 0.00184 / -0.22 / 0.724;
    \`"laplace"\` 0.00110 / 0.00229 / +0.28 / 0.817; \`"mixture"\`
    0.00140 / 0.00296 / -0.33 / 0.668

  - general \`ll()\` – \`"normal"\` 0.00173 / 0.00170 / -0.66 / 0.658;
    \`"laplace"\` 0.00214 / 0.00227 / -0.31 / 0.798; \`"mixture"\`
    0.00144 / 0.00208 / -0.67 / 0.611

  - 8 ETAs (warfarin PK/PD) – \`"normal"\` 0.01459 / 0.01941 / +0.97 /
    0.551; \`"laplace"\` 0.01161 / 0.01430 / +1.07 / 0.576;
    \`"mixture"\` 0.01905 / 0.02003 / +1.19 / 0.527

  \`"auto"\` stays the default – no family wins everywhere – but the
  8-ETA row is the informative one and it is worth reading carefully. It
  is the only fixture here whose sampler actually fails (2.1 subjects
  above k-hat 0.7 under the default), and there \`"laplace"\` gives the
  best theta RMSE of every setting measured, 20 effective sample size.
  On the three easy fixtures no subject exceeds 0.7 under any setting,
  so there is nothing for a bounded-weight proposal to repair and it
  merely costs accuracy. That is the rule of thumb: reach for
  \`"laplace"\` when \`fit\$env\$impPsisK\` shows a tail the \`df\`
  ladder is not fixing, not otherwise.

  One honest caveat on the mechanism. \`"laplace"\` improved ACCURACY at
  8 ETAs without improving k-hat itself (+1.07 against +0.97). Bounded
  weights are an asymptotic guarantee, and 300 samples in 8 dimensions
  is not the asymptotic regime; the tail here is also not the
  non-identified kind the \`autoNonmemSparse\` note describes, since
  every subject has 13 observations for its 8 random effects. So take
  the accuracy gain as measured and the tail argument as theory that
  this fixture does not confirm.

- propMixScale:

  Component variance multipliers for \`proposal="mixture"\`, length 2 or
  3, starting at \`1\` and strictly increasing. Used as given: component
  1 is the Laplace-approximation covariance itself and the rest are the
  wider defensive components, so the mixture's covariance is \`(sum_k
  w_k c_k) \* gamma \* Sigma\` – deliberately wider than \`gamma \*
  Sigma\`, which is the whole mechanism. \`gamma\` still scales the
  entire mixture, so \`iscaleMin\`/\`iscaleMax\` still bound it.

- propMixWeight:

  Component weights for \`proposal="mixture"\`, same length as
  \`propMixScale\`, positive, normalized internally.

- auto:

  NONMEM \`AUTO=1\` equivalent: adapt the proposal degrees of freedom,
  the sample count and the acceptance target \*\*per subject\*\* rather
  than applying one global setting to everybody.

  \* \*\*\`df\`\*\* – any subject gets a heavy-tailed t proposal when
  the model is not transformably normal, and any subject whose Pareto
  k-hat reports tail failure gets one on that evidence. The tutorial's
  other trigger, "fewer data points than there are ETAs", is \*\*not\*\*
  applied on its own – see \`autoNonmemSparse\`. An escalation that
  fails to improve k-hat over two iterations is withdrawn, since a heavy
  tail the data creates is not repairable by proposal shape. \*
  \*\*\`isample\`\*\* – the total sample budget (\`isample \* nsub\`) is
  reallocated toward subjects whose effective-sample fraction is lowest,
  the tutorial's "many ETAs or ... large stochastic fluctuations". It is
  load-balancing, not a cost increase. Note sample count is deliberately
  \*not\* driven by Pareto k-hat: a heavy tail is not repairable by more
  draws (see \`df\`). \* \*\*\`iaccept\`\*\* – lowered to 0.2 for the
  same sparse/categorical subjects, per the tutorial.

  The concrete values (the \`df\` ladder, the k-hat thresholds, the
  budget reallocation rule) are \*\*nlmixr2's choices\*\*. NONMEM does
  not publish what \`AUTO=1\` picks internally; only the triggers and
  \`IACCEPT ~ 0.2\` are documented. Unlike NONMEM's AUTO, which its own
  tutorial warns "may result in lack of stochastic reproducibility",
  this remains seeded and thread-count independent.

  Per-subject values are reported in \`fit\$env\$impDfInd\`,
  \`fit\$env\$impNsampleInd\` and \`fit\$env\$impIacceptInd\`.

  \*\*Measured trade-off.\*\* \`auto = TRUE\` (the default) improves the
  \*tail\* behaviour of the importance weights and the accuracy of
  \`Omega\`, at some cost in Monte-Carlo noise on the objective. On
  theophylline with three ETAs – a one-ETA model has no tail failure to
  fix, so it shows nothing either way – against a reference computed at
  \`isample = 8000\`, over 8 seeds:

  - \`auto = FALSE\`: max Pareto k-hat 0.941, 2.38 subjects above 0.7;
    objective RMSE 0.0113, \`Omega\` RMSE 0.00406

  - \`auto = TRUE\`: max Pareto k-hat 0.571, 0.25 subjects above 0.7;
    objective RMSE 0.0165, \`Omega\` RMSE 0.00301

  So \`auto\` removes the tail failure and estimates \`Omega\` about 26
  accurately, for about 46 on by default because weights with infinite
  variance are a correctness problem – their error is unbounded in the
  worst case – whereas the extra noise is a bounded, measurable cost.
  Setting \`df\` globally gives a better objective still on a model that
  is uniformly heavy-tailed, but costs 75 more objective RMSE on one
  that is not, which is why it is not the default.

  Set \`auto = FALSE\` to recover the previous behaviour exactly (that
  path remains bit-identical to earlier versions), which is worth doing
  when you want the tightest possible objective on a model whose
  \`fit\$env\$impPsisK\` values are already comfortably below 0.7.

- autoNonmemSparse:

  Apply the NONMEM tutorial's \`nobs \< neta\` trigger unconditionally,
  so a subject with fewer observations than random effects always gets a
  t proposal. \`FALSE\` (default) leaves such subjects to the Pareto
  k-hat evidence like any other, and withdraws an escalation that does
  not improve k-hat within two iterations.

  The default diverges from the tutorial deliberately, on measurement.
  On a fixture built to the tutorial's own definition of sparse (2
  observations, 3 etas), 8 seeds against an \`isample = 8000\`
  reference, applying the rule made every number worse – objective RMSE
  0.1517 -\> 0.2059, \`Omega\` RMSE 0.02379 -\> 0.03163, and \*more\*
  failing subjects (3.88 -\> 4.88). With fewer observations than random
  effects the individual posterior is not identified, so the heavy tail
  is structural and no proposal shape repairs it; the reference itself
  still reads max k-hat 0.794.

  Set \`TRUE\` to get the documented NONMEM rule anyway. NONMEM's own
  testing is not published and was not done on these models, so someone
  who measures the opposite on their own problem should be able to have
  it.

- autoDfPatience:

  Number of consecutive iterations an escalated \`df\` may fail to
  improve Pareto k-hat before that escalation is withdrawn and the
  subject returned to the proposal it had without it. \`0\` never
  withdraws. Withdrawal is final, so escalation and withdrawal cannot
  oscillate, and it never goes below the \`df\` the model itself
  requires – a non-normal endpoint keeps its t proposal.

  A rung is judged against what the subject manages WITHOUT any
  escalation, measured while it sits there, so deterioration that
  happens before any escalation is tracked and an escalation is never
  credited for it.

- iscaleMin, iscaleMax:

  Lower/upper bounds for the adapted \`gamma\` (NONMEM ISCALE_MIN /
  ISCALE_MAX). Both bounds are reachable under \`gammaRule="target"\`;
  under \`"floor"\` \`gamma\` only ever moves up, so only \`iscaleMax\`
  can bind.

- iaccept:

  Minimum importance-sampling effective-sample fraction (NONMEM
  IACCEPT). The proposal scale \`gamma\` is kept at its efficient
  starting value while the achieved fraction stays at or above
  \`iaccept\`, and is inflated (toward \`iscaleMax\`) only when it drops
  below this floor.

- ctol:

  Convergence tolerance on the windowed objective-function change;
  \`NULL\` derives it from \`sigdig\`.

- nConvWindow:

  Length of the trailing iteration window used to average the
  objective-function change for convergence (NONMEM-style CTYPE).

- impSeed:

  Base seed for the per-subject thread-safe (threefry) RNG streams;
  results are reproducible and independent of the thread count.

- covMethod:

  Covariance method. \`"imp"\` (default) computes the Monte-Carlo
  importance-sampling observed-information covariance for the estimated
  thetas and Omega parameters (a finite-difference Hessian of the
  importance-sampling objective over fixed common-random-number
  samples), stashed as \`\$impCov\` / \`\$impSe\` and installed as the
  fit covariance; the theta standard errors match the Hessian-based
  FOCEI covariance, though the variance of a tightly-determined random
  effect (an Omega diagonal) can be over-estimated because the fixed
  samples barely span its prior variation. \`"analytic"\`, \`"r,s"\`,
  \`"r"\`, \`"s"\` instead compute the FOCEI covariance post-fit at the
  converged estimates (see \[foceiControl()\]); \`""\` skips the
  covariance step.

- qr:

  When \`TRUE\`, draw quasi-random (Sobol low-discrepancy) importance
  samples instead of pseudo-random Gaussian samples (QRPEM, Leary &
  Dunlavey PAGE 2012); the E-step integrals converge at O(1/N) instead
  of O(1/sqrt(N)).

- qrShift:

  Only used with \`qr=TRUE\`. When \`TRUE\` each (iteration, subject)
  applies a random Cranley-Patterson shift to the Sobol points (seeded,
  thread-count independent); \`FALSE\` reuses one fixed Sobol point set
  everywhere (fully deterministic E-step, no RNG in the draw).

- qrRefresh:

  Only used with \`qr=TRUE\` and \`qrShift=TRUE\`. When \`TRUE\` the
  shift is redrawn each iteration so residual quasi-random error
  averages out over the EM; \`FALSE\` draws one shift per subject at the
  fit start, making each EM iteration a deterministic map (smoothest
  objective trace).

- qrScramble:

  Only used with \`qr=TRUE\`. Scrambling of the Sobol point set:
  \`"none"\` (default) uses the raw sequence, randomized only by the
  Cranley-Patterson shift; \`"owen"\` applies a hash-based nested
  uniform (Owen) scramble; \`"lms"\` applies a linear matrix scramble
  with a digital shift.

  A Cranley-Patterson shift randomizes the point set but leaves the
  correlation structure between the sequence's dimensions intact;
  scrambling permutes the digits and breaks it.

  \*\*Measured trade-off, and it is not the one the paragraph above
  suggests.\*\* The scramble's advantage over the shift is real but does
  NOT grow with the number of random effects: on a smooth test integrand
  at \`isample = 300\` it is 1.79x at one dimension and 1.53x at twelve.
  And it lands on the E-STEP INTEGRAL rather than on the estimates –
  single-iteration importance- sampling \`-2LL\` RMSE 0.056 under
  \`"owen"\` against 0.080 unscrambled (3 ETAs), while converged
  parameter accuracy is unmoved or worse. Against an \`isample = 8000\`
  reference over 8 seeds, theta RMSE:

  - 1 ETA – shift 0.00022, \`"owen"\` 0.00018, \`"lms"\` 0.00022

  - 3 ETAs – shift 0.00031, \`"owen"\` 0.00045, \`"lms"\` 0.00049

  - general \`ll()\` – shift 0.00020, \`"owen"\` 0.00032, \`"lms"\`
    0.00028

  - 8 ETAs – shift 0.01493, \`"owen"\` 0.01787, \`"lms"\` 0.01924

  So \`"none"\` stays the default: at the eta dimension where a scramble
  was expected to pay off most it is the worst of the three on theta,
  and it never improves Pareto k-hat. Set \`"owen"\` when the quantity
  you care about is the reported importance-sampling objective
  (\`fit\$env\$impObj\`, e.g. for comparing models by IS likelihood)
  rather than the parameter estimates – that is the integral it
  measurably improves.

  Scrambling IS the randomization, so it REPLACES the shift rather than
  composing with it: \`qrShift\` is ignored when this is not \`"none"\`,
  while \`qrRefresh\` still decides whether the randomization is redrawn
  each iteration (\`FALSE\` pins one scramble per subject, making each
  EM iteration a deterministic map). The scramble key is derived
  arithmetically from \`impSeed\` and the (iteration, subject,
  dimension) indices rather than drawn from the RNG, so it consumes no
  draws and the fit stays reproducible and independent of the thread
  count.

  \`"lms"\` is named for the method (Matousek/Tezuka linear matrix
  scrambling) and is not a claim to reproduce any particular vendor's
  variant.

- sir:

  When \`TRUE\`, accelerate the non-mu / residual-error M-step by SIR
  (sampling-importance-resampling): the theta-sensitivity Newton step
  uses \`sirSample\` equal-weight resampled points per subject instead
  of all \`isample\` weighted samples.

- sirSample:

  Number of SIR resampled points per subject; \`NULL\` uses \`max(25,
  ceiling(isample/10))\`. Must be at most \`isample\`.

- muModel:

  Mu-referencing variant for the MAP inner problem; for
  \`impmapControl()\` this is always \`"lin"\` and cannot be changed.

- combSens:

  When \`TRUE\` (the default) and there are non-mu (structural or
  residual-error) thetas to estimate, carry their sensitivity columns on
  the INNER model itself (the combined eta+theta sensitivity build,
  \#958) instead of compiling and solving a separate, dedicated model
  for them. With \`sir=FALSE\` (the default) the E-step's own per-sample
  inner solve then supplies the M-step's Newton step directly, with no
  second solve – roughly halving the ODE solving the theta-sensitivity
  Newton step costs. Set \`FALSE\` to use the older two-model path
  instead. (Earlier versions shipped \`FALSE\` as the default: on a
  \`linCmt()\` model with a single non-mu structural theta,
  \`combSens=TRUE\`/\`FALSE\` measurably disagreed – not \`combSens\`
  moving the answer, but an unrelated pre-existing bug in how the shared
  solve pool swapped between peer models, since fixed. The two now agree
  on every fixture measured, so \`TRUE\` is the default.)

## Value

impmapControl object

## Author

Matthew L. Fidler

## Examples

``` r

impmapControl()
#> $maxOuterIterations
#> [1] 5000
#> 
#> $maxInnerIterations
#> [1] 1000
#> 
#> $n1qn1nsim
#> [1] 10001
#> 
#> $iterPrintControl
#> $every
#> [1] 1
#> 
#> $ncol
#> [1] 4
#> 
#> $headerEvery
#> [1] 10
#> 
#> $useColor
#> [1] TRUE
#> 
#> $simple
#> [1] FALSE
#> 
#> attr(,"class")
#> [1] "iterPrintControl" "list"            
#> 
#> $lbfgsLmm
#> [1] 7
#> 
#> $lbfgsPgtol
#> [1] 0
#> 
#> $lbfgsFactr
#> [1] 4.5036e+10
#> 
#> $scaleTo
#> [1] 1
#> 
#> $epsilon
#> [1] 0.001
#> 
#> $derivEps
#> [1] 2.980232e-07 2.980232e-07
#> 
#> $derivMethod
#> [1] 3
#> 
#> $covDerivMethod
#> [1] 1
#> 
#> $covMethod
#> [1] 2
#> 
#> $covType
#> [1] "analytic"
#> 
#> $covMethodDeferred
#> [1] NA
#> 
#> $covSolveTol
#> NULL
#> 
#> $covFull
#> [1] TRUE
#> 
#> $fast
#> [1] FALSE
#> 
#> $priorMethod
#> [1] "auto"
#> 
#> $fdOutlierZ
#> [1] 3.5
#> 
#> $fdOutlierScale
#> [1] 1
#> 
#> $fdRefine
#> [1] "chartrand"
#> 
#> $fdLanczosM
#> [1] 2
#> 
#> $fdRichardsonR
#> [1] 2
#> 
#> $fdRichardsonV
#> [1] 2
#> 
#> $fdChartrandAll
#> [1] 0
#> 
#> $fdOutlierAny
#> [1] 0
#> 
#> $fdIndividualStep
#> [1] 1
#> 
#> $fdChartrand
#> [1] 1
#> 
#> $centralDerivEps
#> [1] 2.980232e-07 2.980232e-07
#> 
#> $eigen
#> [1] 1
#> 
#> $diagXform
#> [1] "sqrt"
#> 
#> $iovXform
#> [1] "sd"
#> 
#> $iovMethod
#> [1] "auto"
#> 
#> $sumProd
#> [1] FALSE
#> 
#> $optExpression
#> [1] TRUE
#> 
#> $literalFix
#> [1] TRUE
#> 
#> $literalFixRes
#> [1] TRUE
#> 
#> $outerOpt
#> [1] -1
#> 
#> $ci
#> [1] 0.95
#> 
#> $sigdig
#> [1] 3
#> 
#> $sigdigTable
#> [1] 3
#> 
#> $scaleObjective
#> [1] 0
#> 
#> $boundTol
#> [1] 0.05
#> 
#> $calcTables
#> [1] TRUE
#> 
#> $noAbort
#> [1] 1
#> 
#> $interaction
#> [1] 1
#> 
#> $foce
#> [1] "nonmem"
#> 
#> $foceType
#> [1] 0
#> 
#> $cholSEtol
#> [1] 6.055454e-06
#> 
#> $foceEbeTol
#> [1] 1e-09
#> 
#> $hessEps
#> [1] 6.055454e-06
#> 
#> $hessEpsLlik
#> [1] 6.055454e-06
#> 
#> $optimHessType
#> [1] 1
#> 
#> $optimHessCovType
#> [1] 1
#> 
#> $hessEtaStepMin
#> [1] 0.05
#> 
#> $censOption
#> [1] 0
#> 
#> $cholAccept
#> [1] 0.001
#> 
#> $resetEtaSize
#> [1] 1.439531
#> 
#> $resetThetaSize
#> [1] Inf
#> 
#> $resetThetaFinalSize
#> [1] Inf
#> 
#> $diagOmegaBoundUpper
#> [1] 5
#> 
#> $diagOmegaBoundLower
#> [1] 100
#> 
#> $cholSEOpt
#> [1] 0
#> 
#> $cholSECov
#> [1] 0
#> 
#> $fo
#> [1] 0
#> 
#> $covTryHarder
#> [1] 0
#> 
#> $outerOptFun
#> function (par, fn, gr, lower = -Inf, upper = Inf, control = list(), 
#>     ...) 
#> {
#>     .ctl <- .controlMaxfun(control)
#>     if (is.null(.ctl$npt)) 
#>         .ctl$npt <- length(par) * 2 + 1
#>     .ctl$iprint <- 0L
#>     .ctl <- .ctl[names(.ctl) %in% c("npt", "rhobeg", "rhoend", 
#>         "iprint", "maxfun")]
#>     .ret <- minqa::bobyqa(par, fn, control = .ctl, lower = lower, 
#>         upper = upper)
#>     .ret <- .bobyqaRetryIfStuck(par, fn, lower, upper, .ctl, 
#>         .ret)
#>     .ret$x <- .ret$par
#>     .ret$message <- .ret$msg
#>     .ret$convergence <- .ret$ierr
#>     .ret$value <- .ret$fval
#>     .ret
#> }
#> <bytecode: 0x55f2f1367cb0>
#> <environment: namespace:nlmixr2est>
#> 
#> $rhobeg
#> [1] 0.2
#> 
#> $rhoend
#> [1] 0.001
#> 
#> $npt
#> NULL
#> 
#> $rel.tol
#> [1] 0.001
#> 
#> $x.tol
#> [1] 0.001
#> 
#> $eval.max
#> [1] 4000
#> 
#> $iter.max
#> [1] 2000
#> 
#> $innerOpt
#> [1] 4
#> 
#> $hessianMethod
#> [1] 1
#> 
#> $trustConf
#> [1] 0.975
#> 
#> $trustRinit
#> NULL
#> 
#> $trustRmax
#> NULL
#> 
#> $trustFterm
#> [1] 1e-05
#> 
#> $trustMterm
#> [1] 1e-05
#> 
#> $abstol
#> [1] 0.001
#> 
#> $reltol
#> [1] 0.001
#> 
#> $derivSwitchTol
#> [1] 0.002
#> 
#> $resetHessianAndEta
#> [1] 0
#> 
#> $muModel
#> [1] "lin"
#> 
#> $muRefCovAlg
#> [1] TRUE
#> 
#> $muModelTol
#> [1] 1e-05
#> 
#> $muModelMaxCycles
#> [1] 20
#> 
#> $muModelClampRetries
#> [1] 10
#> 
#> $stateTrim
#> [1] Inf
#> 
#> $gillK
#> [1] 10
#> 
#> $gillKcov
#> [1] 10
#> 
#> $gillKcovLlik
#> [1] 10
#> 
#> $gillRtol
#> [1] 1.490116e-08
#> 
#> $gillStep
#> [1] 4
#> 
#> $gillStepCov
#> [1] 2
#> 
#> $gillStepCovLlik
#> [1] 4.5
#> 
#> $scaleType
#> [1] 2
#> 
#> $normType
#> [1] 1
#> 
#> $scaleC
#> NULL
#> 
#> $scaleCmin
#> [1] 1e-05
#> 
#> $scaleCband
#> [1]  0.1 10.0
#> 
#> $scaleCmax
#> [1] 1e+05
#> 
#> $scaleC0
#> [1] 1e+05
#> 
#> $outerOptTxt
#> [1] "bobyqa"
#> 
#> $outerOptDefault
#> [1] TRUE
#> 
#> $rmatNorm
#> [1] 1
#> 
#> $rmatNormLlik
#> [1] 1
#> 
#> $smatNorm
#> [1] 1
#> 
#> $smatNormLlik
#> [1] 1
#> 
#> $covGillF
#> [1] 1
#> 
#> $optGillF
#> [1] 1
#> 
#> $gillFtol
#> [1] 0
#> 
#> $gillFtolCov
#> [1] 0
#> 
#> $gillFtolCovLlik
#> [1] 0
#> 
#> $covSmall
#> [1] 1e-05
#> 
#> $adjLik
#> [1] TRUE
#> 
#> $gradTrim
#> [1] Inf
#> 
#> $gradCalcCentralSmall
#> [1] 1e-04
#> 
#> $gradCalcCentralLarge
#> [1] 10000
#> 
#> $etaNudge
#> [1] 1.131586
#> 
#> $etaNudge2
#> [1] 1.518182
#> 
#> $maxOdeRecalc
#> [1] 5
#> 
#> $odeRecalcFactor
#> [1] 3.162278
#> 
#> $nRetries
#> [1] 3
#> 
#> $seed
#> [1] 42
#> 
#> $resetThetaCheckPer
#> [1] 0.1
#> 
#> $etaMat
#> NULL
#> 
#> $repeatGillMax
#> [1] 1
#> 
#> $stickyRecalcN
#> [1] 4
#> 
#> $outerMaxOdeRecalc
#> [1] 5
#> 
#> $outerOdeRecalcFactor
#> [1] 3.162278
#> 
#> $outerStickyRecalcN
#> [1] 4
#> 
#> $indTolRelax
#> [1] TRUE
#> 
#> $eventType
#> [1] 2
#> 
#> $eventSens
#> [1] "jump"
#> 
#> $gradProgressOfvTime
#> [1] 10
#> 
#> $addProp
#> [1] "combined2"
#> 
#> $badSolveObjfAdj
#> [1] 100
#> 
#> $compress
#> [1] FALSE
#> 
#> $rxControl
#> $scale
#> NULL
#> 
#> $method
#> liblsoda 
#>        2 
#> 
#> $atol
#> [1] 1e-06
#> 
#> $rtol
#> [1] 0.001
#> 
#> $maxsteps
#> [1] 500000
#> 
#> $hmin
#> [1] 0
#> 
#> $hmax
#> [1] NA
#> 
#> $hini
#> [1] 0
#> 
#> $maxordn
#> [1] 12
#> 
#> $maxords
#> [1] 5
#> 
#> $covsInterpolation
#> locf 
#>    1 
#> 
#> $addCov
#> [1] TRUE
#> 
#> $returnType
#> rxSolve 
#>       0 
#> 
#> $sigma
#> NULL
#> 
#> $sigmaDf
#> NULL
#> 
#> $nCoresRV
#> [1] 1
#> 
#> $sigmaIsChol
#> [1] FALSE
#> 
#> $sigmaSeparation
#> [1] "auto"
#> 
#> $sigmaXform
#> identity 
#>        4 
#> 
#> $nDisplayProgress
#> [1] 10000
#> 
#> $amountUnits
#> [1] NA
#> 
#> $timeUnits
#> [1] "hours"
#> 
#> $addDosing
#> [1] FALSE
#> 
#> $stateTrim
#> [1] Inf
#> 
#> $updateObject
#> [1] FALSE
#> 
#> $omega
#> NULL
#> 
#> $omegaDf
#> NULL
#> 
#> $omegaIsChol
#> [1] FALSE
#> 
#> $omegaSeparation
#> [1] "auto"
#> 
#> $omegaXform
#> variance 
#>        6 
#> 
#> $nSub
#> [1] 1
#> 
#> $thetaMat
#> NULL
#> 
#> $thetaDf
#> NULL
#> 
#> $thetaIsChol
#> [1] FALSE
#> 
#> $nStud
#> [1] 1
#> 
#> $dfSub
#> [1] 0
#> 
#> $dfObs
#> [1] 0
#> 
#> $seed
#> NULL
#> 
#> $nsim
#> NULL
#> 
#> $minSS
#> [1] 10
#> 
#> $maxSS
#> [1] 10000
#> 
#> $strictSS
#> [1] 1
#> 
#> $infSSstep
#> [1] 12
#> 
#> $istateReset
#> [1] TRUE
#> 
#> $subsetNonmem
#> [1] TRUE
#> 
#> $hmaxSd
#> [1] 0
#> 
#> $maxAtolRtolFactor
#> [1] 0.1
#> 
#> $from
#> NULL
#> 
#> $to
#> NULL
#> 
#> $by
#> NULL
#> 
#> $length.out
#> NULL
#> 
#> $iCov
#> NULL
#> 
#> $keep
#> NULL
#> 
#> $keepF
#> character(0)
#> 
#> $drop
#> NULL
#> 
#> $warnDrop
#> [1] TRUE
#> 
#> $omegaLower
#> [1] -Inf
#> 
#> $omegaUpper
#> [1] Inf
#> 
#> $sigmaLower
#> [1] -Inf
#> 
#> $sigmaUpper
#> [1] Inf
#> 
#> $thetaLower
#> [1] -Inf
#> 
#> $thetaUpper
#> [1] Inf
#> 
#> $indLinPhiM
#> [1] 0
#> 
#> $indLinPhiTol
#> [1] 1e-07
#> 
#> $indLinMatExpType
#> Al-Mohy 
#>       3 
#> 
#> $indLinMatExpOrder
#> [1] 6
#> 
#> $idFactor
#> [1] TRUE
#> 
#> $mxhnil
#> [1] 0
#> 
#> $hmxi
#> [1] 0
#> 
#> $warnIdSort
#> [1] TRUE
#> 
#> $ssAtol
#> [1] 1e-05
#> 
#> $ssRtol
#> [1] 0.01
#> 
#> $safeZero
#> [1] 1
#> 
#> $sumType
#> pairwise 
#>        1 
#> 
#> $prodType
#> long double 
#>           1 
#> 
#> $resample
#> NULL
#> 
#> $resampleID
#> [1] TRUE
#> 
#> $maxwhile
#> [1] 100000
#> 
#> $cores
#> [1] 0
#> 
#> $atolSens
#> [1] 1e-06
#> 
#> $rtolSens
#> [1] 0.001
#> 
#> $ssAtolSens
#> [1] 1e-05
#> 
#> $ssRtolSens
#> [1] 0.01
#> 
#> $simVariability
#> [1] NA
#> 
#> $nLlikAlloc
#> NULL
#> 
#> $useStdPow
#> [1] 0
#> 
#> $naTimeHandle
#> ignore 
#>      1 
#> 
#> $addlKeepsCov
#> [1] FALSE
#> 
#> $addlDropSs
#> [1] TRUE
#> 
#> $ssAtDoseTime
#> [1] TRUE
#> 
#> $ss2cancelAllPending
#> [1] FALSE
#> 
#> $naInterpolation
#> locf 
#>    1 
#> 
#> $keepInterpolation
#> na 
#>  2 
#> 
#> $safeLog
#> [1] 1
#> 
#> $safePow
#> [1] 1
#> 
#> $ssSolved
#> [1] TRUE
#> 
#> $linCmtSensType
#> auto 
#>  100 
#> 
#> $linCmtSensH
#> [1] 1e-04
#> 
#> $linCmtGillFtol
#> [1] 0
#> 
#> $linCmtGillK
#> [1] 20
#> 
#> $linCmtGillStep
#> [1] 4
#> 
#> $linCmtGillRtol
#> [1] 1.490116e-08
#> 
#> $linCmtShiErr
#> [1] 1.490116e-08
#> 
#> $linCmtShiMax
#> [1] 20
#> 
#> $linCmtScale
#> [1] 0 0 0 0 0 0 0
#> 
#> $linCmtHcmt
#> [1] 1
#> 
#> $linCmtHmeanI
#> geometric 
#>         2 
#> 
#> $linCmtHmeanO
#> geometric 
#>         2 
#> 
#> $linCmtSuspect
#> [1] 1e-06
#> 
#> $linCmtForwardMax
#> [1] 2
#> 
#> $indOwnAlloc
#> [1] -1
#> 
#> $maxExtra
#> [1] 1000
#> 
#> $tolFactor
#> NULL
#> 
#> $serializeFile
#> NULL
#> 
#> $dense
#> [1] FALSE
#> 
#> $cvodeLinSolver
#> dense 
#>     1 
#> 
#> $stiff2
#> [1] 0
#> 
#> $autoSwitchMaxStiff
#> [1] 10
#> 
#> $autoSwitchMaxNonstiff
#> [1] 3
#> 
#> $autoSwitchStiffFirst
#> [1] 0
#> 
#> $autoSwitchNonstifftol
#> [1] 0.9
#> 
#> $autoSwitchStifftol
#> [1] 0.9
#> 
#> $autoSwitchDtfac
#> [1] 2
#> 
#> $autoSwitchSwitchMax
#> [1] 5
#> 
#> $useLinCmt
#> [1] TRUE
#> 
#> $file
#> NULL
#> 
#> $chunkSize
#> NULL
#> 
#> $parallel
#> [1] 0
#> 
#> $.zeros
#> NULL
#> 
#> $zeroVarParamHandle
#> [1] "warn"
#> 
#> $indLinStepSearch
#> [1] 1
#> 
#> $indLinMaxIter
#> [1] 20
#> 
#> $indLinRichardson
#> [1] 2
#> 
#> $indLinIteration
#> [1] 3
#> 
#> $indLinJac
#> [1] 0
#> 
#> $indLinForcing
#> [1] 1
#> 
#> $usePrior
#> [1] NA
#> 
#> $priorPdRetry
#> [1] 10
#> 
#> $priorOmega
#> NULL
#> 
#> $priorOmegaEl
#> NULL
#> 
#> $priorSigmaEl
#> NULL
#> 
#> $linCmtSensPhi
#> [1] 2
#> 
#> attr(,"class")
#> [1] "rxControl"
#> 
#> $genRxControl
#> [1] TRUE
#> 
#> $skipCov
#> NULL
#> 
#> $fallbackFD
#> [1] FALSE
#> 
#> $shi21maxOuter
#> [1] 0
#> 
#> $shi21maxInner
#> [1] 20
#> 
#> $shi21maxInnerCov
#> [1] 20
#> 
#> $shi21maxFD
#> [1] 20
#> 
#> $shi21hMax
#> [1] 2
#> 
#> $shi21hMin
#> [1] 1e-04
#> 
#> $smatPer
#> [1] 0.6
#> 
#> $sdLowerFact
#> [1] 0.001
#> 
#> $zeroGradFirstReset
#> [1] TRUE
#> 
#> $zeroGradRunReset
#> [1] TRUE
#> 
#> $zeroGradBobyqa
#> [1] TRUE
#> 
#> $mceta
#> [1] -2
#> 
#> $warm
#> [1] 1
#> 
#> $nAGQ
#> [1] 0
#> 
#> $agqHi
#> [1] Inf
#> 
#> $agqLow
#> [1] -Inf
#> 
#> $sensMethod
#> [1] "default"
#> 
#> $linCmtSensCarry
#> [1] "auto"
#> 
#> $boundedTransform
#> [1] TRUE
#> 
#> $zeroTheta
#> [1] 0.001
#> 
#> $impCov
#> [1] TRUE
#> 
#> $isample
#> [1] 300
#> 
#> $nIter
#> [1] 100
#> 
#> $mapIter
#> [1] 1
#> 
#> $nBurn
#> [1] 0
#> 
#> $burnFreezeOmega
#> [1] FALSE
#> 
#> $gamma
#> [1] 1
#> 
#> $gammaMethod
#> [1] "auto"
#> 
#> $gammaRule
#> [1] "target"
#> 
#> $df
#> [1] 0
#> 
#> $auto
#> [1] TRUE
#> 
#> $autoNonmemSparse
#> [1] FALSE
#> 
#> $autoDfPatience
#> [1] 2
#> 
#> $iscaleMin
#> [1] 0.1
#> 
#> $iscaleMax
#> [1] 10
#> 
#> $iaccept
#> [1] 0.4
#> 
#> $nConvWindow
#> [1] 20
#> 
#> $impSeed
#> [1] 42
#> 
#> $qr
#> [1] FALSE
#> 
#> $qrShift
#> [1] TRUE
#> 
#> $qrRefresh
#> [1] TRUE
#> 
#> $qrScramble
#> [1] "none"
#> 
#> $proposal
#> [1] "auto"
#> 
#> $propMixScale
#> [1] 1 9
#> 
#> $propMixWeight
#> [1] 0.9 0.1
#> 
#> $sir
#> [1] FALSE
#> 
#> $sirSample
#> [1] 30
#> 
#> $combSens
#> [1] TRUE
#> 
#> attr(,"class")
#> [1] "impmapControl"
```
