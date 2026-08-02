# Re-instrument AUTO (imp/impmap) after the pooling fixes

## Why

`impmapControl(auto=TRUE)` escalates the proposal `df` per subject from the Pareto
k-hat diagnostic.  Both the escalation ladder and everything that tests or documents
it were calibrated on plain theophylline, which used to produce genuine tail failures.

The imp/impmap/focei pooling fixes removed those failures.  Measured now on
theophylline (`impmap`, `isample=300`, `nIter=12`):

| quantity                    | documented | measured now |
|-----------------------------|-----------|--------------|
| subjects with k-hat > 0.7   | 2 of 12   | **0**        |
| max k-hat, `auto=FALSE`     | 2.44      | **-1.46**    |
| max k-hat, `auto=TRUE`      | 0.49      | -1.18        |

Consequences, all of which are live:

1. `test-imp-auto.R:60` asserts its own premise (`sum(k0 > 0.7) > 0`) and that premise
   now fails.  The test is not detecting an AUTO regression -- there is nothing left to
   escalate.  It would go green again the moment any subject crossed 0.7 and silently
   tests nothing until then.
2. The ladder's severity mapping (`k>2 -> df 12`, `k>1 -> 20`, `k>0.7 -> 30`,
   `src/imp.cpp`) comes from a df sweep against an `isample=20000` reference **on
   theophylline**.  That sweep is no longer reproducible on that model, so the constants
   are currently unvalidated rather than known-wrong.
3. `~/src/nlmixr2/vignettes/articles/imp-impmap-qrpem.Rmd` publishes the stale numbers
   (k-hat 2.60/1.28 at line ~435, the 2.44 -> 0.49 claim at ~454, and the comparison
   tables around ~921-983).

The goal is the one stated: get the defaults right **where tail failure occurs**, while
not firing -- and not costing accuracy -- where it does not.

## Scope against the NONMEM sources

Both source documents are in `~/src/private`:

- `bauer-nonmem-tutorial-part2-2019.pdf` (the "tutorial" the code comments cite)
- `NONMEM7_Technical_Guide.pdf`

What they actually specify, checked rather than recalled:

- **Triggers for a t proposal** (tutorial, p. ~505-517): "when there are fewer data points
  than there are ETAs to be estimated (sparse data) or data are categorical, then DF
  should be set to a nonzero number to use the t-distribution sampler and/or decrease the
  acceptance rate IACCEPT from the default value of 0.4 to about 0.2".
- **ISAMPLE**: default 300 for IMP, raise "when there are many ETAs or the objective
  function has large stochastic fluctuations".
- **AUTO=1**: "requests NONMEM to make decisions about ISAMPLE for SAEM and IACCEPT, DF,
  and ISAMPLE for IMP for each subject".
- **Reproducibility**: "using this feature may result in lack of stochastic
  reproducibility" -- nlmixr2's AUTO is seeded and thread-count independent, a deliberate
  divergence worth keeping stated.

What they do **not** specify: any DF value, any threshold, any per-subject decision rule,
and nothing resembling a k-hat diagnostic.  The Technical Guide covers IACCEPT only as
the MCMC scaling target and t-sampling only as generation mechanics (Appendix J).

So the ladder constants are nlmixr2's own and cannot be validated against these documents
-- only the triggers and the IACCEPT 0.4/0.2 anchors can.  Two consequences:

1. Keep the documented triggers implemented verbatim and tested separately from the
   k-hat escalation, so a change in our numbers can never be mistaken for a change in
   NONMEM-documented behavior.
2. The ladder must be justified by measurement (Phase 3), because there is no external
   authority to defer to.

## Phase 0 -- is the improvement real?

A max k-hat of -1.46 across every subject is not merely "better", it is extreme.  Settle
this before tuning anything to it, because the remedy is opposite in the two cases.

- Compare in-kernel `impKhatIter` against `R/impPsis.R`'s `impPsisK` (the copy validated
  against `loo::psis`).  `test-imp-auto.R` already asserts these agree -- confirm it still
  passes, since agreement on a degenerate input proves nothing on its own.
- Inspect the raw importance weights per subject: ESS fraction (`Neff/isample`), the
  weight range, and how many draws carry the mass.  Genuinely light tails and a collapsed
  near-uniform weight vector both drive k-hat down, and only one is good news.
- Anchor on accuracy, not the diagnostic: objective and `Omega` RMSE against a
  high-`isample` (20000) reference, `auto` on vs off.  If the sampler really is better,
  RMSE improves or holds; if the weights have collapsed, RMSE degrades while k-hat looks
  excellent.

Exit: a one-line verdict, real or degenerate.  If degenerate, stop -- that is a sampler
bug and this plan is the wrong response to it.

## Phase 1 -- re-measure the k-hat landscape

Sweep the models that plausibly still stress PSIS, `auto=FALSE`, fixed seed, recording
per-subject k-hat, ESS fraction and objective:

- theophylline, `linCmt()` and the ODE form (the current fixture)
- sparse data: few observations per subject, `nobs < neta`
- multi-eta models (the tutorial's "many ETAs" trigger)
- non-normal: the Poisson fixture already in `test-imp-auto.R`
- heavy-tailed / outlier-contaminated data

Deliverable: a table of which model/data combinations produce genuine `k > 0.7` under
current code, with the seed that reproduces it.

## Phase 2 -- a STRESSED fixture, not theo_sd

theo_sd is the wrong instrument now and arguably always was: 11 observations, 1 eta,
transformably normal.  By the tutorial's own criteria none of its triggers fire there --
`test-imp-auto.R` says exactly that -- so it only ever exercised AUTO through the k-hat
path, which is the one the pooling fix removed.  It should not be the fixture that proves
escalation works.

Build the stressed fixture from the tutorial's own stress conditions, so the scenario is
scoped to the source rather than invented:

- **sparse: `nobs < neta`** -- the tutorial's first documented trigger.  Few observations
  per subject with several etas; this is where a Gaussian proposal genuinely fails.
- **categorical / non-normal** -- the second documented trigger.  The Poisson fixture
  already in the file is the seed of this; make it sparse as well.
- optionally heavy-tailed or outlier-contaminated data, as a k-hat stressor beyond what
  the tutorial documents (label it as ours, not NONMEM's).

Requirements: reliably produces `k > 0.7` for **at least one and not all** subjects under
a pinned seed (selectivity is the property under test); cheap enough for the weekly batch;
stable across thread counts -- AUTO is documented as seeded and thread-count independent
where NONMEM's explicitly is not, so verify that rather than assume it.

Note the interaction with sample budget: k-hat is deliberately not a sample-count lever
(`src/imp.cpp` records 10x samples moving k-hat 0.76 -> 3.28), so the fixture must fail on
proposal SHAPE, not on being under-sampled, or Phase 3 will tune the wrong knob.

## Phase 3 -- re-calibrate the ladder

Redo the df sweep that produced `12/20/30` and the `2.0/1.0/0.7` breakpoints, on the
Phase 2 fixture, against an `isample=20000` reference:

- per rung: max k-hat, subjects above 0.7, objective RMSE, `Omega` RMSE, runtime
- confirm or move each breakpoint and each df
- re-check the one-way property.  The comment argues escalation must not reverse because
  a working t proposal drives k-hat down.  That reasoning is sound and independent of the
  pooling fix, so expect it to survive -- but it is now untested, so test it.

Change constants only where the sweep says so, and record the measurement next to them.

## Phase 4 -- reinstrument the tests

Split the two things `test-imp-auto.R` currently conflates:

1. **Escalation fires where it should** -- on the Phase 2 fixture.  Assert the premise
   explicitly and make an absent premise a loud skip with reason, never a silent pass.
   Keep the existing selectivity assertions (some subjects escalate, not all) and that
   the tail is cleared.
2. **Escalation does NOT fire where it should not** -- on theophylline as it now behaves.
   This is the "still allow other failures not to occur" half: with every k-hat well under
   0.7, `auto=TRUE` must leave `impDfInd` at 0, must not degrade ESS, and must not move
   the objective beyond tolerance versus `auto=FALSE`.  That is a real regression guard,
   and it is exactly what the current vacuous premise destroyed.

Both belong in a weekly `.slowBatches` entry (fit-based); keep the cheap control
round-trip checks in the essential subset.

## Phase 5 -- update the article

`~/src/nlmixr2/vignettes/articles/imp-impmap-qrpem.Rmd`: refresh the k-hat numbers
(~435), the `2.44 -> 0.49` claim (~454), and the comparison tables (~921-983) from the
Phase 1/3 measurements.  Note in the text that plain theophylline no longer exhibits the
failure, and point the worked example at the Phase 2 fixture, so a reader running the
article reproduces what it claims.

## Phase 6 -- guard

Add a cheap assertion that the Phase 2 fixture still produces its tail failure under
`auto=FALSE`.  If a future change removes it again -- as the pooling fix did -- that
fires immediately, instead of quietly hollowing out the AUTO tests a second time.

## Phase 7 -- the withdrawal rule is not reachable by the test suite

Eight rounds of independent review found 11 defects in the AUTO withdrawal logic.
The last three (a stale bar, an easy bar, a pinned bar) were all LATENT: no
fixture in the sweep produces the k-hat sequences that trigger them, so the suite
was green throughout.  They were found by reasoning about the state machine, and
confirmed by replaying k-hat sequences through a standalone model of it.

`kAtFloor`, `noImp` and `escDead` are locals inside `impOuter`.  The tests can
only observe `impDfInd` (the last E-step's df) and the objective, so a subtle
regression -- reverting the margin comparison from `<=` to `<`, say -- passes
everything.

The fix is to make the decision testable rather than to add more fits:

* extract the per-subject decision into a pure function of
  `(kh, df, dfFloor, kAtFloor, noImp, escDead, patience, nonmemSparse)`
  returning the new state, and have the loop call it;
* export a thin wrapper for the tests (note `src/init.c` hand-maintains the
  `.Call` table with hardcoded arities -- see CLAUDE.md);
* drive it with the sequences that actually broke it:
  - floor 30, k-hat 0.6,0.6,0.6,1.3 -> escalate -> 0.9        must KEEP
  - floor 0, 0.85,0.85, spike 1.05 -> escalate -> 0.85,0.85   must WITHDRAW
  - floor 0, 0.85, 1.25 -> escalate -> 0.75,0.75              must KEEP

That is a refactor with `.Call` registration consequences, so it is deliberately
NOT part of this release; the rule as it stands is reviewed clean and measured.

## Phase 8 -- test-imp-auto.R is now too slow for its batch

The reinstrumentation took the file from about 6 fits to 14, several of them on
the sparse (nobs < neta) fixture, and it measures ~40 minutes and ~24 GB on one
worker.  `.slowBatches` sizes each batch to stay well under an hour IN TOTAL, so
one file at 40 minutes breaks that budget on its own.

This is a cost the reinstrumentation added, not a pre-existing one.  Options, in
rough order of preference:

* drop `nIter` for the sparse-based tests.  They assert WHICH subjects escalate
  and whether withdrawal fires, not convergence, so the default 12 iterations is
  more than the assertions need.  Verify each assertion still holds at the lower
  count rather than assuming.
* halve `isample` (300 -> 150) for the same tests; k-hat is noisier but the
  assertions are ordinal (some/not-all, fewer-than), not threshold-exact.
* collapse the `autoDfPatience` coverage: the 0-vs-2 comparison needs two sparse
  fits, and one of them could reuse a fit already computed in the same file.

Whatever is chosen, re-time the file afterwards -- the point is the budget, so an
untimed "should be faster" does not close this.

## Phase 9 -- test-imp-xi-gamma.R premises also went stale

Three more failures in the same family as the k-hat ones, and from the same
cause: the pooling fixes made the sampler better behaved, so tests that keyed on
the OLD pathological behaviour no longer hold.

    test-imp-xi-gamma.R:283  individual gamma > 1.2 * global    (now 0.98x)
    test-imp-xi-gamma.R:315  min(impGammaInd) > 1.05            (now 1.015)
    test-imp-xi-gamma.R:388  impXiTrace[1] < 2                  (now FALSE)

Measured, and NOT caused by the AUTO work in this branch: with `auto = FALSE`,
where no subject receives a t proposal at all (`impDfInd` all 0), the ratio is
0.975 -- essentially identical to the 0.981 with `auto = TRUE`.  The df ladder,
the sparsity gate and the withdrawal rule are all irrelevant here.

What actually changed is that `gammaMethod = "individual"` and `"global"` now
converge to nearly the same proposal scale (1.2241 vs 1.2484), where the test
expects them to differ by at least 20%.  That is arguably GOOD -- the two
adaptation laws agreeing is what a well-behaved sampler should do -- but it means
the test no longer demonstrates what it claims.

Treat as the k-hat cases were treated: establish what the current sampler
actually does, then decide whether the assertion should be re-aimed at a fixture
where the two laws genuinely diverge, or the claim itself retired.  Do not
simply widen the tolerance -- the premise, not the threshold, is what moved.
