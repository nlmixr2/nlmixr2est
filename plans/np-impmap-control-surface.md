# Make the npag/npb control surface honest

## Why

`npagControl()` and `npbControl()` are wrappers around `impmapControl()`: whatever
lands in `...` is passed straight through. That buys the shared FOCEI-family
scaffolding, but it also means the entire importance-sampling control surface is
ACCEPTED and then ignored, with no error and no warning.

Measured against the source (every name in `.impmapIsControlNames` grepped through
`R/np*.R` and `src/np*.cpp`, with the false positives -- `gammaOptimize`,
`gridBounds="auto"`, and npag's own `"gamma"` OUTPUT name -- excluded):

**Inert under `npag`/`npb` (22):** `isample`, `nIter`, `mapIter`, `gamma`,
`gammaMethod`, `gammaMethodUser`, `df`, `auto`, `autoNonNormal`,
`autoNonmemSparse`, `autoDfPatience`, `iscaleMin`, `iscaleMax`, `iaccept`,
`ctol`, `nConvWindow`, `impSeed`, `qr`, `qrShift`, `qrRefresh`, `sir`,
`sirSample`.

**Actually consumed:** only `impCov` and the internal M-step index maps
(`impMuThetaIdx`, `impMuEtaIdx`, `impThetaSensIdx`, `impOmegaFixedEta`).

Some of these are merely inapplicable -- there is no importance-sampling proposal
under a nonparametric engine, so `df`/`auto`/`iaccept`/`qr` cannot mean anything.
Others are outright traps:

* `nIter` -- npag counts `cycles`. Someone raising `nIter` to lengthen a fit
  changes nothing, silently.
* `impSeed` -- set for reproducibility; ignored. npb has its own `seed`, npag is
  Sobol-deterministic.
* `ctol` / `nConvWindow` -- convergence controls that do not drive npag's
  convergence (`rhoend`, `cycles`) at all.
* `gamma` -- collides in NAME with npag's `gammaOptimize` while meaning something
  unrelated (proposal-variance inflation / NONMEM `ISCALE`, versus a global
  assay-error multiplier on the residual). Passing `gamma` to `npag` looks like
  it should tune the assay error. It does nothing.

The promise of "npagControl is an impmapControl" is only half kept, and the
unkept half fails silently. That is the deficiency to fix.

## Phase 1 -- decide the policy per control

Three buckets, and every one of the 22 goes in exactly one. Do this first and
write it down, because the implementation is trivial once the classification is
settled and unarguable afterwards.

* **Inapplicable** -- no meaning without an importance-sampling proposal:
  `isample`, `gamma`, `gammaMethod`, `gammaMethodUser`, `df`, `auto`,
  `autoNonNormal`, `autoNonmemSparse`, `autoDfPatience`, `iscaleMin`,
  `iscaleMax`, `iaccept`, `qr`, `qrShift`, `qrRefresh`, `sir`, `sirSample`.
  These should be REJECTED when passed explicitly.
* **Should work and does not** -- a real np analogue exists, so wire it rather
  than reject it: `impSeed` (npb's `seed`), and possibly `mapIter` (npag does run
  a MAP pass through `impMapPass`; check whether it honours an iteration count).
* **Arguably remappable** -- `nIter`, `ctol`, `nConvWindow` have np counterparts
  (`cycles`, `rhoend`). Decide deliberately: alias, or reject with a message
  naming the right control. Aliasing is friendlier but creates two names for one
  thing; rejecting with a pointer is more honest. Prefer rejecting with a pointer.

Exit: a table in this file, control by control, with the decision.

## Phase 2 -- reject explicitly-passed inapplicable controls

The mechanism is small but has one subtlety: `impmapControl()` returns every
field populated with defaults, so the control object cannot tell you what the
caller actually asked for. Capture `...` BEFORE delegating:

    .dots <- list(...)
    .bad <- intersect(names(.dots), .npInertImpControls)
    if (length(.bad)) stop(...)

Decide error vs warning. An error is right here: these are silent no-ops today,
and a warning in a long fit scrolls past. But note the repo convention if a
warning is chosen anywhere -- warnings raised during a run land in `$runInfo`
and must stay under 75 characters and not name the method (CLAUDE.md).

Keep the message actionable: name the control, say it does not apply to a
nonparametric engine, and where relevant name the control that does
(`nIter` -> `cycles`, `ctol` -> `rhoend`, `gamma` -> `gammaOptimize`).

Watch: `npagControl()` is also called internally (e.g. `getValidNlmixrCtl.impmap`
rebuilds a control via `do.call`), so the rejection must not fire on a
round-trip of an already-built control. Test the round-trip explicitly -- this is
exactly the shape of bug that broke `est="fo"/"foi"` before (a control field the
round-trip rejected).

## Phase 3 -- wire what should work

* `impSeed`: make it drive npb's sampler seed, or alias `seed` and `impSeed` to
  one another with a single source of truth. npag is Sobol-deterministic, so
  state that in the docs rather than silently accepting a seed there.
* `mapIter`: determine whether `impMapPass` honours it under np. If it does, it
  belongs in the "consumed" list and the docs; if it does not, either wire it or
  reject it.

## Phase 4 -- resolve the `gamma` / `gammaOptimize` collision

Renaming `gammaOptimize` is a user-visible break and is NOT proposed. Instead:

* reject `gamma` under np with a message pointing at `gammaOptimize`;
* make the two roxygen entries cross-reference each other explicitly, each
  stating what the other one is not.

Revisit only if a future release is already breaking np controls for other
reasons.

## Phase 5 -- tests

The failure being fixed is silence, so the tests must assert the noise:

* every control in the inapplicable bucket raises when passed to `npagControl()`
  and `npbControl()`;
* a control that is legitimately shared (`impCov`) still passes through;
* `do.call(npagControl, npagControl())` round-trips -- the regression risk in
  Phase 2;
* the wired controls from Phase 3 demonstrably change the fit (a seed change
  moves npb draws), rather than merely being accepted.

Cheap control-level checks belong in the push/PR subset; anything needing a fit
goes in a `.slowBatches` entry.

## Phase 6 -- documentation

* `~/src/nlmixr2/vignettes/articles/nonparametric-npag-npb.Rmd` has a
  "Relationship to `impmap`" section that currently lists SEVEN inert controls.
  The real count is 22. Correct it, and replace the prose list with the Phase 1
  table.
* `npagControl()`/`npbControl()` roxygen should state plainly that they extend
  `impmapControl()` for the shared FOCEI-family scaffolding only, and that the
  importance-sampling controls are rejected rather than accepted.

## Note on scope

This is a user-visible behaviour change: code that today passes `isample` to
`npagControl()` and is silently ignored would start erroring. That is the point
-- such code is already not doing what its author intended -- but it belongs in a
release where the change can be called out, not slipped into a patch.
