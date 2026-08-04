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

**Actually consumed:** `impCov`, and -- as of the fix below -- the internal M-step
index maps (`impMuThetaIdx`, `impMuEtaIdx`, `impThetaSensIdx`,
`impOmegaFixedEta`).

An independent review of an earlier draft of this plan corrected that second
half.  The index maps were NOT in fact consumed: `foceiSetup_` loaded them inside
`if (op_focei.isImpmap)`, which is `(impmap || imp || qrpem)` and therefore false
for np, so R computed all four in `.npFamilyControl` and C++ discarded them.
`npagOuter` calls `impMuInterceptStep()` (which iterates `impMuThetaIdx` and does
`fullTheta[th] += delta`) and `impGetOmegaFixedEta()`, so in a clean session the
mu-intercept step ran zero iterations and no eta was reported omega-fixed, while
after an impmap fit in the SAME session the process-global `op_focei` still held
that fit's indices and np indexed `fullTheta` with another model's map.  Fixed in
2d5cffaaf by loading all four under `isNpag || isNpb` and clearing them
otherwise.  Session-order-dependent behaviour of this kind is invisible to the
test suite, since a test file that runs np alone only ever sees the silent
no-op.

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

* **Internal, never user-facing -- must NOT be rejected**: `gammaMethodUser`
  (stamped to record the pre-`"auto"` preference), `impCov`, and the four M-step
  index maps (`impMuThetaIdx`, `impMuEtaIdx`, `impThetaSensIdx`,
  `impOmegaFixedEta`). A caller never types these; a REBUILD always carries them.
  Rejecting them breaks exactly the round-trips Phase 2 is trying to survive, so
  they come off the reject list entirely -- and if the "detect internal fields"
  option is chosen, this is precisely the detection signature.

  Note `sirSample` is computed and stamped by `impmapControl()` when absent, so a
  rebuilt control always carries it too; it is nonetheless a real user-facing
  argument, so it stays in the reject bucket and the rebuild detection must not
  depend on it.
* **Inapplicable** -- no meaning without an importance-sampling proposal:
  `isample`, `gamma`, `gammaMethod`, `df`, `auto`, `autoNonNormal`,
  `autoNonmemSparse`, `autoDfPatience`, `iscaleMin`, `iscaleMax`, `iaccept`,
  `mapIter`, `qr`, `qrShift`, `qrRefresh`, `sir`, `sirSample`.
  These should be REJECTED when passed explicitly. `mapIter` is here on
  evidence, not assumption: `impMapPass()` (`src/inner.cpp`) runs a single
  `foceiOuterFinal()` pass, and `mapIter` is never loaded into C++ by
  `foceiSetup_` nor read by `impMapPass()`, so it cannot mean anything under np.
* **Should work and does not** -- a real np analogue exists, so wire it:
  `impSeed` (npb's `seed`).
* **Arguably remappable** -- `nIter`, `ctol`, `nConvWindow` have np counterparts
  (`cycles`, `rhoend`). Decide deliberately: alias, or reject with a message
  naming the right control. Prefer rejecting with a pointer.


Exit: a table in this file, control by control, with the decision.

## Phase 2 -- reject explicitly-passed inapplicable controls

The obvious mechanism is WRONG, and review caught it:

    .dots <- list(...)                                   # NOT sufficient
    .bad <- intersect(names(.dots), .npInertImpControls)
    if (length(.bad)) stop(...)

`impmapControl()` returns every field populated, so a round-trip --
`do.call(npagControl, npagControl())`, `.foceiFamilyControl()`, or
`getValidNlmixrCtl.impmap()` rebuilding a control -- passes `isample`, `nIter`,
`df` and the rest through `...` as ordinary named arguments. The naive check
rejects them and breaks every internal rebuild. This is the same shape as the
bug that broke `est="fo"/"foi"`: a control field the round-trip could not
survive.

So the check must distinguish an explicit call from a rebuild. Options, in
preference order:

1. **Stamp the control.** Have `npagControl()`/`npbControl()` set a marker
   (e.g. `.ctl$.npBuilt <- TRUE`) on what they return. A rebuild passes that
   marker back through `...`; its presence means "not a fresh user call, skip
   validation". Explicit and cheap to reason about.
2. **Detect internal-only fields.** A rebuild carries names no caller would
   type (`gammaMethodUser`, `impMuThetaIdx`, `impThetaSensIdx`, `impCov`). Any
   of those present implies a rebuild. Works, but relies on that set staying
   internal.

Both options share two failure modes that must be designed against, not
discovered:

* **Validating only in the constructor is not enough.** `getValidNlmixrCtl.npag`
  routes a RAW list through `.npValidCtl` -> `getValidNlmixrCtl.impmap` ->
  `do.call(impmapControl, .)`, never entering `npagControl()`. So
  `getValidNlmixrCtl.npag(list(isample = 500))` bypasses any check that lives
  only in the constructor. Validation has to sit where every path converges.
* **A control edited after construction bypasses it too.**
  `ctl <- npagControl(); ctl$isample <- 5000` carries the marker (or the internal
  fields) and so reads as a rebuild, silently ignoring the edit -- which is the
  original bug wearing a different hat. Either validate at fit time as well, or
  accept and document that post-hoc edits are not checked.

A sequencing detail that decides where the check can live: `npagControl()` calls
`impmapControl(...)` BEFORE stamping `.ctl$est <- "npag"`, so a check placed
inside `impmapControl()` cannot tell it was invoked from the np constructor. The
check therefore belongs in `npagControl()`/`npbControl()` and in `.npValidCtl()`
(which is where the raw-list path converges), not in `impmapControl()`.

Whichever is chosen, the round-trip must be an explicit test (Phase 5), not an
assumption.


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

## Phase 2b -- npbControl's inert FORMAL arguments

`npbControl()` declares `cycles = 100L` and `gammaOptimize = FALSE` as named
FORMALS, documented "Unused for npb". Because they are formals they never reach
`...`, so a `list(...)` check cannot see them: `npbControl(cycles = 500)` stays a
silent no-op even after Phase 2 lands.

Handle them explicitly -- validate the formals as well as `...`, or remove them
from the signature if nothing depends on their presence. Check first whether any
internal caller passes them positionally before changing the signature.

## Phase 3 -- wire what should work

* `impSeed`: `npbControl()` takes its own `seed = 42L`, and `src/npb.cpp` reads
  `control["npSeed"]`, populated from `seed`. `impSeed` reaches
  `impmapControl()` and is then ignored, so `npbControl(impSeed = 123L)`
  silently does nothing where someone plainly meant to set the sampler seed.
  Either alias the two to a single source of truth or reject `impSeed` pointing
  at `seed`. npag is Sobol-deterministic, so say that in the docs rather than
  accepting a seed there at all.

(`mapIter` was investigated and is NOT honoured -- `impMapPass()` runs a single
`foceiOuterFinal()` pass and `mapIter` never reaches C++ -- so it moved to the
reject bucket in Phase 1 rather than being wired.)

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
* a control that is legitimately shared (`impCov`) still passes through, and no
  internal field (`gammaMethodUser`, the index maps) is ever rejected;
* `do.call(npagControl, npagControl())` AND `do.call(npbControl, npbControl())`
  round-trip -- both engines, since they are separate constructors;
* `getValidNlmixrCtl.npag(list(isample = 500))` and the npb equivalent raise,
  covering the raw-list path that bypasses the constructor;
* `getValidNlmixrCtl.npag(fit$control)` on a completed fit still validates;
* a post-construction edit (`ctl$isample <- 5000`) behaves as Phase 2 decided --
  whichever way that goes, pin it so it cannot change silently;
* the wired controls from Phase 3 demonstrably change the fit (a seed change
  moves npb draws), rather than merely being accepted.


Cheap control-level checks belong in the push/PR subset; anything needing a fit
goes in a `.slowBatches` entry.

## Phase 6 -- documentation

* `~/src/nlmixr2/vignettes/articles/nonparametric-npag-npb.Rmd` has a
  "Relationship to `impmap`" section that lists SEVEN inert controls.
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

## Phase 7 -- the general form: op_focei fields that leak across fits

Both bugs found while writing this plan were the same shape: a field loaded into
the process-global `op_focei` only `if (foceiO.containsElementNamed(...))`, with
no clear on absence, so a later fit inherits an earlier one's value. Two were
fixed (np in 2d5cffaaf, advi in 46f00ceb2). A sweep of `foceiSetup_` finds the
pattern is broader:

* the `isImpmap` index-map loads have no `else` clear, so an impmap fit whose
  control lacks them inherits the previous fit's;
* `est="focei"`, `est="saem"` and friends enter none of the three blocks, so the
  index maps are never cleared for them at all;
* the `imp*` scalars (`impCov`, `isample`, `gamma`, `nIter`, `impSeed`, ...) are
  likewise conditional, and `rxOptionsFreeFocei()` does not reset them.

The np and advi cases mattered because those methods READ the fields without
always setting them. For methods that never read them, a stale value is inert
today -- but it is inert by accident, and the next method to reuse one of these
helpers inherits the hazard.

The durable fix is to reset the whole `imp*` group unconditionally at the top of
`foceiSetup_` (or in `rxOptionsFreeFocei()`) and let each block populate what it
needs, rather than clearing field by field at each call site. That is a
self-contained change and is NOT part of this release; it is recorded here so the
pattern is not rediscovered a third time.
