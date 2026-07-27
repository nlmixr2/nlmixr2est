# Limiting the VAE covariate search with `shapes=`

Lets `shapes=` say that a covariate may take **no** shape on a parameter, which
removes that (parameter, covariate) pair from the automatic search. Today
`shapes=` can only choose among parameterizations; there is no way to tell the
search "consider `GA` and `Mage`, nothing else" or "search `GA` only on `lW0`"
without writing the effects into the model.

## Why the current answer is not enough

The only existing way to restrict which pairs are searched is `pinCovariates`,
and it is driven by the model rather than by the control:

```r
neonatalPin <- neonatalModel |>
  model(W0 <- exp(lW0 + beta.lW0.GA.power * log(GA / 40) + eta.W0)) |>
  ini(beta.lW0.GA.power <- 1)
```

That works, but it makes the user commit to three things at once to express one:

1. **A shape.** A pinned pair is allowed only the column matching the shape it
   was written in (`.vaeData.R` `.vaePinColumn`), so pinning `GA` on `lW0` also
   decides it is a `power`. There is no way to say "search `GA` on `lW0`, any
   shape".
2. **A starting value.** The coefficient has to be promoted in `ini()`, or the
   symbol is read as a data covariate and *nothing* is pinned -- a silent
   fall-through to the full search.
3. **All-or-nothing scope.** Declaring one effect pins the search to the
   declared set. There is no way to narrow one parameter while leaving the rest
   of the search alone.

`shapes=` already resolves per (parameter, covariate) with a specificity ladder,
aliases and case-insensitive covariate names. Exclusion is the one thing that
ladder cannot express, and it is the thing users ask for.

## Proposed surface: `fixCov`

Naming a covariate in `shapes=` **is** the statement that it belongs in the
search. `fixCov = TRUE` -- the default -- fixes the searched covariate set to
exactly the ones `shapes=` names; `fixCov = FALSE` restores today's behavior,
where `shapes=` restricts parameterizations only and every other covariate is
still searched with the default shapes.

```r
## search GA and Mage, nothing else -- fixCov = TRUE is the default
vaeControl(shapes = list(GA   = c("power", "hockey"),
                         Mage = "lin"))

## restrict GA's shapes, but keep searching every other covariate
## (this is what the same call means today)
vaeControl(shapes = list(GA = c("power", "hockey"),
                         fixCov = FALSE))

## per-pair: search GA only on lW0, and Mage only on lkin
vaeControl(shapes = list(
  list(var = "lW0",  covar = "GA",   shapes = "power"),
  list(var = "lkin", covar = "Mage", shapes = "lin")
))
```

The point is that the common case costs nothing. A modeler who knows the
covariates they care about lists them and is done -- there is no per-covariate
opt-out to write for the ones they left out, and no deny-then-allow preamble.
The escape hatch (`fixCov = FALSE`) is one element, written once, for the less
common case of "restrict these parameterizations, keep looking at the rest".

`fixCov` is an element of the `shapes` list, so it travels with the rules it
governs rather than as a second control argument that could disagree with them.

### Interaction with the rules

`fixCov` decides which pairs are **eligible**; the existing ladder
(`var`+`covar` > `covar` > `var` > global, ties to the last rule listed) then
decides which shapes an eligible pair may take. The two are separate passes,
which keeps the ladder untouched.

### Why not `"none"`

An earlier draft reserved a `"none"` shape token and had the user deny
per-covariate, or write a global deny rule and then allow back. Both put the
burden on the common case: an explicit covariate list is what modelers
actually write, so that should be the default reading, not an idiom assembled
out of exclusions.

### Backward compatibility

This changes what an existing `shapes = list(...)` call means. Today
`list(GA = "power")` searches every covariate with `GA` restricted to `power`;
under the new default it searches `GA` alone. Any user who was restricting
shapes and relying on the rest of the search continuing has to add
`fixCov = FALSE`. See decision 3.

## Decisions (confirmed with the user before implementation)

| # | Decision | Choice |
|---|----------|--------|
| 1 | Pair rules under `fixCov=TRUE` | a `var`+`covar` rule makes only that PAIR eligible; a `covar`-only rule makes it eligible on every parameter |
| 2 | Naming a categorical | `Sex = TRUE` -- "eligible, shapes not applicable". `TRUE` works for a continuous covariate too, meaning "eligible with the default shapes" |
| 2b | **Mixed list forms must work** | `list(list(var=, covar=, shapes=), Sex = TRUE)` -- named entries and pair entries in ONE list. This is new; today mixing is an error |
| 3 | Backward compatibility | no separate warning: report it through `$runInfo`, following the conventions already used there |
| 4 | Where `fixCov` lives | element of the `shapes` list only; NOT a `vaeControl()` argument, so it cannot disagree with the rules it governs |
| 5 | `shapes` as a character vector | names no covariate, so `fixCov` has nothing to fix: every covariate stays searched. `fixCov` alongside a character vector is an error |
| 6 | Global rule under `fixCov=TRUE` | a rule with neither `var` nor `covar` makes everything eligible, negating `fixCov`; error rather than a silent contradiction |
| 7 | `var`-only rule under `fixCov=TRUE` | every covariate eligible on that parameter, and only that parameter |
| 8 | A data column named `fixCov` | matched exactly and case-sensitively as a control element; covariates are upper-cased, so a real `FIXCOV` column is still reachable. Error on collision |
| 9 | Nothing eligible | error, naming `covariateSelection = FALSE` as the intended spelling |
| 10 | Interaction with pinning | a model that declares covariates already pins; declaration wins, and `$runInfo` says so when `fixCov` would have disagreed |
| 11 | Reporting | the eligible set, and every covariate `fixCov` excluded, in `$runInfo` (this subsumes decision 3) |

### Decision 2b is the invasive one

`.vaeResolveShapes` currently decides the form with
`all(vapply(spec, is.list, logical(1)))` -- every element a list means the pair
form, otherwise the named form, and anything mixed dies on "shapes list must be
named by covariate". Supporting `list(list(...), Sex = TRUE)` means replacing
that whole-list test with a **per-element** one:

- element is a list -> a pair rule (`var`/`covar`/`shapes`)
- element is named -> shorthand for `list(covar = <name>, shapes = <value>)`,
  where `<value>` may be a shape vector or `TRUE`
- element is named `fixCov` -> the control flag, not a rule

A named entry becoming exact shorthand for a `covar`-only pair rule is what
keeps this from adding a second set of semantics: after parsing there is one
rule table, and the existing specificity ladder is untouched.

Note that this makes the "the two forms cannot be mixed" paragraph in the
`vaeNeonatal` vignette wrong once shipped; the article needs updating in the
same PR that lands this.

## Implementation

Two separable pieces: parsing `fixCov` out of the `shapes` list, and turning an
eligibility set into mask zeros. The mask machinery already exists and already
propagates.

1. **`R/vaeCovShapes.R`**
   - `.vaeResolveShapes`: replace the whole-list form detection with the
     per-element dispatch of decision 2b, and pull `fixCov` out before any of
     it. Ordering matters -- a bare `fixCov = FALSE` is not a list, so leaving
     it in would flip the old whole-list test on its own. Return
     `list(rules =, fixCov =)` rather than an attribute; there are three call
     sites and the explicit shape reads better at all of them.
   - `TRUE` as a shape value (decision 2) resolves to "eligible, default
     shapes" -- for a categorical it means eligible with no parameterization to
     choose, which is already how a `cat` column behaves.
   - New `.vaeEligible(rules, fixCov, etaNames, thetaForEta, cov)` returning an
     (eta x raw covariate) logical matrix, following decisions 1, 6 and 7. With
     `fixCov = FALSE` it is all-TRUE and nothing downstream changes.
   - `.vaeShapesFor`: unchanged. Eligibility and parameterization stay separate
     passes so the specificity ladder is untouched.
   - `.vaeShapeAllowMask`: zero every column of an ineligible cell. Both
     existing guards need care -- the `cat` skip and the "requested family
     unavailable -> leave selectable" fallback each protect a covariate from
     being masked out entirely, which is right for a parameterization
     restriction and wrong for an eligibility one.
   - `.vaeAssertContShapes` gains no token, which is the main simplification
     over the `"none"` draft.
2. **`R/vaeData.R`**
   - `.vaeDataPrep` (~line 1007): the mask is built only when
     `!.searchOff && !.pinActive`. Under pinning, resolve `fixCov` far enough
     to emit the decision-10 warning, then discard it.
   - `.vaeCovariateSearch` builds the candidate columns; a covariate ineligible
     on every parameter can be dropped there rather than masked, which is
     cheaper -- but only if that does not disturb `covGroup`/`covBlock`
     numbering. Masking is the safe default; measure before optimizing.
   - Line ~1218 (`.sel <- which(colSums(.covAllow) > 0L)`) already drops a
     column allowed for no dimension, so an excluded covariate leaves the
     design matrix with no extra work.
   - `$runInfo` note (decision 11), next to the existing
     `catDropped`/`logDrop`/`hockeyDrop` warnings; decision 9's error.
3. **`R/vaeFit.R`** -- nothing. `.nCand` is computed from `prepC$covAllow`
   (~line 271), so an exclusion shrinks the bit count that picks
   `bnb`/`l0learn` automatically. Worth a test rather than a change.
4. **`R/vae.R`** -- `shapes=` roxygen: `fixCov`, the default, the
   backward-compatibility note, and that pinning overrides it.
5. **`NEWS.md`** -- decision 3 makes this a user-visible behavior change and it
   needs an entry whatever we choose.

## Tests

- `shapes = list(GA = "power")` searches `GA` and nothing else (assert on the
  `covAllow` mask -- cheap and exact -- and on a fit outcome once).
- The same call with `fixCov = FALSE` reproduces today's mask exactly. This is
  the regression guard for decision 3.
- Pair rules: `list(list(var = "lW0", covar = "GA", shapes = "power"))` makes
  `GA` eligible on `lW0` only (decision 1).
- `fixCov` in a pair-form list does not break form detection -- the specific
  failure the strip-first ordering exists to prevent.
- Mixed forms (decision 2b): `list(list(var=, covar=, shapes=), Sex = TRUE)`
  parses, and a named entry resolves identically to the `covar`-only pair rule
  it is shorthand for.
- `TRUE` as a shape value on a continuous covariate means "eligible, default
  shapes"; on a categorical it means eligible.
- A categorical covariate named in the list is searched; one not named is not
  (decision 2, spelling TBD).
- Nothing eligible raises the decision-9 error.
- A pinned model with a disagreeing `fixCov` warns, and the declaration wins
  (decision 10).
- `covSelectMethod`: excluding covariates moves a dimension from `l0learn` back
  to `bnb`.

## Working agreement

- **Isolated worktree.** Implement on a branch stacked on `feat/vae-hockey`, in
  its own worktree (the repo already keeps one per line of work, e.g.
  `~/src/nlmixr2est-vae-hockey`), so the shapes work never lands on whatever
  branch the main checkout happens to be on.
- **Commit often.** One commit per numbered implementation step above, not one
  commit at the end -- the parse change, the eligibility matrix, the mask
  change, the `runInfo`/error handling, the docs and the NEWS entry are each
  independently reviewable and independently revertable.
- **Review with antigravity.** Run the `antigravity-code-review` skill on the
  branch diff before the work is considered done, and again after any
  non-trivial round of fixes. A Gemini-family reviewer is genuinely independent
  of the reading that produced the code, which is the point.

## Note on branches

The code above is on `feat/vae-hockey` (`.vaeDefaultShapes`, `.vaeHockeyArms`,
the block-aware `.bitsOf`), not on `main`. This plan is written against that
branch and should be implemented on top of it or on a branch stacked on it.
