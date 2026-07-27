# VAE hockey-stick covariate shape

Adds `"hockey"` to the VAE covariate search shape vocabulary: a continuous
covariate may enter a parameter as a two-armed piecewise-linear relationship
with the knot at the covariate's centering value.

## Written form

```r
ka <- exp(tka
          + beta.tka.WT.hockey.low * (WT <  70.5) * (WT - 70.5)
          + beta.tka.WT.hockey.hi  * (WT >= 70.5) * (WT - 70.5)
          + eta.ka)
```

Continuous at the knot (both arms vanish at `WT == 70.5`), so the structural
theta keeps its usual meaning: the parameter value AT the center.  The arms are
disjoint (`<` / `>=`) so they partition subjects exactly and
`arm.low + arm.hi == (WT - ctr)` identically.  Both terms are written with `+`;
the sign of each coefficient comes from the fit.

## Decisions (confirmed with the user before implementation)

| # | Decision | Choice |
|---|----------|--------|
| 1 | Arm expression | `*(cov - center)` -- continuous kink, not `*cov` |
| 2 | Knot boundary | disjoint `<` / `>=` |
| 3 | Joint selection | atomic block in C++ (`covBlock`) |
| 4 | Degenerate arm guard | >= `catCutoff` (5%) of subjects each side, else warn + drop |
| 5 | Coefficient names | `beta.<par>.<COV>.hockey.low` / `.hockey.hi` |
| 6 | Hand-written hockey | no detection needed -- mu2 already handles it (below) |
| 7 | Default `shapes=` | hockey IS in the default set |
| 8 | L0Learn | promoted `Suggests:` -> `Imports:`, own commit in this PR |

## Why no `.vaeDetectShape` / pinning work

A hand-written hockey stick is ALREADY handled end to end.  Each arm's
derivative with respect to its own coefficient is a data-only expression, so
rxode2's mu2 detection recognizes both arms independently.  Verified against
rxode2 5.1.5:

```
  theta       covariate covariateParameter       modelExpression
1   tka WT * (WT >= 70)                hk2 hk2 * (WT >= 70) * WT
2   tka WT * (WT <= 70)                hk1 hk1 * (WT <= 70) * WT
```

`.uiModifyForCovs` materializes one `nlmixrMuDerCov#` data column per arm and
rewrites the term as `nlmixrMuDerCov# * <coef>`.  `R/vaeData.R:238` already
encodes `nlmixrMuDerCov#` columns as ordinary `family = "lin"` search columns,
one coefficient each, each in its own mutual-exclusion group.  So pinning,
`.vaeCoefFactor`, `.vaeDetectShape` and `.vaeUpdateModelPinned` need no change:
the work here is only making the search PROPOSE a hockey for a covariate the
user did not declare, and writing it back.

## Why an atomic block is required

The two arms sum to the linear column:

```
arm.low + arm.hi == (cov - ctr)
```

so `span{1, arm.low, arm.hi}` and `span{1, (cov - ctr), arm.low}` are the SAME
space.  Both cost two columns, so both carry the same BICc penalty and produce
identical RSS -- an exact tie, broken lexicographically on column index.  Left
alone the search would write `beta...lin*(WT-70) + beta...hockeyLow*...`: same
fit, but a model that does not read as a hockey stick.  Making the two arms an
all-or-none block that is mutually exclusive with the covariate's `lin`/`power`
columns removes the tie at the source.

## Phases

Each phase is one commit: merge `origin/main`, run the relevant tests,
antigravity review, push.

0. **plan** (this file).
1. **L0Learn `Suggests:` -> `Imports:`** plus removal of the now-dead
   `requireNamespace` guards and "install it" errors in `R/vaeCovSelectL0.R`.
   Independent of everything below; motivated by phase 6 (hockey raises the
   per-covariate bit cost, so more searches reach the L0 engine).
2. **Shape vocabulary** (`R/vaeCovShapes.R`).  Two internal arm shapes
   `hockeyLo` / `hockeyHi` in a new `"hockey"` family; `.vaeShapeExpr`,
   `.vaeShapeBeta`, `.vaeShapeValue`, `.vaeShapeUsable`, `.vaeShapeFamily`.
   Also renames the four tests that currently use `"hockey"` as their canonical
   INVALID shape name.  Pure functions, no behavior change yet.
3. **Design-matrix emission** (`R/vaeData.R`).  `.vaeCovFamilies` hockey branch,
   the two arm columns sharing one `covBlock` id and the covariate's `covGroup`,
   the 5%-per-side guard with its warning, `covBlock` on the prep list.
4. **Block-atomic branch and bound** (`src/inner.cpp`).  Branch on blocks rather
   than columns; group exclusion compares blocks; `vaeBestSubset_` /
   `vaeScoreSupports_` gain `block=`.  `Rcpp::compileAttributes(".")` AND the
   hand-edited `src/init.c` prototype + arity.
5. **L0 candidate block completion** (`R/vaeCovSelectL0.R`, `src/inner.cpp`).
   L0Learn proposes supports blind to blocks, so partial blocks are completed
   (or dropped) before scoring, mirroring the existing group repair.
6. **Bits accounting** (`R/vaeFit.R:249`).  Count feasible BLOCKS per group, not
   columns, so a hockey covariate costs `log2(4)` and not `log2(5)`.
7. **Write-back** (`R/vaeOutput.R`).  Two terms and two coefficients from one
   selected block; the `beta.<par>.<COV>.hockey.low` / `.hi` name branch.
8. **Default set, docs, baselines**.  Add `"hockey"` to the two literal default
   vectors, regenerate `man/`, `NEWS.md`, re-bless any moved selection baseline,
   place new slow fit tests into a `.slowBatches` batch.

## Evaluation criteria

What "done and correct" means.  Each is a check that can fail, not a sentiment.

1. **No unrequested drift.**  With `shapes=` set to today's five, every VAE
   covariate fit must be numerically identical to `origin/main` -- same
   `selected`, `beta`, objective.  Guaranteed structurally by `covBlock`
   defaulting to singleton blocks (identical BnB tree, strength order and
   tie-break), and checked by running the vae covariate tests against both
   trees.  Drift is only acceptable where it comes from hockey now being in the
   DEFAULT set, and every such case must be explained, not just re-blessed.
2. **The block search is still exact.**  For small designs, enumerate every
   feasible support in R, score with the same BICc objective, and compare to
   `vaeBestSubset_`.  Must cover: hockey wins, `lin` wins, nothing is selected,
   and a design where the tie described above exists.
3. **Round-trip fidelity.**  The single strongest check.  For every subject and
   every latent dim, the WRITTEN model text evaluated on the data must reproduce
   the M-step's `zPopMat[i, k]` to ~1e-10: `theta_written + sum(beta_written *
   expr_written) == zPopMat`.  This catches every reparameterization, centering
   and intercept-adjustment error, and it is what makes the `<` / `>=` boundary
   choice verifiable rather than assumed.
4. **The written model parses and fits.**  Pipe the written model back through
   `rxode2::assertRxUi` and a `maxOuterIterations = 0` focei fit; confirms the
   indicator comparisons survive rxode2 parsing and that the injected terms keep
   the additive `exp(theta + ... + eta)` form `muRefCurEval` needs.
5. **Nesting and penalty behave.**  Hockey's span strictly contains `lin`'s, so
   RSS(hockey) <= RSS(lin) always.  On simulated data with a true kink hockey is
   selected; on data with a true straight line `lin` is selected and hockey is
   rejected -- no false positives bought by the extra column.
6. **Engine parity.**  With `covSelectMethod = "l0learn"` the selected support
   must match the exact path on a design small enough to run both, including a
   design where L0Learn proposes a partial block.
7. **Guard fires correctly.**  A `covCenter=` override putting <5% of subjects
   on one side must drop that covariate's hockey columns, warn once in under 75
   characters with no method-name prefix, and leave `lin`/`power` selectable.
8. **Registration hygiene.**  Signature changes mean `Rcpp::compileAttributes(".")`
   AND a hand-edited `src/init.c` prototype and arg count; `rm -f src/*.o src/*.so`
   before rebuilding.  Verified by a clean rebuild plus a direct `.Call` smoke test,
   since a missed arity survives a clean build and only fails at run time.
9. **House style.**  camelCase R names, ASCII only, warnings under 75 characters
   with no method prefix, short roxygen, `NEWS.md` entry under `## New features`.
10. **Test placement.**  Unit tests (shape expressions, block feasibility,
    brute-force parity, guard) stay in the essential push/PR set; multi-iteration
    hockey FIT tests go into a `.slowBatches` batch and must not use
    `skip_on_ci()`.
11. **Independent review.**  Every commit is reviewed with antigravity; each
    finding is either fixed or dismissed with a stated reason.
