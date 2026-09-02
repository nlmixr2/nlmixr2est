# Two-level (IOV) SAEM after Panhard & Samson (arXiv:0803.4437)

## Context

`nlmixr2est` advertises IOV for SAEM, but **there is no IOV in the SAEM kernel at all**
-- `grep -i "iov\|occ\|kappa\|level" src/saem.cpp` finds nothing relevant, and
`R/saem.R:1544` only sets `attr(nlmixr2Est.saem, "iov") <- TRUE` to opt into a generic
R-level pre-processing hook.  That hook, `.uiApplyIov()` (`R/iov.R:104-349`), rewrites

```r
iov.cl ~ 0.1 | occ
```

into **one population theta carrying the magnitude** plus **K fixed unit-variance etas**,
combined *multiplicatively* (`R/iov.R:306-318`):

```r
rx.iov.cl <- sqrt((iov.cl)^2) * (rx.iov.cl.1*(occ == 1) + rx.iov.cl.2*(occ == 2))
```

so the kernel sees `iov.cl` in the **phi0** block (a fixed effect carried as a degenerate
random effect) and `rx.iov.cl.*` in **phi1** with variance pinned at 1.
`ui$nonMuEtas == c("rx.iov.cl.1","rx.iov.cl.2")` -- **the occasion etas are
non-mu-referenced by construction**, because they enter multiplied by a theta.

### Why this is worth changing

**It measurably does not work.**  Same model, same `theo_md` + `occ` data, fitted IOV
variance:

| method | fitted IOV variance |
|---|---|
| `focei` (default) | **0.003348** |
| `saem` nBurn=200 / nEm=300, seed 42 | **0.000302** |
| `saem` nBurn=15 / nEm=15 | **1.65e-06** |

An order of magnitude low at production iteration counts, and still drifting toward zero.
The mechanism is structural, not a tuning problem:

1. **Mu-referencing is destroyed**, so there is no closed-form M-step.  `.configsaem`
   splits phi purely on `diag(covstruct)` (`R/saem_fit.R:447-450`), so the IOV *variance*
   is estimated as a phi0 parameter -- first by a stochastic sampled mean of MCMC draws
   whose pseudo-variance `Gamma2_phi0` is deliberately decayed (`coef_phi0 = .9638`,
   `src/saem.cpp:3388-3401`), then by a bounded bobyqa in a +-0.75 trust region
   (`refinePhi0Lik`, `src/saem.cpp:968-1170`).  Every other variance component in SAEM
   gets a closed-form update; this one gets simulated annealing.
2. **`w` and the eta variance are jointly unidentified during early burn-in.**  Only
   `w * sigma_eta` is identified, and the unit variance is re-pinned only after
   `nb_fixOmega = round(nBurn * perFixOmega)` (default 10% of burn-in),
   `src/saem.cpp:3374-3383`.
3. **No positivity constraint.**  `R/iov.R:276` sets `.curTheta$lower <- 0` with the
   comment `# doesn't work with saem`; combined with the deliberate `sqrt(theta^2)`
   spelling of `|theta|` (required for a correct symengine derivative, #952) the
   objective is exactly symmetric about zero, so there are two equivalent modes at
   `+-sd`.  (The bound *is* honoured, via `phi0Lower` -- but `warnRxBounded` still emits
   `"...boundaries: iov.cl which are ignored in 'saem'"` into `$runInfo` on every IOV
   fit.)
4. **Correlated IOV is impossible** -- refused at `R/iov.R:128-138`, and again by
   rxode2's `assertRxUiIovNoCor` from `R/saem.R:1518`.
5. **`nphi1` inflates by K per IOV parameter**, and `calc.2LL`'s objective is a full
   `nnodes.gq^nphi1` grid (`R/saem_fit_aux.R:137`).  This is where **issue #1000** lives.

Panhard & Samson (2009, *Biostatistics*; arXiv:0803.4437) is the SAEM extension for
exactly this model.  It uses the collapsed, **strictly mu-referenced** representation

```
phi_ik = mu + beta_k + b_i + c_ik,   b_i ~ N(0, Omega),  c_ik ~ N(0, Psi)
```

and gives **closed-form M-step updates for both Omega and Psi** by integrating `b_i` out
analytically (it enters linearly), leaving MCMC to sample only `phi_ik`.

**Intended outcome:** IOV stops being a nuisance theta and becomes a second variance
component with the same closed-form treatment Omega already gets; correlated IOV and
occasion fixed effects `beta_k` become expressible; SAEM's IOV estimates stop collapsing
toward zero; issue #1000 is closed; and the paper's simulation study is reproduced to
prove the algorithm is right.

## The algorithm being implemented

With `V(theta) = (Omega^-1 + K Psi^-1)^-1` and
`m(phi_i, theta) = V(theta) (Psi^-1 sum_k (phi_ik - beta_k) + Omega^-1 mu)`, the
sufficient statistics (each updated by `s <- s + gamma_l (S - s)`) and M-step are

```
s1_i = sum_k phi_ik                                  mu     = V Psi^-1 (mean_i s1_i - sum_k beta_k) + V Omega^-1 mu_old
s2_k = sum_i phi_ik                                  beta_k = s2_k/n - mu            (k >= 2, beta_1 = 0)
s3   = sum_i m(phi_i)' m(phi_i)                      Omega  = V + s3/n    - mu' mu
s4   = sum_ik (phi_ik - m)'(phi_ik - m)              Psi    = V + s4/(nK) - (1/K) sum_k beta_k' beta_k
s5   = sum_ijk ((y - f)/g)^2                         sigma2 = s5 / sum_ik n_ik
```

MCMC: three kernels on `phi_i` -- independence proposal from the prior
`N(mu + beta, Gamma)`, block random walk `N(phi_i, rho Gamma)`, componentwise random walk
-- with `Gamma` the joint covariance of `phi_i` across the K occasions:
`Omega + Psi` on the diagonal blocks, `Omega` off-diagonal.

**This maps onto the existing kernel almost exactly.**  `do_mcmc`'s three kernels
(`src/saem.cpp:4493-4512`) *are* the paper's three kernels, and kernel 1 already proposes
from `chol(Gamma2_phi1)` -- the full covariance.  Giving `Gamma2_phi1` the
compound-symmetry structure above gets the paper's proposals for free.  The real work is
the structured Omega/Psi M-step replacing the unstructured moment update at
`src/saem.cpp:3276-3386`, and letting one theta feed K phi columns in the `LCOV`/`MCOV`
design.

## Design

### Control gate -- the existing `"iov"` attribute is the switch

`.isIovMethod(est, control)` (`R/iov.R:13-23`) already accepts a
`function(control)` for the `"iov"` attribute.  So the legacy rewrite is turned off, and
SAEM's own formal handling turned on, with no new hook plumbing:

```r
attr(nlmixr2Est.saem, "iov") <- function(control) {
  # TRUE  -> run the legacy .uiApplyIov() magnitude-theta rewrite
  # FALSE -> saem handles `| occ` natively (two-level Panhard & Samson)
  identical(control$iovMethod, "theta")
}
```

`saemControl(iovMethod = c("twoLevel", "theta"))`, defaulting to `"twoLevel"`.  When the
hook is off, `.uiApplyIov()` takes its existing early-out (`R/iov.R:104-113`), clears
`.uiIovEnv`, and the ui reaches the SAEM build with the user's `iov.cl ~ 0.1 | occ` row
intact.  `.uiFinalizeIov()` likewise no-ops (it already guards on
`is.null(.uiIovEnv$ui)`).  Every other method keeps `attr(..., "iov") <- TRUE` and the
legacy path unchanged.

### Column layout

Split the mu-referenced parameters into A (no IOV, one phi column each) and B (IOV, K
columns each, `rx.<v>.<k>`).  The IIV eta of a B parameter is **absorbed** -- its K
columns carry `b + c_k` jointly, exactly as in the paper.  So

```
Gamma[A , A  ] = Omega_AA
Gamma[A , B_k] = Omega_AB                       (the same for every k)
Gamma[B_k,B_l] = Omega_BB + (k == l) * Psi
```

a strict generalization of the paper (which has A empty), and what makes correlated IOV
(a full `Psi`) and IIV/IOV cross-correlation representable at all.

The M-step generalizes the same way.  `b_A` is observed exactly (`b_A = phi_A`), so
condition on it first: `Omega_B|A = Omega_BB - Omega_BA Omega_AA^-1 Omega_AB`,
`V_i = (Omega_B|A^-1 + K_i Psi^-1)^-1` (per-subject `K_i` handles unbalanced occasions),
`m_i` as in the paper using `Omega_B|A` and the conditional mean.  Then

```
Omega_AA += (phi_A - mu_A)(.)'
Omega_AB += (phi_A - mu_A)(m_i - mu_B)'
Omega_BB += V_i + (m_i - mu_B)(.)'
Psi      = ( sum_ik [ (phi_Bik - beta_k - m_i)(.)' + V_i ] ) / sum_i K_i
```

With A empty, balanced `K_i == K` and `beta == 0` these collapse to the paper's formulas
verbatim -- which is the unit test.

### Fixed effects stay on the existing GLS step -- the mu-referencing payoff

`mprior_phi1 = COV1 * MCOV1` with
`Plambda1 = inv(CGamma21) * sum(D1Gamma21 % (COV1' statphi11))` (`src/saem.cpp:3215`) is
already `lambda = (sum_i C_i' Gamma^-1 C_i)^-1 sum_i C_i' Gamma^-1 E[phi_i]`, i.e. the
correct GLS for `phi_i ~ N(D_i lambda, Gamma)`.  Nothing about it changes except that
**one theta must be allowed to feed all K columns of a B parameter**.  Today each
`lambda` maps to exactly one phi column (`R/saem_fit.R:489-530`, `jcov1 <- grep(1,
LCOV1)`).  Add a `lambdaGroup` map so several live `MCOV1` entries share one `Plambda1`
element, and replace the `COV21 % D2Gamma21` Hadamard shortcut (valid only in the
one-to-one case) with the explicit `sum_i R_i Gamma^-1 R_i'`.

`beta_k` then needs **no new syntax**: it is a theta whose `mcov` incidence covers only
the occasion-k columns.

### Fallback conditions

Fall back to `iovMethod = "theta"` (with a `$runInfo` note, not an error) when: more than
one occasion variable is present, an IOV parameter is not mu-referenced, or the user asks
for a non-`sd` `iovXform`.  The legacy code path stays intact and reachable, and remains
the only path for every non-SAEM method.

### Known risk to spike early

With the hook off, `cl <- exp(tcl + eta.cl + iov.cl)` puts **two etas on one theta**.
rxode2's mu-reference detection must accept that and report both in `muRefDataFrame`; if
it instead demotes one to `nonMuEtas`, the whole design is undercut.  Spike this in
Phase 2 before building anything on top of it.  Similarly, rxode2's
`assertRxUiIovNoCor` (called from `R/saem.R:1518`) refuses correlated IOV before SAEM
ever runs -- it has to be made conditional on `iovMethod` (here or upstream) for the
correlated-`Psi` capability to be reachable.

## Phases

Every phase ends with a commit **and a `git push`** -- the branch `iov-refine-saem-1000`
is not yet on origin.

### Phase 0 -- setup
Push the branch; commit this plan to `plans/saem-iov-two-level.md` (matching the existing
`plans/*.md` convention).

### Phase 1 -- diagnose and fix #1000 (independent of the refactor, lands first)
The issue localizes the blowup to the objective, not the fit.  Two mechanisms are already
visible in `calc.2LL` (`R/saem_fit_aux.R:86-190`) and both are plausible:

- **No clamp on the `ll()` path.**  The C++ kernel guards general-likelihood values --
  `_saemGenLikCeiling = 700`, `_saemGenLikBadSolvePenalty = -1e10`, and
  `fk.elem(find(fk >= 1.0e99)).fill(...)` (`src/saem.cpp:3008, 4463-4464`).  `calc.2LL`
  does the raw assignment with no clamp at all:
  `.dyf[.isLL] <- fsave[.isLL]` (line 169).  One 1e99 sentinel or extreme-node
  log-density enters `lQ` ungated, and `-2*sum(lQ)` at that magnitude is the reported
  `1.28e13`.
- **Grid inflation.**  `nx = nnodes.gq^nphi1`, and IOV adds K columns to `nphi1`; corner
  nodes at `cond.mean +- nsd.gq*condsd` -- where two of those "posterior SDs" belong to
  *fixed unit-variance* etas -- push `exp(tcl + w*phi)` far outside the supported region.

Instrument rather than guess: print per-subject `ly`, `lphi1`, `rowSums(log(b))`,
`log(det(Omega))`, and the central-node `fsave[.isLL]`.  Also check the `DV`/`occ` input
wiring, since with IOV `inPars` becomes `c("occ","DV")` -- two entries for the first time
(`R/saemRxUiGet.R:31-52`; `R/saem_fit.R:284-293, 397, 419-441` including the
`.dvCol != 6L` layout assertion; `src/saem.cpp:5173-5180`).

Fix the root cause; add `tests/testthat/test-iov-loglik-1000.R` asserting the `ll()`
twin's `objf` is finite, order-100, and within tolerance of the `add()` twin's.  This
test must keep passing on both `iovMethod` paths.

### Phase 2 -- spike + control surface
Spike the two-etas-on-one-theta mu-ref question above.  Add
`saemControl(iovMethod = c("twoLevel","theta"))` (`R/saemControl.R`), convert
`attr(nlmixr2Est.saem, "iov")` to the `function(control)` form, make
`assertRxUiIovNoCor` conditional, and record the resolved method on the fit so tests can
assert which path ran.  Tests here are ui/ini assertions -- no fit needed.

### Phase 3 -- the SAEM-side additive expansion
Inside the SAEM model build (`R/saemRxUiGetModel.R` / a new `.saemIovExpand()` in
`R/iov.R`, *not* the generic hook), expand each IOV parameter into K **estimated** etas
entering purely additively:

```r
iov.cl <- rx.iov.cl.1*(occ == 1) + rx.iov.cl.2*(occ == 2)
```

dropping the parameter's own IIV eta (absorbed) and attaching the metadata the kernel
needs: `iovGroup` (per phi1 column: base-parameter id), `iovOcc` (occasion index, 0 for
an A column), `iovK`, per-subject `iovKi`, and the initial `Omega`/`Psi` blocks.

### Phase 4 -- `.configsaem` wiring (`R/saem_fit.R`)
Build the initial `Gamma2_phi1` as the compound-symmetry matrix above; teach `covstruct`
the CS pattern; add the `lambdaGroup` generalization of `LCOV1`/`jcov1` and the explicit
`CGamma21` sum; pass `iovGroup`/`iovOcc`/`iovK`/`iovKi` through `ARGS`.  Add the IOV
variances to `parHistOmegaKeep`/`parHistNames` (`R/saemRxUiGet.R:617-684`) -- today they
are excluded, so IOV never appears in the iteration history.

### Phase 5 -- C++ M-step (`src/saem.cpp`)
Replace the unstructured `G1` moment update (lines 3276-3386) with `updateOmegaPsi()`
when IOV metadata is present: per-subject `V_i`/`m_i`, the block accumulators above, then
rebuild `Gamma2_phi1` as CS.  Preserve the existing SA floor (`coef_sa`, `nb_sa`),
`covstruct1` masking, `Gamma2_phi1fixedIx` pinning (so `iov.ka ~ fix(0.05) | occ` still
works) and the `minv` floor, applied to both blocks.  Generalized GLS solve for
`Plambda1`.  Any new random draws go through `rxNormEng`/`rxUnifEng` with a mixed seed
per CLAUDE.md.  If an `[[Rcpp::export]]` arity changes, run `Rcpp::compileAttributes(".")`
**and** hand-edit `src/init.c`; `rm -f src/*.o src/*.so` after header edits.

### Phase 6 -- FIM / standard errors
`d1_loggamma2_phi1` (`src/saem.cpp:3129`) scores only `diag(Gamma2_phi1)`.  Add score
entries for `log diag(Omega)` and `log diag(Psi)`, extend `nb_param`, and update
`.saemFimToCov`'s reordering + delta method (`R/saem.R:673-787`) and `calc.COV`'s `blocB`
(`R/saem_fit_aux.R:459+`).  The paper only claims a *linearized* FIM, so `linFim` (the
saem default) is the reference and `sa`/`fim` are checked against it.

### Phase 7 -- a feasible objective for IOV models
`nnodes.gq^nphi1` is already unusable at `nphi1 = 5` and worsens with K.  Exploit the
conditional independence the new parameterization exposes -- given `b_i` the K occasions
are independent -- or fall back to the Laplace route (`nnodes.gq = 1`).  Add a node-budget
guard that loudly downgrades instead of silently running hundreds of thousands of solves.

### Phase 8 -- output and hand-off
`$omega`'s `$occ` block now reads straight off `Psi`; the `nlmixr2iov*Cv/Sd`
back-transforms become variance-component reporting; `$iov` deviations become the
posterior `c_ik = phi_ik - m_i`.  Fix the three reporting inconsistencies the legacy path
carries, since the new path has no excuse for them: `parFixedDf` and the printed
`parFixed` put the IOV value in *different columns*; an `SE` survives on a row whose
`Estimate` is `NA` and is on the SD scale without saying so; and `fit$cov` still carries a
row/column literally named `iov.cl` for a parameter that is an eta by the time the user
sees it.  Rebuild `$etaMat` (`R/nmObjGetEtaMat.R`) so other methods still inherit values,
keeping the round-trip asserted at `test-iov.R:46-53` working for both paths.

### Phase 9 -- reproduce the paper (`design/iov-panhard/`)
Following the `design/qrpem` + `plans/vae-reference-reproduction.md` convention.
Paper Section 4: one-compartment first-order absorption,
`f = D*Ka/(V*Ka - Cl) * (exp(-Cl/V*t) - exp(-Ka*t))`, `phi = (log V, log Ka, log AUC)`
with `AUC = D/Cl`; `D = 4`; times `0.25, 0.5, 1, 2, 3.5, 5, 7, 9, 12, 24`; `K = 2`
periods, `J = 10` per period; `n = 24` and `n = 40`; **1000 replications each**.
True values `mu = (-0.73, 0.39, 4.61)`, `Omega = diag(0.01, 0.04, 0.04)`,
`Psi = diag(0.0025, 0.01, 0.01)`, `sigma^2 = 0.01`.  SAEM tuning: 500 iterations,
`gamma_l = 1` for the first 200 then `1/(l-200)` -- i.e. `saemControl(nBurn = 200,
nEm = 300)`, already the package default.
*Confirm the residual model (`g = 1` vs `g = sigma*f`) and the H1 `beta` values from the
paper PDF before coding -- the HTML rendering does not state them unambiguously.*

Deliverables:
- `design/iov-panhard/simulate.R`, `run.R`, `report.R`, `results.rds`;
- relative bias % and RMSE % for `V, Ka, AUC, beta_V, beta_Ka, beta_AUC, omega^2 x3,
  psi^2 x3, sigma^2`, plus the Wald type-I error, tabulated **next to the paper's
  Table 1** in `plans/saem-iov-two-level.md`.  Targets: bias/RMSE of the same order as
  the paper's, type-I error near the reported 5.6-6.0%;
- a reduced (~25-replication, loose-tolerance) version as
  `tests/testthat/test-saem-iov-panhard.R`, added to a `.slowBatches` batch in
  `tests/testthat.R` -- fit-heavy, so it must stay out of the essential push/PR subset,
  and must **not** use `skip_on_ci()` (the weekly runner sets `CI=true` too).

### Phase 10 -- close out
`NEWS.md` under the existing `## New features` / `## Bug fixes` headings, including an
explicit note that SAEM IOV estimates change because `twoLevel` is now the default.
Roxygen for `iovMethod` (`\donttest`, American spelling, no Unicode, never "users").
`habit-hooks`; the `codefactor` skill for the real PR lint gate; warm the fixture cache
with one serial run, then the full suite.

## Verification

1. **Unit** -- `updateOmegaPsi()` reduces to the paper's formulas: with A empty, balanced
   `K_i == K` and `beta == 0`, compare the C++ update against a small R reimplementation
   of `V`, `m`, `Omega`, `Psi` on a fixed `phi` sample.
2. **Mechanism, not just equivalence** -- assert the new path actually ran: `nphi0` no
   longer holds an IOV magnitude, the fit records `iovMethod == "twoLevel"`,
   `Gamma2_phi1` is compound-symmetric to machine precision, and `Psi` appears in
   `parHistData`.
3. **The regression that motivated this** -- the `theo_md` + `occ` fit must no longer
   return an IOV variance an order of magnitude below FOCEi's (0.000302 vs 0.003348);
   assert agreement with FOCEi and with the simulated truth, and that it does not shrink
   as iterations increase.
4. **#1000** -- `test-iov-loglik-1000.R`, passing on both `iovMethod` paths.
5. **Legacy path intact** -- `test-iov.R`, `test-iov-template-copy.R`,
   `test-iov-zero-eta.R`, `test-vae-iov.R` all still pass, with `iovMethod = "theta"`
   exercised explicitly and the automatic fallback (two occasion variables) covered.
   Note the current SAEM IOV tests assert only `inherits(fit, "nlmixr2FitCore")` -- there
   is no numeric IOV assertion anywhere today, and there needs to be.
6. **Determinism** -- a repeat fit at a fixed `seed=` reproduces bit-for-bit.
7. **Paper** -- the Phase 9 table.
8. **Full suite** -- serial warm run for the fixture cache, then everything; failures
   reported with their output rather than summarized.

---

# Progress and decisions taken during implementation

## Phase 1 -- #1000 (DONE, commit "fix(saem): stop the IOV magnitude running away")

The issue's premise was wrong: the blowup is not in `calc.2LL`.  `refinePhi0Lik`
froze the ODE for every `distribution == 4` fit on the premise that a
general-likelihood phi0 is always a likelihood SD the solve never sees.  With
the solve frozen the objective was measurably constant in the IOV magnitude
(`sum = -2523.23` for every candidate from 382 to 1000), and `localTrust` was
off for `distribution == 4`, so bobyqa got `(0, Inf)` and walked to its bound.
The magnitude then grew geometrically from `niter_phi0` onward
(0.0608 -> 6.47 -> 42 -> ... -> 1.24e17) and the objf followed to 3.6e33.
Fixed by measuring whether the freeze is valid (`phi0NeedsLiveSolve()`) and by
clamping an ODE-driven general-lik phi0 to the existing absolute trust region.
`objf` 3.6e33 -> 683 on the issue's model.

## Spike result (changes the plan)

`rxode2` REFUSES `theta + eta1 + eta2` in a user model:
`mu-ref err: currently do not theta + eta1 + eta2`.  But `iov.x ~ v | occ`
parses fine and already yields a two-level `$omega` (`$id` + `$occ`), and
`theta + eta + <variable>` parses.  So the collapsed representation must be
built in the GENERATED saem model, and the per-occasion combination has to sit
on a line of its own.

## Phase 3 -- the (b, c) parameterization first (DONE)

The plan's collapsed `phi_ik` layout needs one theta to feed K phi columns,
which the `LCOV`/`MCOV` design cannot express today (each live `MCOV1` entry is
its own lambda) and which ripples into the FIM, `parHist` and the covariance
splice.  The non-collapsed `(b_i, c_ik)` form of the SAME model needs none of
that -- it drops straight into the existing machinery:

```r
rx.iov.cl <- (occ == 1)*rx.iov.cl.1 + (occ == 2)*rx.iov.cl.2
cl <- exp(tcl + eta.cl + rx.iov.cl)
```

gives `nonMuEtas = rx.iov.cl.1, rx.iov.cl.2` (phi means pinned at 0), five phi1
columns, no phi0 column at all, and the occasion variances as ordinary omega
entries.  The single new constraint is that the K per-occasion variances
estimate ONE `Psi`; under that equality the maximizer of the complete-data
likelihood is the plain mean of the K per-occasion moments, so
`poolOmegaGroups()` (src/saem.cpp) stays a closed-form M-step.

Measured on `theo_md` + `occ` (nBurn=200, nEm=300, seed 42):

| path | IOV variance | objf |
|---|---|---|
| `iovMethod="theta"` (shared rewrite) | 0.000302 | 351.19 |
| `iovMethod="twoLevel"` (pooled) | **0.000999** | 346.03 |
| `focei` | 0.002965 | 350.68 |

The two occasion variances come out exactly equal, so the constraint is real.

**Phase 4 therefore becomes**: collapse those columns to the paper's `phi_ik`
(one theta feeding K columns via a new `lambdaGroup`, compound-symmetry
`Gamma`, and the Rao-Blackwellized `V(theta)`/`m(phi_i, theta)` M-step), with
Phase 3 kept as the fallback and as a cross-check for the paper reproduction.

## Blocker found while reproducing the paper: the residual transform cache

The first 20-replicate pilot of the paper's design came back with `add.sd`
biased +400% (0.5 against a truth of 0.1) and `Psi` biased about -55%.  Four
experiments isolated it, and none of it is about IOV:

1. **FOCEi on the same replicate recovers the truth** (`add.sd` 0.112,
   `Omega` 0.0099/0.0348/0.0359, `Psi` 0.0029/0.0272/0.0102), so the simulator
   and the model are right.  The simulator's standardized residuals have SD
   0.0975 against a nominal 0.1.
2. **Fixing the residual at the truth** brings `saem`'s variance components
   back: `iovMethod="twoLevel"` then gives `Psi` 0.0014/0.0276/0.0098 against
   FOCEi's 0.0029/0.0272/0.0102, while `iovMethod="theta"` gives
   0.00067/0.00047/0.0055 -- so the two-level path tracks FOCEi and the shared
   rewrite does not.
3. **A model with no IOV at all** shows the same residual bias (`add.sd` 0.31
   against 0.1 while FOCEi gets 0.112), so it is not the occasion handling.
4. **The same fit, same data, same `seed=`, three times in one session gives
   three different answers**: `add.sd` 0.0859, 0.0877, 0.1807.

Cause: `ensureSaemFixedTransformCache()` (`src/saem.cpp`) keys the residual
step's transformed-prediction cache on the ADDRESS of `ysb`/`fsb`.  Those are
per-iteration locals (`vec ysb, fsb;`), freed and reallocated every M-step, so
the allocator routinely returns the same address holding different data.  Every
other field in the guard (`_saemLen`, `_saemYj`, `_saemPropT`, `_saemLambda`,
`_saemLow`, `_saemHi`) is constant across iterations, so nothing invalidated the
cache and the optimizer scored later iterations against the FIRST iteration's
predictions -- which still carry all the between-subject variability as error.
The residual stayed near its warm start and the variance components shrank to
compensate.  Only `obj`/`objC`/`objD` use the cache, i.e. the residual models
that need the optimizer (combined1/combined2, add+pow and their TBS variants);
pure `add()`/`prop()` have closed-form M-steps, which is why the test corpus
never caught it.

Fixed by invalidating on a generation counter bumped wherever the cache's source
data is set.  After the fix the three repeated fits are bit-identical
(`add.sd` 0.0859 every time, `prop.sd` 0.0963, `psi.lKa` 0.02750).

## Pilot reproduction (20 replicates, n = 24), both paths, cache fix applied

Relative bias %, `design/iov-panhard/run.R 24 20 <method>`:

| | true | `iovMethod="theta"` | `iovMethod="twoLevel"` |
|---|---|---|---|
| mu (log V) | -0.73 | 1.09 | 1.17 |
| mu (log Ka) | 0.39 | 3.33 | 2.49 |
| mu (log AUC) | 4.61 | -0.10 | -0.11 |
| omega (log V) | 0.01 | 9.33 | 2.16 |
| omega (log Ka) | 0.04 | 8.02 | -3.40 |
| omega (log AUC) | 0.04 | 1.46 | -0.07 |
| **psi (log V)** | 0.0025 | **-95.28** | **-18.25** |
| **psi (log Ka)** | 0.01 | **-92.68** | **8.09** |
| **psi (log AUC)** | 0.01 | **-38.76** | **-4.64** |
| sigma (add) | 0.1 | 18.99 | -5.57 |
| sigma (prop) | 0.1 | 1.48 | 1.69 |

The two fixes are independent and both needed: the residual-cache fix repairs
`sigma` and, through it, `Omega` on BOTH paths; the two-level parameterization is
what makes `Psi` recoverable at all.  The shared rewrite still loses 39-95% of
the inter-occasion variance even with the cache repaired.

The legacy path is also about twice as slow here (403 s vs 214 s for 20 fits),
because the IOV magnitude goes through the bounded phi0 optimizer every
iteration past `niter_phi0`.

## Phase 8 -- output and hand-off (DONE)

`.uiFinalizeIov()` is shared rather than duplicated: `.saemIovExpandUi()` leaves
the same state behind that the shared rewrite does (`$ui`, `$iovVars`,
`$iovDrop`, `$lines`, `$iovRename`), plus `$iovTwoLevel` to say the variance
comes off the pooled occasion etas instead of a magnitude theta.  Three branches
were needed:

- `getEstimateDf()` reads the variance from the pooled eta rows.  It runs over
  the final frame AND over `iniDf0`, and `iniDf0` is the user's own frame -- it
  never went through the expansion, so its `iov.x ~ v | occ` row is already the
  thing being rebuilt and is left alone (this was a real crash: "replacement has
  0 rows, data has 1").
- `$iov` is not rescaled: the shared rewrite's occasion etas are unit-variance
  and have to be multiplied by the fitted SD, while the two-level ones already
  carry `c_ik`.
- the `$fixef`/`parFixed` blocks that strip the magnitude theta are skipped --
  there is none, and `x[-integer(0)]` empties the vector rather than leaving it
  alone.

Both paths now produce the same shape: `$omega` split into `$id`/`$occ`, the
user's `iov.cl` row restored with `condition = "occ"`, an `$iov` deviation
table, and no `rx.*` leaking into the fit or its `iniDf`.


## Default flip to "twoLevel" (DONE) -- and the two bugs it caught

Flipping the default immediately broke `test-iov.R`, which is exactly what the
flip was for.  Both were real:

1. **Pipeline position.**  The expansion ran at dispatch, after every hook.
   `.preProcessBoundedTransform` runs LAST and rewrites a bounded theta
   (`tcl <- log(c(0, 2.7, 100))`) into `tcl <- 4.6 - exp(rxBoundedTr.tcl)`,
   which takes `tcl` out of `muRefDataFrame`.  So `.saemIovInfo()` said "in
   scope" at hook time (letting the shared rewrite stand down) and "needs a
   mu-referenced parameter" at dispatch -- and a decline at dispatch was
   silently ignored, leaving the occasion term unhandled ("subscript out of
   bounds").  The expansion is now a pre-processing hook,
   `.uiApplyIovTwoLevel`, which sorts immediately after `.uiApplyIov` and so
   runs at the same point in the pipeline as the rewrite it replaces.
2. **Missing per-record alias.**  The shared rewrite emits a second line per
   parameter, `iov.x.rx <- rx.iov.x`, which `.uiFinalizeIov()`'s data-frame
   branch renames back to `iov.x` -- that is how the fit's data frame carries
   the realized per-record IOV value.  The two-level expansion did not, so
   `any(names(fit) == "iov.cl")` was FALSE.

All eight IOV test files pass under the new default (156 assertions).
## Pre-existing failures in the essential suite (NOT from this branch)

Running the 191-file essential subset surfaced two failing files.  Both were
checked against a clean `origin/main` worktree and fail there identically, so
neither comes from this work -- they arrive with the merged upstream PR #1032
(`lincmt-origin-sens-1119`):

Running the SAME sequential suite on a clean `origin/main` worktree (188 files
there; this branch adds 3) gives:

| file | origin/main | this branch |
|---|---|---|
| `test-focei-char.R` | error=1, passed=0 | error=1, passed=0 |
| `test-focei-lincmt-alag-sens.R` | failed=8, passed=26 | failed=8, passed=26 |
| `test-simModelCache.R` | failed=2, passed=29 | failed=2, passed=29 |
| **total** | failed 10, error 1, passed 4937 | failed 10, error 1, passed 5004 |

Same failures, +67 passing assertions from the new tests: no regression.

The alag failures are analytic-vs-finite-difference gradient mismatches in
"a MIXED-ROUTE regimen gets the right alag()/f() eta gradient
(rxode2/rxode2#1119)" -- that PR's own subject, with no saem or IOV involvement.

`test-simModelCache.R` needed the like-for-like run to attribute: it passes
31/31 in isolation on BOTH branches and only fails inside a long sequential
session, so an isolated check would have wrongly blamed this branch.

Harness note: several test files `assign(..., globalenv())`, which clobbers a
loop variable of the same name in a driver script run through `Rscript`.  That
killed the suite driver at the same file twice before the loop was moved inside
a function.


# Phase 4 -- the collapsed (Panhard & Samson) sampler, opt in

`saemControl(iovMethod = "collapsed")`.  Where `"twoLevel"` samples `b_i` and
`c_ik` as separate columns, this carries `phi_ik = mu + b_i + c_ik` in ONE
column per occasion, which is the paper's own parameterization:

```r
rx.cl.1 <- exp(rx.tcl.1 + rx.eta.cl.1)
rx.cl.2 <- exp(rx.tcl.2 + rx.eta.cl.2)
cl <- rx.cl.1*(occ == 1) + rx.cl.2*(occ == 2)
```

**What the paper's algorithm actually buys here.**  Its two advantages are the
collapsed sampler and Rao-Blackwellizing `b_i` out of the M-step -- but the
second exists only BECAUSE the paper samples `phi_ik` and must then infer `b_i`.
`"twoLevel"` samples `b_i` and `c_ik` directly, so its Omega/Psi updates are
already the exact complete-data maximizers and there is nothing to
Rao-Blackwellize.  Phase 4's real deliverable is therefore the mixing, which is
an efficiency claim to be measured rather than asserted.

Design notes, each of which cost a debugging cycle:

- **One mu-reference per LINE.**  rxode2 detects only the first additive
  `theta + eta` group in a line, so both occasions on one line leaves the second
  eta non-mu-referenced.  Split across lines, every occasion column is an
  ordinary mu-referenced parameter -- its own phi column, its own omega entry,
  mean estimated rather than pinned -- and the whole existing machinery applies.
- **The equality constraints are exact, not projections.**  The group's block is
  compound-symmetric, so `1` is an eigenvector of it and `Gamma^-1 1` is
  proportional to `1`; the constrained GLS for the shared `mu` is therefore the
  equal-weight average of the group's unconstrained solutions.  Verified: the
  two thetas come out bit-identical (1.0439118 each), and the block exactly CS
  (diagonals 0.05152835, off-diagonals 0.04412665).
- **`rxUiDecompress` is not a copy.**  The two-level path escapes this only
  because its `rxRename` returns a fresh ui; the collapsed path has no rename,
  so mutating the decompressed ui also mutated the object stashed as "the
  original".  The ini and model are snapshotted as VALUES before mutating.
- **`.getUiFunFromIniAndModel()` returns a model FUNCTION, not a ui.**  The
  shared finalizer only ends up with a ui because its `rxRename` call evaluates
  one.
- **`.uiIovEnv` is process-global**, and each path installs its own post-fit
  restoration, so every path has to clear the others' markers or a later fit
  runs through the wrong finalizer.

Known limitation: the Gaussian-quadrature objective is not comparable across the
two parameterizations when `Psi << Omega` -- the CS prior is then nearly
degenerate along `b` and the axis-aligned grid covers it badly (measured: objf
478.9 collapsed against 351.9 two-level on theo_md, from near-identical
parameter estimates).  Compare on estimates, not `objf`.

## Is the FOCEi objective under IOV a valid yardstick?

The `iovMethod` comparison used `est="focei"` on each fit's restored ui as a common
scale.  That is worth questioning, because FOCEi always takes the legacy `.uiApplyIov()`
rewrite: the occasion effect becomes a magnitude theta times etas whose variance is
FIXED at 1, so `Psi` never enters `log|Omega|` or the eta penalty -- it acts only by
scaling the prediction, and FOCEi then re-optimizes those unit-variance etas freely.  If
that left the objective flat in `Psi`, the comparison would have been meaningless.

Profiling the objective over `Psi` alone, on data simulated with a known
`psi(log Ka) = 0.01` (`panhardSim(60, 4242)`, everything else at the truth):

| `Psi` | FOCEi objf |
|---|---|
| 0.00001 | 236.6897 |
| 0.00010 | 236.2073 |
| 0.00100 | 231.8497 |
| 0.00300 | 224.6282 |
| **0.01000** | **213.8754** |
| 0.03000 | 216.6466 |
| 0.10000 | 252.7838 |

The minimum sits at the true value and the curvature is large (23 units out to
`Psi = 1e-5`).  The unit-variance-eta reparameterization is a change of variables that
the marginal likelihood is invariant to; only the Laplace approximation to it differs,
and not enough to move the optimum off the truth here.  So no IOV-specific adjustment of
the ui is needed for the FOCEi objective, and the comparison is on a sound common scale.

Kept as a regression test ("the FOCEi objective under IOV is minimized at the true Psi",
`test-saem-iov-twolevel-fit.R`).

**But it does invalidate the reading of the `theo_md` comparison.**  Given this much
sensitivity to `Psi` on data that has real IOV, the three `iovMethod` settings landing
within 0.2 objective units on `theo_md` means `theo_md` carries almost no IOV
information -- its occasions are just `TIME >= 144` of one continuous profile.  That
comparison cannot discriminate between the methods, and the earlier claim that
`twoLevel` "picks the best etas/iovs" is not supported by it.  The default has to be
settled on the Phase 9 reproduction, which has known simulated IOV.

## Phase 9 result -- the paper is reproduced

1000 replicates per arm, both complete, `iovMethod = "twoLevel"`, SAEM tuning as
in the paper (`nBurn = 200`, `nEm = 300`).  Against Panhard & Samson (2009)
Table 1, extended-SAEM columns (arXiv:0803.4437v1):

Relative bias (%):

| | ours n=24 | ours n=40 | paper n=24 | paper n=40 |
|---|---|---|---|---|
| mu log V | 0.01 | -0.20 | 0.01 | -0.06 |
| mu log Ka | -0.12 | -0.10 | 0.48 | 0.02 |
| mu log AUC | -0.03 | -0.02 | -0.08 | -0.11 |
| omega^2 log V | -7.46 | -4.58 | -5.13 | -3.45 |
| omega^2 log Ka | -4.98 | -3.06 | -3.99 | -3.23 |
| omega^2 log AUC | -4.25 | -2.37 | -4.88 | -1.51 |
| psi^2 log V | **-2.19** | **-3.44** | -8.67 | -5.93 |
| psi^2 log Ka | **-6.83** | **-4.22** | -10.94 | -7.06 |
| psi^2 log AUC | **-3.60** | **-3.63** | -5.37 | -4.92 |

Relative RMSE (%):

| | ours n=24 | ours n=40 | paper n=24 | paper n=40 |
|---|---|---|---|---|
| mu log V | 3.25 | 2.56 | 3.9 | 2.91 |
| mu log Ka | 12.41 | 9.61 | 14.4 | 10.79 |
| mu log AUC | 0.93 | 0.73 | 1.0 | 0.79 |
| omega^2 log V | 40.86 | 32.92 | 38.7 | 30.30 |
| omega^2 log Ka | 42.07 | 33.13 | 42.4 | 33.49 |
| omega^2 log AUC | 34.46 | 26.59 | 34.5 | 27.41 |
| psi^2 log V | 63.70 | 55.88 | 69.4 | 58.78 |
| psi^2 log Ka | 74.86 | 63.93 | 73.5 | 62.00 |
| psi^2 log AUC | 47.32 | 36.11 | 43.6 | 33.31 |

Every row moves the right way from n=24 to n=40 -- bias and RMSE both shrink --
and every figure lands within a few points of the paper's.  The mu RMSEs are
consistently a little BETTER than the published ones, the omega^2 RMSEs are
essentially identical, and on `psi^2` -- the inter-occasion variance this whole
change is about -- our estimates are markedly LESS downward-biased than the
paper's own implementation.

Two rows of the paper's table are deliberately not reproduced, and the harness
says so rather than inventing them:

- `beta_V`, `beta_ka`, `beta_AUC`: the arXiv HTML never states the beta values
  used to generate the data.  The harness simulates under H0 (`beta = 0`), where
  a relative bias in beta is undefined.  (The paper's reported beta biases are
  finite, so its simulation clearly used nonzero betas.)
- `sigma^2`: the paper ties a single sigma to `g = 1 + f`, while this harness
  estimates `add.sd` and `prop.sd` separately with `combined1()`.  Estimating two
  parameters where the paper estimates one makes each individually noisier, so
  the sigma rows are not comparable.  (Ours: bias +2.93/-0.63% at n=24,
  +1.94/-0.30% at n=40.)

`design/iov-panhard/report.R` carries the paper's numbers inline and prints the
comparison, so the claim is reproducible rather than transcribed into prose.

## The default: twoLevel, decided on the reproduction corpus

`theo_md` could not settle this (see the FOCEi section above -- its "occasions"
are `TIME >= 144` of one continuous profile, so it carries almost no IOV
information and all three settings scored within 0.2 objective units).  The
1000-replicate corpus can.  Both arms at n = 24, 1000/1000 usable:

| | twoLevel bias | collapsed bias | twoLevel RMSE | collapsed RMSE |
|---|---|---|---|---|
| mu log V | 0.01 | 0.03 | 3.25 | 3.24 |
| mu log Ka | -0.12 | 0.07 | 12.41 | 12.36 |
| mu log AUC | -0.03 | -0.01 | 0.93 | 0.93 |
| omega^2 log V | -7.46 | -7.29 | 40.86 | 40.21 |
| omega^2 log Ka | -4.98 | -4.00 | 42.07 | 41.67 |
| omega^2 log AUC | -4.25 | -4.63 | 34.46 | 34.36 |
| **psi^2 log V** | **-2.19** | 6.40 | **63.70** | 69.76 |
| **psi^2 log Ka** | **-6.83** | **-16.08** | **74.86** | 76.04 |
| **psi^2 log AUC** | -3.60 | **0.09** | 47.32 | **46.17** |

`mu` and `omega^2` are indistinguishable -- every figure agrees to well inside
Monte Carlo error.  The two separate only on `psi^2`, the inter-occasion
variance, and there `twoLevel` wins 2 of 3 components on both bias and RMSE.
Where `collapsed` loses it loses badly: -16.08% on log Ka is more than double
`twoLevel`'s -6.83% and worse than the -10.94% Panhard and Samson report for
their own implementation.  Mean absolute `psi^2` bias is 4.21% for `twoLevel`
against 7.52% for `collapsed`.

So `twoLevel` stays the default and `collapsed` stays opt-in.

`collapsed` was not slower -- 8308 s against 10938 s for the same 1000
replicates -- and in fact shared the machine with the n = 40 arm for part of its
run while `twoLevel` did not, so its real margin is probably wider.  That is
worth knowing but does not change the default: the setting exists to estimate an
inter-occasion variance, and it is the less accurate of the two at doing so.
These arms are single-seed-stream runs of one design; the ordering on `psi^2` is
consistent and large, but it is not a claim about every model shape.

## The legacy path, measured at 1000 replicates

The default flip now rests on matched arms rather than a pilot on one side.
n = 24, 1000/1000 usable in each:

| relative bias % | `theta` (legacy) | `twoLevel` |
|---|---|---|
| mu log V / Ka / AUC | -0.03 / 0.57 / -0.02 | 0.01 / -0.12 / -0.03 |
| omega^2 log V / Ka / AUC | 2.34 / 2.17 / -0.65 | -7.46 / -4.98 / -4.25 |
| **psi^2 log V / Ka / AUC** | **-94.70 / -89.00 / -43.56** | -2.19 / -6.83 / -3.60 |
| sigma add / prop | **25.01** / -0.02 | 2.93 / -0.63 |

| relative RMSE % | `theta` | `twoLevel` |
|---|---|---|
| psi^2 log V / Ka / AUC | 94.98 / 91.58 / 61.53 | 63.70 / 74.86 / 47.32 |
| sigma add | 47.33 | 38.40 |

The legacy path collapses the inter-occasion variance to near zero, and the
variance it fails to capture reappears as a 25% inflation of the additive
residual error -- the two findings corroborate each other, which is stronger
evidence than the psi^2 bias alone.  `mu` is unaffected in both (the fixed
effects are recovered either way); `omega^2` is slightly LESS biased under the
legacy path, which is consistent with the same leak, since a fit that cannot put
occasion variability in `psi` has to put it somewhere.

The 20-replicate pilot that previously backed this claim said -95/-93/-39, so it
was directionally right, but it is now superseded.
