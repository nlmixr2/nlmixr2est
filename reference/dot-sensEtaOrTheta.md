# Calculate d(state)/d(eta) or d(state)/d(theta) sensitivities

Calculate d(state)/d(eta) or d(state)/d(theta) sensitivities

## Usage

``` r
.sensEtaOrTheta(
  s,
  theta = FALSE,
  extraThetaVars = NULL,
  rxui = NULL,
  matExpForcing = TRUE
)
```

## Arguments

- s:

  symengine environment (from \`.loadSymengine()\`)

- theta:

  when \`TRUE\` calculate the sensitivities with respect to
  \`THETA\[#\]\`; otherwise with respect to \`ETA\[#\]\`

- rxui:

  rxode2 UI object backing \`s\`; only needed for a matExp() model
  (\`.sensMatExpNative()\` needs the raw natural-name model text),
  \`NULL\` otherwise

- matExpForcing:

  when \`FALSE\`, a matExp() model with an \`indLin()\` forcing term
  (Michaelis-Menten) falls back to the ordinary variational-ODE flatten
  instead of \`.sensMatExpNative()\`, even when \`rxui\` is supplied.
  nlm's population log-likelihood gradient (\`rxUiGet.nlmHdTheta\`) does
  not yet agree with a finite-difference reference for a forcing model's
  explicit (non-state-mediated) theta sensitivity – FOCEi's EBE gradient
  is unaffected and keeps the native path (issue \#860; the nlm gap is
  tracked separately, matching how \`.rxKeepMatExpNative()\` already
  excludes SAEM from forcing models for its own reason).

## Value

the symengine environment \`s\` augmented with the sensitivity equations
(\`..sens\`, \`..ddt\`, \`..stateInfo\`, ...)

## Author

Matthew L. Fidler
