# Set up an nlm-family objective for repeated hook-firing evaluation

Preprocesses the data and LOADS the nlm population (predOnly) problem
into the C++ engine, so that repeated `nlmSolveR(theta)` calls evaluate
the population objective – firing any registered likelihood-contribution
hook (e.g. a plugin's per-observation cotangent capture) – WITHOUT
re-running the optimizer. One compiled setup is reused across
evaluations. Intended for extension packages (e.g. nlmixr2nn) that
optimize an out-of-band parameter block (network weights injected via a
par-loader) and need the exact error-model cotangent from the nlm C++
solve at each iterate. Free the loaded problem with
[`.nlmFreeEnv()`](https://nlmixr2.github.io/nlmixr2est/reference/dot-nlmFreeEnv.md)
when done.

## Usage

``` r
nlmObjectiveSetup(ui, data, control = NULL)
```

## Arguments

- ui:

  rxode2/nlmixr2 model. Uses `ui$control` when present.

- data:

  event data.

- control:

  optional nlm-family control; defaults to
  [`nlmControl()`](https://nlmixr2.github.io/nlmixr2est/reference/nlmControl.md)
  (or `ui$control` if that is an nlm-family control).

## Value

(invisibly) the scaled starting parameter vector to hand to
`nlmSolveR()`; the C++ problem is left loaded.

## Author

Matthew L. Fidler
