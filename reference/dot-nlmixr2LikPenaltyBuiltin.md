# Built-in likelihood penalty pointer

Built-in likelihood penalty pointer

## Usage

``` r
.nlmixr2LikPenaltyBuiltin(name)
```

## Arguments

- name:

  \`"thetaL2"\` (population: \`par\[1\] \* sum((theta - par\[2\])^2)\`
  on the -2LL scale), \`"etaL2"\` (subject: \`par\[2\] - par\[1\]/2 \*
  sum(eta^2)\` added to each subject's log-likelihood) or \`"omegaL2"\`
  (population: \`par\[1\] \* sum((diag(omega)\[k\] - par\[2\])^2)\` over
  the 1-based etas \`k = par\[-(1:2)\]\`, every eta when none)

## Value

external pointer to the C \`nlmixrLikPenalty\`
