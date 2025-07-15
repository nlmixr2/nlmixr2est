# nlmixr2est 4.0.0

- When using a model to start a new focei model, the ETAs from the
  last fit are used as the starting point.  Now you can use
  `foceiControl(etaMat=NA)` to skip this and use `eta=0` for all
  items.

- When using `foceiControl(etaMat=fit)`, this will extract the ETAs
  from a fit for use in the next optimization.

- When using a `foceiControl(etaMat=)` option nlmixr2 no longer only
  evaluates the inner problem with the `etaMat` value.

- Add `mceta` option to `"focei"`.

  - `mceta=-1` is the default; the eta restarts at the best eta from
    the last step to start the inner optimization.
  - `mceta=0` the eta starts at `0` to start the inner optimization.
  - `mceta=1` the eta starts at either `0` or the best `eta`, which
    ever gives the lowest objective function to start the inner
    optimization.
  - `mceta=n` under the assumption of `omega` sample `n-1` `eta`
     values and use the lowest objective function of eta sampled, last
     best eta and eta=0 to start the inner optimization.

- Fix Rstudio print (issue #536)

- Support rxode2's new `+var()` definition in `saem`

- Support literal fixing of residuals (#524).  All methods that
  support a literal fix of residuals have an option `literalFixRes`
  which defaults to `TRUE`.  To get the behavior from older models you can use
  `literalFixRes=FALSE`
