# nlmixr2 4.1.0

- Updated inferring the estimation method from the control
  object. Requires the control object to have a class of length one
  and match the estimation method.  For example `foceiControl()` would
  assume that the estimation method is related to `focei`.

- Changed Rstudio completion to not evaluate (in case it gets turned
  on for data.frames) (See #568)

- Turned on data completion for items like `$fitMergeInner`

- **Breaking change:** Changed the estimation method `posthoc` to add
  tables and calculate the covariance by default.  It is now a method
  with it's own control, `posthocControl()`.  As previously the
  default is not to include the interaction term (but you can turn it
  on with `posthocControl(interaction=TRUE)`).

- Added `foceControl()`, `foControl()` and `foiControl()` for the
  `foce`, `fo` and `foi` methods, respectively.  They try to convert
  the related control structures to the correct control structure for
  the estimation method.

- Added iov support for `focei`,  `foce`, and `saem` (#614)

- Added new estimation method `agq` which uses adaptive Gauss-Hermite
  Quadrature to fit a nonlinear-mixed effect model. In this method,
  you can choose the number of quadrature points to estimate the
  likelihood, with higher numbers giving more accurate likelihoods.
  The AGQ implementation in nlmixr2est allows you to specify the
  number of quadrature points via the `agqControl()` function, and
  supports both single and multiple subject models. This method is
  particularly useful for models where accurate likelihood estimation
  is critical.  Also added a `laplace` method which is the same as
  `agq` with 1 node (and is numerically the same as `focei`, `foce` or
  log-likelihood `focei`/`laplace`, etc)

- Fixed saem mu-reference display by not compressing the internal item
  `saem0`.
