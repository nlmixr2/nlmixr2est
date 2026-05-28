# Adding an estimation routine to nlmixr2est

## Overview

`nlmixr2est` adds estimation methods by S3 dispatch on
[`nlmixr2Est()`](https://nlmixr2.github.io/nlmixr2est/reference/nlmixr2Est.md).
At a high level,
[`nlmixr2()`](https://nlmixr2.github.io/nlmixr2est/reference/nlmixr2.md):

1.  Builds an environment containing `ui`, `data`, `control`, and
    `table`.

2.  Normalizes the control object with
    `getValidNlmixrControl(control, est)`.

3.  Runs the registered preprocessing hooks.

4.  Sets `class(env) <- c(est, "nlmixr2Est")`.

5.  Calls `nlmixr2Est(env, ...)`, which dispatches to
    `nlmixr2Est.<est>()`.

That means a new routine is usually added by providing one new S3
estimation method plus the control helpers that let the rest of the
package treat the routine like the built-in methods.

The most compact built-in example is the `newuoa` implementation in
`R/newuoa.R`. Wrappers such as `foce`, `laplace`, and `agq` are useful
when your new routine is really a specialization of an existing family.

## The minimum pieces

In practice, a new routine named `myest` needs these pieces.

### 1. A control constructor

The control constructor should return a named list with class
`"myestControl"`. Reuse `sharedControl()` where possible so standard
options such as `literalFix`, `boundedTransform`, and `rxControl` work
the same way as they do for the existing methods.

``` r

myestControl <- function(...) {
  .ctl <- sharedControl(...) # like from foceiControl(...)
  class(.ctl) <- "myestControl"
  .ctl
}
```

### 2. A control validator

[`nlmixr2()`](https://nlmixr2.github.io/nlmixr2est/reference/nlmixr2.md)
never hands your method the raw `control=` input. Instead it calls
`getValidNlmixrCtl.<est>()`, which should:

- create defaults when `control` is `NULL`

- coerce plain lists with `do.call(myestControl, .ctl)`

- reject incompatible classes and fall back to defaults

- return a fresh `myestControl` object

``` r

getValidNlmixrCtl.myest <- function(control) {
  .ctl <- control[[1]]
  if (is.null(.ctl)) .ctl <- myestControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list")) {
    .ctl <- do.call("myestControl", .ctl)
  }
  if (!inherits(.ctl, "myestControl")) {
    .minfo("invalid control for `est=\"myest\"`, using default")
    .ctl <- myestControl()
  } else {
    .ctl <- do.call(myestControl, .ctl)
  }
  .ctl
}
```

### 3. A way to store the control on the fit object

Most estimation methods implement both
`nmObjHandleControlObject.<controlClass>` and `nmObjGetControl.<est>()`.
This keeps the original method-specific control available after fitting.

``` r

nmObjHandleControlObject.myestControl <- function(control, env) {
  assign("myestControl", control, envir = env)
}

nmObjGetControl.myest <- function(x, ...) {
  .env <- x[[1]]
  if (exists("myestControl", .env)) {
    .control <- get("myestControl", .env)
    if (inherits(.control, "myestControl")) return(.control)
  }
  if (exists("control", .env)) {
    .control <- get("control", .env)
    if (inherits(.control, "myestControl")) return(.control)
  }
  stop("cannot find myest related control object", call. = FALSE)
}
```

All nlmixr2 fit output is adapted through a FOCEI-like post-processing
path; if you want a `nlmixr2` fit object you also want
`nmObjGetFoceiControl.myest()` so downstream code can recover a FOCEI
compatible control object.

### 4. The estimation method itself

Your method receives the preprocessed environment. At minimum it should:

- inspect `env$ui`, `env$data`, `env$control`, and `env$table`

- run any model checks specific to the routine

- fit the model

- return a regular `nlmixr2` fit object or another supported object

``` r

nlmixr2Est.myest <- function(env, ...) {
  .ui <- env$ui
  rxode2::assertRxUiPopulationOnly(
    .ui,
    " for the estimation routine 'myest'",
    .var.name = .ui$modelName
  )

  .control <- env$control
  assign("control", .control, envir = .ui)
  on.exit({
    if (exists("control", envir = .ui)) rm("control", envir = .ui)
  }, add = TRUE)

  ## fit model here and return a standard nlmixr2 object
}
```

For routines that create a standard fit object, the easiest pattern is
usually:

1.  preprocess data into the same structures used by existing methods

2.  compute final parameter values and any optimizer diagnostics

3.  call
    [`nmObjHandleControlObject()`](https://nlmixr2.github.io/nlmixr2est/reference/nmObjHandleControlObject.md)
    to attach the control

4.  hand the result to one of the existing output constructors such as
    [`nlmixr2CreateOutputFromUi()`](https://nlmixr2.github.io/nlmixr2est/reference/nlmixr2CreateOutputFromUi.md)

Again, `R/newuoa.R`, `R/nlminb.R`, and `R/optim.R` are good templates
for this.

## What is in `env` when your routine runs?

[`nlmixr2Est()`](https://nlmixr2.github.io/nlmixr2est/reference/nlmixr2Est.md)
checks the basic contract before dispatch:

- `env$ui` is an `rxUi`
- `env$data` is a `data.frame`
- `env$control` is already validated for your method
- `env$table` is already normalized as a `tableControl`

By the time `nlmixr2Est.<est>()` runs, preprocessing hooks may already
have rewritten the UI, the dataset, or even the routine name. The method
should use what is in `env` rather than re-reading the original call.

## The method attributes

The attributes attached to `nlmixr2Est.<est>()` are important because
they are read before your method executes. They do not just document
capabilities; they change the preprocessing path that builds the
estimation problem.

### `covPresent`

`covPresent` is consulted by the covariate-presence preprocessing hook.

``` r

attr(nlmixr2Est.myest, "covPresent") <- TRUE
```

When `TRUE`,
[`.nlmixr0preProcessCovariatesPresent()`](https://nlmixr2.github.io/nlmixr2est/reference/dot-nlmixr0preProcessCovariatesPresent.md)
checks that the dataset contains `TIME` plus every model covariate, with
a case-normalization pass that uppercases non-covariate columns. Set
this for essentially all fitting and simulation routines that consume
model covariates directly.

### `unbounded`

`unbounded` controls the bounded-parameter transform hook.

``` r

attr(nlmixr2Est.myest, "unbounded") <- TRUE
```

This attribute can be:

- `TRUE`: your optimizer works on an unconstrained scale, so bounded
  parameters are rewritten to internal unconstrained parameters before
  fitting

- `FALSE`: your method handles bounds natively, so the hook is skipped

- `function(control)`: the choice depends on control settings

This is one of the most important attributes because it changes the
*built* estimation problem. For unbounded methods, the preprocessing
hook rewrites the model, injects back-transforms, and later maps results
back to the natural parameter scale.

### `mu`

`mu` tells `nlmixr2est` whether the routine supports mu-referencing
helpers.

``` r

attr(nlmixr2Est.myest, "mu") <- TRUE
```

This attribute is read by `.isMuMethod()`. If it resolves to `TRUE`,
mu-related preprocessing can replace `mu2`/`mu3`/`mu4` covariate
expressions with derived data columns before estimation. `saem` uses a
simple `TRUE`; `nlme` uses a function so support depends on
`control$muRefCovAlg`.

If your method does **not** understand the transformed mu-reference
path, leave this attribute unset or `FALSE`.

### `iov`

`iov` declares support for inter-occasion variability.

``` r

attr(nlmixr2Est.myest, "iov") <- TRUE
```

This attribute is read through `.isIovMethod()`. Use it when the routine
can work with IOV models after preprocessing. If the routine cannot
support IOV, do not set the attribute; the rest of the pipeline will
then avoid advertising IOV support for that method.

### `random`

`random` is currently used as a capability marker on `rxSolve`-style
routines such as `rxSolve`, `simulate`, `simulation`, and `predict` and
is used in `nlmixr2save` so that when using the `:=` it will
save/restore the seed state for long running simulations that are
cached.

``` r

attr(nlmixr2Est.myest, "random") <- TRUE
```

Unlike `covPresent`, `unbounded`, `mu`, and `iov`, this is not part of
the main estimation preprocessing path.

## How the attributes affect the routine build

The easiest way to think about the attributes is that they steer how
`nlmixr2est` builds the final problem seen by your method:

1.  `covPresent` says whether the incoming data must be checked against
    model covariates.

2.  `mu` says whether the UI and data may be rewritten for mu-referenced
    covariate handling.

3.  `unbounded` says whether bounded parameters should be converted to
    an unconstrained internal optimization problem.

4.  `iov` says whether the method is treated as IOV-capable.

5.  `random` says the method uses random numbers and the output will
    depend on the seed state.

So if `nlmixr2Est.myest()` sees an unexpected UI or dataset, the first
thing to check is usually the attributes on the method and the
preprocessing hooks they activate.

## A compact example

The following is the smallest realistic registration skeleton:

``` r

myestControl <- function(...) {
  .ctl <- sharedControl(...)
  class(.ctl) <- "myestControl"
  .ctl
}

getValidNlmixrCtl.myest <- function(control) {
  .ctl <- control[[1]]
  if (is.null(.ctl)) .ctl <- myestControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list")) {
    .ctl <- do.call("myestControl", .ctl)
  }
  if (!inherits(.ctl, "myestControl")) {
    .ctl <- myestControl()
  } else {
    .ctl <- do.call(myestControl, .ctl)
  }
  .ctl
}

nmObjHandleControlObject.myestControl <- function(control, env) {
  assign("myestControl", control, envir = env)
}

nmObjGetControl.myest <- function(x, ...) {
  .env <- x[[1]]
  get("myestControl", envir = .env)
}

nlmixr2Est.myest <- function(env, ...) {
  .ui <- env$ui
  .control <- env$control
  assign("control", .control, envir = .ui)
  on.exit({
    if (exists("control", envir = .ui)) rm("control", envir = .ui)
  }, add = TRUE)

  ## fit model and return result
}

attr(nlmixr2Est.myest, "covPresent") <- TRUE
attr(nlmixr2Est.myest, "unbounded") <- TRUE
attr(nlmixr2Est.myest, "mu") <- FALSE
attr(nlmixr2Est.myest, "iov") <- FALSE
```

Once that skeleton exists, roxygen will add the S3 registrations to
`NAMESPACE`, and the routine becomes available through:

``` r

fit <- nlmixr2(model, data, est = "myest", control = myestControl())
```

## Suggested development workflow

When adding a new routine, the usual order is:

1.  copy a nearby method family (`newuoa`, `nlminb`, `optim`, or a FOCEI
    wrapper)
2.  add the control constructor and validator
3.  add `nlmixr2Est.<est>()`
4.  set the method attributes deliberately
5.  add `nmObjHandleControlObject.<controlClass>` and
    `nmObjGetControl.<est>()`
6.  run
    [`devtools::document()`](https://devtools.r-lib.org/reference/document.html)
    so S3 registration and the vignette are rebuilt

If the new routine behaves oddly before it reaches the optimizer,
inspect the method attributes first; they often explain why the UI or
data no longer look like the original model call.
