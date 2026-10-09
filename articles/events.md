# nlmixr2est events for loggers

`nlmixr2est` tells the `rxode2` *event bus* when a user-facing piece of
work finishes: a fit, a simulation, or an update to an existing fit.
This vignette is for people who want to **listen** to those events, for
example to keep a log of every fit in a project. The main consumer is
the [nlmixr2log](https://github.com/nlmixr2/nlmixr2log) package; if you
simply want a full log, use that package instead of writing your own
listener.

## The event bus in brief

The bus lives in `rxode2` (see
[`?rxode2::rxEventListen`](https://nlmixr2.github.io/rxode2/reference/rxEventListen.html)
for the full API). A listener is a function `function(event, ...)`: it
receives the event name and the payload fields as named arguments. It is
registered with `rxode2::rxEventListen(id, fun)` and removed with
`rxode2::rxEventUnlisten(id)`.

Events are delivered only from the **outermost** operation. Every
emitting function enters a shared operation scope, and an event is
delivered only when no scope is active. So the many
[`rxSolve()`](https://nlmixr2.github.io/rxode2/reference/rxSolve.html)
calls, refits and covariance steps done inside
[`nlmixr2()`](https://nlmixr2.github.io/nlmixr2est/reference/nlmixr2.md)
are silent, and a bootstrap or a covariate search that runs many fits
reports only its own result (from the package that runs it). Worker
processes started inside a scope stay silent too.

A listener that errors gives a warning; it never stops the fit or the
other listeners. Events triggered while a listener runs are dropped.

nlmixr2est checks whether the installed `rxode2` has the bus; without it
every hook is a no-op.

## Events emitted by nlmixr2est

| Function | Event | Notes |
|----|----|----|
| [`nlmixr2()`](https://nlmixr2.github.io/nlmixr2est/reference/nlmixr2.md) that returns a fit | `fitComplete` | any `est`, including `nlmixr2(fit, est = ...)` |
| [`nlmixr2()`](https://nlmixr2.github.io/nlmixr2est/reference/nlmixr2.md) with `est = "rxSolve"`, `"simulate"` or `"predict"` | `solveComplete` | also reached by `simulate(fit)` (`kind = "rxSolve"`) and `predict(fit)` |
| [`vpcSim()`](https://nlmixr2.github.io/nlmixr2est/reference/vpcSim.md) | `solveComplete` | `kind = "vpcSim"`; the simulations inside are silent |
| [`augPred()`](https://rdrr.io/pkg/nlme/man/augPred.html) on a fit | `solveComplete` | `kind = "augPred"` |
| [`addCwres()`](https://nlmixr2.github.io/nlmixr2est/reference/addCwres.md) | `fitUpdate` | `what = "cwres"` |
| [`addNpde()`](https://nlmixr2.github.io/nlmixr2est/reference/addNpde.md) | `fitUpdate` | `what = "npde"` |
| [`addTable()`](https://nlmixr2.github.io/nlmixr2est/reference/addTable.md) | `fitUpdate` | `what = "table"` |
| [`setOfv()`](https://nlmixr2.github.io/nlmixr2est/reference/setOfv.md) | `fitUpdate` | `what = "ofv"` |
| reading a deferred objective (`fit$objf`, `fit$AIC`, …) | `fitUpdate` | `what = "ofv"`; only the first read, when it computes the value |

`predict(fit, level = "individual")` calls
[`rxode2::rxSolve()`](https://nlmixr2.github.io/rxode2/reference/rxSolve.html)
directly, so its `solveComplete` comes from `rxode2` (with
`kind = "rxSolve"` and the fit as `object`).

### `fitComplete`

| Field | Content |
|----|----|
| `fit` | the new fit (an `nlmixr2FitCore` object, with its final timing) |
| `object` | what was passed to [`nlmixr2()`](https://nlmixr2.github.io/nlmixr2est/reference/nlmixr2.md): a model function, a `rxUi`, or a previous fit |
| `call` | the [`nlmixr2()`](https://nlmixr2.github.io/nlmixr2est/reference/nlmixr2.md) call, normalized (see below) |
| `objName` | the variable name of `object` when it was a plain name (e.g. `"one.cmt"`), otherwise `NULL`; the pipe placeholders `.` and `.x` give `NULL` |
| `source` | always `"fit"` for nlmixr2est |

### `solveComplete`

| Field | Content |
|----|----|
| `result` | the simulation or prediction (an `rxSolve` object for [`nlmixr2()`](https://nlmixr2.github.io/nlmixr2est/reference/nlmixr2.md); the returned value for [`vpcSim()`](https://nlmixr2.github.io/nlmixr2est/reference/vpcSim.md) and [`augPred()`](https://rdrr.io/pkg/nlme/man/augPred.html)) |
| `object` | the fit (or model) it was made from, so the result can be linked to its fit |
| `call` | the call, normalized |
| `kind` | the `est` used by [`nlmixr2()`](https://nlmixr2.github.io/nlmixr2est/reference/nlmixr2.md) (`"rxSolve"`, `"simulate"` or `"predict"`), or `"vpcSim"` or `"augPred"` |

### `fitUpdate`

| Field | Content |
|----|----|
| `fit` | the updated fit |
| `original` | the fit that was passed in |
| `name` | the variable name of the user’s fit when the update is *in place*, otherwise `NULL` |
| `what` | `"cwres"`, `"npde"`, `"table"` or `"ofv"` |
| `inPlace` | `TRUE` when the user’s object now holds the result |

`inPlace` follows what actually happened to the user’s object:

- [`addCwres()`](https://nlmixr2.github.io/nlmixr2est/reference/addCwres.md)
  and
  [`addNpde()`](https://nlmixr2.github.io/nlmixr2est/reference/addNpde.md)
  (with the default `updateObject = TRUE`) rebind the variable that was
  passed in. `inPlace` is `TRUE` only when that rebinding happened,
  which needs a plain variable name, bound in the calling environment to
  the same fit. `addCwres(fits[[1]])`, `addCwres(env$fit)` or
  `updateObject = FALSE` give `inPlace = FALSE` and `name = NULL`.
- [`addTable()`](https://nlmixr2.github.io/nlmixr2est/reference/addTable.md)
  works on a copy by default (`updateObject = FALSE`), giving
  `inPlace = FALSE`. `addTable(updateObject = TRUE)` modifies the fit’s
  shared environment, so `inPlace` is `TRUE` whenever the result shares
  the original’s environment, even for
  `addTable(fits[[1]], updateObject = TRUE)` (where `name` is then
  `NULL` because there is no plain name).
- [`setOfv()`](https://nlmixr2.github.io/nlmixr2est/reference/setOfv.md)
  always changes the fit in place: `inPlace = TRUE`.
- Reading a deferred objective computes it inside the fit:
  `inPlace = TRUE`, `name = NULL`, and `original` is the fit itself.

### The `call` field

`rxode2` normalizes every `call` field before delivery so a recorded
call never embeds a large object: the function name is restored
(e.g. `nlmixr2` or `vpcSim`), and arguments that are values rather than
code (for example a data frame inlined by
[`do.call()`](https://rdrr.io/r/base/do.call.html)) become
`` `<value>` `` (or a single `` `<...>` `` marker when there are many).

## What is not emitted

- **Nested work.** The fits, solves and covariance steps inside
  [`nlmixr2()`](https://nlmixr2.github.io/nlmixr2est/reference/nlmixr2.md);
  the internal FOCEi refit of
  [`addCwres()`](https://nlmixr2.github.io/nlmixr2est/reference/addCwres.md);
  the simulations inside
  [`vpcSim()`](https://nlmixr2.github.io/nlmixr2est/reference/vpcSim.md);
  and every internal re-estimation done while a prior is applied.
- **No-op updates.** An update whose result is identical to the fit
  passed in emits nothing:
  [`addCwres()`](https://nlmixr2.github.io/nlmixr2est/reference/addCwres.md)
  on a fit that already has CWRES,
  [`addNpde()`](https://nlmixr2.github.io/nlmixr2est/reference/addNpde.md)
  on a fit that already has NPDE, or `addTable(updateObject = TRUE)`
  that recomputes the same table. A second read of an objective that is
  already computed emits nothing either.
  ([`setOfv()`](https://nlmixr2.github.io/nlmixr2est/reference/setOfv.md)
  is the exception: it always emits, since it changes the fit in place
  and returns the same object.)
- **Results that are not fits or solves.** `nlmixr2(model)` without data
  (which returns the model) and
  [`nlmixr2()`](https://nlmixr2.github.io/nlmixr2est/reference/nlmixr2.md)
  with no arguments emit nothing.
- **Errors.** A function that stops with an error emits nothing.
- **Accessors.** Asking a fit for data that needs CWRES (through
  [`nmObjGetData()`](https://nlmixr2.github.io/nlmixr2est/reference/nmObjGetData.md))
  adds them silently: it is not a user update.
- A fit passed *inline*, as in
  `nlmixr2(nlmixr2(one.cmt, data, "focei"), est = "posthoc")`, is
  evaluated before the outer call starts, so it is reported on its own,
  followed by the outer fit.

## Listening

A small recording listener keeps the event names and a few payload
fields:

``` r

library(nlmixr2est)
#> Loading required package: nlmixr2data

events <- list()
rxode2::rxEventListen("demo", function(event, ...) {
  p <- list(...)
  fields <- p[c("kind", "what", "objName", "name", "inPlace")]
  fields <- fields[!vapply(fields, is.null, logical(1))]
  events[[length(events) + 1L]] <<- paste0(
    event, "(",
    paste(names(fields), vapply(fields, deparse, character(1)),
          sep = " = ", collapse = ", "),
    ")"
  )
})
show <- function() {
  if (length(events) == 0L) cat("(no events)\n") else writeLines(unlist(events))
  events <<- list()
}
```

A small model fit with `posthoc`:

``` r

one.cmt <- function() {
  ini({
    tka <- log(1.57)
    tcl <- log(2.72)
    tv <- log(31.5)
    eta.ka ~ 0.6
    eta.cl ~ 0.3
    eta.v ~ 0.1
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    linCmt() ~ add(add.sd)
  })
}

fit <- suppressMessages(
  nlmixr2(one.cmt, nlmixr2data::theo_sd, est = "posthoc")
)
```

``` r

show()
#> fitComplete(objName = "one.cmt")
```

Only one `fitComplete` is delivered, although the fit solved the model
many times. Simulations and predictions from the fit each give one
`solveComplete`:

``` r

s <- suppressMessages(simulate(fit))
p <- suppressMessages(predict(fit, nlmixr2data::theo_sd))
v <- suppressMessages(vpcSim(fit, n = 3))
show()
#> solveComplete(kind = "rxSolve")
#> solveComplete(kind = "predict")
#> solveComplete(kind = "vpcSim")
```

Updates report whether the user’s object now holds the result.
[`addTable()`](https://nlmixr2.github.io/nlmixr2est/reference/addTable.md)
works on a copy by default, while `updateObject = TRUE` changes the
fit’s shared environment. The second call below recomputes the table the
fit already has, so its result is identical to its input and nothing is
emitted; the third adds NPDE in place:

``` r

f2 <- suppressMessages(addTable(fit))
fit <- suppressMessages(addTable(fit, updateObject = TRUE))
fit <- suppressMessages(addTable(fit, updateObject = TRUE,
                                 table = tableControl(npde = TRUE)))
show()
#> fitUpdate(what = "table", inPlace = FALSE)
#> fitUpdate(what = "table", name = "fit", inPlace = TRUE)
```

A `saem` fit defers its objective function. The first read computes it
inside the fit and emits a `fitUpdate`; later reads do not:

``` r

fit2 <- suppressMessages(nlmixr2(one.cmt, nlmixr2data::theo_sd, est = "saem",
                                 control = saemControl(print = 0)))
show()
#> fitComplete(objName = "one.cmt")
o1 <- suppressMessages(fit2$objf)
o2 <- fit2$objf
show()
#> fitUpdate(what = "ofv", inPlace = TRUE)
```

[`addCwres()`](https://nlmixr2.github.io/nlmixr2est/reference/addCwres.md)
rebinds the variable it was given. Calling it again has nothing to do,
so the second call emits nothing:

``` r

fit2 <- suppressMessages(addCwres(fit2))
show()
#> fitUpdate(what = "cwres", name = "fit2", inPlace = TRUE)
fit2 <- suppressMessages(addCwres(fit2))
show()
#> (no events)
```

Remove the listener when you are done:

``` r

rxode2::rxEventUnlisten("demo")
```

## Full logging

[nlmixr2log](https://github.com/nlmixr2/nlmixr2log) listens to these
events (and to those of the other nlmixr2 packages) and keeps a project
log of fits, simulations, updates and saved results, so you rarely need
to write a listener yourself.

## For developers: `nlmixrUpdateObject()`

`nlmixrUpdateObject(fit, objName, envir, origFitEnv)` rebinds the
variable `objName` in `envir` to the updated `fit` when that variable
still holds the original fit. It now returns `TRUE` (invisibly) when it
rebound the variable and `FALSE` otherwise, for instance when `objName`
is not a single name (as for `addCwres(fits[[1]])`, where it used to
error) or the variable is not bound to the original fit. A package that
writes its own updating function can use this value to set `inPlace` on
the `fitUpdate` event it emits.
