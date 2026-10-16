# The decline chain of covMethod="analytic", for the always-run static gate in
# test-cov-analytic-units.R and the weekly enumeration in
# test-cov-analytic-declines.R.

# Helpers the enumeration makes fail, one entry at a time (a pair fails both).
.covAnalyticDeclineLeaves <- list(
  ".foceiEtaThetaMap",
  ".foceiAnalyticErrFull",
  ".foceiAnalyticDirections",
  ".foceiOmegaPairs",
  ".omegaVarCovDeriv",
  ".foceiAnalyticAugModelDirs",
  ".nmSubjectObs",
  ".foceiAnalyticSolveAll",
  ".foceiAnalyticRxSolve",
  ".foceiAnalyticSolveAllFD3",
  c(".foceiAnalyticSolveAllFD3", ".foceiAnalyticSolveSubjectFD3"),
  c(".foceiAnalyticSolveAllFD3", ".foceiAnalyticSolveFA"),
  ".foceiAnalyticExpandSigma",
  ".foceiAnalyticSubjectR",
  ".foceiAnalyticSubjectRfoce",
  "foceiRAllFR_",
  "foceiRAllFoceFR_",
  ".foceiAnalyticSgNames",
  ".foceiAnalyticSolve"
)

# Functions the walk reaches, outside the leaves, with a NULL or noted exit that
# are not leaves, and why.
.covAnalyticDeclineExempt <- c(
  .foceiCalcRanalyticAssemble = "an entry point: every return() is checked by the return gate",
  .foceiAnalyticAssembleR = "an assembler: every return() is checked by the return gate",
  .foceiAnalyticAssembleRFR = "an assembler: every return() is checked by the return gate",
  .foceiAnalyticFallback = "the decline itself"
)

# Entry points of the chain
.covAnalyticDeclineEntries <- c(".foceiCalcRanalytic", ".foceiCalcRanalyticAssemble", ".foceiCovAnalyticCalc")

# Collect, into `acc`, the symbols an expression names and the ui properties it
# reads: `x$name` and `rxUiGet(x, "name")` add "name" to acc$fields (rxode2
# dispatches them to rxUiGet.name()), and the field after `$` or `@` is not a
# symbol it names.  A string counts as a name too, for a function called by
# name (do.call(".f", ...), match.fun(), lapply(x, ".f")); any string equal to
# a package function's name is included, so labels such as
# ini$err %in% c("boxCox", "yeoJohnson") add those functions as well.
.covAnalyticWalk <- function(e, acc) {
  if (is.symbol(e)) {
    acc$syms <- c(acc$syms, as.character(e))
    return(invisible())
  }
  if (is.character(e)) {
    acc$syms <- c(acc$syms, e)
    return(invisible())
  }
  if (!is.call(e)) {
    return(invisible())
  }
  .h <- e[[1]]
  if (is.symbol(.h) && as.character(.h) %in% c("$", "@")) {
    .covAnalyticWalk(e[[2]], acc)
    .f <- e[[3]]
    if (identical(as.character(.h), "$") && (is.symbol(.f) || is.character(.f))) {
      acc$fields <- c(acc$fields, as.character(.f))
    }
    return(invisible())
  }
  .uiGet <- identical(.h, as.name("rxUiGet")) ||
    (is.call(.h) && identical(.h[[1]], as.name("::")) && identical(.h[[3]], as.name("rxUiGet")))
  if (.uiGet && length(e) >= 3L && is.character(e[[3]])) {
    acc$fields <- c(acc$fields, e[[3]])
  }
  for (.a in as.list(e)) {
    # an empty argument (x[, 1]) is the missing symbol
    if (!missing(.a)) .covAnalyticWalk(.a, acc)
  }
  invisible()
}

# The package functions a function calls or passes: the symbols in its body
# that are package functions, and the rxUiGet.<name>() methods behind the ui
# properties it reads
.covAnalyticNamedFns <- function(fn, ns) {
  .acc <- new.env(parent = emptyenv())
  .acc$syms <- character(0)
  .acc$fields <- character(0)
  .covAnalyticWalk(body(get(fn, envir = ns)), .acc)
  .nm <- unique(c(.acc$syms, paste0("rxUiGet.", .acc$fields)))
  # exists() rejects an empty name and one over 10000 bytes
  .nm <- .nm[!is.na(.nm) & nzchar(.nm) & nchar(.nm, type = "bytes") <= 10000L]
  .nm[vapply(
    .nm,
    function(.n) exists(.n, envir = ns, inherits = FALSE) && is.function(get(.n, envir = ns)),
    logical(1)
  )]
}

# Every package function the entry points can reach.  With `stopAtLeaves`, the
# walk does not go below a leaf: a function only reached through leaves can make
# the covariance decline only by making a leaf give up, which the enumeration does
# for every leaf.  The model build behind .foceiAnalyticAugModelDirs (the
# rxUiGet.* methods it reads) is such a subtree.
.covAnalyticDeclineClosure <- function(
  ns = asNamespace("nlmixr2est"),
  stopAtLeaves = TRUE,
  entries = .covAnalyticDeclineEntries,
  leaves = unique(unlist(.covAnalyticDeclineLeaves))
) {
  .seen <- character(0)
  .todo <- entries
  while (length(.todo) > 0L) {
    .f <- .todo[1L]
    .todo <- .todo[-1L]
    if (.f %in% .seen) {
      next
    }
    .seen <- c(.seen, .f)
    if (!(stopAtLeaves && .f %in% leaves)) {
      .todo <- c(.todo, setdiff(.covAnalyticNamedFns(.f, ns), .seen))
    }
  }
  .seen
}

# Whether a function can give up: it notes a detail or returns a bare NULL
.covAnalyticGivesUp <- function(fn, ns = asNamespace("nlmixr2est")) {
  .b <- paste(deparse(body(get(fn, envir = ns))), collapse = "\n")
  grepl(".foceiAnalyticNote(", .b, fixed = TRUE) || grepl("return(NULL)", .b, fixed = TRUE)
}

# The return() calls in an expression: those of the function itself
# (`nested = FALSE`, not looking inside the functions it defines), or those of the
# functions it defines (`nested = TRUE`).  The gate covers explicit return() only:
# the error handler of an assembly (tryCatch(error = .foceiAnalyticErrWarn(...)))
# sets the reason by assignment, which the "an error (" pins in
# test-cov-analytic-declines.R and -reason.R cover.
.covAnalyticReturns <- function(e, nested = FALSE, inside = FALSE) {
  if (!is.call(e)) {
    return(list())
  }
  if (identical(e[[1]], as.name("function"))) {
    return(if (nested || inside) .covAnalyticReturns(e[[3]], nested, inside = TRUE) else list())
  }
  .r <- if (identical(e[[1]], as.name("return")) && (inside || !nested)) list(e) else list()
  # an empty argument (x[, 1]) is the missing symbol
  .sub <- lapply(as.list(e)[-1], function(.a) if (missing(.a)) list() else .covAnalyticReturns(.a, nested, inside))
  c(.r, unlist(.sub, recursive = FALSE))
}

# What each return() returns: the called function's name, or the value deparsed
.covAnalyticReturnKinds <- function(fn, ns = asNamespace("nlmixr2est"), nested = FALSE) {
  vapply(
    .covAnalyticReturns(body(get(fn, envir = ns)), nested = nested),
    function(.r) {
      .a <- .r[[2]]
      if (is.call(.a)) deparse(.a[[1]])[1] else deparse(.a)[1]
    },
    character(1)
  )
}
