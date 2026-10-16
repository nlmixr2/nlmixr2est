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

# Functions on the chain with a NULL or noted exit that are not leaves, and why.
.covAnalyticDeclineExempt <- c(
  .foceiCalcRanalyticAssemble = "an entry point: every return() is checked by the return gate",
  .foceiAnalyticAssembleR = "an assembler: every return() is checked by the return gate",
  .foceiAnalyticAssembleRFR = "an assembler: every return() is checked by the return gate",
  .foceiAnalyticFallback = "the decline itself",
  .foceiAnalyticObsAligned = "returns FALSE with a note; its callers are the solve leaves",
  .foceiAnalyticCovDirMap = "NULL means no eta-scaled covariate directions, not a failure",
  .vaeOuterCols = "NULL selects the rxSolve route of .foceiAnalyticSolveAll, not a decline",
  .toRx = "part of the model build, whose failure is the .foceiAnalyticAugModelDirs leaf"
)

# Entry points of the chain, and the functions a closure walk stops at (the
# control constructor only reads settings).
.covAnalyticDeclineEntries <- c(".foceiCalcRanalytic", ".foceiCalcRanalyticAssemble", ".foceiCovAnalyticCalc")
.covAnalyticDeclineStop <- "foceiControl"

# The package functions a function names, by the symbols in its body
.covAnalyticNamedFns <- function(fn, ns) {
  .nm <- unique(all.names(body(get(fn, envir = ns))))
  .nm[vapply(
    .nm,
    function(.n) exists(.n, envir = ns, inherits = FALSE) && is.function(get(.n, envir = ns)),
    logical(1)
  )]
}

# Every package function the entry points can reach
.covAnalyticDeclineClosure <- function(ns = asNamespace("nlmixr2est")) {
  .seen <- character(0)
  .todo <- .covAnalyticDeclineEntries
  while (length(.todo) > 0L) {
    .f <- .todo[1L]
    .todo <- .todo[-1L]
    if (.f %in% .seen) {
      next
    }
    .seen <- c(.seen, .f)
    if (!(.f %in% .covAnalyticDeclineStop)) {
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

# The return() calls anywhere in an expression
.covAnalyticReturns <- function(e) {
  if (!is.call(e)) {
    return(list())
  }
  .r <- if (identical(e[[1]], as.name("return"))) list(e) else list()
  # an empty argument (x[, 1]) is the missing symbol
  .sub <- lapply(as.list(e)[-1], function(.a) if (missing(.a)) list() else .covAnalyticReturns(.a))
  c(.r, unlist(.sub, recursive = FALSE))
}

# What each return() returns: the called function's name, or the value deparsed
.covAnalyticReturnKinds <- function(fn, ns = asNamespace("nlmixr2est")) {
  vapply(
    .covAnalyticReturns(body(get(fn, envir = ns))),
    function(.r) {
      .a <- .r[[2]]
      if (is.call(.a)) deparse(.a[[1]])[1] else deparse(.a)[1]
    },
    character(1)
  )
}
