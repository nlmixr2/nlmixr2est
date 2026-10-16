# A warning() that reports on a covariance method goes through .covWarn(), which
# clips it to the width of a $runInfo line.  The package's own functions are
# walked, so a warning added later is judged without being listed here.

# a warning whose text is about a covariance method has one of these in it
.covGatePattern <- paste(
  "covariance",
  "\\bFIM\\b",
  "sandwich",
  "\\b[RS] matrix",
  "\\$cov[A-Z]",
  "\\bcov(Full|Method)\\b",
  "\\bSE\\b",
  sep = "|"
)

.covGateWarners <- c("warning", "warningCondition", "simpleWarning", "warn", "cli_warn")

# the function a call names, "" when it is not a name or pkg::name
.covGateHead <- function(e) {
  .h <- e[[1L]]
  if (is.name(.h)) {
    return(as.character(.h))
  }
  if (is.call(.h) && length(.h) == 3L && as.character(.h[[1L]]) %in% c("::", ":::")) {
    return(as.character(.h[[3L]]))
  }
  ""
}

# every call to one of `heads` anywhere in `e`, argument defaults included
.covGateWalk <- function(e, heads) {
  if (!is.call(e) && !is.pairlist(e)) {
    return(list())
  }
  .own <- if (is.call(e) && .covGateHead(e) %in% heads) list(e) else list()
  .sub <- lapply(seq_along(e), function(i) {
    if (identical(e[[i]], quote(expr = ))) list() else .covGateWalk(e[[i]], heads)
  })
  c(.own, unlist(.sub, recursive = FALSE))
}

# every string constant in `e`
.covGateStrings <- function(e) {
  if (is.character(e)) {
    return(e)
  }
  if (!is.call(e) && !is.pairlist(e)) {
    return(character())
  }
  .sub <- lapply(seq_along(e), function(i) {
    if (identical(e[[i]], quote(expr = ))) character() else .covGateStrings(e[[i]])
  })
  as.character(unlist(.sub))
}

# The functions of `fns` that raise a covariance warning without .covWarn().  A
# warning() with no text of its own passes a condition on, so the text of its
# function stands for it.
.covGateBypasses <- function(fns, pattern = .covGatePattern, funnel = ".covWarn") {
  .hit <- character()
  for (.n in setdiff(names(fns), funnel)) {
    for (.call in .covGateFnCalls(fns[[.n]], .covGateWarners)) {
      .own <- .covGateStrings(.call)
      .txt <- if (length(.own) > 0L) .own else .covGateFnStrings(fns[[.n]])
      if (any(grepl(pattern, .txt))) {
        .hit <- c(.hit, .n)
      }
    }
  }
  unique(.hit)
}

# The functions of `fns` with a .covWarn() text whose fixed part alone is wider
# than a $runInfo line, labels and counts left out.
.covGateTooWide <- function(fns, width = 74L) {
  .hit <- character()
  for (.n in names(fns)) {
    for (.call in .covGateFnCalls(fns[[.n]], ".covWarn")) {
      .fixed <- gsub("%%", "%", gsub("%[-+ 0#]*[0-9]*(\\.[0-9]+)?[sdif]", "", .covGateStrings(.call)))
      if (any(cli::ansi_nchar(.fixed, type = "width") > width)) {
        .hit <- c(.hit, .n)
      }
    }
  }
  unique(.hit)
}

# every call to one of `heads` in a function, its argument defaults included
.covGateFnCalls <- function(f, heads) {
  c(.covGateWalk(body(f), heads), .covGateWalk(formals(f), heads))
}

# every string constant in a function, its argument defaults included
.covGateFnStrings <- function(f) {
  c(.covGateStrings(body(f)), .covGateStrings(formals(f)))
}

.covGateFunctions <- function() {
  .ns <- asNamespace("nlmixr2est")
  Filter(function(x) is.function(x) && !is.primitive(x), mget(ls(.ns, all.names = TRUE), envir = .ns))
}

test_that("the gate finds a covariance warning that bypasses .covWarn()", {
  .bad <- list(
    plain = function() warning("covariance step failed"),
    sandwich = function() warning("since sandwich matrix is corrected", call. = FALSE),
    fim = function(x) warning(sprintf("linearized FIM failed; using %s", x)),
    rMatrix = function() warning(paste0("R matrix ", "is not finite")),
    sMatrix = function() warning("S matrix had problems"),
    covR = function() warning("compare to $covR"),
    covFull = function() warning("covFull variance block dropped"),
    covMethod = function() warning("switching to covMethod='r'"),
    se = function() warning("no SE for a theta"),
    cond = function() warningCondition("covariance step failed"),
    qualified = function() rlang::warn("covariance step failed"),
    nested = function() lapply(1, function(i) if (i) warning("covariance step failed")),
    default = function(a = warning("covariance step failed")) a,
    # no text of its own: the function's other text says what it is about
    passOn = function(w) {
      .why <- "covariance step failed"
      warning(w, call. = FALSE)
    }
  )
  expect_identical(.covGateBypasses(.bad), names(.bad))
})

test_that("the gate lets a covariance warning through .covWarn() and other warnings alone", {
  .ok <- list(
    viaFunnel = function() .covWarn("covariance step failed"),
    omega = function() warning("omega not positive definite; tables use its nearest PD matrix"),
    covariate = function() warning("covariate selection pinned to model-specified covariates"),
    covariates = function() warning("time-varying covariate(s) were excluded"),
    other = function() warning("outer bobyqa search stalled"),
    passOn = function(w) {
      .why <- "omega block zero"
      warning(w, call. = FALSE)
    },
    none = function() 1
  )
  expect_identical(.covGateBypasses(.ok), character())
  # the funnel is the one function allowed a warning() about a covariance
  .funnel <- list(.covWarn = function(w) {
    .why <- "covariance step failed"
    warning(w, call. = FALSE)
  })
  expect_identical(.covGateBypasses(.funnel), character())
  expect_identical(.covGateBypasses(.funnel, funnel = "other"), ".covWarn")
})

test_that("no warning in the package reports on a covariance method except through .covWarn()", {
  .fns <- .covGateFunctions()
  # the funnel itself is the one warning() call it makes, and the package uses it
  expect_length(.covGateFnCalls(.covWarn, "warning"), 1L)
  expect_true(any(vapply(.fns, function(f) length(.covGateFnCalls(f, ".covWarn")) > 0L, logical(1))))
  expect_identical(.covGateBypasses(.fns), character())
})

# a function whose body is `call`, with the text inside it as a literal
.covGateFn <- function(call) {
  eval(bquote(function(a = 1) .(call)))
}

test_that("the gate finds a .covWarn() text too wide before any label is added", {
  .x <- function(n) strrep("x", n)
  .wide <- list(
    long = .covGateFn(call(".covWarn", .x(75))),
    format = .covGateFn(call(".covWarn", call("sprintf", paste0("\"%s\" ", .x(72), " %d"), quote(a), 1L))),
    branch = .covGateFn(call(".covWarn", call("if", quote(a), "short", .x(80)))),
    nested = .covGateFn(call("lapply", 1, call("function", as.pairlist(alist(i = )), call(".covWarn", .x(75)))))
  )
  expect_identical(.covGateTooWide(.wide), names(.wide))
  .fit <- list(
    exact = .covGateFn(call(".covWarn", .x(74))),
    format = .covGateFn(call(".covWarn", call("sprintf", paste0("\"%s\" ", .x(70), " %d"), quote(a), 1L))),
    percent = .covGateFn(call(".covWarn", call("sprintf", paste0("100%% ", .x(69))))),
    dynamic = function(w) .covWarn(w)
  )
  expect_identical(.covGateTooWide(.fit), character())
})

test_that("every .covWarn() text in the package fits a $runInfo line before any label is added", {
  expect_identical(.covGateTooWide(.covGateFunctions()), character())
})
