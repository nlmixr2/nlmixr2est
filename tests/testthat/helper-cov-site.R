# Helpers for the tests that drive a covariance-warning site (a function that
# calls .covWarn() or .covClip()): the labels a fit really carries, and a wrapper
# that checks the site was reached.  test-cov-warning-gate.R fails for a site that
# no test drives.

# Covariance-method labels as a fit carries them: every covMethod a control
# accepts, the correction decorations of the finite-difference methods (see
# .covFdType()), and the sqrtm-repaired forms.  `decorated = FALSE` keeps the
# names a user requests, `full = FALSE` leaves out the " (full)" scope suffix.
.covRealLabels <- function(decorated = TRUE, full = TRUE) {
  .ns <- asNamespace("nlmixr2est")
  .controls <- mget(grep("Control$", ls(.ns), value = TRUE), envir = .ns)
  .choices <- unlist(lapply(.controls, function(f) {
    if (!is.function(f) || !("covMethod" %in% names(formals(f)))) {
      return(character())
    }
    tryCatch(as.character(eval(formals(f)$covMethod)), error = function(e) character())
  }))
  .base <- unique(c(.choices, "Ha", "linFim"))
  if (decorated) {
    .base <- unique(c(
      .base,
      outer(c("r", "r+", "|r|"), c("s", "s+", "|s|"), paste, sep = ","),
      c("r+", "|r|", "s+", "|s|", "|linFim|", "|Ha|")
    ))
  }
  .base <- .base[vapply(.base, .covIsName, logical(1))]
  if (full) {
    .base <- unique(c(.base, paste0(.base, .covFullSuffix)))
  }
  .base
}

# a fit environment holding a covariance
.covSitePd <- function(nm = c("tka", "tcl")) {
  matrix(c(0.04, 0.01, 0.01, 0.09), 2, dimnames = list(nm, nm))
}

.covSiteEnv <- function(cov = .covSitePd(), covMethod = "r,s") {
  .e <- new.env(parent = emptyenv())
  .e$cov <- cov
  .e$covMethod <- covMethod
  .e$objDf <- data.frame(OBJF = 1, "Condition#(Cov)" = 99, "Condition#(Cor)" = 99, check.names = FALSE)
  .e
}

# does `e` contain `x` as a sub-expression?
.covSiteHas <- function(e, x) {
  if (identical(e, x)) {
    return(TRUE)
  }
  if (!is.call(e) && !is.pairlist(e)) {
    return(FALSE)
  }
  any(vapply(
    seq_along(e),
    function(i) !identical(e[[i]], quote(expr = )) && .covSiteHas(e[[i]], x),
    logical(1)
  ))
}

# is `f` the function `target`, or a function defined inside it?
.covSiteIn <- function(f, target) {
  identical(body(f), body(target)) || .covSiteHas(body(target), body(f))
}

# Evaluate `code`, which has to reach .covWarn() or .covClip() from the package
# function named `caller`, and return its value.  The name is what the gate in
# test-cov-warning-gate.R reads.
.covSite <- function(caller, code) {
  .ns <- asNamespace("nlmixr2est")
  .target <- get(caller, envir = .ns)
  .rec <- new.env(parent = emptyenv())
  .rec$hit <- FALSE
  .warn <- get(".covWarn", envir = .ns)
  .clip <- get(".covClip", envir = .ns)
  # the frame the call was written in: a promise (code given to a timer, say)
  # runs inside the function that wrote it
  .note <- function(frame) {
    if (frame > 0L && .covSiteIn(sys.function(frame), .target)) {
      .rec$hit <- TRUE
    }
  }
  .out <- testthat::with_mocked_bindings(
    code,
    .covWarn = function(w) {
      .note(sys.parent())
      .warn(w)
    },
    .covClip = function(w, width = 74L) {
      .note(sys.parent())
      .clip(w, width)
    }
  )
  testthat::expect_true(.rec$hit, label = paste0("a covariance warning from ", caller))
  .out
}

# the width a $runInfo line renders in, and whether `w` is one whole line
.covSiteWhole <- function(w) {
  all(cli::ansi_nchar(w, type = "width") <= 74L) && !any(endsWith(w, "..."))
}
