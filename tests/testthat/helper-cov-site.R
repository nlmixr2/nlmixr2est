# Helpers for the tests that drive a covariance-warning site (a function that
# calls .covWarn() or .covClip()): the labels a fit really carries, a wrapper that
# checks the site was reached, and the walkers over the package's functions that
# test-cov-warning-gate.R uses to find the sites.

# --- the labels a fit carries ---------------------------------------------------

# the alternatives of one regex group of .covFdType(): "r\\+?" is "r" or "r+",
# "\\|r\\|" is "|r|"
.covFdAlternatives <- function(group) {
  .alt <- strsplit(group, "(?<!\\\\)\\|", perl = TRUE)[[1L]]
  unlist(lapply(.alt, function(a) {
    if (endsWith(a, "\\+?")) {
      .b <- gsub("\\", "", substr(a, 1L, nchar(a) - 3L), fixed = TRUE)
      c(.b, paste0(.b, "+"))
    } else {
      gsub("\\", "", a, fixed = TRUE)
    }
  }))
}

# The finite-difference labels, expanded from the anchored patterns .covFdType()
# matches them with: "r", "r+", "|r|", "s", "s+", "|s|" and the sandwiches.
.covFdLabels <- function() {
  .pats <- grep("^\\^", .covGateStrings(body(.covFdType)), value = TRUE)
  unique(unlist(lapply(.pats, function(p) {
    .core <- sub("\\$$", "", sub("^\\^", "", p))
    .groups <- regmatches(.core, gregexpr("\\(([^()]*)\\)", .core))[[1L]]
    .alt <- lapply(sub("^\\((.*)\\)$", "\\1", .groups), .covFdAlternatives)
    if (length(.alt) == 1L) .alt[[1L]] else as.vector(outer(.alt[[1L]], .alt[[2L]], paste, sep = ","))
  })))
}

# The covariances that are repaired by sqrtm(x %*% x), as (requested, installed)
# pairs.  SAEM's linFim and Ha come in the theta shape only; the decoupled engines
# ("sa", "imp") also have a full shape.
.covRealRepaired <- function() {
  .req <- c("linFim", .saemHaThetaName, "sa", "imp", paste0(c("sa", "imp"), .covFullSuffix))
  data.frame(request = .req, label = sub("^([^ ]+)", "|\\1|", .req))
}

# Covariance-method labels as a fit carries them: every covMethod a control
# accepts, SAEM's "Ha", the finite-difference decorations and the sqrtm-repaired
# forms, each with and without " (full)".  `decorated = FALSE` keeps the names a
# user requests, `full = FALSE` leaves out the " (full)" scope suffix.
.covRealLabels <- function(decorated = TRUE, full = TRUE) {
  .ns <- asNamespace("nlmixr2est")
  .controls <- mget(grep("Control$", ls(.ns), value = TRUE), envir = .ns)
  .choices <- unlist(lapply(.controls, function(f) {
    if (!is.function(f) || !("covMethod" %in% names(formals(f)))) {
      return(character())
    }
    tryCatch(as.character(eval(formals(f)$covMethod)), error = function(e) character())
  }))
  .base <- unique(c(.choices, .saemHaThetaName, "linFim"))
  if (decorated) {
    .rep <- .covRealRepaired()$label
    .base <- unique(c(.base, .covFdLabels(), .rep[!endsWith(.rep, .covFullSuffix)]))
  }
  .base <- .base[vapply(.base, .covIsName, logical(1))]
  if (full) {
    .base <- unique(c(.base, paste0(.base, .covFullSuffix)))
  }
  .base
}

# The label an "r" covariance keeps when the analytic step declined it (see
# .foceiCalcRanalytic())
.covRealAnalyticKept <- function() {
  .r <- Filter(function(l) identical(.covFdType(l), "r"), .covFdLabels())
  paste0(c(.r, paste0(.r, .covFullSuffix)), " (analytic)")
}

# --- a fit environment holding a covariance --------------------------------------

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

# --- driving a site --------------------------------------------------------------

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

# Evaluate `code`, which has to reach .covWarn(), .covClip() or .covNotUsableMsg()
# from each of the package functions named in `callers`, and return its value.
# The names are what the gate in test-cov-warning-gate.R reads.
.covSite <- function(callers, code) {
  .ns <- asNamespace("nlmixr2est")
  .targets <- stats::setNames(lapply(callers, get, envir = .ns), callers)
  .rec <- new.env(parent = emptyenv())
  .rec$hit <- stats::setNames(rep(FALSE, length(callers)), callers)
  .warn <- get(".covWarn", envir = .ns)
  .clip <- get(".covClip", envir = .ns)
  .msg <- get(".covNotUsableMsg", envir = .ns)
  # the frame the call was written in: code handed to a timer runs inside the
  # function that wrote it
  .note <- function(frame) {
    if (frame > 0L && !all(.rec$hit)) {
      .f <- sys.function(frame)
      for (.n in callers[!.rec$hit]) {
        if (.covSiteIn(.f, .targets[[.n]])) {
          .rec$hit[[.n]] <- TRUE
        }
      }
    }
  }
  .out <- testthat::with_mocked_bindings(
    code,
    .covWarn = function(w) {
      .note(sys.parent())
      .warn(w)
    },
    .covClip = function(w, width = .covWarnWidth) {
      .note(sys.parent())
      .clip(w, width)
    },
    .covNotUsableMsg = function(request, why, tail) {
      .note(sys.parent())
      .msg(request, why, tail)
    }
  )
  for (.n in callers) {
    testthat::expect_true(.rec$hit[[.n]], label = paste0("a covariance warning from ", .n))
  }
  .out
}

# the first warning `expr` gives, as text; NA when it gives none
.covSiteWarning <- function(expr) {
  tryCatch(
    {
      expr
      NA_character_
    },
    warning = conditionMessage
  )
}

# is every line of `w` a whole line of a $runInfo line's width?
.covSiteWhole <- function(w) {
  all(cli::ansi_nchar(w, type = "width") <= .covWarnWidth) && !any(endsWith(w, "..."))
}

# --- walking the package's functions ---------------------------------------------

# the function a call names, "" when it is not a name or pkg::name
.covGateHead <- function(e) {
  .h <- e[[1L]]
  if (is.name(.h)) {
    return(as.character(.h))
  }
  if (is.call(.h) && length(.h) == 3L && is.name(.h[[1L]]) && as.character(.h[[1L]]) %in% c("::", ":::")) {
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

