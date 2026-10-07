#' Deparse a value so that evaluating it gives back the identical value
#' @param value value to deparse
#' @return deparsed string
#' @noRd
.deparseValue <- function(value) {
  .d <- deparse1(value)
  if (identical(tryCatch(eval(str2lang(.d)), error = function(e) NULL), value)) {
    return(.d)
  }
  # the default 15 significant digits do not always read back the same double
  deparse1(value, control = c("keepNA", "keepInteger", "niceNames", "showAttributes", "digits17"))
}

#' The `print = ` argument that rebuilds an `iterPrintControl()`
#' @param value the control's `iterPrintControl`
#' @return `print = <every>`, or `print = iterPrintControl(...)` when more than
#'   `every` differs from this session's default
#' @noRd
.deparseIterPrint <- function(value) {
  .def <- iterPrintControl()
  .w <- names(.def)[!vapply(names(.def), function(n) identical(.def[[n]], value[[n]]), logical(1))]
  if (identical(.w, "every")) {
    return(paste0("print = ", deparse1(value$every)))
  }
  paste0(
    "print = iterPrintControl(",
    paste(paste0(.w, " = ", vapply(.w, function(n) deparse1(value[[n]]), character(1))), collapse = ", "),
    ")"
  )
}

#' The default a control is deparsed against, built at the control's `sigdig`
#'
#' Tolerances derived from `sigdig` then match and are not written out.  A
#' control that does not keep `sigdig` (`saemControl()`) is taken to have been
#' built with the `sigdig` its `sigdigTable` or `tol` implies when that rebuilds
#' its `rxControl`; the default then carries that `sigdig` as an attribute.
#' @param default the constructor's default control
#' @param object the control being deparsed
#' @param ctor name of the constructor
#' @return the default control at the control's `sigdig`
#' @noRd
.deparseSigdigDefault <- function(default, object, ctor = class(default)[1]) {
  if (!any(c("sigdig", "...") %in% names(formals(ctor)))) {
    return(default)
  }
  .rebuild <- function(sig) tryCatch(do.call(ctor, list(sigdig = sig)), error = function(e) NULL)
  .sig <- object[["sigdig"]]
  if (!is.null(.sig)) {
    if (identical(.sig, default[["sigdig"]])) {
      return(default)
    }
    .ret <- .rebuild(.sig)
    return(if (is.null(.ret)) default else .ret)
  }
  if (!isTRUE(object[["genRxControl"]]) || identical(object[["rxControl"]], default[["rxControl"]])) {
    return(default)
  }
  # not kept: sigdigTable and tol (10^-sigdig) are what sigdig sets by default
  .tol <- object[["tol"]]
  .cand <- c(object[["sigdigTable"]], if (is.numeric(.tol) && length(.tol) == 1L && .tol > 0) -log10(.tol))
  for (.sig in unique(.cand[is.finite(.cand)])) {
    if (abs(.sig - round(.sig)) < 1e-8) {
      .sig <- round(.sig)
    }
    .ret <- .rebuild(.sig)
    if (!is.null(.ret) && identical(.ret[["rxControl"]], object[["rxControl"]])) {
      attr(.ret, "sigdig") <- .sig
      return(.ret)
    }
  }
  default
}

.deparseShared <- function(x, value) {
  if (x == "iterPrintControl") {
    .deparseIterPrint(value)
  } else if (x == "rxControl") {
    .rx <- rxUiDeparse(value, "a")
    .rx <- .rx[[3]]
    # rxode2 writes 15 significant digits; give a value that does not read back
    # the same its full precision
    .args <- vapply(
      seq_along(.rx)[-1L],
      function(i) {
        .n <- names(.rx)[i]
        .v <- value[[.n]]
        .d <- if (is.numeric(.v) && !identical(tryCatch(eval(.rx[[i]]), error = function(e) NULL), .v)) {
          .deparseValue(.v)
        } else {
          deparse1(.rx[[i]])
        }
        if (is.null(.n) || .n == "") .d else paste0(.n, " = ", .d)
      },
      character(1)
    )
    paste0("rxControl = ", deparse1(.rx[[1]]), "(", paste(.args, collapse = ", "), ")")
  } else if (x == "scaleType") {
    if (is.integer(value)) {
      paste0("scaleType =", deparse1(names(.scaleTypeIdx[which(value == .scaleTypeIdx)])))
    } else {
      paste0("scaleType =", deparse1(value))
    }
  } else if (x == "normType") {
    if (is.integer(value)) {
      paste0("normType =", deparse1(names(.normTypeIdx[which(value == .normTypeIdx)])))
    } else {
      paste0("normType =", deparse1(value))
    }
  } else if (x == "solveType") {
    if (is.integer(value)) {
      .solveTypeIdx <- c("hessian" = 3L, "grad" = 2L, "fun" = 1L)
      paste0("solveType =", deparse1(names(.solveTypeIdx[which(value == .solveTypeIdx)])))
    } else {
      paste0("solveType =", deparse1(value))
    }
  } else if (x == "eventType") {
    if (is.integer(value)) {
      .eventTypeIdx <- c("central" = 2L, "forward" = 1L, "forward" = 3L)
      paste0("eventType = ", deparse1(names(.eventTypeIdx[which(value == .eventTypeIdx)])))
    } else {
      paste0("eventType = ", deparse1(value))
    }
  } else if (x == "censMethod") {
    if (is.integer(value)) {
      .censMethodIdx <- c("truncated-normal" = 3L, "cdf" = 2L, "omit" = 1L, "pred" = 5L, "ipred" = 4L, "epred" = 6L)
      paste0("censMethod = ", deparse1(names(.censMethodIdx[which(value == .censMethodIdx)])))
    } else {
      paste0("censMethod = ", deparse1(value))
    }
  } else {
    NA_character_
  }
}

#' Identify Differences Between Standard and New Objects but used in rxUiDeparse
#'
#' Compares elements of a standard object against a new one so
#' `rxUiDeparse` only shows values that differ from the default.
#'
#' @param standard The standard object used for comparison. (for example `foceiControl()`)
#'
#' @param new The new object to be compared against the standard. This
#'   would be what the user supplide like
#'   `foceiControl(outerOpt="bobyqa")`
#' @param internal A character vector of element names to be ignored
#'   during the comparison. Default is an empty character
#'   vector. These are for internal items of the list that flag
#'   certain properties like if the `rxControl()` was generated by the
#'   `foceiControl()` procedure or not.
#' @return A vector of indices indicating which elements of the
#'   standard object differ from the new object.
#' @examples
#' standard <- list(a = 1, b = 2, c = 3)
#' new <- list(a = 1, b = 3, c = 3)
#' .deparseDifferent(standard, new)
#' @export
#' @keywords internal
#' @author Matthew L. Fidler
.deparseDifferent <- function(standard, new, internal = character(0)) {
  which(vapply(
    names(standard),
    function(x) {
      if (x %in% internal) {
        FALSE
      } else if (is.function(standard[[x]])) {
        warning(paste0("'", x, "' as a function not supported in ", class(standard), "() deparsing"), call. = FALSE)
        FALSE
      } else {
        !identical(standard[[x]], new[[x]])
      }
    },
    logical(1),
    USE.NAMES = FALSE
  ))
}

#' Deparse finalize a control or related object into a language object
#'
#' Deparses an object into a language expression, optionally using a
#' custom function for specific elements.
#'
#' @param default A default object used for comparison; This is the
#'   estimation control procedure.  It should have a class matching
#'   the function that created it.
#' @param object The object to be deparsed into a language expression
#' @param w A vector of indices indicating which elements are
#'   different and need to be deparsed. This likely comes from
#'   `.deparseDifferent()`
#' @param var A string representing the variable name to be assigned
#'   in the deparsed expression.
#' @param fun An optional custom function to handle specific elements
#'   during deparsing. Default is NULL. This handles things that are
#'   specific to an estimation control and is used by functions like
#'   `rxUiDeparse.saemControl()`
#' @return A language object representing the deparsed expression.
#' @keywords internal
#' @author Matthew L. Fidler
#' @export
.deparseFinal <- function(default, object, w, var, fun = NULL) {
  .cls <- class(object)
  if (length(w) == 0) {
    return(str2lang(paste0(var, " <- ", .cls, "()")))
  }
  .retD <- vapply(
    names(default)[w],
    function(x) {
      .val <- .deparseShared(x, object[[x]])
      if (!is.na(.val)) {
        return(.val)
      }
      if (is.function(fun)) {
        .val <- fun(default, x, object[[x]])
        if (!is.na(.val)) {
          return(.val)
        }
      }
      paste0(x, "=", .deparseValue(object[[x]]))
    },
    character(1),
    USE.NAMES = FALSE
  )
  str2lang(paste(var, " <- ", .cls, "(", paste(.retD, collapse = ","), ")"))
}

#' Write out the arguments a deparsed control does not rebuild
#'
#' A default that follows another argument (`n1qn1nsim` follows
#' `maxInnerIterations`) is not seen by comparing with the constructor's own
#' default; evaluating the call is.  Only constructor arguments are added, and
#' only while the call still evaluates.
#' @param ret `var <- ctor(...)` language object
#' @param object the control
#' @param allowed names that may be added
#' @param token function giving the `name = value` string of a name
#' @return `ret`, with any missing arguments added
#' @noRd
.deparseFixup <- function(ret, object, allowed, token) {
  .eval <- function(r) tryCatch(suppressWarnings(suppressMessages(eval(r[[3]]))), error = function(e) NULL)
  .r <- .eval(ret)
  for (.k in 1:3) {
    if (is.null(.r)) {
      return(ret)
    }
    .d <- intersect(allowed, names(object))
    .d <- .d[!vapply(.d, function(n) identical(object[[n]], .r[[n]]), logical(1))]
    if (length(.d) == 0L) {
      return(ret)
    }
    .args <- as.list(str2lang(paste0("f(", paste(vapply(.d, token, character(1)), collapse = ", "), ")")))[-1]
    .cur <- as.list(ret[[3]])
    for (.nm in names(.args)) {
      .cur[[.nm]] <- .args[[.nm]]
    }
    .new <- ret
    .new[[3]] <- as.call(.cur)
    .r <- .eval(.new)
    if (is.null(.r)) {
      return(ret)
    }
    ret <- .new
  }
  ret
}

#' Deparse a control as a call to its constructor
#'
#' The body of the `rxUiDeparse()` methods of the estimation controls: the call
#' sets only what differs from the constructor's default.
#' @param object the control
#' @param var name the call is assigned to
#' @param default the control the constructor builds by default
#' @param internal names that are never deparsed
#' @param fun see [.deparseFinal()]
#' @return the language object `var <- <constructor>(...)`
#' @noRd
.deparseControl <- function(object, var, default, internal = "genRxControl", fun = NULL) {
  .default <- .deparseSigdigDefault(default, object)
  .w <- .deparseDifferent(.default, object, internal)
  # an rxControl that was supplied, even one equal to the generated one
  if (!identical(object[["genRxControl"]], .default[["genRxControl"]])) {
    .w <- sort(union(.w, which(names(.default) == "rxControl")))
  }
  if (!identical(object[["sigdig"]], default[["sigdig"]])) {
    .w <- sort(union(.w, which(names(.default) == "sigdig")))
  }
  .ret <- .deparseFinal(.default, object, .w, var, fun = fun)
  .sig <- attr(.default, "sigdig")
  if (!is.null(.sig)) {
    # sigdig recovered for a control that does not keep it goes first
    .ret[[3]] <- as.call(c(as.list(.ret[[3]])[1], list(sigdig = .sig), as.list(.ret[[3]])[-1]))
  }
  .tok <- function(x) {
    .val <- .deparseShared(x, object[[x]])
    if (is.na(.val) && is.function(fun)) {
      .val <- fun(.default, x, object[[x]])
    }
    if (is.na(.val)) paste0(x, " = ", .deparseValue(object[[x]])) else .val
  }
  .deparseFixup(.ret, object, setdiff(names(formals(class(default)[1])), c("...", internal, "rxControl")), .tok)
}
