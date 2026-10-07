# Sensitivities through history functions (lag(), diff(), ...) of a calculated
# variable (#1176).  rxode2's symengine binds a lagged calculated
# variable to a bare symbol, so its derivatives stop there.  A history function
# is linear, so d(lag(v))/d(p) = lag(d(v)/d(p)): each lagged variable gets a
# sensitivity lhs and every derivative is chained through it.

.foceiHistFn <- c("lag", "lead", "diff", "first", "last", "lag0", "lead0", "diff0")

#' The definitions of the lagged calculated variables
#'
#' The AR(1) residual's own lagged variables (`rx_ar*`) are excluded; the
#' AR(1) gradient correction handles those.
#' @param s symengine environment
#' @return the `var=expr` lines of `s$..lhs` that define them
#' @noRd
.foceiLagDefs <- function(s) {
  .v <- s$..laggedVars
  if (length(.v) == 0L || is.null(s$..lhs)) {
    return(character(0))
  }
  .v <- .v[!grepl("^rx_ar", .v)]
  s$..lhs[sub("=.*$", "", s$..lhs) %in% .v]
}

#' Whether rxode2 text uses any of the given variables
#'
#' @param txt rxode2 text
#' @param vars variable names
#' @param hist when `TRUE`, only count a use as a history function argument
#' @return logical
#' @noRd
.foceiLagRefs <- function(txt, vars, hist = FALSE) {
  .pre <- if (hist) paste0("\\b(", paste(.foceiHistFn, collapse = "|"), ")\\(\\s*") else "(?<![A-Za-z0-9_.])"
  .post <- if (hist) "\\s*[,)]" else "(?![A-Za-z0-9_.])"
  any(vapply(
    vars,
    function(v) any(grepl(paste0(.pre, gsub(".", "\\.", v, fixed = TRUE), .post), txt, perl = TRUE)),
    logical(1)
  ))
}

#' Whether a model uses a history function of a calculated variable
#'
#' A dosing `lag(cmt) <-` is not a history function, and a covariate's
#' history needs no sensitivity.
#' @param ui rxode2 UI
#' @return logical
#' @noRd
.foceiUsesLagVar <- function(ui) {
  .lhs <- tryCatch(rxode2::rxModelVars(ui)$lhs, error = function(e) character(0))
  .lhs <- .lhs[!grepl("^rx_ar", .lhs)]
  if (length(.lhs) == 0L) {
    return(FALSE)
  }
  .found <- FALSE
  .walk <- function(e) {
    if (.found || !is.call(e)) {
      return(invisible())
    }
    .f <- e[[1]]
    if (is.name(.f) && as.character(.f) %in% c("<-", "=", "~")) {
      .walk(e[[3]])
      return(invisible())
    }
    if (
      is.name(.f) &&
        as.character(.f) %in% .foceiHistFn &&
        length(e) >= 2L &&
        is.name(e[[2]]) &&
        as.character(e[[2]]) %in% .lhs
    ) {
      .found <<- TRUE
      return(invisible())
    }
    lapply(as.list(e)[-1], .walk)
    invisible()
  }
  .lst <- tryCatch(ui$lstExpr, error = function(e) NULL)
  lapply(.lst, .walk)
  .found
}

#' Sensitivities through the lagged variables
#'
#' Each lagged variable `v` gets the lhs `rx_lsens_<i>_<p>` = d(v)/d(p) for
#' each parameter `p`, and a history call `h(v, ...)` differentiates to
#' `h(rx_lsens_<i>_<p>, ...)`.
#' @param s symengine environment holding the state sensitivities
#' @param stateVars the model states
#' @param pars the parameter symbols, e.g. `ETA_1_`
#' @param exprs names of the symengine expressions that will be differentiated
#' @param statePars the parameters with state sensitivities
#' @return list with `lines` (the
#'   sensitivity lhs, in model order), `dfe(e, p, lagOnly)`, the symengine
#'   total derivative of `e` by `p` (only its terms through the lagged
#'   variables when `lagOnly`), and `txt(d, p)`, its rxode2 text
#' @noRd
.foceiLagSens <- function(s, stateVars, pars, exprs = "rx_pred_", statePars = pars) {
  .defs <- .foceiLagDefs(s)
  .defVar <- sub("=.*$", "", .defs)
  .defRhs <- sub("^[^=]*=", "", .defs)
  .vars <- unique(.defVar)
  .sym <- function(x) symengine::S(x)
  .exprs <- lapply(exprs[vapply(exprs, exists, logical(1), envir = s)], get, envir = s)
  .rhsSE <- lapply(.defRhs, function(x) .sym(rxode2::rxToSE(x)))
  # history calls of a lagged variable, taken from the symengine trees: a
  # text round trip does not keep them identical (`lag(c0, 1)` vs `1.0`)
  .histSE <- list()
  .walk <- function(e) {
    if (symengine::get_type(e) == "FunctionSymbol") {
      .a <- as.list(symengine::get_args(e))
      if (
        sub("\\(.*$", "", as.character(e)) %in% .foceiHistFn && length(.a) >= 1L && as.character(.a[[1]]) %in% .vars
      ) {
        .histSE[[as.character(e)]] <<- e
        return(invisible())
      }
    }
    lapply(as.list(symengine::get_args(e)), .walk)
    invisible()
  }
  lapply(c(.exprs, .rhsSE), function(e) .walk(.sym(e)))
  .histSE <- unname(.histSE)
  .calls <- vapply(
    .histSE,
    function(e) {
      # rxFromSE() is NSE: hand it a plain variable
      rxode2::rxFromSE(e)
    },
    character(1)
  )
  .lsens <- function(i, p) paste0("rx_lsens_", i, "_", p)
  # the history call with its variable replaced by that variable's sensitivity
  .histSens <- vapply(
    .calls,
    function(x) {
      .e <- str2lang(x)
      .e[[2]] <- as.name(.lsens(match(as.character(.e[[2]]), .vars), "rxLagParPh"))
      deparse1(.e)
    },
    character(1),
    USE.NAMES = FALSE
  )
  .subsHist <- function(e) {
    for (j in seq_along(.histSE)) {
      e <- symengine::subs(e, .histSE[[j]], .sym(paste0("rx_hist_", j, "_")))
    }
    e
  }
  .chain <- function(e, p, lagOnly = FALSE) {
    .ret <- symengine::S(0)
    if (!lagOnly) {
      .ret <- symengine::D(e, .sym(p))
      for (.st in if (p %in% statePars) stateVars) {
        .ret <- .ret + symengine::D(e, .sym(.st)) * .sym(paste0("rx__sens_", .st, "_BY_", p, "__"))
      }
    }
    for (i in seq_along(.vars)) {
      .ret <- .ret + symengine::D(e, .sym(.vars[i])) * .sym(.lsens(i, p))
    }
    for (j in seq_along(.histSE)) {
      .ret <- .ret + symengine::D(e, .sym(paste0("rx_hist_", j, "_"))) * .sym(paste0("rx_hsens_", j, "_"))
    }
    .ret
  }
  .unsub <- function(e, p) {
    .txt <- rxode2::rxFromSE(e)
    for (j in seq_along(.histSE)) {
      .txt <- gsub(paste0("\\brx_hist_", j, "_(?![A-Za-z0-9_.])"), .calls[j], .txt, perl = TRUE)
      .txt <- gsub(
        paste0("\\brx_hsens_", j, "_(?![A-Za-z0-9_.])"),
        gsub("rxLagParPh", p, .histSens[j], fixed = TRUE),
        .txt,
        perl = TRUE
      )
    }
    .txt
  }
  .rhsSE <- lapply(.rhsSE, .subsHist)
  .lines <- unlist(lapply(seq_along(.defRhs), function(d) {
    vapply(
      pars,
      function(p) paste0(.lsens(match(.defVar[d], .vars), p), "=", .unsub(.chain(.rhsSE[[d]], p), p)),
      character(1),
      USE.NAMES = FALSE
    )
  }))
  list(
    lines = .lines,
    dfe = function(e, p, lagOnly = FALSE) .chain(.subsHist(.sym(e)), p, lagOnly),
    txt = .unsub
  )
}

#' The lagged-variable eta sensitivities of a FOCEi inner model
#'
#' @param x list holding the rxode2 UI
#' @param s symengine environment holding the eta sensitivities
#' @param stateVars the model states
#' @return `NULL` when the model uses no history function of a variable; else
#'   the `.foceiLagSens()` list, with its lines stored on `s$..lagSens`
#' @noRd
.foceiLagEtaSens <- function(x, s, stateVars) {
  s$..lagSens <- NULL
  s$..lagEta <- NULL
  if (!.foceiUsesLagVar(x[[1]])) {
    return(NULL)
  }
  if (.foceiLagInOde(s)) {
    stop("lag() of a variable inside an ODE is not supported", call. = FALSE)
  }
  .eta <- paste0("ETA_", seq_len(s$..maxEta), "_")
  # the combined eta+theta build (#958) carries theta columns too
  .theta <- if (!is.null(s$..combThetaIdx)) paste0("THETA_", s$..combThetaIdx, "_")
  .lag <- .foceiLagSens(
    s,
    stateVars,
    c(.eta, .theta),
    c("rx_pred_", "rx_r_", "rx_pred_f_", "rx_lambda_"),
    statePars = c(.eta, if (length(s$..combThetaStruct)) paste0("THETA_", s$..combThetaStruct, "_"))
  )
  s$..lagSens <- .lag$lines
  s$..lagEta <- .lag
  .lag
}

#' rxode2 text of a derivative that may hold lagged-variable sensitivities
#'
#' @param lag `.foceiLagSens()` list, or `NULL`
#' @param d symengine derivative
#' @param dfe name of the sensitivity it defines (`..._BY_ETA_<n>___`)
#' @return rxode2 text
#' @noRd
.foceiLagTxt <- function(lag, d, dfe) {
  if (is.null(lag)) {
    return(rxode2::rxFromSE(d))
  }
  lag$txt(d, sub("^.*_BY_((ETA|THETA)_[0-9]+_)__$", "\\1", dfe))
}
