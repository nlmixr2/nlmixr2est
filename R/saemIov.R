# Two-level (inter-occasion) variability for saem.
#
# saemControl(iovMethod = "twoLevel") turns the shared IOV pre-processing
# rewrite off (see the "iov" attribute on nlmixr2Est.saem) and hands the
# occasion term to the code here instead.
#
# The shared rewrite (.uiApplyIov(), R/iov.R) carries the occasion magnitude as
# a population parameter multiplying per-occasion unit-variance etas.  That
# makes the magnitude non-mu-referenced, so saem estimates a VARIANCE through
# its fixed-effect-only (phi0) path -- a stochastic sampled mean over draws
# whose pseudo-variance is deliberately annealed, then a bounded direct
# optimization -- while every other variance component in the algorithm gets a
# closed-form M-step.  On theo_md that estimate collapses toward zero.
#
# Panhard X, Samson A (2009), "Extension of the SAEM algorithm for nonlinear
# mixed effects models with two levels of random effects", Biostatistics 10(1),
# 121-135, write the same model with the occasion term as a second variance
# component,
#
#   phi_ik = mu + b_i + c_ik,  b_i ~ N(0, Omega),  c_ik ~ N(0, Psi)
#
# which keeps the parameter mu-referenced and makes Psi an ordinary variance.
# The expansion here writes that model out for the kernel: one zero-mean eta
# per occasion level, entering ADDITIVELY behind occasion indicators, with the
# K variances constrained equal in the M-step (poolOmegaGroups, src/saem.cpp).

#' Flatten a `+` chain into its terms
#'
#' @param x language object
#' @return list of the chain's terms (the object itself when it is not a `+`)
#' @noRd
.saemIovFlattenPlus <- function(x) {
  if (is.call(x) && identical(x[[1]], quote(`+`)) && length(x) == 3L) {
    c(.saemIovFlattenPlus(x[[2]]), .saemIovFlattenPlus(x[[3]]))
  } else {
    list(x)
  }
}

#' Every maximal `+` chain in an expression
#'
#' @param x language object
#' @return list of chains, each a list of terms
#' @noRd
.saemIovPlusChains <- function(x) {
  if (!is.call(x)) return(list())
  if (identical(x[[1]], quote(`+`)) && length(x) == 3L) {
    return(list(.saemIovFlattenPlus(x)))
  }
  .sub <- as.list(x)[-1]
  if (length(.sub) == 0L) return(list())
  do.call(`c`, c(list(list()), lapply(.sub, .saemIovPlusChains)))
}

#' The mu-referenced theta an occasion eta rides on
#'
#' The occasion eta has to sit in the same additive position as a
#' mu-referenced eta -- that additive position is what makes the parameter
#' `mu + b_i + c_ik` rather than something the closed-form M-step cannot reach.
#'
#' @param ui rxode2 ui
#' @param iovName name of the occasion eta
#' @return one-row data frame of `theta`/`eta`, or `NULL` when the occasion eta
#'   does not share an additive position with exactly one mu-referenced eta
#' @noRd
.saemIovThetaFor <- function(ui, iovName) {
  .muRef <- ui$muRefDataFrame
  .chains <- do.call(`c`, c(list(list()),
                            lapply(ui$lstExpr, .saemIovPlusChains)))
  .hit <- NULL
  for (.ch in .chains) {
    .nm <- vapply(.ch, function(t) {
      if (is.name(t)) as.character(t) else ""
    }, character(1), USE.NAMES=FALSE)
    if (!(iovName %in% .nm)) next
    .w <- which(.nm %in% .muRef$eta)
    if (length(.w) != 1L) next
    .e <- .nm[.w]
    .hit <- unique(rbind(.hit, .muRef[.muRef$eta == .e, c("theta", "eta")]))
  }
  if (is.null(.hit) || nrow(.hit) != 1L) return(NULL)
  .hit
}

#' Describe a model's two-level (IOV) structure for saem
#'
#' @param ui rxode2 ui, with the `iov.x ~ v | occ` rows still in `iniDf`
#' @param data data set the fit will use; the occasion levels come from it
#' @return `NULL` when the model has no IOV.  Otherwise a list with `occVar`,
#'   `levels`, `pars` (a data frame of `iov`/`theta`/`eta`/`est`/`fix`) and
#'   `etaNames` (a list, per IOV parameter, of the per-occasion eta names).  A
#'   model outside the scope of this handling returns a character string saying
#'   why, so the caller can fall back to `iovMethod = "theta"`.
#' @noRd
.saemIovInfo <- function(ui, data) {
  .iniDf <- ui$iniDf
  .w <- which(!is.na(.iniDf$condition) &
                .iniDf$condition != "id" &
                is.na(.iniDf$err))
  if (length(.w) == 0L) return(NULL)
  .occ <- unique(.iniDf$condition[.w])
  if (length(.occ) != 1L) {
    return("two-level IOV needs one occasion variable")
  }
  .off <- .w[which(!is.na(.iniDf$neta1[.w]) &
                     .iniDf$neta1[.w] != .iniDf$neta2[.w])]
  if (length(.off) > 0L) {
    return("two-level IOV cannot do correlated occasions")
  }
  if (is.null(data[[.occ]])) {
    stop("IOV variable '", .occ, "' is not present in the data ", call.=FALSE)
  }
  .lvl <- sort(unique(data[[.occ]]))
  if (!is.numeric(.lvl)) {
    stop("IOV variable '", .occ, "' must be numeric", call.=FALSE)
  }
  .pars <- NULL
  for (.i in .w) {
    .nm <- .iniDf$name[.i]
    .th <- .saemIovThetaFor(ui, .nm)
    if (is.null(.th)) {
      return("two-level IOV needs a mu-referenced parameter")
    }
    .pars <- rbind(.pars,
                   data.frame(iov=.nm, theta=.th$theta, eta=.th$eta,
                              est=.iniDf$est[.i], fix=.iniDf$fix[.i],
                              stringsAsFactors=FALSE))
  }
  list(occVar=.occ, levels=.lvl, pars=.pars,
       etaNames=setNames(lapply(.pars$iov,
                                function(v) paste0("rx.", v, ".", .lvl)),
                         .pars$iov))
}

#' Rewrite a ui so saem sees the occasion term as a second variance component
#'
#' Each `iov.x ~ v | occ` row becomes one zero-mean eta per observed occasion
#' level, all sharing the declared variance, combined behind occasion
#' indicators on a line of their own:
#'
#' \preformatted{  rx.iov.cl <- (occ == 1)*rx.iov.cl.1 + (occ == 2)*rx.iov.cl.2}
#'
#' Keeping the combination on its own line matters: rxode2 refuses
#' `theta + eta1 + eta2` in one additive position ("mu-ref err: currently do
#' not theta + eta1 + eta2"), but accepts `theta + eta + <variable>`, which
#' leaves `eta.cl` mu-referenced to `tcl` and puts the occasion etas in
#' `nonMuEtas` -- exactly the phi layout the kernel needs (their phi means are
#' pinned at 0 and their variances are ordinary omega entries).
#'
#' @param ui rxode2 ui
#' @param info the result of [.saemIovInfo()]
#' @return the rewritten ui
#' @noRd
.saemIovExpandUi <- function(ui, info) {
  .ui <- rxode2::rxUiDecompress(ui)
  .nm <- info$pars$iov
  # free the user's symbol: it names an eta today and has to name the combined
  # per-record value instead
  .ui <- suppressWarnings(
    eval(str2lang(paste0("rxode2::rxRename(.ui, ",
                         paste(paste0("rx.", .nm, "=", .nm), collapse=", "),
                         ")"))))
  .iniDf <- .ui$iniDf
  .thetas <- .iniDf[is.na(.iniDf$neta1), , drop=FALSE]
  .etas <- .iniDf[is.na(.iniDf$ntheta), , drop=FALSE]
  .template <- .etas[1, ]
  if (any(names(.template) == "prior")) .template$prior <- NA_character_
  .template$label <- NA_character_
  # drop the occasion rows (renamed above) and renumber what is left
  .drop <- paste0("rx.", .nm)
  .etas <- .etas[!(.etas$name %in% .drop), , drop=FALSE]
  .maxEta <- 0L
  if (nrow(.etas) > 0L) {
    .etas$neta1 <- as.integer(factor(.etas$neta1, levels=sort(unique(.etas$neta1))))
    .etas$neta2 <- as.integer(factor(.etas$neta2, levels=sort(unique(.etas$neta2))))
    .maxEta <- max(.etas$neta1)
  }
  for (.i in seq_along(.nm)) {
    for (.l in info$levels) {
      .cur <- .template
      .cur$name <- paste0("rx.", .nm[.i], ".", .l)
      .cur$label <- paste0(.nm[.i], "(", info$occVar, "==", .l, ")")
      .cur$est <- info$pars$est[.i]
      .cur$fix <- info$pars$fix[.i]
      .cur$condition <- "id"
      .maxEta <- .cur$neta1 <- .cur$neta2 <- .maxEta + 1L
      .etas <- rbind(.etas, .cur)
    }
  }
  .lines <- c(
    lapply(seq_along(.nm), function(.i) {
      str2lang(paste0("rx.", .nm[.i], " <- ",
                      paste(paste0("(", info$occVar, " == ", info$levels, ")*",
                                   info$etaNames[[.i]]),
                            collapse=" + ")))
    }),
    # the realized per-record value, kept as an output column so the fit's data
    # frame carries `iov.x` the way the shared rewrite's does
    # (.uiFinalizeIov() renames `iov.x.rx` back to `iov.x`)
    lapply(.nm, function(v) str2lang(paste0(v, ".rx <- rx.", v))))
  assign("iniDf", rbind(.thetas, .etas), envir=.ui)
  assign("lstExpr", c(.lines, .ui$lstExpr), envir=.ui)
  # .uiFinalizeIov() (R/iov.R) undoes the rewrite after the fit -- restoring the
  # `iov.x ~ v | occ` row, splitting $omega into $id/$occ, building $iov and the
  # shrinkage table.  All of that is shared; hand it the same state the shared
  # rewrite leaves behind, plus `iovTwoLevel` to say the variance comes off the
  # pooled occasion etas rather than a magnitude theta.
  .uiIovEnv$ui <- ui
  .uiIovEnv$iovVars <- .nm
  .uiIovEnv$iovDrop <- unlist(info$etaNames, use.names=FALSE)
  .uiIovEnv$lines <- .lines
  .uiIovEnv$iovTwoLevel <- info$etaNames
  .uiIovEnv$muModel <- NULL
  .uiIovEnv$iovRename <-
    str2lang(paste0("rxode2::rxRename(.ui, ",
                    paste(paste0(.nm, "=", "rx.", .nm), collapse=", "), ")"))
  rxode2::rxUiDecompress(suppressWarnings(suppressMessages(.ui$fun())))
}

#' Strip `(` wrappers from an expression
#'
#' `(occ == 1)` parses as a call to `(` around the comparison, so a structural
#' match has to look through it.
#'
#' @param x language object
#' @return `x` with any enclosing `(` calls removed
#' @noRd
.saemIovUnparen <- function(x) {
  while (is.call(x) && identical(x[[1]], quote(`(`)) && length(x) == 2L) {
    x <- x[[2]]
  }
  x
}

#' Is this expression a `==` comparison?
#'
#' @param o language object
#' @return logical
#' @noRd
.saemIovIsCmp <- function(o) is.call(o) && identical(o[[1]], quote(`==`))

#' Is this expression one of the model's diagonal etas?
#'
#' @param o language object
#' @param etas names of the model's diagonal etas
#' @return logical
#' @noRd
.saemIovIsEta <- function(o, etas) is.name(o) && as.character(o) %in% etas

#' The eta of an `(occ == level) * eta` term
#'
#' @param term language object, one term of a `+` chain
#' @param etas names of the model's diagonal etas
#' @return the eta name, or `NA_character_` when the term is not of that shape
#' @noRd
.saemIovIndicatorEta <- function(term, etas) {
  .t <- .saemIovUnparen(term)
  if (!is.call(.t)) return(NA_character_)
  if (!identical(.t[[1]], quote(`*`))) return(NA_character_)
  if (length(.t) != 3L) return(NA_character_)
  .a <- .saemIovUnparen(.t[[2]])
  .b <- .saemIovUnparen(.t[[3]])
  if (.saemIovIsCmp(.a) && .saemIovIsEta(.b, etas)) return(as.character(.b))
  if (.saemIovIsCmp(.b) && .saemIovIsEta(.a, etas)) return(as.character(.a))
  NA_character_
}

#' The right-hand side of an assignment, or NULL
#'
#' @param line language object, one model line
#' @return the right-hand side, or `NULL` when the line is not an assignment
#' @noRd
.saemIovRhs <- function(line) {
  if (!is.call(line)) return(NULL)
  if (!(identical(line[[1]], quote(`<-`)) || identical(line[[1]], quote(`=`)))) {
    return(NULL)
  }
  .saemIovUnparen(line[[3]])
}

#' The etas of a line that is exactly a sum of `(occ == level) * eta` terms
#'
#' @param line one model line
#' @param etas names of the model's diagonal etas
#' @return character vector of the etas, or `NULL` when the line is not that shape
#' @noRd
.saemIovIndicatorSumEtas <- function(line, etas) {
  .rhs <- .saemIovRhs(line)
  if (is.null(.rhs) || !is.call(.rhs)) return(NULL)
  if (!identical(.rhs[[1]], quote(`+`))) return(NULL)
  .terms <- .saemIovFlattenPlus(.rhs)
  if (length(.terms) < 2L) return(NULL)
  .e <- vapply(.terms, .saemIovIndicatorEta, character(1), etas=etas,
               USE.NAMES=FALSE)
  if (anyNA(.e) || anyDuplicated(.e) > 0L) return(NULL)
  .e
}

#' Occasion-eta pool groups read back out of the model
#'
#' The pooling constraint is recovered from the model text rather than carried
#' alongside it: a line whose right-hand side is exactly a sum of
#' `(occ == level) * eta` terms says those etas are one occasion parameter
#' observed at different levels, so they estimate ONE variance (`Psi`).  The
#' legacy rewrite's line is a PRODUCT (`magnitude * (sum)`), so it does not
#' match and is left alone.
#'
#' @param ui rxode2 ui, already expanded by [.saemIovExpandUi()]
#' @return list of character vectors, one per pool group; empty when there is
#'   nothing to pool
#' @noRd
.saemIovPoolFromModel <- function(ui) {
  .idf <- ui$iniDf
  .etas <- .idf[!is.na(.idf$neta1) & .idf$neta1 == .idf$neta2, "name"]
  .groups <- lapply(ui$lstExpr, .saemIovIndicatorSumEtas, etas=.etas)
  .groups[!vapply(.groups, is.null, logical(1))]
}

# Pool-group id for every phi1 column, in Gamma2_phi1 (saemEtaNames) order.
# Columns sharing a non-zero id estimate one variance; 0 means the column has
# its own.
#' @export
rxUiGet.saemOmegaPool <- function(x, ...) {
  .ui <- x[[1]]
  .names <- rxUiGet.saemEtaNames(x, ...)
  .ret <- rep(0L, length(.names))
  .groups <- .saemIovPoolFromModel(.ui)
  .g <- 0L
  for (.grp in .groups) {
    .w <- which(.names %in% .grp)
    if (length(.w) < 2L) next
    .g <- .g + 1L
    .ret[.w] <- .g
  }
  .ret
}
attr(rxUiGet.saemOmegaPool, "rstudio") <- c(0L, 0L)

#' Two-level IOV expansion, as a pre-processing hook
#'
#' Runs at the same point in the pipeline as the shared rewrite it replaces --
#' which matters: `.preProcessBoundedTransform` runs LAST and rewrites a bounded
#' theta into an lhs expression (`tcl <- 4.6 - exp(rxBoundedTr.tcl)`), which
#' takes that theta out of `muRefDataFrame`.  Expanding after that would see a
#' parameter that is no longer mu-referenced and decline, while the shared
#' rewrite had already stood down -- leaving the occasion term unhandled.
#'
#' @param ui rxode2 ui
#' @param est estimation routine name
#' @param data data set for the fit
#' @param control control object
#' @return `NULL`, or a list with the rewritten `ui`
#' @noRd
#' @author Matthew L. Fidler
.uiApplyIovTwoLevel <- function(ui, est, data, control) {
  if (!identical(est, "saem")) return(NULL)
  if (!identical(control$iovMethod, "twoLevel")) return(NULL)
  .info <- .saemIovInfo(ui, data)
  # a character is a decline; .uiApplyIov() has already fallen back to the
  # shared rewrite for it, so there is nothing left to do here
  if (!is.list(.info)) return(NULL)
  list(ui = .saemIovExpandUi(ui, .info))
}

preProcessHooksAdd(".uiApplyIovTwoLevel", .uiApplyIovTwoLevel)
