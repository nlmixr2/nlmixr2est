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
    .we <- which(.iniDf$name == .th$eta & is.na(.iniDf$ntheta))
    .wt <- which(.iniDf$name == .th$theta & is.na(.iniDf$neta1))
    .pars <- rbind(.pars,
                   data.frame(iov=.nm, theta=.th$theta, eta=.th$eta,
                              est=.iniDf$est[.i], fix=.iniDf$fix[.i],
                              # the IIV variance and the theta's own value, which
                              # the collapsed form needs to build its blocks
                              iiv=.iniDf$est[.we], thetaEst=.iniDf$est[.wt],
                              thetaFix=.iniDf$fix[.wt],
                              stringsAsFactors=FALSE))
  }
  list(occVar=.occ, levels=.lvl, pars=.pars,
       etaNames=setNames(lapply(.pars$iov,
                                function(v) paste0("rx.", v, ".", .lvl)),
                         .pars$iov))
}

#' Reason the collapsed sampler cannot take this model
#'
#' On top of [.saemIovInfo()]'s own scope: the collapsed form shares ONE theta
#' across the occasion columns, and the exact constrained solve for it is the
#' equal-weight average only because the block is compound-symmetric (`1` is an
#' eigenvector of a CS matrix, so `Gamma^-1 1` is proportional to `1`).  That
#' argument needs the intercept-only design, so a mu-referenced covariate on the
#' occasion parameter's theta is out of scope.
#'
#' @param ui rxode2 ui
#' @param info the result of [.saemIovInfo()]
#' @return `NULL`, or a short character reason
#' @noRd
.saemIovCollapsedDecline <- function(ui, info) {
  if (!is.list(info)) return(NULL)
  .cov <- ui$muRefCovariateDataFrame
  if (is.data.frame(.cov) && any(.cov$theta %in% info$pars$theta)) {
    return("collapsed IOV cannot do a covariate on that theta")
  }
  NULL
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
  .uiIovEnv$iovCollapsed <- NULL
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
.saemIovIndicatorSym <- function(term) {
  .t <- .saemIovUnparen(term)
  if (!is.call(.t)) return(NA_character_)
  if (!identical(.t[[1]], quote(`*`))) return(NA_character_)
  if (length(.t) != 3L) return(NA_character_)
  .a <- .saemIovUnparen(.t[[2]])
  .b <- .saemIovUnparen(.t[[3]])
  if (.saemIovIsCmp(.a) && is.name(.b)) return(as.character(.b))
  if (.saemIovIsCmp(.b) && is.name(.a)) return(as.character(.a))
  NA_character_
}

#' The eta an `(occ == level) * <sym>` term ultimately names
#'
#' @param term language object, one term of a `+` chain
#' @param etas names of the model's diagonal etas
#' @param ui rxode2 ui, needed only to follow the collapsed form's variable
#' @return the eta name, or `NA_character_`
#' @noRd
.saemIovIndicatorEta <- function(term, etas, ui = NULL) {
  .sym <- .saemIovIndicatorSym(term)
  if (is.na(.sym)) return(NA_character_)
  if (.sym %in% etas) return(.sym)
  # collapsed form: the indicator multiplies a VARIABLE whose own line carries
  # the occasion's mu-referenced eta (rx.cl.1 <- exp(rx.tcl.1 + rx.eta.cl.1)),
  # so follow it through to that eta
  if (is.null(ui)) return(NA_character_)
  .li <- .saemIovLineFor2(ui, .sym, etas)
  if (is.na(.li)) return(NA_character_)
  .li
}

#' The single diagonal eta on the line defining `sym`, if there is exactly one
#'
#' @param ui rxode2 ui
#' @param sym name of an assigned variable
#' @param etas names of the model's diagonal etas
#' @return the eta name, or `NA_character_`
#' @noRd
.saemIovLineFor2 <- function(ui, sym, etas) {
  for (.l in ui$lstExpr) {
    if (!is.call(.l)) next
    if (!(identical(.l[[1]], quote(`<-`)) || identical(.l[[1]], quote(`=`)))) next
    if (!is.name(.l[[2]]) || as.character(.l[[2]]) != sym) next
    .v <- intersect(all.vars(.l[[3]]), etas)
    if (length(.v) == 1L) return(.v)
    return(NA_character_)
  }
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
.saemIovIndicatorSumEtas <- function(line, etas, ui = NULL) {
  .rhs <- .saemIovRhs(line)
  if (is.null(.rhs) || !is.call(.rhs)) return(NULL)
  if (!identical(.rhs[[1]], quote(`+`))) return(NULL)
  .terms <- .saemIovFlattenPlus(.rhs)
  if (length(.terms) < 2L) return(NULL)
  .e <- vapply(.terms, .saemIovIndicatorEta, character(1), etas=etas, ui=ui,
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
  .groups <- lapply(ui$lstExpr, .saemIovIndicatorSumEtas, etas=.etas, ui=ui)
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
  .m <- control$iovMethod
  if (!(identical(.m, "twoLevel") || identical(.m, "collapsed"))) return(NULL)
  .info <- .saemIovInfo(ui, data)
  # a character is a decline; .uiApplyIov() has already fallen back to the
  # shared rewrite for it, so there is nothing left to do here
  if (!is.list(.info)) return(NULL)
  if (identical(.m, "collapsed")) {
    return(list(ui = .saemIovExpandUiCollapsed(ui, .info)))
  }
  list(ui = .saemIovExpandUi(ui, .info))
}

preProcessHooksAdd(".uiApplyIovTwoLevel", .uiApplyIovTwoLevel)

#' Contract the pooled occasion columns of a covariance matrix
#'
#' `$cov` carries one row/column per phi1 column, so the K per-occasion columns
#' of a two-level parameter appear K times under their internal `om.rx.<iov>.<k>`
#' names.  They estimate ONE variance -- the M-step pins them equal
#' (`poolOmegaGroups`, src/saem.cpp) -- so contract them by averaging, which is
#' the delta method for `Psi = mean(v_1, ..., v_K)`.  Note this contracts the
#' LINEARIZED (unconstrained) covariance rather than computing a constrained
#' Fisher information, so it is an approximation, but a named one: the reported
#' variance is `Var(mean(v_k))`, not `Var(v_1)`.
#'
#' The collapsed row keeps the position of the group's first member so the rest
#' of the matrix keeps its order.
#'
#' @param cv covariance matrix, or `NULL`
#' @param groups named list, user IOV parameter -> its per-occasion eta names
#' @return the contracted matrix
#' @noRd
.saemIovCollapseCov <- function(cv, groups) {
  if (!is.matrix(cv) || is.null(rownames(cv)) || length(groups) == 0L) return(cv)
  .nm <- rownames(cv)
  .grp <- lapply(groups, function(g) {
    .w <- match(paste0("om.", g), .nm)
    .w[!is.na(.w)]
  })
  .grp <- .grp[vapply(.grp, length, integer(1)) >= 2L]
  if (length(.grp) == 0L) return(cv)
  .rows <- .saemIovCovRows(.nm, .grp)
  .a <- matrix(0, nrow = length(.rows), ncol = length(.nm))
  for (.k in seq_along(.rows)) {
    .a[.k, .rows[[.k]]$idx] <- 1 / length(.rows[[.k]]$idx)
  }
  .out <- .a %*% cv %*% t(.a)
  .rn <- vapply(.rows, function(r) r$nm, character(1))
  dimnames(.out) <- list(.rn, .rn)
  .out
}

#' The rows a contracted covariance matrix will have
#'
#' One entry per output row, in input order: the input indices that feed it and
#' the name it takes.  A pooled group contributes a single entry at the position
#' of its first member, so the rest of the matrix keeps its order.
#'
#' @param nm rownames of the input matrix
#' @param grp named list of pooled index vectors, each of length >= 2
#' @return list of `list(idx, nm)`
#' @noRd
.saemIovCovRows <- function(nm, grp) {
  # group id of every row; 0 where the row is not pooled
  .of <- integer(length(nm))
  for (.j in seq_along(grp)) .of[grp[[.j]]] <- .j
  .rows <- list()
  .seen <- integer(0)
  for (.i in seq_along(nm)) {
    .j <- .of[.i]
    if (.j == 0L) {
      .rows[[length(.rows) + 1L]] <- list(idx = .i, nm = nm[.i])
    } else if (!(.j %in% .seen)) {
      .seen <- c(.seen, .j)
      .rows[[length(.rows) + 1L]] <-
        list(idx = grp[[.j]], nm = paste0("om.", names(grp)[.j]))
    }
  }
  .rows
}

#' Substitute symbols in an expression, dropping some from `+` chains
#'
#' @param x language object
#' @param sub named character vector, `old = new`
#' @param drop character vector of symbols to remove from any `+` chain
#' @return the rewritten expression
#' @noRd
.saemIovSubst <- function(x, sub, drop = character(0)) {
  if (is.name(x)) {
    .n <- as.character(x)
    if (.n %in% names(sub)) return(as.name(sub[[.n]]))
    return(x)
  }
  if (!is.call(x)) return(x)
  if (identical(x[[1]], quote(`+`)) && length(x) == 3L) {
    .terms <- .saemIovFlattenPlus(x)
    .keep <- Filter(function(t) !(is.name(t) && as.character(t) %in% drop), .terms)
    .keep <- lapply(.keep, .saemIovSubst, sub = sub, drop = drop)
    if (length(.keep) == 0L) return(0)
    .out <- .keep[[1]]
    for (.i in seq_along(.keep)[-1]) .out <- call("+", .out, .keep[[.i]])
    return(.out)
  }
  as.call(lapply(as.list(x), .saemIovSubst, sub = sub, drop = drop))
}

#' The model line that defines a parameter from its mu-referenced eta
#'
#' @param ui rxode2 ui
#' @param eta name of the mu-referenced eta
#' @return the index of the line, or `NA_integer_`
#' @noRd
.saemIovLineFor <- function(ui, eta) {
  .w <- which(vapply(ui$lstExpr, function(l) {
    .rhs <- .saemIovRhs(l)
    !is.null(.rhs) && eta %in% all.vars(.rhs)
  }, logical(1)))
  if (length(.w) != 1L) return(NA_integer_)
  .w
}

#' Build one IOV parameter's collapsed lines, thetas and eta block
#'
#' @param ui rxode2 ui being rewritten
#' @param lst its model lines
#' @param info the `.saemIovInfo()` result
#' @param i row of `info$pars` to build
#' @param lvl occasion levels
#' @param thetaTpl,etaTpl template `iniDf` rows
#' @return list of `lines`, `dropTheta`, `dropEta`, `theta` and `etaBlock`
#' @noRd
.saemIovCollapsedOne <- function(ui, lst, info, i, lvl, thetaTpl, etaTpl) {
  .th <- info$pars$theta[i]
  .et <- info$pars$eta[i]
  .iv <- info$pars$iov[i]
  .li <- .saemIovLineFor(ui, .et)
  if (is.na(.li)) {
    stop("cannot find the model line for '", .et, "'", call. = FALSE)
  }
  .line <- lst[[.li]]
  .lhs <- as.character(.line[[2]])
  .rhs <- .line[[3]]
  .colNames <- paste0("rx.", .lhs, ".", lvl)
  .thNames <- paste0("rx.", .th, ".", lvl)
  .etNames <- paste0("rx.", .et, ".", lvl)
  .lines <- lapply(seq_along(lvl), function(.k) {
    call("<-", as.name(.colNames[.k]),
         .saemIovSubst(.rhs,
                       sub = stats::setNames(c(.thNames[.k], .etNames[.k]),
                                             c(.th, .et)),
                       drop = .iv))
  })
  .lines[[length(.lines) + 1L]] <-
    str2lang(paste0(.lhs, " <- ",
                    paste(paste0(.colNames, "*(", info$occVar, " == ", lvl, ")"),
                          collapse = " + ")))
  .theta <- do.call(rbind, lapply(seq_along(lvl), function(.k) {
    .cur <- thetaTpl
    .cur$name <- .thNames[.k]
    .cur$est <- info$pars$thetaEst[i]
    .cur$fix <- info$pars$thetaFix[i]
    .cur$label <- NA_character_
    .cur$lower <- -Inf
    .cur$upper <- Inf
    .cur$condition <- NA_character_
    .cur$err <- NA_character_
    .cur
  }))
  list(lines = .lines, dropTheta = .th, dropEta = c(.et, .iv), theta = .theta,
       # the joint covariance of phi_i over occasions: Omega + Psi on the
       # diagonal, Omega off it
       etaBlock = list(names = .etNames, tpl = etaTpl,
                       diag = info$pars$iiv[i] + info$pars$est[i],
                       off = info$pars$iiv[i], fix = info$pars$fix[i]))
}

#' Rewrite a ui for the collapsed (Panhard & Samson) sampler
#'
#' Where [.saemIovExpandUi()] keeps `b_i` and `c_ik` as separate columns, this
#' carries `phi_ik = mu + b_i + c_ik` in ONE column per occasion:
#'
#' \preformatted{  rx.cl.1 <- exp(rx.tcl.1 + rx.eta.cl.1)
#'   rx.cl.2 <- exp(rx.tcl.2 + rx.eta.cl.2)
#'   cl <- rx.cl.1*(occ == 1) + rx.cl.2*(occ == 2)}
#'
#' One mu-reference per LINE is load bearing: rxode2 detects only the first
#' additive `theta + eta` group in a line, so putting both occasions on one line
#' leaves the second eta non-mu-referenced.  Split this way, every occasion
#' column is an ordinary mu-referenced parameter -- its own phi column, its own
#' omega entry, mean estimated rather than pinned.
#'
#' The etas are declared as one CORRELATED block so `covstruct` carries the
#' off-diagonals: the block is `Omega + Psi` on the diagonal and `Omega` off it,
#' which is the joint covariance of `phi_i` across occasions.  The equality
#' constraints that make it compound-symmetric (and tie the K thetas to the
#' single `mu` the user declared) are imposed in the M-step.
#'
#' @param ui rxode2 ui
#' @param info the result of [.saemIovInfo()]
#' @return the rewritten ui
#' @noRd
.saemIovExpandUiCollapsed <- function(ui, info) {
  .ui <- rxode2::rxUiDecompress(ui)
  # rxUiDecompress is NOT a copy, and this path (unlike the two-level one) does
  # no rxRename to force a fresh object -- so the assign()s below mutate `ui`
  # itself.  Snapshot what the finalizer needs as VALUES first.
  info$origIniDf <- .ui$iniDf
  info$origLstExpr <- .ui$lstExpr
  .iniDf <- .ui$iniDf
  .lst <- .ui$lstExpr
  .lvl <- info$levels
  .thetas <- .iniDf[is.na(.iniDf$neta1), , drop = FALSE]
  .etas <- .iniDf[is.na(.iniDf$ntheta), , drop = FALSE]
  .thetaTpl <- .thetas[1, ]
  .etaTpl <- .etas[1, ]
  if (any(names(.thetaTpl) == "prior")) .thetaTpl$prior <- NA_character_
  if (any(names(.etaTpl) == "prior")) .etaTpl$prior <- NA_character_

  .built <- lapply(seq_along(info$pars$iov), function(.i) {
    .saemIovCollapsedOne(.ui, .lst, info, .i, .lvl, .thetaTpl, .etaTpl)
  })
  .newLines <- do.call(`c`, lapply(.built, function(b) b$lines))
  .dropTheta <- vapply(.built, function(b) b$dropTheta, character(1))
  .dropEta <- unlist(lapply(.built, function(b) b$dropEta), use.names = FALSE)
  .addTheta <- do.call(rbind, lapply(.built, function(b) b$theta))
  .addEtaBlocks <- lapply(.built, function(b) b$etaBlock)

  .lst <- .lst[-vapply(info$pars$eta, function(e) .saemIovLineFor(.ui, e),
                       integer(1))]
  .thetas <- .thetas[!(.thetas$name %in% .dropTheta), , drop = FALSE]
  .etas <- .etas[!(.etas$name %in% .dropEta), , drop = FALSE]
  .thetas <- rbind(.thetas, .addTheta)
  .thetas$ntheta <- seq_len(nrow(.thetas))
  .maxEta <- 0L
  if (nrow(.etas) > 0L) {
    .etas$neta1 <- as.integer(factor(.etas$neta1, levels = sort(unique(.etas$neta1))))
    .etas$neta2 <- as.integer(factor(.etas$neta2, levels = sort(unique(.etas$neta2))))
    .maxEta <- max(.etas$neta1)
  }
  for (.b in .addEtaBlocks) {
    .idx <- .maxEta + seq_along(.b$names)
    .maxEta <- max(.idx)
    for (.a in seq_along(.b$names)) {
      for (.c in seq_len(.a)) {
        .cur <- .b$tpl
        .cur$neta1 <- .idx[.a]
        .cur$neta2 <- .idx[.c]
        .cur$est <- if (.a == .c) .b$diag else .b$off
        .cur$fix <- .b$fix
        .cur$condition <- "id"
        .cur$label <- NA_character_
        .cur$name <- if (.a == .c) .b$names[.a] else {
          paste0("(", .b$names[.c], ",", .b$names[.a], ")")
        }
        .etas <- rbind(.etas, .cur)
      }
    }
  }
  assign("iniDf", rbind(.thetas, .etas), envir = .ui)
  assign("lstExpr", c(.newLines, .lst), envir = .ui)
  .uiIovEnv$ui <- ui
  .uiIovEnv$iovVars <- info$pars$iov
  .uiIovEnv$iovDrop <- unlist(lapply(.addEtaBlocks, function(b) b$names),
                              use.names = FALSE)
  .uiIovEnv$lines <- .newLines
  .uiIovEnv$iovTwoLevel <- NULL
  .uiIovEnv$iovCollapsed <- info
  .uiIovEnv$muModel <- NULL
  .uiIovEnv$iovRename <- NULL
  rxode2::rxUiDecompress(suppressWarnings(suppressMessages(.ui$fun())))
}

#' Recover Omega, Psi and mu from a collapsed fit's compound-symmetric block
#'
#' The fitted block is `Omega + Psi` on the diagonal and `Omega` off it, and the
#' M-step holds it compound-symmetric, so reading any one entry of each kind is
#' enough; the mean over the group is used so a fit stopped mid-annealing still
#' gives a sensible answer.
#'
#' @param om fitted omega matrix over the collapsed etas
#' @param th named fitted theta vector
#' @param info the `.saemIovInfo()` result stashed at expansion time
#' @return list of `theta`, `omega` and `psi`, each named by IOV parameter
#' @noRd
.saemIovCollapsedParts <- function(om, th, info) {
  .lvl <- info$levels
  .theta <- .omega <- .psi <- stats::setNames(rep(NA_real_, nrow(info$pars)),
                                              info$pars$iov)
  for (.i in seq_len(nrow(info$pars))) {
    .en <- paste0("rx.", info$pars$eta[.i], ".", .lvl)
    .tn <- paste0("rx.", info$pars$theta[.i], ".", .lvl)
    .en <- .en[.en %in% rownames(om)]
    .tn <- .tn[.tn %in% names(th)]
    if (length(.en) < 2L || length(.tn) < 1L) next
    .d <- mean(diag(om[.en, .en, drop = FALSE]))
    .o <- om[.en, .en, drop = FALSE]
    .o <- mean(.o[upper.tri(.o)])
    .theta[.i] <- mean(th[.tn])
    .omega[.i] <- .o
    .psi[.i] <- .d - .o
  }
  list(theta = .theta, omega = .omega, psi = .psi)
}

#' Restore the user's model after a collapsed (Panhard & Samson) fit
#'
#' Registered as a post-final hook.  The collapsed expansion replaced the user's
#' line outright, so this rebuilds from the ORIGINAL ui rather than unpicking the
#' rewritten one: it takes the pre-rewrite `iniDf`/`model` and writes the fitted
#' values back into it -- the shared `mu` from the pooled thetas, the
#' between-subject variance from the block's off-diagonal, and the
#' inter-occasion variance from diagonal minus off-diagonal.
#'
#' @param ret fit object
#' @return the fit, with the user's parameterization restored
#' @noRd
.saemIovFinalizeCollapsed <- function(ret) {
  .info <- .uiIovEnv$iovCollapsed
  if (is.null(.info) || is.null(.uiIovEnv$ui)) return(ret)
  if (is.environment(ret$env) && !is.null(ret$ui)) {
    .fit <- ret$env$ui
    .om <- .fit$omega
    if (is.list(.om)) .om <- .om$id
    .parts <- .saemIovCollapsedParts(.om, ret$env$fixef, .info)
    .orig <- rxode2::rxUiDecompress(.uiIovEnv$ui)
    .ini <- .info$origIniDf
    # carry every estimate that survived the rewrite unchanged
    .fitIni <- .fit$iniDf
    .keep <- match(.ini$name, .fitIni$name)
    .ok <- !is.na(.keep)
    .ini$est[.ok] <- .fitIni$est[.keep[.ok]]
    # then the three the collapsed block owns
    for (.i in seq_len(nrow(.info$pars))) {
      .p <- .info$pars[.i, ]
      .ini$est[.ini$name == .p$theta & is.na(.ini$neta1)] <- .parts$theta[[.i]]
      .ini$est[.ini$name == .p$eta & is.na(.ini$ntheta)] <- .parts$omega[[.i]]
      .ini$est[.ini$name == .p$iov & is.na(.ini$ntheta)] <- .parts$psi[[.i]]
    }
    .newIni <- as.expression(lotri::as.lotri(.ini))
    .newIni[[1]] <- quote(`ini`)
    # .getUiFunFromIniAndModel() hands back a model FUNCTION, not a ui -- the
    # shared finalizer only ends up with a ui because its rxRename() call
    # evaluates one.  There is no rename here, so build the ui explicitly.
    .uiFun <- .getUiFunFromIniAndModel(.orig, .newIni,
                                       rxode2::as.model(.info$origLstExpr))
    .ui <- rxode2::rxUiDecompress(
      suppressWarnings(suppressMessages(.uiFun())))
    assign("ui", .ui, envir = ret$env)
    assign("iniDf0", .info$origIniDf, envir = ret$env)
    assign("omega", .ui$omega, envir = ret$env)
    .fx <- ret$env$fixef
    .drop <- unlist(lapply(seq_len(nrow(.info$pars)), function(i) {
      paste0("rx.", .info$pars$theta[i], ".", .info$levels)
    }), use.names = FALSE)
    .fx <- .fx[!(names(.fx) %in% .drop)]
    for (.i in seq_len(nrow(.info$pars))) {
      .fx[[.info$pars$theta[.i]]] <- .parts$theta[[.i]]
    }
    assign("fixef", .fx, envir = ret$env)

    # The collapsed columns hold phi_ik = mu + b_i + c_ik jointly, so the
    # between-subject and inter-occasion deviations come out of them by the
    # obvious decomposition: b_i is the subject's mean over occasions and c_ik
    # what is left.  That is also the split the CS block assumes.
    .split <- .saemIovSplitRanef(ret$env$ranef, .info)
    if (!is.null(.split)) {
      assign("ranef", .split$ranef, envir = ret$env)
      if (!is.null(.split$iov)) assign("iov", .split$iov, envir = ret$env)
    }
  }
  if (inherits(ret, "data.frame")) {
    .w <- which(grepl("^rx[.]", names(ret)))
    if (length(.w) > 0L) {
      .cls <- class(ret)
      class(ret) <- "data.frame"
      ret <- ret[, -.w]
      class(ret) <- .cls
    }
  }
  ret
}

postFinalObjectHooksAdd(".saemIovFinalizeCollapsed", .saemIovFinalizeCollapsed)

#' Split the collapsed etas into between-subject and inter-occasion parts
#'
#' The collapsed columns hold `phi_ik = mu + b_i + c_ik` jointly, so `b_i` is the
#' subject's mean over occasions and `c_ik` is what is left -- which is also the
#' split the compound-symmetric block assumes.
#'
#' @param re the fit's `ranef` data frame
#' @param info the `.saemIovInfo()` result
#' @return `NULL`, or a list of the rewritten `ranef` and the `iov` tables
#' @noRd
.saemIovSplitRanef <- function(re, info) {
  if (!is.data.frame(re)) return(NULL)
  .tab <- list()
  for (.i in seq_len(nrow(info$pars))) {
    .en <- paste0("rx.", info$pars$eta[.i], ".", info$levels)
    .en <- .en[.en %in% names(re)]
    if (length(.en) < 2L) next
    .m <- as.matrix(re[, .en, drop = FALSE])
    .b <- rowMeans(.m)
    re[[info$pars$eta[.i]]] <- .b
    .id <- if ("ID" %in% names(re)) re$ID else seq_len(nrow(re))
    .one <- data.frame(ID = rep(.id, times = length(info$levels)),
                       occ = rep(info$levels, each = nrow(.m)),
                       dev = as.vector(.m - .b))
    names(.one) <- c("ID", info$occVar, info$pars$iov[.i])
    .tab[[info$pars$iov[.i]]] <- .one
    re <- re[, !(names(re) %in% .en), drop = FALSE]
  }
  if (length(.tab) == 0L) return(list(ranef = re, iov = NULL))
  # one table per occasion variable, matching the shared rewrite's shape
  .one <- .tab[[1]]
  for (.n in names(.tab)[-1]) .one[[.n]] <- .tab[[.n]][[.n]]
  .one <- .one[order(.one[[1]], .one[[2]]), , drop = FALSE]
  rownames(.one) <- NULL
  list(ranef = re, iov = stats::setNames(list(.one), info$occVar))
}
