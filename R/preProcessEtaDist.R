## Declared non-Gaussian random effect (eta) distributions.
##
## `lotri` parses `dist(eta.cl) ~ dgamma(...)` and `rxode2` turns it into
## a model (`rxEtaDistExpand()`): a latent standard normal with a FIXED
## identity omega, plus `phiU()` + the family's inverse CDF, plus
## unconstrained `rxCor.*` thetas carrying the Gaussian copula's
## correlation.  Two things are left for nlmixr2est.
##
## 1. Run that expansion before estimation, which is all the support most
##    methods need: what they are handed afterwards is an ordinary model
##    with a fixed identity omega.
##
## 2. Refuse the methods for which it is NOT ordinary.  A declared
##    distribution that a method quietly ignored would fit a different
##    model than the one written, with nothing to say so -- the same
##    reasoning, and the same attribute-on-the-S3-method mechanism, as the
##    prior gate in R/priors.R.
##
##      attr(nlmixr2Est.myMethod, "etaDist") <- TRUE
##
## The methods that support it are the FOCEi family (whose inner problem
## needs only d(eta)/d(latent), which rxode2 differentiates exactly
## through `phiU()` and the inverse CDF), SAEM (a declared eta has no
## `theta + eta` form, so it lands in the already-exercised `nonMuEtas`
## path and is still Gibbs/Metropolis sampled with the same sample
## covariance update), and simulation.  Refused: `npag`/`npb`, which model
## the random effect distribution nonparametrically, so a declared one
## contradicts them outright; `nlme` and `nls`, which are Gaussian by
## construction; and `vae`/`emvi`/`fbvi`, whose ELBO hardcodes the normal
## family.

#' The `"etaDist"` attribute of the dispatched estimation method
#'
#' Read from the `nlmixr2Est.<method>` S3 method, so a method registered
#' by another package can declare support without editing this file.  The
#' attribute may be `TRUE`/`FALSE`, the string `"native"`, or a
#' `function(control)` returning one of those.
#'
#' `"native"` means the method translates the declaration ITSELF and must
#' see it unexpanded -- babelmixr2's `est="nonmem"` writes Bauer's own
#' `$ABBR FUNCTION GAMMACDFINV` control stream, which reads nothing like
#' the expansion and could not be recovered from it.
#'
#' @param est estimation routine name
#' @param control control object
#' @return `TRUE`, `FALSE` or `"native"`
#' @noRd
#' @author Matthew L. Fidler
.etaDistMethodAttr <- function(est, control=NULL) {
  if (!is.character(est) || length(est) != 1L) return(FALSE)
  .v <- as.character(utils::methods("nlmixr2Est"))
  if (!(paste0("nlmixr2Est.", est) %in% .v)) return(FALSE)
  .a <- attr(utils::getS3method("nlmixr2Est", est), "etaDist")
  if (is.null(.a)) return(FALSE)
  if (is.function(.a)) .a <- .a(control)
  if (identical(.a, "native")) return("native")
  isTRUE(.a)
}

#' Does the dispatched estimation method support a declared eta distribution?
#'
#' @param est estimation routine name
#' @param control control object
#' @return boolean
#' @noRd
#' @author Matthew L. Fidler
.isEtaDistMethod <- function(est, control=NULL) {
  !isFALSE(.etaDistMethodAttr(est, control))
}

#' The error a method that cannot use a declared distribution gets
#'
#' One message, raised from two places: the pre-processing hook (early,
#' before any work is done) and the gate in `nlmixr2Est()` (the backstop).
#'
#' @param d declaring random effects, as `rxUiEtaDists()` returns them
#' @param est estimation routine name
#' @param control control object
#' @return nothing, called for the error
#' @noRd
#' @author Matthew L. Fidler
.etaDistRefuse <- function(d, est, control) {
  if (nrow(d) == 0L) return(invisible())
  if (!is.character(est) || length(est) != 1L) return(invisible())
  if (.isEtaDistMethod(est, control)) return(invisible())
  stop("est=\"", est, "\" cannot use the declared non-normal random effect ",
       "distribution(s) on '", paste(d$name, collapse="', '"), "'",
       call.=FALSE)
}

#' Refuse a declared eta distribution the dispatched method cannot use
#'
#' @param env nlmixr2 estimation environment
#' @return nothing, called for the error
#' @noRd
#' @author Matthew L. Fidler
.nlmixr2AssertEtaDist <- function(env) {
  .ui <- try(get("ui", envir=env), silent=TRUE)
  if (inherits(.ui, "try-error") || is.null(.ui)) return(invisible())
  .d <- .rxUiEtaDists(.ui)
  if (nrow(.d) == 0L) return(invisible())
  .est <- class(env)[1]
  .control <- if (exists("control", envir=env)) get("control", envir=env) else NULL
  .etaDistRefuse(.d, .est, .control)
}

#' `rxode2::rxUiEtaDists()` when the installed rxode2 has it
#'
#' Looked up rather than called directly so that an older rxode2 -- which
#' cannot have produced a declaration in the first place -- degrades to
#' "no declarations" instead of erroring.
#'
#' @param ui rxode2 ui
#' @return the declaring random effects, zero rows when there are none
#' @noRd
#' @author Matthew L. Fidler
.rxUiEtaDists <- function(ui) {
  .ns <- asNamespace("rxode2")
  if (!exists("rxUiEtaDists", envir=.ns, inherits=FALSE)) {
    return(data.frame(name=character(0), etaDist=character(0),
                      stringsAsFactors=FALSE))
  }
  get("rxUiEtaDists", envir=.ns)(ui)
}

#' Pre-processing hook: expand declared eta distributions
#'
#' @param ui rxode2 ui object
#' @param est estimation routine name
#' @param data data
#' @param control control object
#' @return list with the expanded `ui`, or NULL when there is nothing to do
#' @noRd
#' @author Matthew L. Fidler
.preProcessEtaDist <- function(ui, est, data, control) {
  if (is.null(ui)) return(NULL)
  .d <- .rxUiEtaDists(ui)
  if (nrow(.d) == 0L) return(NULL)
  ## Refuse HERE rather than leaving it to the gate in nlmixr2Est().  The
  ## hooks run first, so a method that cannot use a declared distribution
  ## would otherwise pay for the expansion and its own pre-processing
  ## before being told no -- for `est="npag"` that is the whole
  ## nonparametric mu-expansion, which is not a wait to impose on someone
  ## who is about to get an error.  The gate stays as the backstop for
  ## paths that reach nlmixr2Est() without running hooks.
  .etaDistRefuse(.d, est, control)
  ## a method that translates the declaration itself has to see it
  ## unexpanded (see `.etaDistMethodAttr()`)
  if (identical(.etaDistMethodAttr(est, control), "native")) return(NULL)
  ## rxEtaDistExpand() clears the iniDf's `etaDist` column (rxode2
  ## R/etaDist.R), and the copula block it replaces with independent latents
  ## plus rxCor.* thetas -- so after expansion rxUiEtaDists() reports nothing
  ## and the declarations are unrecoverable.  Everything downstream that needs
  ## to know a random effect WAS declared (the ODE-free M-step in particular)
  ## reads this stash instead.
  .decl <- .etaDistDeclStash(ui, .d)
  ## Decompress BEFORE stashing: rxUiDecompress() on a compressed ui returns a
  ## new object, so assigning into it would write to a temporary and the stash
  ## would never reach the ui that is returned.
  .ui <- rxode2::rxUiDecompress(rxode2::rxEtaDistExpand(ui))
  ## In `meta`, which is the ONLY container that survives to the estimators.
  ## Measured, on a real saem fit, by planting a probe in each candidate and
  ## seeing which arrived: the ui environment, the control and an extra iniDf
  ## column were all gone by the time the M-step asked, because the ui is
  ## rebuilt and the est method installs its own freshly-built control after
  ## the hooks have run.  `meta` is rxode2's own metadata environment and is
  ## deliberately carried across model rewrites, so it is the one that holds.
  if (!is.null(.decl)) .etaDistDeclSet(.ui, .decl)
  list(ui=.ui)
}

#' Preserve the declarations across `rxEtaDistExpand()`
#'
#' The expansion is lossy by design: it rewrites the model into one with
#' ordinary standard-normal random effects, so the declaration it consumed is
#' no longer anywhere in the ui.  This records the parts that cannot be
#' reconstructed afterwards -- which etas were declared, with what family, and
#' how the copula paired them -- keyed so the expanded model's own thetas and
#' etas can be found from it.
#'
#' Copula pairing is read here rather than later because the block that carries
#' it (a unit-diagonal omega whose off-diagonal IS the correlation) is exactly
#' what the expansion removes.  Only a PAIR is recorded: the M-step drivers
#' reconstruct a single partner, so a larger block returns `NULL` and takes the
#' general path instead of being silently wrong.
#'
#' @param ui the UNEXPANDED rxode2 ui
#' @param d declared random effects, as `rxUiEtaDists()` returns them
#' @return a list, or `NULL` when the block is not one this can describe
#' @noRd
.etaDistDeclStash <- function(ui, d) {
  .ini <- rxode2::rxUiDecompress(ui)$iniDf
  .n <- nrow(d)
  .netaOf <- function(.nm) {
    .w <- which(.ini$name == .nm & .ini$neta1 == .ini$neta2)
    if (length(.w) == 1L) as.integer(.ini$neta1[.w]) else NA_integer_
  }
  .id <- vapply(d$name, .netaOf, integer(1))
  .cw <- rep(-1L, .n)                       # 0-based partner, or -1
  .ct <- rep(NA_character_, .n)             # the rxCor.* theta carrying it
  .off <- .ini[!is.na(.ini$neta1) & !is.na(.ini$neta2) &
                 .ini$neta1 != .ini$neta2, , drop = FALSE]
  if (nrow(.off) > 0L) {
    for (.r in seq_len(nrow(.off))) {
      .a1 <- which(.id == .off$neta1[.r]); .a2 <- which(.id == .off$neta2[.r])
      if (length(.a1) != 1L || length(.a2) != 1L) next  # not a declared pair
      .hi <- max(.a1, .a2); .lo <- min(.a1, .a2)
      if (.cw[.hi] >= 0L) return(NULL)                  # >2 declared partners
      .cw[.hi] <- .lo - 1L
      ## rxEtaDistExpand() names the Cholesky theta rxCor.<later>.<earlier>
      .ct[.hi] <- paste0("rxCor.", d$name[.hi], ".", d$name[.lo])
    }
  }
  list(name = as.character(d$name), etaDist = as.character(d$etaDist),
       corWith = .cw, corTheta = .ct)
}


preProcessHooksAdd(".preProcessEtaDist", .preProcessEtaDist)

#' The correlation matrix a fit's `rxCor.*` thetas encode
#'
#' Inverts the row-normalized Cholesky parameterization
#' `rxEtaDistExpand()` writes:
#'
#'   L[i, j] = tanh(y[i, j]) * s[i, j - 1],  L[i, i] = s[i, i - 1]
#'
#' with `s[i, 0] = 1` and `s[i, j] = s[i, j-1]*sqrt(1 - tanh(y[i,j])^2)`,
#' then returns `R = L L'`.
#'
#' @param nms the block's random effect names, in block order
#' @param y named numeric vector of the `rxCor.<i>.<j>` estimates
#' @return the correlation matrix, with `nms` as dimnames
#' @noRd
#' @author Matthew L. Fidler
.etaDistCorFromY <- function(nms, y) {
  .k <- length(nms)
  .L <- diag(.k)
  for (.i in seq_len(.k)) {
    .s <- 1.0
    for (.j in seq_len(.i - 1L)) {
      .t <- tanh(y[[paste0("rxCor.", nms[.i], ".", nms[.j])]])
      .L[.i, .j] <- .t * .s
      .s <- .s * sqrt(1 - .t * .t)
    }
    .L[.i, .i] <- .s
  }
  .R <- .L %*% t(.L)
  dimnames(.R) <- list(nms, nms)
  .R
}

#' Report the declared distributions on a fit
#'
#' Computed on demand through the `nmObjGet` accessors rather than stored
#' by a post-estimation hook, so they are there for every method that can
#' fit such a model -- the post-final hooks only run on the FOCEi path.
#'
#' `$etaDist` is the declarations as the model wrote them; `$etaDistCor`
#' is the copula correlation matrix of each declared block, rebuilt from
#' the `rxCor.*` estimates.
#'
#' The `rxCor.*` rows already read as correlations in the fit's
#' back-transformed column: `rxEtaDistExpand()` gives them
#' `backTransform("tanh")`, and `tanh()` of one is the partial correlation
#' between its two random effects given the ones before them -- which for
#' a 2x2 block, the usual case, is simply the correlation.
#'
#' @param x list of the fit environment and the exact flag, as
#'   `nmObjGet()` dispatches it
#' @param ... ignored
#' @return the declarations, or NULL when the model declared none
#' @export
#' @keywords internal
#' @author Matthew L. Fidler
nmObjGet.etaDist <- function(x, ...) {
  .info <- .etaDistInfo(x[[1]])
  if (is.null(.info)) return(NULL)
  .info$etaDist
}
attr(nmObjGet.etaDist, "desc") <-
  "The non-normal random effect distributions the model declared"

#' @rdname nmObjGet.etaDist
#' @export
nmObjGet.etaDistCor <- function(x, ...) {
  .env <- x[[1]]
  .fix <- try(get("fixef", envir=.env), silent=TRUE)
  if (inherits(.fix, "try-error") || is.null(.fix)) return(NULL)
  .info <- .etaDistInfo(.env)
  .blocks <- if (is.null(.info)) NULL else .info$blocks
  if (is.null(.blocks) || length(.blocks) == 0L) {
    ## `etaDistInfo` does not survive onto the fit's ui, so the blocks are
    ## recovered from the fit itself -- see .etaDistBlocksFromFit()
    .blocks <- .etaDistBlocksFromFit(x[[1]])
  }
  if (length(.blocks) == 0L) return(NULL)
  .cor <- lapply(.blocks, function(.nms) {
    .need <- unlist(lapply(seq_along(.nms), function(.i) {
      if (.i == 1L) return(NULL)
      paste0("rxCor.", .nms[.i], ".", .nms[seq_len(.i - 1L)])
    }), use.names=FALSE)
    if (length(.need) == 0L || !all(.need %in% names(.fix))) return(NULL)
    .etaDistCorFromY(.nms, as.list(.fix[.need]))
  })
  names(.cor) <- vapply(.blocks, function(.n) .n[1], character(1),
                        USE.NAMES=FALSE)
  .cor <- .cor[!vapply(.cor, is.null, logical(1))]
  if (length(.cor) == 0L) return(NULL)
  .cor
}
attr(nmObjGet.etaDistCor, "desc") <-
  "The Gaussian copula correlation of each declared random effect block"

#' What `rxEtaDistExpand()` recorded on the fit's model
#'
#' @param env fit environment
#' @return the `etaDistInfo` list, or NULL
#' @noRd
#' @author Matthew L. Fidler
.etaDistInfo <- function(env) {
  .ui <- try(get("ui", envir=env), silent=TRUE)
  if (inherits(.ui, "try-error") || is.null(.ui)) return(NULL)
  .ui <- try(rxode2::rxUiDecompress(.ui), silent=TRUE)
  if (inherits(.ui, "try-error")) return(NULL)
  .info <- try(get("etaDistInfo", envir=.ui), silent=TRUE)
  if (inherits(.info, "try-error")) return(NULL)
  .info
}

#' Drop the latent random effects from the reported parameter table
#'
#' `rxEtaDistExpand()` leaves the latent standard normals (`rxz.<eta>`) in
#' `parFixed`, where they print as `NA` in every column.  Their variance is
#' fixed at one by construction -- that is what makes the copula a copula --
#' so they are not estimates and there is nothing to report for them.
#'
#' The copula correlations (`rxCor.<i>.<j>`) are NOT touched here: the
#' expansion already gives them `backTransform = "tanh"`, so their
#' back-transformed column is the correlation and has always been correct.
#' (An earlier version of this hook recomputed it, on the strength of my
#' having misread the raw Estimate column as a correlation.  It was the
#' reading that was wrong, not the table.)
#'
#' @param ret fit object
#' @return `ret`, without the latent rows
#' @noRd
#' @author Matthew L. Fidler
.postFinalEtaDistParFixed <- function(ret) {
  .env <- try(ret$env, silent=TRUE)
  if (inherits(.env, "try-error") || is.null(.env)) return(ret)
  .pfd <- try(get("parFixedDf", envir=.env), silent=TRUE)
  if (inherits(.pfd, "try-error") || is.null(.pfd)) return(ret)
  .nm <- rownames(.pfd)
  if (is.null(.nm) || !any(grepl("^rxz[.]", .nm))) return(ret)
  assign("parFixedDf", .pfd[!grepl("^rxz[.]", .nm), , drop=FALSE], envir=.env)
  .pf <- try(get("parFixed", envir=.env), silent=TRUE)
  if (!inherits(.pf, "try-error") && !is.null(.pf)) {
    .nm2 <- rownames(.pf)
    if (!is.null(.nm2) && any(grepl("^rxz[.]", .nm2))) {
      assign("parFixed", .pf[!grepl("^rxz[.]", .nm2), , drop=FALSE], envir=.env)
    }
  }
  ret
}

#' Report the declared block in `$omega`, under the names the model used
#'
#' The user writes `eta.cl + eta.v1 ~ c(1, 0.5, 1)` -- a covariance block with
#' the correlation in it.  What comes back is the expansion's internals: a 2x2
#' identity named `rxz.eta.cl` / `rxz.eta.v1`, with the fitted correlation
#' living in a `rxCor.*` theta instead.  The block the user wrote is nowhere in
#' `$omega`.
#'
#' The latent random effects are standard normals, so their covariance matrix
#' IS the correlation matrix -- unit diagonal is what the declaration requires
#' -- and the fitted block goes back into `$omega` on the covariance scale,
#' named as the model named it.  `$omegaR` then derives the correlation view
#' through the machinery every other model uses.
#'
#' The `rxCor.*` rows stay in `parFixed`: that is where their standard error
#' is, on the estimated scale like every other row, and moving the value into
#' `$omega` must not take the uncertainty out of the output.
#'
#' @param ret fit object
#' @return `ret`, with `$omega` carrying the declared block
#' @noRd
#' @author Matthew L. Fidler
.postFinalEtaDistOmega <- function(ret) {
  .env <- try(ret$env, silent=TRUE)
  if (inherits(.env, "try-error") || is.null(.env)) return(ret)
  .om <- try(get("omega", envir=.env), silent=TRUE)
  if (inherits(.om, "try-error") || is.null(.om) || !is.matrix(.om)) return(ret)
  .blocks <- try(.etaDistBlocksFromFit(ret), silent=TRUE)
  if (inherits(.blocks, "try-error") || length(.blocks) == 0L) return(ret)
  .fix <- try(get("fixef", envir=.env), silent=TRUE)
  if (inherits(.fix, "try-error") || is.null(.fix)) return(ret)
  .nm <- rownames(.om)
  if (is.null(.nm)) return(ret)
  for (.b in .blocks) {
    .need <- unlist(lapply(seq_along(.b), function(.i) {
      if (.i == 1L) return(NULL)
      paste0("rxCor.", .b[.i], ".", .b[seq_len(.i - 1L)])
    }), use.names=FALSE)
    if (length(.need) == 0L || !all(.need %in% names(.fix))) next
    .R <- try(.etaDistCorFromY(.b, as.list(.fix[.need])), silent=TRUE)
    if (inherits(.R, "try-error")) next
    .row <- match(paste0("rxz.", .b), .nm)
    if (anyNA(.row)) next
    .om[.row, .row] <- .R
    .nm[.row] <- .b          # report them as the model named them
  }
  dimnames(.om) <- list(.nm, .nm)
  assign("omega", .om, envir=.env)
  ret
}

postFinalObjectHooksAdd(".postFinalEtaDistOmega", .postFinalEtaDistOmega)

postFinalObjectHooksAdd(".postFinalEtaDistParFixed", .postFinalEtaDistParFixed)

#' Recover the declared correlation blocks from the fit itself
#'
#' `rxEtaDistExpand()` records what it did in `etaDistInfo` on the ui it
#' returns, but that does not survive to the fit object (measured: present on
#' the expanded ui, absent on `fit$ui`), which is why `fit$etaDistCor` came
#' back NULL for a model that plainly has a declared block.
#'
#' Everything needed is still in the fit, so it is read from there instead of
#' carried: the latent random effects are `rxz.<declared eta>` and appear in
#' the omega in their block order, and the copula parameters are
#' `rxCor.<i>.<j>` thetas naming the pair they connect.  Two random effects
#' are in the same block exactly when such a theta joins them.
#'
#' @param ret fit object
#' @return list of character vectors, one per block, in omega order; empty
#'   when the model declares nothing
#' @noRd
#' @author Matthew L. Fidler
.etaDistBlocksFromFit <- function(ret) {
  .ui <- try(rxode2::rxUiDecompress(ret$ui), silent=TRUE)
  if (inherits(.ui, "try-error") || is.null(.ui)) return(list())
  .ini <- .ui$iniDf
  if (is.null(.ini) || !any(names(.ini) == "neta1")) return(list())
  .e <- .ini[!is.na(.ini$neta1) & .ini$neta1 == .ini$neta2, ]
  if (nrow(.e) == 0L) return(list())
  .e <- .e[order(.e$neta1), ]
  .lat <- .e$name[grepl("^rxz[.]", .e$name)]
  if (length(.lat) == 0L) return(list())
  .dec <- sub("^rxz[.]", "", .lat)
  .th <- .ini$name[!is.na(.ini$ntheta)]
  .cor <- .th[grepl("^rxCor[.]", .th)]
  ## adjacency from the copula thetas; a lone declared random effect is its
  ## own block and simply has no correlation to report
  .grp <- seq_along(.dec)
  for (.c in .cor) {
    .p <- sub("^rxCor[.]", "", .c)
    .i <- which(vapply(.dec, function(.d) startsWith(.p, paste0(.d, ".")),
                       logical(1)))
    for (.ii in .i) {
      .j <- which(.dec == sub(paste0("^", .dec[.ii], "[.]"), "", .p))
      if (length(.j) == 1L) {
        .keep <- min(.grp[.ii], .grp[.j])
        .drop <- max(.grp[.ii], .grp[.j])
        .grp[.grp == .drop] <- .keep
      }
    }
  }
  .out <- lapply(sort(unique(.grp)), function(.g) .dec[.grp == .g])
  .out[vapply(.out, length, integer(1)) > 1L]
}
