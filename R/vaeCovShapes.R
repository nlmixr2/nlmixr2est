## Covariate shapes (parameterizations) for the VAE covariate search.
##
## The shape vocabulary matches nlmixr2scm's `shapes=` argument.  The VAE
## covariate M-step is an OLS fit of the latent mean on [1 | X_S], so its
## objective depends only on the column SPAN.  With a free intercept
## log(cov/ctr) and log(cov) span the same space, as do (cov - ctr), cov and
## cov/ctr.  Shapes therefore collapse to two searchable FAMILIES; the shape
## itself only chooses how the accepted relationship is written back.

## user-selectable shapes for a continuous covariate, in canonical order
.vaeContShapes <- c("power", "lin", "log", "identity", "center")
## every shape name, including the categorical one (never user-selectable)
.vaeAllShapes <- c(.vaeContShapes, "cat")

#' Searchable family a covariate shape belongs to
#'
#' Shapes in one family span the same column space and are indistinguishable to
#' the selection objective.
#' @param shape character vector of shape names
#' @return character vector of families, one of `"log"`, `"lin"`, `"cat"`
#' @noRd
.vaeShapeFamily <- function(shape) {
  .f <- c(power = "log", log = "log",
          lin = "lin", identity = "lin", center = "lin",
          cat = "cat")[shape]
  if (anyNA(.f)) {
    stop("unknown covariate shape: ",
         paste(unique(shape[is.na(.f)]), collapse = ", "),
         "\navailable: ", paste(.vaeAllShapes, collapse = ", "),
         call. = FALSE)
  }
  unname(.f)
}

#' Validate a vector of user-supplied continuous shape names
#' @param shape character vector
#' @return the vector, unchanged, in the order given
#' @noRd
.vaeAssertContShapes <- function(shape) {
  checkmate::assertCharacter(shape, min.len = 1, any.missing = FALSE)
  .bad <- setdiff(shape, .vaeContShapes)
  if (length(.bad) > 0L) {
    stop("unknown covariate shape: ", paste(.bad, collapse = ", "),
         "\navailable: ", paste(.vaeContShapes, collapse = ", "),
         call. = FALSE)
  }
  if (anyDuplicated(shape)) stop("duplicate covariate shape", call. = FALSE)
  shape
}

#' Model text for one covariate shape
#'
#' The returned string is what multiplies the coefficient in the written-back
#' model, e.g. `log(WT/70.5)`.
#' @param shape shape name
#' @param col data column name
#' @param center centering value (used by `power`, `lin`, `center`)
#' @param level factor level, for `shape = "cat"`
#' @param raw TRUE when a `cat` column is a bare 0/1 indicator already in its
#'   natural parameterization -- it then enters as the column itself
#' @return length-one character string
#' @noRd
.vaeShapeExpr <- function(shape, col, center = NA_real_, level = NULL,
                          raw = FALSE) {
  .c <- if (is.na(center)) NA_character_ else as.character(signif(center, 12))
  switch(shape,
         power = paste0("log(", col, "/", .c, ")"),
         log = paste0("log(", col, ")"),
         lin = paste0("(", col, " - ", .c, ")"),
         identity = col,
         center = paste0("(", col, "/", .c, ")"),
         cat = if (raw) col else paste0("(", col, " == ", .vaeLevelLit(level), ")"),
         stop("unknown covariate shape: ", shape, call. = FALSE))
}

#' Literal for a factor level inside generated model text
#'
#' A numeric level code compares numerically; anything else is quoted so
#' rxode2's string-comparison path handles it.  Whether the level IS numeric is
#' decided upstream by `.vaeCovLevelValue` from the covariate's own type -- do
#' not re-coerce here, or a character level that merely looks numeric (`"01"`)
#' would be emitted as `== 1`, comparing against the wrong thing and no longer
#' matching its own column when the model is piped back.
#' @param level the level, already typed by `.vaeCovLevelValue`
#' @return length-one character string
#' @noRd
.vaeLevelLit <- function(level) {
  if (is.numeric(level)) return(as.character(signif(level, 12)))
  paste0("\"", as.character(level), "\"")
}

#' Re-express a family coefficient in a chosen shape's parameterization
#'
#' The search fits `mu = ic + b*x` where `x` is the family column.  The written
#' model uses the shape expression `e` instead, so the coefficient and the
#' intercept both move such that `ic2 + b2*e == ic + b*x` for every covariate
#' value.
#' @param shape shape being written
#' @param center centering value of the family column
#' @param beta fitted family coefficient
#' @return list(beta = written coefficient, interceptAdj = amount to ADD to the
#'   structural theta)
#' @noRd
.vaeShapeBeta <- function(shape, center, beta) {
  switch(shape,
         power = list(beta = beta, interceptAdj = 0),
         log = list(beta = beta, interceptAdj = -beta * log(center)),
         lin = list(beta = beta, interceptAdj = 0),
         identity = list(beta = beta, interceptAdj = -beta * center),
         center = list(beta = beta * center, interceptAdj = -beta * center),
         cat = list(beta = beta, interceptAdj = 0),
         stop("unknown covariate shape: ", shape, call. = FALSE))
}

#' Normalize a user `shapes=` specification into match rules
#'
#' Accepts a character vector (one global rule), a list named by covariate, or
#' a list of `list(var=, covar=, shapes=)` items in nlmixr2scm's `pairsVec`
#' form.  Covariate names match case-insensitively because the VAE upper-cases
#' data columns.
#' @param spec the `shapes` control value
#' @return data.frame with columns `var`, `cov` (NA meaning "any") and a
#'   `shapes` list-column
#' @noRd
.vaeResolveShapes <- function(spec) {
  .mk <- function(var, cov, shapes) {
    .d <- data.frame(var = as.character(var), cov = as.character(cov),
                     stringsAsFactors = FALSE)
    .d$shapes <- list(.vaeAssertContShapes(shapes))
    .d
  }
  if (is.null(spec)) spec <- .vaeContShapes
  if (is.character(spec)) return(.mk(NA_character_, NA_character_, spec))
  if (!is.list(spec)) stop("shapes must be a character vector or a list", call. = FALSE)
  if (length(spec) == 0L) return(.mk(NA_character_, NA_character_, .vaeContShapes))
  .isPair <- all(vapply(spec, function(e) is.list(e), logical(1)))
  .out <- list()
  if (.isPair) {
    for (.e in spec) {
      .cov <- if (is.null(.e$covar)) .e$cov else .e$covar
      .sh <- if (is.null(.e$shapes)) .e$shape else .e$shapes
      if (is.null(.sh)) .sh <- .vaeContShapes
      .out[[length(.out) + 1L]] <-
        .mk(if (is.null(.e$var)) NA_character_ else as.character(.e$var),
            if (is.null(.cov)) NA_character_ else toupper(as.character(.cov)),
            as.character(.sh))
    }
  } else {
    .nm <- names(spec)
    if (is.null(.nm) || any(!nzchar(.nm))) {
      stop("shapes list must be named by covariate, or be list(var=, covar=, shapes=) items",
           call. = FALSE)
    }
    for (.i in seq_along(spec)) {
      .out[[length(.out) + 1L]] <-
        .mk(NA_character_, toupper(.nm[.i]), as.character(spec[[.i]]))
    }
  }
  do.call(rbind, .out)
}

#' Numeric column a shape's model text evaluates to
#'
#' The inverse of `.vaeShapeExpr`: given the raw per-subject covariate values it
#' returns exactly what the written expression computes.  Used to rebuild a
#' pinned column as the model's OWN expression, so an estimated slope transfers
#' back with no correction whatever shape the user wrote.
#' @param shape shape name
#' @param v raw per-subject covariate values
#' @param center centering value written in the model
#' @return numeric vector, same length as `v`
#' @noRd
.vaeShapeValue <- function(shape, v, center) {
  switch(shape,
         power = log(v / center),
         log = log(v),
         lin = v - center,
         identity = v,
         center = v / center,
         stop("unknown covariate shape: ", shape, call. = FALSE))
}

#' The sub-expression a coefficient multiplies
#'
#' A model line may carry several covariate effects, so the shape has to be read
#' off the factor THIS coefficient multiplies rather than off the whole line.
#'
#' The walk is restricted to ADDITIVE context -- the assignment, one
#' mu-referencing transform, then only `+`/`-` and parentheses.  It deliberately
#' does not descend into an arbitrary call or a nested product, because the term
#' must be a standalone additive `coef * <expr>` for its slope to transfer.
#' Searching everywhere would return `WT` for `sqrt(coef * WT)` or for the
#' interaction `coef * WT * AGE`, and the pinned fit would then silently replace
#' the user's term with a plain linear WT effect.
#' @param e expression to walk (a model line)
#' @param coef coefficient (theta) name
#' @return the multiplied expression, or `NULL` when there is no such term
#' @noRd
.vaeCoefFactor <- function(e, coef) {
  .unwrap <- function(x) {
    while (is.call(x) && identical(x[[1L]], as.name("(")) && length(x) == 2L) {
      x <- x[[2L]]
    }
    x
  }
  ## a standalone `coef * <expr>` (either order), else NULL
  .term <- function(x) {
    x <- .unwrap(x)
    if (is.call(x) && identical(x[[1L]], as.name("*")) && length(x) == 3L) {
      if (is.name(x[[2L]]) && identical(as.character(x[[2L]]), coef)) return(x[[3L]])
      if (is.name(x[[3L]]) && identical(as.character(x[[3L]]), coef)) return(x[[2L]])
    }
    NULL
  }
  ## Only POSITIVE additive position: the right operand of a binary `-`, and
  ## anything under a unary `-`, carries a negation the pinned slope would not
  ## see -- writing the fitted beta back into `theta - beta*cov` would flip the
  ## sign of the effect.  Those are left unmatched so the coefficient regresses.
  .walk <- function(x) {
    x <- .unwrap(x)
    if (!is.call(x)) return(NULL)
    .op <- if (is.name(x[[1L]])) as.character(x[[1L]]) else ""
    if (.op == "+") {
      if (length(x) == 2L) return(.walk(x[[2L]]))          # unary plus
      if (length(x) == 3L) {
        .r <- .walk(x[[2L]])
        if (!is.null(.r)) return(.r)
        return(.walk(x[[3L]]))
      }
    }
    if (.op == "-" && length(x) == 3L) return(.walk(x[[2L]]))
    .term(x)
  }
  ## The whole line must mention the coefficient exactly once.  `b*x1 + b*x2`
  ## would otherwise be fitted on x1 alone and then written back to both terms,
  ## and a coefficient reused inside another call would escape the walk entirely.
  .count <- function(x) {
    if (is.name(x)) return(as.integer(identical(as.character(x), coef)))
    if (is.call(x)) {
      return(sum(vapply(as.list(x)[-1L], .count, integer(1))))
    }
    0L
  }
  if (.count(e) != 1L) return(NULL)
  .e <- e
  if (is.call(.e) && is.name(.e[[1L]]) &&
        as.character(.e[[1L]]) %in% c("<-", "=", "~") && length(.e) == 3L) {
    .e <- .e[[3L]]
  }
  ## a single mu-referencing transform wrapper, e.g. exp(theta + beta*cov + eta)
  .e <- .unwrap(.e)
  if (is.call(.e) && is.name(.e[[1L]]) && length(.e) == 2L &&
        as.character(.e[[1L]]) %in% c("exp", "log", "logit", "expit",
                                      "probit", "probitInv")) {
    .e <- .e[[2L]]
  }
  .walk(.e)
}

#' Which shape a written covariate expression is in
#'
#' Recognizes every form `.vaeShapeExpr` emits, so a model written (or written
#' BACK) by the VAE can be piped into another fit and have each covariate
#' effect pinned to the shape it was written in.  Anything else is unrecognized
#' and the caller routes that coefficient to the regress M-step.
#' Classification is STRICT: the expression must be one of the emitted forms
#' (bare, or wrapped in parentheses).  It deliberately does not search inside an
#' unrecognized call -- `sqrt(WT)` must not be read as the identity shape just
#' because `WT` appears in it, or a slope that does not transfer would be pinned
#' as though it did.
#' @param e expression the coefficient multiplies (see `.vaeCoefFactor`)
#' @param cov raw covariate (data column) name
#' @return list(shape, center, level); `shape` is `NA` when unrecognized
#' @noRd
.vaeDetectShape <- function(e, cov) {
  .no <- list(shape = NA_character_, center = NA_real_, level = NULL)
  if (is.null(e)) return(.no)
  .isCov <- function(x) is.name(x) && identical(as.character(x), cov)
  .num <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  ## strip redundant parentheses -- "(WT - 70)" is a call to `(`
  while (is.call(e) && identical(e[[1L]], as.name("(")) && length(e) == 2L) {
    e <- e[[2L]]
  }
  if (!is.call(e)) {
    ## a bare covariate multiplied by the coefficient is the identity shape
    if (.isCov(e)) return(list(shape = "identity", center = 0, level = NULL))
    return(.no)
  }
  .op <- if (is.name(e[[1L]])) as.character(e[[1L]]) else ""
  if (.op == "log" && length(e) == 2L) {
    .a <- e[[2L]]
    while (is.call(.a) && identical(.a[[1L]], as.name("(")) && length(.a) == 2L) {
      .a <- .a[[2L]]
    }
    if (.isCov(.a)) return(list(shape = "log", center = 1, level = NULL))
    if (is.call(.a) && identical(.a[[1L]], as.name("/")) && length(.a) == 3L &&
          .isCov(.a[[2L]]) && .num(.a[[3L]])) {
      return(list(shape = "power", center = as.numeric(.a[[3L]]), level = NULL))
    }
    return(.no)          # log() of something not transferable
  }
  if (.op == "/" && length(e) == 3L && .isCov(e[[2L]]) && .num(e[[3L]])) {
    return(list(shape = "center", center = as.numeric(e[[3L]]), level = NULL))
  }
  if (.op == "-" && length(e) == 3L && .isCov(e[[2L]]) && .num(e[[3L]])) {
    return(list(shape = "lin", center = as.numeric(e[[3L]]), level = NULL))
  }
  if (.op == "==" && length(e) == 3L && .isCov(e[[2L]])) {
    return(list(shape = "cat", center = 0, level = e[[3L]]))
  }
  .no
}

#' Is a shape usable at this centering value?
#'
#' `center` divides by the centering value and `log` takes its logarithm, so
#' neither is expressible when the center is zero (or negative).  Without this
#' guard a zero-centered covariate would be written as `beta*(COV/0)` with the
#' coefficient rescaled to exactly 0 -- silently erasing a selected effect.
#' @param shape shape name(s)
#' @param center centering value
#' @return logical, one per shape
#' @noRd
.vaeShapeUsable <- function(shape, center) {
  vapply(shape, function(.s) {
    switch(.s,
           center = is.finite(center) && center != 0,
           log = is.finite(center) && center > 0,
           power = is.finite(center) && center > 0,
           TRUE)
  }, logical(1), USE.NAMES = FALSE)
}

#' Per-(latent dim, column) mask of shapes the user allows
#'
#' `shapes=` may restrict a single (parameter, covariate) pair, but the design
#' matrix is shared across latent dimensions, so the restriction has to be
#' enforced as a mask rather than by omitting columns.  Categorical columns are
#' never restricted -- `shapes=` governs continuous parameterizations only.
#' @param cov output of `.vaeCovariateSearch`
#' @param rules output of `.vaeResolveShapes`
#' @param etaNames per-latent-dim random-effect names
#' @keywords internal
#' @param thetaForEta per-latent-dim mu-referenced theta names (may be NA)
#' @return integer 0/1 matrix, `length(etaNames)` by `ncol(cov$covMat)`
#' @noRd
.vaeShapeAllowMask <- function(cov, rules, etaNames, thetaForEta) {
  .nCov <- length(cov$covNames)
  .m <- matrix(1L, length(etaNames), .nCov)
  if (.nCov == 0L || is.null(rules)) return(.m)
  for (.k in seq_along(etaNames)) {
    .al <- c(etaNames[.k], thetaForEta[.k], sub("^eta\\.", "", etaNames[.k]))
    .al <- unique(.al[!is.na(.al)])
    for (.j in seq_len(.nCov)) {
      if (identical(cov$covFamily[.j], "cat")) next
      .ok <- .vaeShapesFor(rules, .al, cov$covRaw[.j])
      .want <- .vaeShapeFamily(.ok)
      ## A covariate whose requested family is unavailable (log shapes on
      ## non-positive data) is carried by a FALLBACK column of another family.
      ## Masking that column would block the covariate from the search entirely
      ## and defeat the fallback, so when none of the requested families exist
      ## for this covariate, leave whatever does exist selectable.
      .have <- setdiff(unique(cov$covFamily[cov$covRaw == cov$covRaw[.j]]), "cat")
      if (length(intersect(.want, .have)) == 0L) next
      if (!(cov$covFamily[.j] %in% .want)) .m[.k, .j] <- 0L
    }
  }
  .m
}

#' Shapes allowed for one (parameter, covariate) pair
#'
#' The most specific matching rule wins: (var, cov) beats cov-only, which beats
#' var-only, which beats the global rule.  `par` is matched against any of the
#' parameter's aliases (eta name, mu-referenced theta name, bare name).  A pair
#' no rule mentions keeps every shape -- `shapes=` restricts parameterizations,
#' never which covariates are searched (that is `pinCovariates`).
#' @param rules output of `.vaeResolveShapes`
#' @param parAliases character vector of names this parameter answers to
#' @param cov covariate (raw data column) name
#' @return character vector of shapes, in the order the user listed them
#' @noRd
.vaeShapesFor <- function(rules, parAliases, cov) {
  .cov <- toupper(cov)
  .mVar <- vapply(rules$var, function(v) is.na(v) || v %in% parAliases, logical(1))
  .mCov <- !is.na(rules$cov) & rules$cov == .cov
  .anyCov <- is.na(rules$cov)
  .spec <- ifelse(.mCov, 2L, 0L) + ifelse(!is.na(rules$var), 1L, 0L)
  .ok <- .mVar & (.mCov | .anyCov)
  if (!any(.ok)) return(.vaeContShapes)
  .w <- which(.ok)
  .w <- .w[.spec[.w] == max(.spec[.w])]
  rules$shapes[[.w[length(.w)]]]
}
