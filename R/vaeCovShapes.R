## Covariate shapes (parameterizations) for the VAE covariate search.
##
## The shape vocabulary matches nlmixr2scm's `shapes=` argument.  The VAE
## covariate M-step is an OLS fit of the latent mean on [1 | X_S], so its
## objective depends only on the column SPAN.  With a free intercept
## log(cov/ctr) and log(cov) span the same space, as do (cov - ctr), cov and
## cov/ctr.  Shapes therefore collapse to searchable FAMILIES; within a family the
## shape only chooses how the accepted relationship is written back.
##
## "hockey" is the exception: a two-armed piecewise-linear relationship knotted at
## the centering value.  Its span strictly CONTAINS the lin family's -- the arms
## sum to the lin column -- so it is a family of its own, not another
## parameterization of "lin", and it is the only shape contributing more than one
## design column (one per arm).  The arms are not user-selectable: `shapes=` takes
## "hockey" and both arms come with it, all-or-none, which is what keeps the
## written form well defined (see plans/vae-hockey-shape.md).

## every shape a user may NAME in `shapes=`, in canonical order
.vaeContShapes <- c("power", "lin", "log", "identity", "center", "hockey")
## the shapes actually TRIED when `shapes=` is not given.  Deliberately a
## separate vector even while it matches the one above: the two roles diverge,
## and conflating them makes any newly-nameable shape a default in the same
## stroke -- which is a behavior change nobody asked for.
.vaeDefaultShapes <- c("power", "lin", "log", "identity", "center", "hockey")
## the two design columns a requested "hockey" expands into, low arm first
.vaeHockeyArms <- c("hockeyLow", "hockeyHi")
## every shape name, including the categorical one (never user-selectable)
.vaeAllShapes <- c(.vaeContShapes, "cat")

#' Searchable family a covariate shape belongs to
#'
#' Shapes in one family span the same column space and are indistinguishable to
#' the selection objective.
#' Accepts the internal hockey arm names as well as the user-facing `"hockey"`,
#' since a selected column carries the arm name rather than the request.
#' @param shape character vector of shape names
#' @return character vector of families, one of `"log"`, `"lin"`, `"hockey"`,
#'   `"cat"`
#' @noRd
.vaeShapeFamily <- function(shape) {
  .f <- c(power = "log", log = "log",
          lin = "lin", identity = "lin", center = "lin",
          hockey = "hockey", hockeyLow = "hockey", hockeyHi = "hockey",
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
         ## the arms are disjoint (`<` against `>=`), so they partition subjects
         ## and sum to the lin column exactly -- a subject sitting ON the knot
         ## belongs to the high arm alone.  Both vanish at the knot, so the
         ## structural theta stays the parameter value AT the center.
         hockeyLow = paste0("(", col, " < ", .c, ")*(", col, " - ", .c, ")"),
         hockeyHi = paste0("(", col, " >= ", .c, ")*(", col, " - ", .c, ")"),
         cat = if (raw) col else paste0("(", col, " == ", .vaeLevelLit(level), ")"),
         stop("unknown covariate shape: ", shape, call. = FALSE))
}

#' Name component a shape contributes to its coefficient name
#'
#' Coefficients are named `beta.<par>.<COV>.<tag>`.  The tag is the shape itself
#' except for the hockey arms, which read as `hockey.low` / `hockey.hi` so the
#' pair is obvious in `ini()` and in the parameter table.
#' @param shape shape name
#' @return length-one character string
#' @noRd
.vaeShapeCoefTag <- function(shape) {
  switch(shape, hockeyLow = "hockey.low", hockeyHi = "hockey.hi", shape)
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
  ## encodeString escapes embedded quotes and backslashes; pasting raw quotes
  ## around a level like `A"B` emits model text that does not even parse
  encodeString(as.character(level), quote = "\"")
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
         ## an arm's written expression IS its design column, and both vanish at
         ## the knot, so no part of the effect moves into the intercept
         hockeyLow = list(beta = beta, interceptAdj = 0),
         hockeyHi = list(beta = beta, interceptAdj = 0),
         cat = list(beta = beta, interceptAdj = 0),
         stop("unknown covariate shape: ", shape, call. = FALSE))
}

## The element name that carries the eligibility flag rather than a rule.  Matched
## EXACTLY and case-sensitively; data columns are upper-cased by the search, so a
## real `FIXCOV` column is still reachable as a covariate name (the collision is
## caught in .vaeEligible, where the covariate names are known).
.vaeFixCovName <- "fixCov"

#' Normalize a user `shapes=` specification into match rules
#'
#' Accepts a character vector (one global rule) or a list.  A list element is
#' dispatched on ITSELF rather than on the whole list, so the two historic forms
#' may be mixed freely:
#'
#'   * a list element -> a `list(var=, covar=, shapes=)` rule (nlmixr2scm's
#'     `pairsVec` form);
#'   * a named element -> shorthand for `list(covar = <name>, shapes = <value>)`,
#'     which is what keeps the named form from needing semantics of its own;
#'   * the element named `fixCov` -> the eligibility flag, not a rule.
#'
#' A shape value of `TRUE` means "eligible, with the default shapes" -- the way
#' to name a covariate whose parameterization is not up for discussion, i.e. a
#' categorical, which takes no shape at all.
#'
#' Covariate names match case-insensitively because the VAE upper-cases data
#' columns.
#' @param spec the `shapes` control value
#' @return list with `rules` (data.frame of `var`, `cov` -- NA meaning "any" --
#'   and a `shapes` list-column) and `fixCov` (length-one logical)
#' @noRd
.vaeResolveShapes <- function(spec) {
  .mk <- function(var, cov, shapes) {
    .d <- data.frame(var = as.character(var), cov = as.character(cov),
                     stringsAsFactors = FALSE)
    ## TRUE == "eligible, default shapes"; anything else must name real shapes
    if (isTRUE(shapes)) shapes <- .vaeDefaultShapes
    if (is.logical(shapes)) {
      stop("shapes value must be TRUE or a shape vector, not ",
           deparse(shapes), call. = FALSE)
    }
    .d$shapes <- list(.vaeAssertContShapes(as.character(shapes)))
    .d
  }
  .ret <- function(rules, fixCov) list(rules = rules, fixCov = fixCov)
  if (is.null(spec)) spec <- .vaeDefaultShapes
  ## a character vector names no covariate, so there is nothing for fixCov to fix
  if (is.character(spec)) {
    return(.ret(.mk(NA_character_, NA_character_, spec), FALSE))
  }
  if (!is.list(spec)) stop("shapes must be a character vector or a list", call. = FALSE)
  if (length(spec) == 0L) {
    return(.ret(.mk(NA_character_, NA_character_, .vaeDefaultShapes), FALSE))
  }
  .nm <- names(spec)
  if (is.null(.nm)) .nm <- rep("", length(spec))
  ## Pull fixCov out BEFORE any rule parsing.  It is not a rule, and a bare
  ## logical among the elements must not be read as one.
  .fixCov <- TRUE
  .fx <- which(.nm == .vaeFixCovName)
  if (length(.fx) > 1L) stop("fixCov given more than once in shapes", call. = FALSE)
  if (length(.fx) == 1L) {
    .fixCov <- spec[[.fx]]
    checkmate::assertLogical(.fixCov, len = 1, any.missing = FALSE,
                             .var.name = "fixCov")
    spec <- spec[-.fx]
    .nm <- .nm[-.fx]
  }
  if (length(spec) == 0L) {
    ## `shapes = list(fixCov = TRUE)` asks to fix the searched set to the
    ## covariates named, and names none.  Silently reading that as FALSE would
    ## override an explicit flag; silently reading it as TRUE would search
    ## nothing.  Neither is what was meant, so say so.
    if (isTRUE(.fixCov) && length(.fx) == 1L) {
      stop("shapes: fixCov=TRUE but no covariate is named\n",
           "  name the covariates to search, or use covariateSelection=FALSE ",
           "to turn the search off", call. = FALSE)
    }
    return(.ret(.mk(NA_character_, NA_character_, .vaeDefaultShapes), FALSE))
  }
  .out <- vector("list", length(spec))
  for (.i in seq_along(spec)) {
    .e <- spec[[.i]]
    if (is.list(.e)) {
      ## pair rule; `$` partial-matches, so cov/covar and shape/shapes both work
      .cov <- if (is.null(.e$covar)) .e$cov else .e$covar
      .sh <- if (is.null(.e$shapes)) .e$shape else .e$shapes
      if (is.null(.sh)) .sh <- .vaeDefaultShapes
      .out[[.i]] <-
        .mk(if (is.null(.e$var)) NA_character_ else as.character(.e$var),
            if (is.null(.cov)) NA_character_ else toupper(as.character(.cov)),
            .sh)
    } else if (nzchar(.nm[.i])) {
      ## named element: exactly a covar-only pair rule
      .out[[.i]] <- .mk(NA_character_, toupper(.nm[.i]), .e)
    } else {
      stop("shapes list element ", .i,
           " must be named by covariate, or be a list(var=, covar=, shapes=) item",
           call. = FALSE)
    }
  }
  .ret(do.call(rbind, .out), .fixCov)
}

#' Which (parameter, covariate) pairs the search may consider
#'
#' `fixCov = TRUE` (the default whenever `shapes=` is given as a list) fixes the
#' searched covariate set to exactly the covariates the rules NAME.  Naming a
#' covariate is the statement that it belongs in the search, so the common case
#' -- "search these, and only these" -- costs nothing beyond listing them, with
#' no per-covariate opt-out for everything left out.
#'
#' Eligibility and parameterization are separate passes: this decides WHICH
#' cells may be searched, `.vaeShapesFor` decides what an eligible cell may look
#' like.  The specificity ladder is therefore untouched.
#'
#' @param rules `$rules` from `.vaeResolveShapes`
#' @param fixCov `$fixCov` from `.vaeResolveShapes`
#' @param etaNames per-latent-dim random-effect names
#' @param thetaForEta per-latent-dim mu-referenced theta names (may be NA)
#' @param covRaw raw (data column) name per search column
#' @return logical matrix, `length(etaNames)` by `length(unique(covRaw))`, over
#'   the RAW covariates in `unique(covRaw)` order
#' @noRd
.vaeEligible <- function(rules, fixCov, etaNames, thetaForEta, covRaw) {
  .raw <- unique(covRaw)
  .m <- matrix(TRUE, length(etaNames), length(.raw),
               dimnames = list(NULL, .raw))
  if (!isTRUE(fixCov) || length(.raw) == 0L) return(.m)
  ## A rule naming neither a parameter nor a covariate makes everything
  ## eligible, which is a direct contradiction of fixCov rather than a
  ## restriction to honor.  Say so instead of silently ignoring one of the two.
  if (any(is.na(rules$var) & is.na(rules$cov))) {
    stop("shapes: a rule with neither var= nor covar= makes every covariate ",
         "eligible, which contradicts fixCov=TRUE\n",
         "  use fixCov=FALSE to restrict shapes without restricting the search",
         call. = FALSE)
  }
  if (.vaeFixCovName %in% .raw || toupper(.vaeFixCovName) %in% toupper(.raw)) {
    stop("shapes: a data covariate is named `", .vaeFixCovName,
         "`, which collides with the eligibility flag", call. = FALSE)
  }
  .m[] <- FALSE
  for (.k in seq_along(etaNames)) {
    .al <- c(etaNames[.k], thetaForEta[.k], sub("^eta\\.", "", etaNames[.k]))
    .al <- unique(.al[!is.na(.al)])
    for (.r in seq_len(nrow(rules))) {
      .vOk <- is.na(rules$var[.r]) || rules$var[.r] %in% .al
      if (!.vOk) next
      if (is.na(rules$cov[.r])) {
        ## var-only: every covariate, on this parameter alone
        .m[.k, ] <- TRUE
      } else {
        ## `which`, not `match`: match() takes only the FIRST hit, which would
        ## leave a second same-named-up-to-case covariate ineligible.  Both entry
        ## points upper-case the data columns before the search, so that cannot
        ## arise today -- this keeps the function correct without relying on it.
        .j <- which(toupper(.raw) == rules$cov[.r])
        if (length(.j) > 0L) .m[.k, .j] <- TRUE
      }
    }
  }
  .m
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
         hockeyLow = (v < center) * (v - center),
         hockeyHi = (v >= center) * (v - center),
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
  ## A single mu-referencing transform wrapper, e.g. exp(theta + beta*cov + eta).
  ## Matched by NAME and always walking the first argument: the bounded forms
  ## take limits too (expit(x, 0, 1), logit(x, lo, hi)), and requiring a
  ## single-argument call would demote a perfectly ordinary bounded parameter to
  ## the regress M-step.
  .e <- .unwrap(.e)
  if (is.call(.e) && is.name(.e[[1L]]) && length(.e) >= 2L &&
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
           ## the knot is the centering value, so it only has to be finite -- a
           ## zero or negative knot is a perfectly ordinary place to bend
           hockey = is.finite(center),
           hockeyLow = is.finite(center),
           hockeyHi = is.finite(center),
           TRUE)
  }, logical(1), USE.NAMES = FALSE)
}

#' Per-(latent dim, column) mask of what the user allows
#'
#' Two independent restrictions land in the same mask, because the design matrix
#' is shared across latent dimensions and neither can be enforced by omitting
#' columns:
#'
#'   * ELIGIBILITY (`fixCov`) -- may this parameter carry this covariate at all?
#'     An ineligible cell has EVERY column of that covariate zeroed, categorical
#'     columns included.
#'   * PARAMETERIZATION (`shapes`) -- given that it may, which forms may it take?
#'     Categorical columns are never restricted this way: `shapes=` governs
#'     continuous parameterizations only.
#'
#' The order matters.  Eligibility is absolute, so it is applied last and
#' overrides the parameterization pass, including that pass's fallback for a
#' covariate whose requested family the data cannot support.
#' @param cov output of `.vaeCovariateSearch`
#' @param resolved output of `.vaeResolveShapes` (`$rules` + `$fixCov`); a bare
#'   rules data.frame is accepted as `fixCov = FALSE` for callers that only ever
#'   restricted parameterizations
#' @param etaNames per-latent-dim random-effect names
#' @keywords internal
#' @param thetaForEta per-latent-dim mu-referenced theta names (may be NA)
#' @return integer 0/1 matrix, `length(etaNames)` by `ncol(cov$covMat)`
#' @noRd
.vaeShapeAllowMask <- function(cov, resolved, etaNames, thetaForEta) {
  .nCov <- length(cov$covNames)
  .m <- matrix(1L, length(etaNames), .nCov)
  if (.nCov == 0L || is.null(resolved)) return(.m)
  if (is.data.frame(resolved)) resolved <- list(rules = resolved, fixCov = FALSE)
  rules <- resolved$rules
  if (is.null(rules)) return(.m)
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
  ## Eligibility last, and unconditionally: an ineligible covariate is out of the
  ## search whatever the loop above decided, including via the fallback `next`.
  .el <- .vaeEligible(rules, resolved$fixCov, etaNames, thetaForEta, cov$covRaw)
  if (!all(.el)) {
    .col <- match(cov$covRaw, colnames(.el))
    for (.k in seq_along(etaNames)) {
      .m[.k, !.el[.k, .col]] <- 0L
    }
  }
  .m
}

#' Shapes allowed for one (parameter, covariate) pair
#'
#' The most specific matching rule wins: (var, cov) beats cov-only, which beats
#' var-only, which beats the global rule.  `par` is matched against any of the
#' parameter's aliases (eta name, mu-referenced theta name, bare name).  A pair
#' no rule mentions keeps the DEFAULT shapes -- `shapes=` restricts parameterizations,
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
  if (!any(.ok)) return(.vaeDefaultShapes)
  .w <- which(.ok)
  .w <- .w[.spec[.w] == max(.spec[.w])]
  rules$shapes[[.w[length(.w)]]]
}
