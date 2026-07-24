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
#' rxode2's string-comparison path handles it.
#' @param level the level, as stored in the data
#' @return length-one character string
#' @noRd
.vaeLevelLit <- function(level) {
  if (is.numeric(level)) return(as.character(signif(level, 12)))
  .n <- suppressWarnings(as.numeric(as.character(level)))
  if (!is.na(.n)) return(as.character(signif(.n, 12)))
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
