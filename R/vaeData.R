# vaeData.R -- VAE data preparation. Builds, from the rxode2 ui + data:
#  - encoder inputs: standardized (time, DV) sequences padded to Tmax + lengths
#  - per-subject decoder inputs: event table, obs times, observed DV
#  - eta<->theta (z_pop) mapping and initial z_pop / omega / a from ini()
# Encoder-input standardization follows the reference (time/max, (DV-mean)/sd).

# Data columns never treated as covariate candidates by the VAE search
.vaeReservedCols <- c("ID", "TIME", "DV", "EVID", "AMT", "CMT", "MDV", "SS", "II",
                      "ADDL", "RATE", "DUR", "DVID", "CENS", "LIMIT", "OCC")

# Finite fallback bound (+/-) for an unbounded model-declared covariate coefficient
# regressed by the VAE M-step; keeps the 1-D optimize() step from running away.
# Generous on a transformed/log scale (an effect past this is implausible); a user
# ini() bound overrides it.
.vaeCovCoefBound <- 10

# Fallback bound half-width for an UNBOUNDED structural non-mu theta estimated by
# the M-step: `ini() estimate +/- max(.vaeNonMuThetaBound, |est| * .vaeNonMuThetaRel)`.
# Purely a divergence guard -- with +-Inf bounds a flat direction runs the estimate
# away (~1e68 on an unbounded theo_sd `tv`).  Deliberately generous so it does not
# bind at a sane optimum; a user `ini()` bound always wins.
.vaeNonMuThetaBound <- 10
.vaeNonMuThetaRel <- 3

# Max plausible log-scale effect used to derive a SCALE-AWARE fallback bound for a
# raw linear covariate coefficient (beta*COV): |beta*max|COV|| <= this, so the
# bound shrinks as the covariate magnitude grows (a raw WT coefficient is ~1/WT).
.vaeCovCoefEffect <- 5

#' Scale-aware finite fallback bounds for unbounded covariate coefficients.
#'
#' For a linear mu-ref covariate coefficient (in `muRefCovariateDataFrame`) the
#' sane magnitude scales like `1/max|COV|`, so a fixed wide bound makes the 1-D
#' optimize() overshoot a shallow interior optimum.  Scale that coefficient's
#' bound by the covariate's data magnitude; leave transformed/unrecognized
#' coefficients (already O(1) regressors) at the flat `.vaeCovCoefBound`.
#' @param ui rxode2 ui
#' @param data raw (un-normalized) modeling data
#' @param names covariate-coefficient theta names needing a bound
#' @return named numeric of positive half-widths (one per `names`)
#' @noRd
.vaeCovCoefBoundVec <- function(ui, data, names) {
  .b <- setNames(rep(.vaeCovCoefBound, length(names)), names)
  .mrc <- ui$muRefCovariateDataFrame
  if (is.null(.mrc) || nrow(.mrc) == 0L) return(.b)
  .dd <- as.data.frame(data); colnames(.dd) <- toupper(colnames(.dd))
  for (.nm in names) {
    .w <- which(.mrc$covariateParameter == .nm)
    if (length(.w) == 0L) next
    .cov <- toupper(.mrc$covariate[.w[1]])
    if (is.null(.dd[[.cov]])) next
    .mx <- max(abs(as.numeric(.dd[[.cov]])), na.rm = TRUE)
    if (is.finite(.mx) && .mx > 0) .b[.nm] <- min(.vaeCovCoefBound, .vaeCovCoefEffect / .mx)
  }
  .b
}

#' Centering value for one continuous covariate
#'
#' A `covCenter` entry (matched case-insensitively) overrides the statistic.
#' Both statistics are computed over SUBJECTS, so a subject with many
#' observations does not pull the center.
#' @noRd
.vaeCovCenterValue <- function(nm, v, covCenterType, covCenter) {
  if (!is.null(covCenter)) {
    .i <- match(toupper(nm), toupper(names(covCenter)))
    if (!is.na(.i)) return(unname(covCenter[[.i]]))
  }
  if (covCenterType == "median") stats::median(v) else mean(v)
}

#' Reference level and testable levels of a categorical covariate
#'
#' The reference is the most frequent level counted per subject (matching
#' nlmixr2scm); levels held by fewer than `catCutoff` of subjects are lumped
#' into it rather than getting a near-singular indicator of their own.
#' @return list(ref, levels, dropped)
#' @noRd
.vaeCovLevels <- function(v, catCutoff) {
  .s <- as.character(v)
  .s <- .s[!is.na(.s)]
  if (length(.s) == 0L) return(list(ref = NA_character_, levels = character(0),
                                    dropped = character(0)))
  ## table() orders alphabetically and sort() is stable, so a frequency tie
  ## resolves to the alphabetically first level -- deterministic either way
  .t <- sort(table(.s), decreasing = TRUE)
  .p <- .t / sum(.t)
  .ref <- names(.t)[1L]
  .keep <- names(.p)[.p >= catCutoff & names(.p) != .ref]
  list(ref = .ref, levels = .keep,
       dropped = setdiff(names(.p), c(.ref, .keep)))
}

#' Families a covariate's columns must span, given the shape rules
#'
#' Over-approximates on purpose: a family requested by ANY rule gets a column,
#' and the per-parameter restriction is applied later through the covAllow mask.
#' @return character vector of families, canonical order, plus the representative
#'   shape naming each one
#' @noRd
.vaeCovFamilies <- function(shapeRules, cov, center = NA_real_) {
  .sh <- unique(unlist(lapply(seq_len(nrow(shapeRules)), function(.i) {
    if (is.na(shapeRules$cov[.i]) || identical(shapeRules$cov[.i], toupper(cov))) {
      shapeRules$shapes[[.i]]
    } else character(0)
  })))
  if (length(.sh) == 0L) .sh <- .vaeDefaultShapes
  ## Drop shapes that cannot be written at this center before choosing which one
  ## names each family, so a usable sibling (e.g. "lin" for "center") wins.  When
  ## NOTHING requested is writable, substitute each family's plain form rather
  ## than dropping the covariate from the search -- the caller reports it.
  .sub <- FALSE
  if (!is.na(center)) {
    .u <- .sh[.vaeShapeUsable(.sh, center)]
    if (length(.u) > 0L) {
      .sh <- .u
    } else {
      .sh <- unname(c(log = "power", lin = "lin",
                      hockey = "hockey")[unique(.vaeShapeFamily(.sh))])
      .sub <- TRUE
    }
  }
  .fam <- .vaeShapeFamily(.sh)
  .u <- unique(.fam)
  list(family = .u, substituted = .sub,
       repShape = vapply(.u, function(.f) .sh[which(.fam == .f)[1L]],
                         character(1), USE.NAMES = FALSE))
}

#' Discover + encode subject-level covariates for the VAE search
#'
#' Returns one row per SEARCH COLUMN, not one per covariate: a continuous
#' covariate contributes one column per eligible shape family, and a categorical
#' one contributes an indicator per testable level.  `covGroup` marks columns
#' that compete for a single slot -- the alternate shapes of one covariate --
#' so at most one of them can be selected for a given parameter.
#' @param d normalized data.frame (upper-case column names)
#' @param ids unique subject ids in estimation order
#' @param shapeRules resolved `shapes=` rules (see `.vaeResolveShapes`)
#' @param covCenterType `"median"` or `"mean"`
#' @param covCenter named numeric of per-covariate centering overrides
#' @param catCutoff minimum subject proportion for a level to be testable
#' @return list(covNames, covRaw, covShape, covFamily, covLevel, covGroup,
#'   covBlock, covMat, covType, covPop, covExpr, covCanon, tvExcl, catDropped,
#'   logDrop, shapeSub, naExcl, hockeyDrop)
#' @noRd
.vaeCovariateSearch <- function(d, ids, shapeRules = NULL,
                                covCenterType = c("median", "mean"),
                                covCenter = NULL, catCutoff = 0.05) {
  covCenterType <- match.arg(covCenterType)
  if (is.null(shapeRules)) shapeRules <- .vaeResolveShapes(NULL)$rules
  N <- length(ids)
  ## auto-discover subject-level covariate candidates (constant within ID),
  ## excluding reserved data columns
  .cand <- setdiff(names(d), .vaeReservedCols)
  .usable <- function(nm) {
    v <- d[[nm]]
    is.numeric(v) || is.character(v) || is.factor(v)
  }
  .allCand <- .cand[vapply(.cand, .usable, logical(1))]
  ## First row of each subject, in estimation order.  A subject's value is read
  ## from its first NON-missing row: a covariate is often blank on dose records,
  ## and taking the literal first row would call it missing (or, before,
  ## time-varying) when it is neither.
  .first <- match(ids, d$ID)
  .subjVal <- function(nm) {
    .v <- d[[nm]][.first]
    if (anyNA(.v)) {
      ## one match() over the non-missing rows rather than a scan per subject,
      ## so this stays linear in the data even when many first rows are blank
      .all <- d[[nm]]
      .ok <- which(!is.na(.all))
      if (length(.ok) > 0L) {
        .hit <- .ok[match(ids, d$ID[.ok])]
        .w <- which(is.na(.v) & !is.na(.hit))
        if (length(.w) > 0L) .v[.w] <- .all[.hit[.w]]
      }
    }
    .v
  }
  ## Constant WITHIN a subject, ignoring missing entries -- an NA on one row does
  ## not make an otherwise fixed covariate time-varying.
  .isSubjConst <- vapply(.allCand, function(nm) {
    all(vapply(ids, function(id) {
      .u <- unique(d[[nm]][d$ID == id])
      length(.u[!is.na(.u)]) <= 1L
    }, logical(1)))
  }, logical(1))
  # The VAE covariate search absorbs covariates as subject-level (time-invariant)
  # effects, so a covariate that varies within a subject (time-varying) cannot be
  # searched.  Unlike saem/mu-focei (whose covariates are declared in the model,
  # detectable via .nlmixrTimeVaryingCovariates), the VAE search scans every
  # eligible data column, so time-varying ones are those that are not
  # subject-constant.  They are reported back (tvExcl) so callers can warn.
  .raw <- .allCand[.isSubjConst]
  .tvExcl <- setdiff(.allCand, .raw)
  ## The design matrix feeds an ordinary least squares M-step, so every column
  ## must be complete and finite.  A covariate missing for even one subject
  ## cannot be searched: silently imputing it would invent data, and the written
  ## model has no ifelse() guard to carry an imputation to solve time.  (An NA
  ## also used to reach `all(v > 0)` and abort the whole fit with "missing value
  ## where TRUE/FALSE needed".)
  .complete <- vapply(.raw, function(nm) {
    v <- .subjVal(nm)
    !anyNA(v) && (!is.numeric(v) || all(is.finite(v)))
  }, logical(1))
  .naExcl <- .raw[!.complete]
  .raw <- .raw[.complete]
  ## A covariate with a single distinct value carries no information, and as a
  ## constant 0/1 indicator it would enter the design as a column of ones --
  ## duplicating the intercept and making the least-squares M-step singular.
  .raw <- .raw[vapply(.raw, function(nm) length(unique(.subjVal(nm))) > 1L,
                      logical(1))]

  .cols <- list()
  .rawVals <- list()
  .catDropped <- character(0)
  .logDrop <- character(0)
  .shapeSub <- character(0)
  .hockeyDrop <- character(0)
  ## `group` and `block` are both KEYS here; each becomes an integer id once
  ## every column is known.  Columns sharing a GROUP are alternate forms of one
  ## relationship and compete for a single slot; columns sharing a BLOCK are the
  ## arms of ONE relationship and are selected all-or-none.  `block` defaults to
  ## the column's own name, so every column is its own block unless it says
  ## otherwise -- which reproduces the historic per-column search exactly.
  .add <- function(name, raw, shape, family, level, group, values, type,
                   center, expr, block = NULL) {
    ## a level containing "_" can collide with another covariate's column name
    ## (covariate WT level "A_B" vs covariate WT_A level "B"), and a duplicate
    ## name makes the second column unreachable through match()
    name <- .vaeUniqueName(name, vapply(.cols, `[[`, character(1), "name"))
    .cols[[length(.cols) + 1L]] <<-
      list(name = name, raw = raw, shape = shape, family = family,
           level = level, group = group, values = values, type = type,
           center = center, expr = expr,
           block = if (is.null(block)) name else block)
  }

  for (nm in .raw) {
    v <- .subjVal(nm)
    ## kept so a pinned column can be rebuilt as the model's OWN expression
    .rawVals[[nm]] <- v
    ## mu2/mu3 covariates are pre-transformed by the hook into a linear
    ## nlmixrMuDerCov# data column (the centering/transform is already baked in),
    ## so encode them LINEARLY -- never re-apply a log transform, and never split
    ## them into shape families.
    .isMuDer <- grepl("^NLMIXRMUDERCOV[0-9]+$", nm, ignore.case = TRUE)
    .isNum <- is.numeric(v)
    ## a 0/1 indicator column (e.g. SEXF) is already in its natural
    ## parameterization: leave it RAW so its coefficient is the level-1 shift and
    ## the structural theta stays the reference (0) value.  (`%in%` yields FALSE
    ## for NA, so this is already a strict TRUE/FALSE; isTRUE makes that explicit.)
    .isInd <- .isNum && isTRUE(all(v %in% c(0, 1)))
    if (.isMuDer) {
      .ctr <- .vaeCovCenterValue(nm, v, covCenterType, covCenter)
      .add(nm, nm, "lin", "lin", NA_character_, nm, v - .ctr, "categorical",
           .ctr, .vaeShapeExpr("lin", nm, .ctr))
    } else if (.isInd) {
      .add(nm, nm, "cat", "cat", NA_character_, nm, v, "categorical", 0,
           .vaeShapeExpr("cat", nm, raw = TRUE))
    } else if (.isNum && length(unique(v)) > 2L) {
      ## continuous: one column per eligible shape family
      .ctr <- .vaeCovCenterValue(nm, v, covCenterType, covCenter)
      ## a shape that cannot be expressed at this center (e.g. "center" when the
      ## center is 0) must not name a column -- it would be written back as a
      ## division by zero with the coefficient rescaled to nothing
      .fam <- .vaeCovFamilies(shapeRules, nm, .ctr)
      if (isTRUE(.fam$substituted)) .shapeSub <- c(.shapeSub, nm)
      .canLog <- all(v > 0) && .ctr > 0
      ## A hockey knot with (almost) nothing on one side leaves that arm a column
      ## of zeros, which makes the least-squares M-step singular.  The median
      ## splits the subjects in half by construction, so this only bites a
      ## covCenter= override; either way the arms are dropped rather than fitted.
      ## the floor of 1 is not the same rule as catCutoff: `catCutoff = 0` is a
      ## documented setting ("test every level"), and against a bare proportion a
      ## knot outside the data range would pass with an EMPTY side
      .canHockey <- min(sum(v < .ctr), sum(v >= .ctr)) >= max(1, catCutoff * N)
      .f <- .fam$family; .r <- .fam$repShape
      .keep <- (.f != "log" | .canLog) & (.f != "hockey" | .canHockey)
      if (any(.f == "hockey" & !.canHockey)) .hockeyDrop <- c(.hockeyDrop, nm)
      if (!any(.keep)) {
        ## every requested family is undefined here (log shapes on a covariate
        ## with non-positive values, or a hockey knot with an empty side); fall
        ## back to the linear family rather than silently dropping the covariate
        ## from the search
        if (any(.f == "log" & !.canLog)) .logDrop <- c(.logDrop, nm)
        .f <- "lin"; .r <- "lin"
      } else {
        if (any(.f == "log" & !.canLog)) .logDrop <- c(.logDrop, nm)
        .f <- .f[.keep]; .r <- .r[.keep]
      }
      for (.i in seq_along(.f)) {
        ## every shape family of this covariate shares the covariate's key, so
        ## the search may take at most one of them
        if (.f[.i] == "hockey") {
          ## one column per ARM, both keyed to a single block so the search takes
          ## them together or not at all.  Selecting one arm beside the lin column
          ## would span the same space for the same penalty -- an exact tie whose
          ## winner is an arbitrary parameterization that does not read as a
          ## hockey stick.
          .blk <- paste0(nm, "|hockey")
          for (.arm in .vaeHockeyArms) {
            .add(paste0(nm, "_", .arm), nm, .arm, "hockey", NA_character_, nm,
                 .vaeShapeValue(.arm, v, .ctr), "continuous", .ctr,
                 .vaeShapeExpr(.arm, nm, .ctr), block = .blk)
          }
        } else {
          .val <- if (.f[.i] == "log") log(v / .ctr) else v - .ctr
          .add(paste0(nm, "_", .r[.i]), nm, .r[.i], .f[.i], NA_character_, nm,
               .val, "continuous", .ctr, .vaeShapeExpr(.r[.i], nm, .ctr))
        }
      }
    } else {
      ## categorical: an indicator per testable level, reference = modal level.
      ## Each level is its OWN group -- several levels of one factor may enter a
      ## parameter together, they are not alternate forms of one relationship.
      .lv <- .vaeCovLevels(v, catCutoff)
      if (length(.lv$dropped) > 0L) {
        .catDropped <- c(.catDropped, paste0(nm, "=", .lv$dropped))
      }
      .s <- as.character(v)
      for (.l in .lv$levels) {
        ## a distinct key per level, so levels never exclude one another
        .add(paste0(nm, "_", .l), nm, "cat", "cat", .l, paste0(nm, "|", .l),
             as.numeric(!is.na(.s) & .s == .l), "categorical", 0,
             .vaeShapeExpr("cat", nm, level = .vaeCovLevelValue(v, .l)))
      }
    }
  }

  .nc <- length(.cols)
  .covNames <- vapply(.cols, `[[`, character(1), "name")
  .covMat <- matrix(0, N, .nc, dimnames = list(NULL, .covNames))
  for (.i in seq_len(.nc)) .covMat[, .i] <- .cols[[.i]]$values
  ## group ids are derived from the keys rather than counted as columns are
  ## emitted, so a covariate that contributes NO column cannot shift them
  .key <- vapply(.cols, `[[`, character(1), "group")
  .covGroup <- match(.key, unique(.key))
  .bkey <- vapply(.cols, `[[`, character(1), "block")
  .covBlock <- match(.bkey, unique(.bkey))
  list(covNames = .covNames,
       covRaw = vapply(.cols, `[[`, character(1), "raw"),
       covShape = vapply(.cols, `[[`, character(1), "shape"),
       covFamily = vapply(.cols, `[[`, character(1), "family"),
       covLevel = vapply(.cols, `[[`, character(1), "level"),
       covGroup = .covGroup,
       covBlock = .covBlock,
       covMat = .covMat,
       covType = vapply(.cols, `[[`, character(1), "type"),
       covPop = vapply(.cols, `[[`, numeric(1), "center"),
       covExpr = vapply(.cols, `[[`, character(1), "expr"),
       ## the encoder head takes one column per GROUP: alternate shapes are
       ## near-collinear copies, but distinct factor levels are not
       covCanon = !duplicated(.covGroup), covRawVal = .rawVals,
       tvExcl = .tvExcl, catDropped = .catDropped, logDrop = unique(.logDrop),
       shapeSub = unique(.shapeSub), naExcl = .naExcl,
       hockeyDrop = unique(.hockeyDrop))
}

#' Search column a model-declared covariate pair pins to
#'
#' Pinning keeps the user's model text verbatim, so a declared pair may only
#' occupy the column whose family matches the form it was WRITTEN in: a
#' `log(cov/center)` effect takes the log-family column, anything else (a linear
#' effect, a mu2/mu3 pre-transformed column, an indicator) takes the linear or
#' indicator column.  Returns `NA` when the covariate has no such column.
#' @param cov output of `.vaeCovariateSearch`
#' @param pair one row of the declared-pair table
#' @noRd
.vaePinColumn <- function(cov, pair) {
  .w <- which(cov$covRaw == pair$covName)
  if (length(.w) == 0L) return(NA_integer_)
  ## the family the coefficient was WRITTEN in, so `beta*(WT - 70)` pins to the
  ## linear column just as `beta*log(WT/70)` pins to the log one
  ## A bare linear multiplier fits either an indicator column or a mu2/mu3
  ## pre-transformed one, which the search stores as the linear family.
  .want <- if (!is.null(pair$family) && !is.na(pair$family)) {
    if (identical(pair$family, "cat")) c("cat", "lin") else pair$family
  } else if (identical(pair$covType, "continuous")) "log" else c("lin", "cat")
  .m <- .w[cov$covFamily[.w] %in% .want]
  ## A written level comparison pins to THAT level's indicator.  If the level has
  ## no column -- it was lumped into the reference by catCutoff, or is absent from
  ## the data -- there is nothing to pin to, so fail rather than fall back to some
  ## other level's column and fit the wrong indicator.
  if (!is.null(pair$level) && !is.na(pair$level)) {
    .m <- .m[!is.na(cov$covLevel[.m]) & cov$covLevel[.m] == pair$level]
  }
  if (length(.m) == 0L) NA_integer_ else .m[1L]
}

#' The level value to compare against in generated model text
#'
#' Keeps a numeric code numeric so the emitted comparison is numeric, while a
#' character/factor level is compared as a string.
#' @noRd
.vaeCovLevelValue <- function(v, level) {
  if (is.numeric(v)) {
    .n <- suppressWarnings(as.numeric(level))
    if (!is.na(.n)) return(.n)
  }
  level
}

#' Covariates explored by the VAE covariate search
#'
#' Returns the candidate columns that `nlmixr2(..., est = "vae")` would explore
#' during automated covariate selection, using the same discovery rules as the
#' fit: every non-reserved data column that is constant within each subject is a
#' candidate.  A numeric candidate with more than two unique values is
#' continuous and contributes one column per eligible shape; anything else is
#' categorical and contributes an indicator per testable level.  Columns sharing
#' a `group` are alternate shapes of one covariate, so at most one of them can
#' enter a given parameter.  Time-varying columns cannot be searched and are
#' excluded with a warning.
#'
#' @param data estimation dataset containing at least an `ID` column; column
#'   names are matched case-insensitively, as in the VAE fit
#' @param warn when `TRUE` (default) warn about time-varying columns excluded
#'   from the search; when `FALSE` exclude them silently
#' @param shapes,covCenterType,covCenter,catCutoff as in [vaeControl()]; control
#'   which shapes are explored and how covariates are centered
#' @return a data frame with one row per candidate search column and columns
#'   `covariate` (the column name), `raw` (upper-cased data column it comes
#'   from), `shape`, `level` (for categorical indicators), `group` (mutual
#'   exclusion group), `block` (columns selected all-or-none, i.e. the two arms
#'   of a `"hockey"` relationship), `type` and `center`; zero rows when nothing
#'   qualifies
#' @export
#' @author Matthew L. Fidler
#' @examples
#' d <- data.frame(id = rep(1:3, each = 2), time = rep(0:1, 3), dv = rnorm(6),
#'                 wt = rep(c(70, 80, 60), each = 2),
#'                 sex = rep(c(0, 1, 0), each = 2))
#' vaeCovariates(d)
#'
#' # restrict the explored shapes
#' vaeCovariates(d, shapes = "power")
vaeCovariates <- function(data, warn = TRUE,
                          shapes = c("power", "lin", "log", "identity", "center", "hockey"),
                          covCenterType = c("median", "mean"),
                          covCenter = NULL, catCutoff = 0.05) {
  checkmate::assertLogical(warn, len = 1, any.missing = FALSE)
  d <- as.data.frame(data)
  names(d) <- toupper(names(d))
  if (is.null(d$ID)) {
    stop("'data' must contain an ID column", call. = FALSE)
  }
  .cov <- .vaeCovariateSearch(d, unique(d$ID), .vaeResolveShapes(shapes)$rules,
                              match.arg(covCenterType), covCenter, catCutoff)
  if (warn && length(.cov$tvExcl) > 0L) {
    warning("time-varying covariate(s) were excluded from automatic covariate search: ",
            paste(.cov$tvExcl, collapse = ", "), call. = FALSE)
  }
  if (warn && length(.cov$naExcl) > 0L) {
    warning("covariate(s) with missing values were excluded from automatic covariate search: ",
            paste(.cov$naExcl, collapse = ", "), call. = FALSE)
  }
  if (warn && length(.cov$hockeyDrop) > 0L) {
    warning("hockey skipped, <5% of subjects one side of the knot: ",
            paste(.cov$hockeyDrop, collapse = ", "), call. = FALSE)
  }
  data.frame(covariate = .cov$covNames, raw = .cov$covRaw, shape = .cov$covShape,
             level = .cov$covLevel, group = .cov$covGroup,
             block = .cov$covBlock, type = .cov$covType,
             center = .cov$covPop, row.names = NULL)
}

#' Detect a clean `log(cov/center)` form for `cov` inside an expression.
#'
#' Walks the parse tree; returns the divisor `center` when `cov` occurs as
#' `log(cov)` (center 1) or `log(cov/<numeric literal>)`.  Any other form
#' (raw `cov`, `log(cov/expr)`, `log(a*cov)`) yields `inLog=FALSE` or a `NA`
#' center so the caller can fall back to the in-place regress M-step.
#' @noRd
.vaeLogCenter <- function(e, cov) {
  if (is.call(e)) {
    if (identical(e[[1L]], as.name("log")) && length(e) == 2L &&
          cov %in% all.vars(e[[2L]])) {
      .a <- e[[2L]]
      if (is.name(.a) && identical(as.character(.a), cov)) {
        return(list(inLog = TRUE, center = 1))
      }
      if (is.call(.a) && identical(.a[[1L]], as.name("/")) && length(.a) == 3L &&
            is.name(.a[[2L]]) && identical(as.character(.a[[2L]]), cov) &&
            is.numeric(.a[[3L]]) && length(.a[[3L]]) == 1L && is.finite(.a[[3L]])) {
        return(list(inLog = TRUE, center = as.numeric(.a[[3L]])))
      }
      return(list(inLog = TRUE, center = NA_real_))
    }
    for (.i in seq_along(e)[-1L]) {
      .r <- .vaeLogCenter(e[[.i]], cov)
      if (isTRUE(.r$inLog)) return(.r)
    }
  }
  list(inLog = FALSE, center = NA_real_)
}

#' Covariate that a coefficient multiplies, within an expression.
#'
#' Disambiguates which data covariate `coef` pairs with when a model line carries
#' several covariate effects (e.g. `wt.cl*log(WT/70) + sex.cl*SEX`): walks to the
#' `*` term containing `coef` and returns the single covariate on the other side.
#' `NULL` when it cannot be resolved to exactly one covariate.
#' @noRd
.vaeCoefCov <- function(e, coef, covs) {
  if (is.call(e)) {
    if (identical(e[[1L]], as.name("*")) && length(e) == 3L) {
      .lv <- all.vars(e[[2L]]); .rv <- all.vars(e[[3L]])
      if (coef %in% .lv) { .c <- intersect(.rv, covs); if (length(.c) == 1L) return(.c) }
      if (coef %in% .rv) { .c <- intersect(.lv, covs); if (length(.c) == 1L) return(.c) }
    }
    for (.i in seq_along(e)[-1L]) {
      .r <- .vaeCoefCov(e[[.i]], coef, covs)
      if (!is.null(.r)) return(.r)
    }
  }
  NULL
}

#' Model-declared covariate/parameter pairs for pinned VAE selection.
#'
#' One row per model-written covariate coefficient, resolving which latent dim
#' (eta) `k` it modifies, the data covariate, the user's coefficient theta, and
#' whether the pair can be handled by the restricted branch-and-bound search
#' (`inPool`) -- i.e. the covariate is in the subject-level search pool AND its
#' written functional form matches the VAE encoding so the estimated slope
#' transfers directly (continuous->`log(cov/center)`, categorical->linear).
#' Pairs that are not `inPool` (out-of-pool covariate, or a form whose slope
#' would not transfer) are estimated in place by the regress M-step instead.
#' @param ui rxode2 ui
#' @param cov output of `.vaeCovariateSearch`; `covName` below is the RAW data
#'   column, resolved to a specific search column later by `.vaePinColumn`
#' @return data frame (k, covName, coefName, thetaName, covType, userCenter,
#'   inPool), or `NULL` when the model declares no covariate effects
#' @noRd
.vaeModelCovariatePairs <- function(ui, cov) {
  .coefThetas <- .vaeCovariateCoefThetas(ui)
  if (length(.coefThetas) == 0L) return(NULL)
  ## A declared pair names a RAW data column, while the search pool holds one
  ## column per shape/level.  Resolve against the raw names here and let
  ## .vaePinColumn pick the column matching the written form.
  covNames <- unique(cov$covRaw)
  covType <- cov$covType[match(covNames, cov$covRaw)]
  .thetaForEta <- .foceiEtaThetaMap(ui)$thetaForEta
  .thetaPool <- .thetaForEta[!is.na(.thetaForEta)]
  .allCov <- ui$allCovs
  if (is.null(.allCov)) .allCov <- character(0)
  .mrc <- ui$muRefCovariateDataFrame
  .lst <- ui$lstExpr
  .rows <- vector("list", 0L)
  for (.coef in .coefThetas) {
    .thName <- NA_character_; .covTok <- NA_character_; .linear <- FALSE
    if (!is.null(.mrc) && nrow(.mrc) > 0L && .coef %in% .mrc$covariateParameter) {
      .r <- .mrc[.mrc$covariateParameter == .coef, , drop = FALSE][1L, ]
      .thName <- as.character(.r$theta)
      ## rxode2 may record an ALGEBRAIC covariate expression here (mu2-style,
      ## e.g. "log(0.0142857 * WT)") rather than a bare data column.  Only take
      ## it as a plain linear effect when it names a pool covariate directly;
      ## otherwise fall through to the model-line scan below.
      .cand <- as.character(.r$covariate)
      if (!is.na(match(toupper(.cand), covNames))) {
        .covTok <- .cand
        .linear <- TRUE
      }
    }
    if (is.na(.covTok)) {
      .lines <- Filter(function(e) .coef %in% all.vars(e), .lst)
      if (length(.lines) == 0L) next
      .vars <- all.vars(.lines[[1L]])
      if (is.na(.thName)) {
        .thHit <- intersect(.thetaPool, .vars)
        if (length(.thHit) == 1L) .thName <- .thHit
      }
      ## a line may carry several covariate effects (e.g.
      ## wt.cl*log(WT/70) + sex.cl*SEX): pick the covariate THIS coefficient
      ## multiplies rather than skipping the coefficient (skipping could drop
      ## pinning to the unrestricted full search).  Unresolved -> not pinnable.
      .cc <- .vaeCoefCov(.lines[[1L]], .coef, .allCov)
      if (!is.null(.cc)) .covTok <- .cc
    }
    ## Always emit a row for a detected coefficient so pinning stays restrictive;
    ## a pair that cannot be resolved/transferred is marked not `inPool` and
    ## estimated in place by the regress M-step.
    .k <- match(.thName, .thetaForEta)
    .j <- if (!is.na(.covTok)) match(toupper(.covTok), covNames) else NA_integer_
    .inPool <- !is.na(.k) && !is.na(.j)
    .ct <- if (.inPool) covType[.j] else NA_character_
    .userCenter <- NA_real_
    .shape <- NA_character_
    .level <- NA_character_
    if (.inPool) {
      if (identical(.ct, "categorical")) {
        ## A plain linear beta*cov transfers (the slope is invariant to the
        ## shift).  An explicit level comparison -- which is what the VAE itself
        ## writes for a factor -- transfers too, and must pin to THAT level's
        ## indicator column.  Anything else is not slope-transferable.
        .cl <- Filter(function(e) .coef %in% all.vars(e), .lst)
        .dc <- if (length(.cl)) {
          .vaeDetectShape(.vaeCoefFactor(.cl[[1L]], .coef), .covTok)
        } else list(shape = NA_character_, level = NULL)
        if (.linear) {
          .userCenter <- 0; .shape <- "cat"
        } else if (identical(.dc$shape, "cat")) {
          .userCenter <- 0; .shape <- "cat"
          if (!is.null(.dc$level)) {
            .level <- as.character(if (is.character(.dc$level)) .dc$level
                                   else deparse(.dc$level))
          }
        } else {
          .inPool <- FALSE
        }
      } else {
        ## Continuous: read the SHAPE the coefficient was written in, so a model
        ## the VAE itself wrote (any of power/lin/log/identity/center) pipes back
        ## into another fit and pins to that same shape.  Anything unrecognized
        ## is not slope-transferable and goes to the regress M-step.
        .cl <- Filter(function(e) .coef %in% all.vars(e), .lst)
        ## no `coef * <expr>` factor -> unresolved -> not pinnable (the safe
        ## direction: the coefficient is estimated in place instead)
        .ds <- if (length(.cl)) {
          .vaeDetectShape(.vaeCoefFactor(.cl[[1L]], .coef), .covTok)
        } else list(shape = NA_character_, center = NA_real_)
        if (!is.na(.ds$shape) && is.finite(.ds$center)) {
          .shape <- .ds$shape
          .userCenter <- .ds$center
        } else {
          .inPool <- FALSE
        }
      }
    }
    .rows[[length(.rows) + 1L]] <- data.frame(
      k = if (is.na(.k)) NA_integer_ else as.integer(.k),
      covName = if (.inPool) covNames[.j] else if (is.na(.covTok)) NA_character_ else toupper(.covTok),
      coefName = .coef, thetaName = if (is.na(.thName)) NA_character_ else .thName,
      shape = .shape, level = .level,
      family = if (is.na(.shape)) NA_character_ else .vaeShapeFamily(.shape),
      covType = if (is.na(.ct)) NA_character_ else .ct,
      userCenter = .userCenter, inPool = .inPool,
      stringsAsFactors = FALSE)
  }
  if (length(.rows) == 0L) return(NULL)
  do.call(rbind, .rows)
}

#' Which regressed thetas may be optimized in `residOptimize="twoStage"` stage 2.
#'
#' Stage 2 pins the ODE states from stage 1 and re-optimizes with only the
#' candidate moving, so a parameter is eligible exactly when the state trajectory
#' cannot depend on it.  The rule is PER PARAMETER:
#'
#'   eligible = it is an `err` parameter, OR no solve-defining expression reaches it
#'
#' The `err` half is the historic rule and keeps every error parameter in stage 2
#' exactly as before.  The second half is what an `ll()`/generalized endpoint
#' needs: its residual-like parameters are plain thetas with no `err` row, so on
#' the old rule stage 2 was empty and `"twoStage"` silently degraded to the joint
#' `"optimize"` solve.
#'
#' A model-level "does this model have error parameters" short-circuit would NOT
#' do -- a multi-endpoint model with one Gaussian and one `ll()` endpoint has
#' `err` rows AND log-density-only thetas, and both belong in stage 2.
#'
#' Deliberately conservative in the risky direction: an unparsable model, or one
#' with no ODE states to pin, contributes nothing beyond the `err` set.
#' @param ui rxode2 ui object
#' @param regressNames regressed theta names, in `regressThetaIdx0` order
#' @param regressErrIdx0 0-based slot in `a` per regressed name, -1 when not one
#' @return integer 0/1 per regressed name
#' @noRd
.vaeRegressStage2 <- function(ui, regressNames, regressErrIdx0) {
  if (length(regressNames) == 0L) return(integer(0))
  ## Recycling here would silently mis-mask: a structural theta labelled stage 2
  ## gets optimized against a frozen ODE, which is wrong rather than merely slow.
  ## The two are built together in .vaeDataPrep, so a mismatch is a caller bug --
  ## fail loudly at prep time instead.
  if (length(regressErrIdx0) != length(regressNames)) {
    stop("vae: regressErrIdx0 (", length(regressErrIdx0), ") must match ",
         "regressNames (", length(regressNames), ")", call. = FALSE)
  }
  .isErr <- regressErrIdx0 >= 0L
  .odeFree <- tryCatch(.vaeOdeFreeThetas(ui, regressNames),
                       error = function(e) rep(FALSE, length(regressNames)))
  as.integer(.isErr | .odeFree)
}

#' Names the ODE state trajectory provably cannot depend on.
#'
#' Fixpoint over the model assignments, seeded from every expression that feeds
#' the solve.  Assignment ORDER is ignored, which can only over-collect symbols
#' and therefore only ever answers FALSE where the truth is TRUE -- the safe
#' direction (the theta stays in stage 1, as it is today).
#' @param ui rxode2 ui object
#' @param thetaNames candidate names
#' @return logical vector, `TRUE` when no solve-defining expression reaches the name
#' @noRd
.vaeOdeFreeThetas <- function(ui, thetaNames) {
  .no <- rep(FALSE, length(thetaNames))
  .lst <- tryCatch(ui$lstExpr, error = function(e) NULL)
  if (is.null(.lst) || length(.lst) == 0L) return(.no)
  .states <- tryCatch(rxode2::rxState(ui), error = function(e) character(0))
  if (length(.states) == 0L) return(.no)
  ## A linCmt() model solves compartments from parameters read by NAME (cl, v,
  ## ka, ...) that are not syntactically connected to the linCmt() call, so the
  ## assignment-graph scan cannot trace them.  When a linCmt() appears alongside
  ## ODE states (e.g. linCmt PK + a PD compartment), classify NOTHING as ODE-free
  ## -- everything stays in stage 1, and an error parameter is still stage-2
  ## eligible via the err rule.  (predDf$linCmt is FALSE here because the endpoint
  ## is the ODE state, so scan the expressions.)
  .hasLinCmt <- any(vapply(.lst, function(.e)
    grepl("(^|[^A-Za-z0-9._])linCmt[BAC]? *\\(",
          paste(deparse(.e), collapse = " ")), logical(1)))
  if (.hasLinCmt) return(.no)
  .isAssign <- function(.ex) is.call(.ex) && length(.ex) == 3L &&
    (identical(.ex[[1]], as.name("<-")) || identical(.ex[[1]], as.name("=")) ||
       identical(.ex[[1]], as.name("~")))
  .syms <- function(.e) {
    if (is.name(.e)) return(as.character(.e))
    if (is.call(.e)) return(unlist(lapply(as.list(.e)[-1L], .syms), use.names = FALSE))
    character(0)
  }
  ## a solve-defining left-hand side: d/dt(x), x(0), and the dosing modifiers
  .solveLhs <- function(.txt) {
    grepl("^d */ *dt *\\(", .txt) ||
      grepl("^[A-Za-z._][A-Za-z0-9._]* *\\( *0 *\\)$", .txt) ||
      grepl("^(f|alag|lag|rate|dur) *\\(", .txt)
  }
  ## `x_0 <- ...` is rxode2's other spelling of the `x(0) <- ...` initial
  ## condition.  It is a plain NAME, so without this it would look like an
  ## ordinary intermediate and its rhs would never seed the solve.
  .initNames <- paste0(.states, "_0")
  .seed <- character(0); .map <- list()
  .add <- function(.ex) {
    .rhs <- .syms(.ex[[3]])
    if (is.name(.ex[[2]])) {
      ## Key by as.character(), the SAME way .syms() renders a reference, so the
      ## fixpoint does not depend on deparse quoting.
      .nm <- as.character(.ex[[2]])
      if (.nm %in% .initNames) { .seed <<- c(.seed, .rhs); return(invisible()) }
      ## Every other name assignment -- including a `~` endpoint line -- becomes a
      ## map edge.  Do NOT special-case the endpoint variable: `conc ~ central / v`
      ## can define a variable that a `d/dt()` also reads, and dropping it would
      ## hide the dependency.  Keeping the edge is at worst conservative (an error
      ## parameter reached this way stays in stage 1, and it is stage-2 eligible
      ## via the err rule anyway).
      .map[[.nm]] <<- unique(c(.map[[.nm]], .rhs))
      return(invisible())
    }
    .txt <- paste(deparse(.ex[[2]]), collapse = "")
    ## `ll(x) ~ <density>` is the likelihood; seeding it would make every
    ## log-density symbol look solve-reachable and defeat the whole scan.
    if (identical(.ex[[1]], as.name("~")) && grepl("^ll *\\(", .txt)) return(invisible())
    ## solve-defining, or an lhs shape not recognized here: treat the rhs as
    ## reachable rather than guess (the safe direction: stage 1)
    .seed <<- c(.seed, .rhs)
    invisible()
  }
  ## Walk into control flow: rxode2 models may wrap assignments in `if`/`else`
  ## blocks, and an assignment feeding a d/dt inside one must still be collected.
  ## A gating CONDITION is seeded too -- it decides whether a d/dt runs, so the
  ## solve depends on it.
  .walk <- function(.ex) {
    if (!is.call(.ex)) return(invisible())
    if (.isAssign(.ex)) return(.add(.ex))
    if (identical(.ex[[1]], as.name("if")) || identical(.ex[[1]], as.name("while"))) {
      .seed <<- c(.seed, .syms(.ex[[2]]))
      for (.k in seq_along(.ex)[-(1:2)]) .walk(.ex[[.k]])
      return(invisible())
    }
    for (.k in seq_along(.ex)[-1L]) .walk(.ex[[.k]])
    invisible()
  }
  for (.ex in .lst) .walk(.ex)
  .seen <- unique(.seed); .todo <- .seen
  while (length(.todo) > 0L) {
    .nxt <- unique(unlist(.map[intersect(.todo, names(.map))], use.names = FALSE))
    .todo <- setdiff(.nxt, .seen)
    .seen <- c(.seen, .todo)
  }
  !(as.character(thetaNames) %in% .seen)
}

#' Prepare VAE inputs from a ui + data
#' @param ui rxode2 ui object
#' @param data estimation data (ID/TIME/DV/EVID/... columns)
#' @return list of prepared VAE inputs
#' @noRd
.vaeDataPrep <- function(ui, data, control = NULL) {
  .inputScale <- if (is.null(control$inputScale)) "reference" else control$inputScale
  .idf <- ui$iniDf
  .map <- .foceiEtaThetaMap(ui)
  .etaNames <- .map$etaNames
  .neta <- length(.etaNames)
  if (.neta == 0L) stop("est=\"vae\" requires at least one random effect", call. = FALSE)

  ## full theta vector (THETA_i_ in ntheta order), from ini estimates
  .thRows <- .idf[!is.na(.idf$ntheta), , drop = FALSE]
  .thRows <- .thRows[order(.thRows$ntheta), , drop = FALSE]
  .th <- setNames(as.numeric(.thRows$est), paste0("THETA_", seq_len(nrow(.thRows)), "_"))
  ## structural theta index (in the full theta vector) paired with each eta.
  ## A random effect that is not mu-referenced to a single theta -- a mixture eta
  ## (mix(exp(lke1+eta.ke),p,exp(lke2+eta.ke))), an eta on a fixed (literalFix-ed)
  ## theta, or a genuinely free eta -- is modeled as theta+eta with theta forced
  ## to 0: it centers at 0, is held there by the M-step, and is excluded from
  ## covariate selection; the rest of the model (literal / component thetas /
  ## covariate expression) carries its structure.
  .zPopThetaIdx <- match(.map$thetaForEta, .thRows$name)
  .isFree <- is.na(.zPopThetaIdx)
  .zPop <- numeric(.neta)                                      # structural population means (transformed)
  .zPop[!.isFree] <- as.numeric(.th[.zPopThetaIdx[!.isFree]])
  ## latent dims whose backing structural theta is FIXED (e.g. nonMuTheta="fix", or
  ## a user-fixed theta carrying an eta): the M-step holds their typical value at
  ## the ini() value, and they are dropped from the iteration print.
  .zPopFix <- logical(.neta)
  .zPopFix[!.isFree] <- isTRUE2(.thRows$fix[.zPopThetaIdx[!.isFree]])

  ## omega init (diagonal) for the etas + which variances are FIXED (held by the
  ## M-step, not estimated)
  .omega <- vapply(.etaNames, function(nm) {
    .r <- .idf[!is.na(.idf$neta1) & .idf$neta1 == .idf$neta2 & .idf$name == nm, , drop = FALSE]
    as.numeric(.r$est[1])
  }, numeric(1))
  .omegaFix <- vapply(.etaNames, function(nm) {
    .r <- .idf[!is.na(.idf$neta1) & .idf$neta1 == .idf$neta2 & .idf$name == nm, , drop = FALSE]
    isTRUE(as.logical(.r$fix[1]))
  }, logical(1))
  ## full ini omega block (declared off-diagonals included) + per-entry fix; the
  ## M-step estimates every nonzero entry of this structure
  .omBlock <- .omegaBlockFromIniDf(.idf, .etaNames)
  ## structural-theta bounds per eta (Inf/-Inf when unbounded or free): the M-step
  ## clamps the population estimate to [lower, upper], giving the constrained
  ## estimate (at the bound when the unconstrained optimum is outside).
  .zPopLower <- rep(-Inf, .neta); .zPopUpper <- rep(Inf, .neta)
  .zPopLower[!.isFree] <- as.numeric(.thRows$lower[.zPopThetaIdx[!.isFree]])
  .zPopUpper[!.isFree] <- as.numeric(.thRows$upper[.zPopThetaIdx[!.isFree]])

  ## normalize data columns (needed early for covariate discovery + pinning)
  d <- as.data.frame(data)
  names(d) <- toupper(names(d))
  ## no EVID: derive from AMT; with no AMT column either (dose-free data), all
  ## rows are observations (d$AMT is NULL -> the ifelse would yield length 0).
  if (is.null(d$EVID)) {
    d$EVID <- if (is.null(d$AMT)) rep(0L, nrow(d)) else ifelse(is.na(d$AMT) | d$AMT == 0, 0L, 1L)
  }
  .ids <- unique(d$ID)
  N <- length(.ids)

  ## subject-level covariate discovery + encoding (shared with vaeCovariates())
  ## tolerate control lists that predate the shape settings (older serialized
  ## objects), which then reproduce the historic single-shape search
  .cct <- if (is.null(control$covCenterType)) "mean" else control$covCenterType
  .cco <- if (is.null(control$catCutoff)) 0.05 else control$catCutoff
  .csh <- if (is.null(control$shapes)) "power" else control$shapes
  .resolvedShapes <- .vaeResolveShapes(.csh)
  .cov <- .vaeCovariateSearch(d, .ids, .resolvedShapes$rules, .cct,
                              control$covCenter, .cco)
  if (length(.cov$tvExcl) > 0L) {
    ## keep the $runInfo note single-line even with many covariates
    .tvPre <- "time-varying covariate(s) not searched: "
    warning(.tvPre, .vaeTruncList(.cov$tvExcl, prefix = .tvPre), call. = FALSE)
  }
  if (length(.cov$catDropped) > 0L) {
    .cdPre <- "level(s) below catCutoff lumped with reference: "
    warning(.cdPre, .vaeTruncList(.cov$catDropped, prefix = .cdPre), call. = FALSE)
  }
  if (length(.cov$logDrop) > 0L) {
    .ldPre <- "non-positive covariate(s), log shapes skipped: "
    warning(.ldPre, .vaeTruncList(.cov$logDrop, prefix = .ldPre), call. = FALSE)
  }
  if (length(.cov$shapeSub) > 0L) {
    .ssPre <- "shape unwritable at center, plain form used: "
    warning(.ssPre, .vaeTruncList(.cov$shapeSub, prefix = .ssPre), call. = FALSE)
  }
  if (length(.cov$hockeyDrop) > 0L) {
    .hkPre <- "<5% of subjects one side of knot, hockey skipped: "
    warning(.hkPre, .vaeTruncList(.cov$hockeyDrop, prefix = .hkPre),
            call. = FALSE)
  }
  if (length(.cov$naExcl) > 0L) {
    .naPre <- "covariate(s) with missing values not searched: "
    warning(.naPre, .vaeTruncList(.cov$naExcl, prefix = .naPre), call. = FALSE)
  }

  ## pinCovariates=FALSE with a model that declares covariates: turn OFF the
  ## automatic search and estimate every declared covariate in place by the
  ## regress M-step (the covariateSelection=FALSE treatment).  Emptying the
  ## search pool makes the C++ M-step skip covariate selection entirely.  With no
  ## model-declared covariates there is nothing to switch off -- the full search
  ## runs.  (Explicit covariateSelection=FALSE keeps its own path below.)
  .declaredCoefs <- .vaeCovariateCoefThetas(ui)
  .searchOff <- isFALSE(control$pinCovariates) && length(.declaredCoefs) > 0L &&
    !isFALSE(control$covariateSelection)
  if (.searchOff) {
    warning("pinCovariates=FALSE: model covariates estimated in place", call. = FALSE)
    .cov$covNames <- character(0)
    .cov$covMat <- matrix(0, N, 0L)
    .cov$covType <- character(0)
    .cov$covPop <- numeric(0)
    .cov$covRaw <- character(0)
    .cov$covShape <- character(0)
    .cov$covFamily <- character(0)
    .cov$covLevel <- character(0)
    .cov$covGroup <- integer(0)
    .cov$covBlock <- integer(0)
    .cov$covExpr <- character(0)
    .cov$covCanon <- logical(0)
  }

  ## pinned covariate selection: restrict the search to model-declared covariate
  ## /parameter pairs.  Build the per-(eta k x covariate j) allow-mask from the
  ## `inPool` declared pairs; zero those coefficients in the training theta so the
  ## decoder stays covariate-free (the effect is recovered by the M-step prior
  ## regression, exactly as in the unconstrained search) and injected back into
  ## the model afterward.  Declared pairs that cannot be searched are routed to
  ## the regress M-step below (`.pinCovCoef`).
  .pinActive <- FALSE
  .covAllow <- NULL
  .pinPairs <- NULL
  .pinCovCoef <- character(0)
  if (isTRUE(control$pinCovariates) && !isFALSE(control$covariateSelection)) {
    .pinPairs <- .vaeModelCovariatePairs(ui, .cov)
    if (!is.null(.pinPairs) && nrow(.pinPairs) > 0L) {
      .pinActive <- TRUE
      .nCov <- length(.cov$covNames)
      if (.nCov > 0L) {
        ## A covariate column carries ONE encoding.  Claim each column for the
        ## first declared pair's center; if the same covariate is declared again
        ## with a DIFFERENT center (e.g. log(WT/70) on CL and log(WT/80) on KA)
        ## that pair cannot share the column, so demote it to the regress M-step.
        .claim <- rep(NA_real_, .nCov)
        .claimShape <- rep(NA_character_, .nCov)
        for (.r in seq_len(nrow(.pinPairs))) {
          if (!.pinPairs$inPool[.r]) next
          .j <- .vaePinColumn(.cov, .pinPairs[.r, , drop = FALSE])
          if (is.na(.j)) {
            .pinPairs$inPool[.r] <- FALSE
          } else if (is.na(.claim[.j])) {
            .claim[.j] <- .pinPairs$userCenter[.r]
            .claimShape[.j] <- .pinPairs$shape[.r]
          } else if (!isTRUE(all.equal(.claim[.j], .pinPairs$userCenter[.r])) ||
                       !identical(.claimShape[.j], .pinPairs$shape[.r])) {
            ## same column, different center OR different written shape: only one
            ## of them can own the column, so the rest go to the regress M-step
            .pinPairs$inPool[.r] <- FALSE
          }
        }
        .inRows <- .pinPairs[.pinPairs$inPool, , drop = FALSE]
        ## restrict the search to the declared in-pool cells.  An all-zero row
        ## means "no covariate may be selected on this dim" -- crucial when every
        ## declared pair is out-of-pool, so a non-declared (or the out-of-pool)
        ## covariate is never auto-selected under pinning.  A pinned pair allows
        ## ONLY the column matching the shape the user wrote, so pinning never
        ## rewrites the model line it promised to keep verbatim.
        .covAllow <- matrix(0L, .neta, .nCov)
        for (.r in seq_len(nrow(.inRows))) {
          .j <- .vaePinColumn(.cov, .inRows[.r, , drop = FALSE])
          if (!is.na(.j)) .covAllow[.inRows$k[.r], .j] <- 1L
        }
        ## decoder covariate-free during training: hold the declared (in-pool)
        ## coefficient at 0 so the covariate enters only through the prior.
        for (.cn in unique(.inRows$coefName)) {
          .ti <- match(.cn, .thRows$name)
          if (!is.na(.ti)) .th[.ti] <- 0
        }
        ## Retain ONLY the model's own centering (already carried by the mu2/mu3
        ## nlmixrMuDerCov# column, or by the written log(cov/center)); do NOT add
        ## the VAE's mean-centering.  Use each pinned covariate at its MODEL value
        ## so zPop is the model intercept and no post-hoc correction is needed.
        ## Centering a predictor in a regression WITH an intercept leaves the slope
        ## and the selection unchanged -- this only relocates the intercept.
        ## Each claimed column is adjusted EXACTLY once, off the original covPop.
        for (.j in which(!is.na(.claim))) {
          .sh <- .claimShape[.j]
          .rv <- .cov$covRawVal[[.cov$covRaw[.j]]]
          if (!is.na(.sh) && !identical(.sh, "cat") && !is.null(.rv) &&
                is.numeric(.rv)) {
            ## Rebuild the column as the model's OWN expression, so the estimated
            ## slope transfers back with no correction whatever shape was written.
            .cov$covMat[, .j] <- .vaeShapeValue(.sh, .rv, .claim[.j])
          } else {
            ## categorical / mu2-derived: strip the VAE centering, leaving the
            ## column the model's linear term already multiplies
            .cov$covMat[, .j] <- .cov$covMat[, .j] + .cov$covPop[.j]
          }
          .cov$covPop[.j] <- 0
        }
      }
      .pinCovCoef <- unique(.pinPairs$coefName[!.pinPairs$inPool])
      warning("covariate selection pinned to model-specified covariates", call. = FALSE)
      if (length(.pinCovCoef) > 0L) {
        warning("pinned covariate(s) outside search pool estimated in place", call. = FALSE)
      }
    }
  }

  ## `shapes=` may restrict a single (parameter, covariate) pair, but the design
  ## matrix is shared across latent dimensions, so a per-pair restriction has to
  ## be enforced as a mask -- omitting the column would remove that shape from
  ## every parameter.  When nothing is restricted the mask stays NULL and C++
  ## runs the unrestricted search.
  ##
  ## Deliberately NOT applied while pinning is active: `shapes=` governs the
  ## automatic search, whereas a pinned cell is an effect the user wrote in the
  ## model.  Intersecting the two can empty a pinned row -- a declared effect
  ## that is then neither searched nor regressed, and so written back as exactly
  ## 0, silently deleting it.  Under pinning the declaration wins; a declared
  ## shape with no column to pin to already falls back to the regress M-step.
  ## The same holds for fixCov: the declaration is the more specific statement,
  ## so it wins, and the disagreement is reported rather than acted on.
  if (!.searchOff && !.pinActive && length(.cov$covNames) > 0L) {
    .shapeMask <- .vaeShapeAllowMask(.cov, .resolvedShapes, .etaNames,
                                     .foceiEtaThetaMap(ui)$thetaForEta)
    ## Only when a search is actually going to run.  With
    ## covariateSelection=FALSE there is nothing for fixCov to narrow, so the
    ## error below would fire against a user who had ALREADY done what it tells
    ## them to do, and the note would blame fixCov for a search that was off
    ## regardless.
    if (isTRUE(.resolvedShapes$fixCov) && isTRUE(control$covariateSelection)) {
      ## fixCov=TRUE with nothing left to search is a contradiction the user
      ## almost certainly did not intend; covariateSelection=FALSE is the way to
      ## ask for no search at all.
      if (all(.shapeMask == 0L)) {
        stop("shapes: fixCov=TRUE leaves no covariate searchable on any parameter\n",
             "  use covariateSelection=FALSE to turn the search off outright",
             call. = FALSE)
      }
      ## Grouped by RAW covariate, not by column.  A covariate restricted to one
      ## shape has its other columns masked to zero while the covariate itself is
      ## still searched through the column that survived, so a per-column test
      ## could name it as excluded when it was not.
      .fxBy <- tapply(colSums(.shapeMask), .cov$covRaw, sum)
      .fxDrop <- names(.fxBy)[.fxBy == 0]
      if (length(.fxDrop) > 0L) {
        .fxPre <- "fixCov=TRUE, covariate(s) not searched: "
        warning(.fxPre, .vaeTruncList(.fxDrop, prefix = .fxPre),
                call. = FALSE)
      }
    }
    if (any(.shapeMask == 0L)) .covAllow <- .shapeMask
  } else if (isTRUE(.resolvedShapes$fixCov) && .pinActive) {
    warning("fixCov=TRUE ignored: the model declares covariates, which pins the search",
            call. = FALSE)
  }

  ## Fixed-effect thetas estimated directly by a bounded bobyqa regression in the
  ## M-step (vs. the latent-space zPop update).  Two sources, unioned:
  ##  * nonMuTheta="regress"/"grad": non-mu-referenced structural thetas (no eta); and
  ##  * covariateSelection=FALSE: model-declared covariate coefficients -- always
  ##    estimated in place, independent of nonMuTheta, so the mu-referenced
  ##    covariate expression the user wrote is fit rather than held fixed.
  ## Carry each 0-based index into the full theta vector (`.th`, ntheta order)
  ## plus the ini() bounds (NA -> +-Inf).
  .regressNames <- character(0)
  .regressThetaIdx0 <- integer(0)
  .regressLower <- numeric(0); .regressUpper <- numeric(0)
  if (.vaeNonMuIsRegress(control$nonMuTheta)) {
    .regressNames <- .vaeNonMuThetas(ui)
  }
  .covCoefNames <- character(0)
  if (isFALSE(control$covariateSelection)) {
    .covCoefNames <- .declaredCoefs
  } else if (.searchOff) {
    ## pinCovariates=FALSE with model-declared covariates: all of them regress.
    .covCoefNames <- .declaredCoefs
  } else if (.pinActive && length(.pinCovCoef) > 0L) {
    ## pinned selection: a declared covariate that cannot be handled by the
    ## restricted search (out-of-pool, or a form whose slope will not transfer)
    ## is estimated in place by the regress M-step, like covariateSelection=FALSE.
    .covCoefNames <- .pinCovCoef
  }
  .regressNames <- c(.regressNames, .covCoefNames)
  .regressNames <- unique(.regressNames)
  ## residOptimize="optimize": the residual-error thetas join the SAME optimizer
  ## as the non-mu structural thetas, against the same full outer objective.  A
  ## FIXED error parameter is excluded (nothing to estimate), as is one already
  ## regressed for another reason.
  .errRegressNames <- character(0)
  if (!identical(control$residOptimize, "moment")) {
    .errAll <- .idf[!is.na(.idf$err) & !is.na(.idf$ntheta), , drop = FALSE]
    if (nrow(.errAll) > 0L) {
      .errFree <- .errAll[!(!is.na(.errAll$fix) & .errAll$fix), , drop = FALSE]
      ## Every free residual parameter enters the optimizer.  The stage-2
      ## objective evaluates the ordinary likelihood with the ODE frozen, so `r`
      ## comes from the model's own rx_r_ and ANY error form is scored correctly
      ## -- there is nothing left for a form-specific filter to protect against.
      ## (It was needed only while the objective recomputed `r` itself from a
      ## hardcoded per-form expression, which could silently ignore a parameter
      ## and let the optimizer move it on noise.)
      .errRegressNames <- setdiff(as.character(.errFree$name), .regressNames)
    }
  }
  .regressNames <- c(.regressNames, .errRegressNames)
  if (length(.regressNames) > 0L) {
    .ri <- match(.regressNames, .thRows$name)
    .regressThetaIdx0 <- as.integer(.ri - 1L)
    .lo <- as.numeric(.thRows$lower[.ri]); .hi <- as.numeric(.thRows$upper[.ri])
    .regressLower <- ifelse(is.na(.lo), -Inf, .lo)
    .regressUpper <- ifelse(is.na(.hi), Inf, .hi)
    ## A residual SCALE parameter must not be allowed to reach zero.  The
    ## likelihood floors a zero variance (r == 0 -> r = 1) to stay finite, which
    ## makes a collapsed residual look attractive rather than forbidden -- a
    ## boxCox fit converged to add.err = 0 exactly this way.  Floor the scale
    ## parameters strictly above zero: absolute forms (add, lnorm) relative to
    ## the spread of the data, relative forms (prop, pow) at a small constant.
    ## Exponents and lambdas are NOT scales and are left alone.
    .errScaleAbs <- as.character(.idf$name[!is.na(.idf$err) & .idf$err %in% c("add", "lnorm")])
    .errScaleRel <- as.character(.idf$name[!is.na(.idf$err) & .idf$err %in% c("prop", "pow")])
    .dvObs <- suppressWarnings(as.numeric(d$DV[d$EVID == 0]))
    .dvSpread <- stats::sd(.dvObs[is.finite(.dvObs)])
    if (!is.finite(.dvSpread) || .dvSpread <= 0) .dvSpread <- 1
    .isAbs <- .regressNames %in% .errScaleAbs
    .isRel <- .regressNames %in% .errScaleRel
    .regressLower[.isAbs] <- pmax(.regressLower[.isAbs], 1e-4 * .dvSpread)
    .regressLower[.isRel] <- pmax(.regressLower[.isRel], 1e-6)
    ## A transform-both-sides lambda (Box-Cox / Yeo-Johnson) is only meaningful
    ## on a narrow interval, and it is unbounded in the ini() block, so the
    ## optimizer would otherwise search a meaningless range.  Constrain it to
    ## (-2, 2); SAEM does the same thing by mapping lambda through a bounded
    ## transform (`toLambda`).  A tighter user bound still wins.
    .lamNames <- as.character(.idf$name[!is.na(.idf$err) &
                                        .idf$err %in% c("boxCox", "yeoJohnson")])
    if (length(.lamNames) > 0L) {
      .isLam <- .regressNames %in% .lamNames
      .regressLower[.isLam] <- pmax(.regressLower[.isLam], -2)
      .regressUpper[.isLam] <- pmin(.regressUpper[.isLam], 2)
    }
    ## An UNBOUNDED covariate coefficient regressed alone routes through the 1-D
    ## optimize() branch of .boundedResidOpt, which searches the whole interval and
    ## overshoots a shallow interior optimum on a too-wide interval.  Give an
    ## unbounded coefficient a finite, scale-aware fallback interval (user ini()
    ## bounds still win).
    .isCov <- .regressNames %in% .covCoefNames
    .noLo <- .isCov & !is.finite(.regressLower)
    .noHi <- .isCov & !is.finite(.regressUpper)
    if (any(.noLo | .noHi)) {
      .bnd <- .vaeCovCoefBoundVec(ui, data, .regressNames[.isCov])
      .regressLower[.noLo] <- -.bnd[.regressNames[.noLo]]
      .regressUpper[.noHi] <- .bnd[.regressNames[.noHi]]
    }
    ## A STRUCTURAL non-mu theta needs the same guard.  With +-Inf bounds nothing
    ## constrains the M-step (bobyqa's interval, or the "grad" Adam projection), and
    ## a theta whose likelihood is flat in one direction runs away: an unbounded
    ## `tv <- 3.45` on theo_sd reaches ~1e68 (an lnorm fit there reports an OFV of
    ## 359315 against focei's 686), while the same model with `tv <- c(2, 3.45, 5)`
    ## converges.  Fall back to a generous window around the ini() ESTIMATE, wide
    ## enough not to bind at a sane optimum but finite so the search cannot diverge.
    .isStruct <- !.isCov
    .sLo <- .isStruct & !is.finite(.regressLower)
    .sHi <- .isStruct & !is.finite(.regressUpper)
    if (any(.sLo | .sHi)) {
      .init <- as.numeric(.thRows$est[.ri])
      .init[!is.finite(.init)] <- 0
      ## scale-aware half-width: the absolute floor covers a log-scale parameter
      ## (init ~ 3.45 -> +-10 is exp(+-10), ample), the relative term keeps a
      ## large-magnitude natural-scale init (say 1000) from being over-constrained
      .hw <- pmax(.vaeNonMuThetaBound, abs(.init) * .vaeNonMuThetaRel)
      .regressLower[.sLo] <- (.init - .hw)[.sLo]
      .regressUpper[.sHi] <- (.init + .hw)[.sHi]
    }
  }

  ## residual error params (all of them, in theta order): value, theta index,
  ## type (add/prop/...), and bounds. Combined models have >1 row; log-likelihood
  ## models may have none. `a` is the (named) error-param vector.
  .errRow <- .idf[!is.na(.idf$err) & !is.na(.idf$ntheta), , drop = FALSE]
  .errRow <- .errRow[order(.errRow$ntheta), , drop = FALSE]
  .a <- if (nrow(.errRow) > 0) setNames(as.numeric(.errRow$est), .errRow$name) else numeric(0)
  ## For each regressed parameter, its 0-based slot in `a` (the error-parameter
  ## vector), or -1 when it is not an error parameter.  vaeBuildTh writes `a`
  ## over the error theta positions, so the objective must substitute the
  ## CANDIDATE value into `a`; without this map it would be flat in every error
  ## parameter and the optimizer would never move one.
  .regressErrIdx0 <- if (length(.regressNames) > 0L) {
    .m <- match(.regressNames, names(.a))
    as.integer(ifelse(is.na(.m), 0L, .m) - 1L)
  } else integer(0)
  .errThetaIdx <- as.integer(.errRow$ntheta)
  .errType <- as.character(.errRow$err)
  .errLower <- as.numeric(.errRow$lower); .errUpper <- as.numeric(.errRow$upper)
  ## residOptimize="twoStage" stage-2 eligibility (see .vaeRegressStage2)
  .regressStage2 <- .vaeRegressStage2(ui, .regressNames, .regressErrIdx0)

  ## per-subject decoder inputs + gather all obs for standardization
  subj <- vector("list", N)
  .allTime <- numeric(0); .allDv <- numeric(0)
  for (i in seq_len(N)) {
    .di <- d[d$ID == .ids[i], , drop = FALSE]
    .obs <- .di[.di$EVID == 0, , drop = FALSE]
    .times <- .obs$TIME
    .y <- .obs$DV
    ## M2/M3/M4 censoring columns (0 / NA when absent)
    .cens <- if (is.null(.obs$CENS)) integer(length(.y)) else as.integer(.obs$CENS)
    .limit <- if (is.null(.obs$LIMIT)) rep(NA_real_, length(.y)) else as.numeric(.obs$LIMIT)
    subj[[i]] <- list(ev = .di, times = .times, y = .y, n = length(.times),
                      cens = .cens, limit = .limit)
    .allTime <- c(.allTime, .times); .allDv <- c(.allDv, .y)
  }
  Tmax <- max(vapply(subj, function(s) s$n, integer(1)))
  .tMax <- max(.allTime)
  ## inputScale: which DV values the encoder-input centering/scaling is computed
  ## over.  The reference takes mean/sd across the WHOLE padded [N, Tmax] matrix,
  ## so the zero padding of short subjects enters both -- on a ragged dataset
  ## that is a materially different scale from the observed-only statistics
  ## (neonatal: sd 1582 vs 506).  "reference" reproduces it; "observed" uses the
  ## observed values only.
  if (identical(.inputScale, "reference")) {
    .padded <- c(.allDv, rep(0, N * Tmax - length(.allDv)))
    .dvMean <- mean(.padded); .dvSd <- stats::sd(.padded)
  } else {
    .dvMean <- mean(.allDv); .dvSd <- stats::sd(.allDv)
  }
  if (!is.finite(.dvSd) || .dvSd <= 0) .dvSd <- 1

  ## encoder inputs: [N, Tmax, 2] standardized (time, DV), padded; lengths
  dataIn <- array(0, c(N, Tmax, 2L))
  lengths <- integer(N)
  for (i in seq_len(N)) {
    s <- subj[[i]]; ni <- s$n; lengths[i] <- ni
    dataIn[i, seq_len(ni), 1L] <- s$times / .tMax
    dataIn[i, seq_len(ni), 2L] <- (s$y - .dvMean) / .dvSd
  }
  ## Encoder-head covariates.  The reference concatenates them to the LSTM's
  ## FINAL HIDDEN STATE before the linear head that emits (mu, logSigma, L)
  ## -- `torch.cat((hidden[-1], covariates), dim=1)` in its encoder -- so the
  ## approximate posterior q(z|x) is conditioned on the covariates.  That is
  ## central to the method: it is how the encoder can express a covariate
  ## relationship at all, and the M-step then only has to read it off the
  ## posterior means.  Feeding zero columns here leaves the posterior
  ## unconditioned, which shows up as less between-subject spread aligned with
  ## the covariates (smaller omega, more variance pushed into residual error) and
  ## weaker covariate effects.  Use the same encoded matrix the selection step
  ## uses (continuous -> log(cov/mean), categorical -> linear).
  ## Only the covariates actually in play are supplied.  With no covariate
  ## search there is nothing for the encoder to condition on (the declared
  ## coefficients are estimated in place by the regress M-step instead), and
  ## under `pinCovariates` only the pinned candidates are; conditioning on a
  ## covariate the search cannot select would let the posterior encode a
  ## relationship the model never reports.
  ## One column per exclusion GROUP: alternate shapes of a covariate are
  ## near-collinear copies, so feeding all of them would widen the encoder head
  ## for no information.  Distinct factor levels are separate groups and all go in.
  ## Dedupe among the SELECTABLE columns, not against a fixed canonical column:
  ## the allowed column of a group need not be the group's first one (a pinned
  ## log(cov/center) pair takes the log column even when `shapes=` lists a linear
  ## shape first), and intersecting the two would drop the covariate entirely.
  covIn <- if (isFALSE(control$covariateSelection) || .searchOff ||
                 ncol(.cov$covMat) == 0L) {
    matrix(0, N, 0L)
  } else {
    .sel <- if (is.null(.covAllow)) seq_len(ncol(.cov$covMat)) else
      which(colSums(.covAllow) > 0L)
    ## Dedupe to one BLOCK per group, not one column: the arms of a hockey block
    ## are complementary halves of one relationship, so keeping only the first
    ## would hand the encoder a covariate truncated at the knot.  A block never
    ## spans groups, and with every column its own block this is exactly the
    ## historic one-column-per-group rule.
    .g <- .cov$covGroup[.sel]
    .b <- .cov$covBlock[.sel]
    .keep <- .sel[.b %in% .b[!duplicated(.g)]]
    .cov$covMat[, .keep, drop = FALSE]
  }
  if (!is.matrix(covIn) || nrow(covIn) != N) covIn <- matrix(0, N, 0L)

  list(N = N, neta = .neta, zDim = .neta, etaNames = .etaNames,
       th = .th, zPopThetaIdx = .zPopThetaIdx, isFree = .isFree, omegaFix = .omegaFix,
       zPopFix = .zPopFix,
       zPopLower = .zPopLower, zPopUpper = .zPopUpper,
       errThetaIdx = .errThetaIdx, errType = .errType,
       errLower = .errLower, errUpper = .errUpper,
       regressNames = .regressNames, regressThetaIdx0 = .regressThetaIdx0,
       regressErrIdx0 = .regressErrIdx0, regressStage2 = .regressStage2,
       regressLower = .regressLower, regressUpper = .regressUpper,
       zPop = .zPop, omega = .omega, a = .a,
       omegaMat = .omBlock$mat, omegaFixMat = .omBlock$fixMat,
       subj = subj, dataIn = dataIn, lengths = lengths, covIn = covIn,
       covNames = .cov$covNames, covMat = .cov$covMat, covType = .cov$covType,
       covPop = .cov$covPop, covRaw = .cov$covRaw, covShape = .cov$covShape,
       covFamily = .cov$covFamily, covLevel = .cov$covLevel,
       covGroup = .cov$covGroup, covBlock = .cov$covBlock,
       covExpr = .cov$covExpr,
       covCanon = .cov$covCanon, shapeRules = .resolvedShapes$rules,
       pinActive = .pinActive, pinPairs = .pinPairs, covAllow = .covAllow,
       tMax = .tMax, dvMean = .dvMean, dvSd = .dvSd, Nobs = length(.allDv))
}
