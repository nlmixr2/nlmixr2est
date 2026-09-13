# Stepwise covariate search over declared-distribution ARGUMENT ROLES.
#
# The unit of search is a (declared eta, family argument ROLE) pair -- gamma's
# `shape` against its `rate` -- because that is what a covariate can enter and
# what a coefficient's sign means.  `.etaDistAddCovariate()` is the applier;
# everything here decides WHICH relation to hand it and whether to keep it.
#
# Deliberately an OUTER search over refits rather than an in-loop selection.
# vae's L0/branch-and-bound selects on `zPop` and explicitly skips free etas
# (`src/inner.cpp`), which is every declared latent, so there is nothing to
# extend there -- see the plan's phase 3.4.

#' Candidate (eta, role) pairs a covariate may enter
#'
#' A role that names a SUPPORT endpoint is refused (`.etaDistRoleNoCovariate`):
#' a subject-varying endpoint makes the density discontinuous in the parameter,
#' so there is nothing for the M-step to follow.
#'
#' A family whose roles are unknown -- an older lotri, or a family missing from
#' the table -- yields NO candidates rather than every argument grouped
#' together.  Guessing the grouping would attach a covariate to whichever
#' argument happened to come first, which is worse than declining.
#'
#' @param ui rxode2 ui carrying `dist()` declarations
#' @return data.frame of eta/role/argument-position, empty when there is nothing
#'   to search
#' @noRd
.etaDistCovTargets <- function(ui) {
  .empty <- data.frame(eta = character(0), role = character(0),
                       pos = integer(0), stringsAsFactors = FALSE)
  .d <- tryCatch(rxode2::rxUiEtaDists(rxode2::rxUiDecompress(ui)),
                 error = function(e) NULL)
  if (is.null(.d) || NROW(.d) == 0L) return(.empty)
  .rows <- lapply(seq_len(NROW(.d)), function(.i) {
    .g <- .etaDistRoleGroups(.d$etaDist[.i])
    if (is.null(.g) || length(.g) == 0L) return(NULL)
    .g <- .g[!(names(.g) %in% .etaDistRoleNoCovariate)]
    if (length(.g) == 0L) return(NULL)
    data.frame(eta = .d$name[.i], role = names(.g),
               pos = vapply(.g, function(.p) as.integer(.p[1]), integer(1)),
               stringsAsFactors = FALSE)
  })
  .rows <- .rows[!vapply(.rows, is.null, logical(1))]
  if (length(.rows) == 0L) return(.empty)
  .out <- do.call(rbind, .rows)
  rownames(.out) <- NULL
  .out
}

#' Every relation the search may add, before any is fitted
#'
#' The candidate set is `target x covariate x shape`, with one shape per
#' covariate per target: alternate shapes of one covariate are mutually
#' exclusive, which is `covGroup`'s existing meaning applied within a role
#' rather than a second rule (`.vaeCovariateSearch`).
#'
#' @param ui rxode2 ui
#' @param cov the `.vaeCovariateSearch()` result for this data
#' @param already data.frame of relations already in the model, excluded here so
#'   a forward step never re-proposes one
#' @return data.frame with eta/role/cov/shape/center, one row per relation
#' @noRd
.etaDistCovCandidates <- function(ui, cov, already = NULL) {
  .t <- .etaDistCovTargets(ui)
  .empty <- data.frame(eta = character(0), role = character(0),
                       cov = character(0), shape = character(0),
                       center = numeric(0), stringsAsFactors = FALSE)
  if (NROW(.t) == 0L || is.null(cov) || length(cov$covNames) == 0L) return(.empty)
  ## covCanon keeps ONE column per covGroup: the alternate shapes of a covariate
  ## are the group, and proposing all of them at once would let two shapes of
  ## the same covariate be selected together.
  .keep <- if (!is.null(cov$covCanon)) which(cov$covCanon) else seq_along(cov$covNames)
  if (length(.keep) == 0L) return(.empty)
  .g <- expand.grid(ti = seq_len(NROW(.t)), ci = .keep,
                    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  .out <- data.frame(
    eta = .t$eta[.g$ti], role = .t$role[.g$ti],
    cov = if (is.null(cov$covRaw)) cov$covNames[.g$ci] else cov$covRaw[.g$ci],
    shape = if (is.null(cov$covShape)) rep("power", NROW(.g)) else cov$covShape[.g$ci],
    center = if (is.null(cov$covPop)) rep(NA_real_, NROW(.g)) else cov$covPop[.g$ci],
    stringsAsFactors = FALSE)
  if (!is.null(already) && NROW(already) > 0L) {
    .key <- function(d) paste(d$eta, d$role, d$cov, d$shape, sep = "\r")
    .out <- .out[!(.key(.out) %in% .key(already)), , drop = FALSE]
  }
  rownames(.out) <- NULL
  .out
}

#' Build a ui with a whole SET of relations applied to the base
#'
#' Every model the search scores is built this way -- base plus a set -- rather
#' than by editing the previous step's ui.  Adding is then exactly the inverse
#' of dropping, so forward and backward compare like with like, and no
#' incremental edit can drift.  It also means backward elimination needs no
#' "remove a relation" operation, which would have to guess what the applier
#' did to the model text.
#'
#' `NULL` when any relation cannot be applied -- callers treat that as a model
#' that cannot be scored, never as one that scored badly.
#' @noRd
.etaDistApplyAll <- function(baseUi, rels) {
  .u <- baseUi
  if (is.null(rels) || NROW(rels) == 0L) return(.u)
  for (.j in seq_len(NROW(rels))) {
    .u <- tryCatch(.etaDistAddCovariate(.u, rels$eta[.j], rels$role[.j],
                                        rels$cov[.j], shape = rels$shape[.j],
                                        center = rels$center[.j]),
                   error = function(e) NULL)
    if (is.null(.u)) return(NULL)
  }
  .u
}

#' Fit one relation SET and report its objective
#'
#' A model that cannot be built or cannot be fitted is one that could not be
#' SCORED, which is reported as such.  Returning `Inf` instead would let a
#' broken candidate lose quietly and look considered.
#' @noRd
.etaDistCovFit <- function(baseUi, data, est, control, rels) {
  .u <- .etaDistApplyAll(baseUi, rels)
  if (is.null(.u)) return(list(ok = FALSE, why = "relation could not be applied"))
  .f <- tryCatch(suppressMessages(suppressWarnings(
    nlmixr2(.u, data, est = est, control = control))), error = function(e) e)
  if (inherits(.f, "error")) {
    return(list(ok = FALSE, why = paste0("fit failed: ", conditionMessage(.f))))
  }
  .o <- tryCatch(as.numeric(.f$objf), error = function(e) NA_real_)
  if (!is.finite(.o)) return(list(ok = FALSE, why = "fit returned no objective"))
  list(ok = TRUE, objf = .o, ui = .u, fit = .f)
}

.etaDistPathRow <- function(step, rels, res, ref, note = NA_character_) {
  data.frame(step = step, rels,
             objf = if (isTRUE(res$ok)) res$objf else NA_real_,
             dOfv = if (isTRUE(res$ok)) res$objf - ref else NA_real_,
             why = if (isTRUE(res$ok)) note else res$why,
             stringsAsFactors = FALSE)
}

#' Stepwise covariate search on declared-distribution arguments
#'
#' Forward selection then backward elimination over (declared eta, family
#' argument role, covariate, shape), scored on the change in objective.
#'
#' The search is over REFITS. A declared random effect enters as
#' `eta = Q(phiU(z); args(theta))` and has no `theta + eta` mu reference, which
#' is what `nlmixr2extra::covarSearchAuto()` and `nlmixr2scm::runSCM()` rewrite
#' -- their stepwise drivers are the reusable idea, their relation appliers are
#' not. Candidate columns, shapes and mutual exclusion come from
#' `.vaeCovariateSearch()` rather than a second discovery pass.
#'
#' Subject-level covariates only, which is `.vaeCovariateSearch()`'s own rule
#' (it reports what it dropped in `tvExcl`). That is not merely inherited: a
#' covariate varying WITHIN a subject currently freezes its declaration's
#' thetas, so searching over one would be comparing fits that never moved.
#'
#' @param object model function or ui carrying `dist()` declarations
#' @param data data to fit
#' @param est estimation method. `"focei"` is the default because it is the one
#'   measured to recover a coefficient on a declaration.
#' @param control control for `est`
#' @param dOfvAdd objective drop a relation must produce to be ADDED
#' @param dOfvDrop objective rise a relation must cost to be KEPT in backward
#'   elimination; conventionally stricter than `dOfvAdd`
#' @param maxSteps maximum forward steps
#' @return list of `ui` (the selected model), `relations` kept, final `objf`,
#'   and `path`: every candidate scored, including the ones that could not be
#'   scored and why
#' @export
etaDistCovarSearch <- function(object, data, est = "focei", control = NULL,
                               dOfvAdd = 3.84, dOfvDrop = 6.63, maxSteps = 10L) {
  .base <- rxode2::rxUiDecompress(rxode2::assertRxUi(object))
  if (NROW(.etaDistCovTargets(.base)) == 0L) {
    stop("etaDistCovarSearch: no dist() declaration with a covariate-eligible ",
         "argument role -- nothing to search", call. = FALSE)
  }
  .idCol <- grep("^id$", names(data), ignore.case = TRUE)
  if (length(.idCol) == 0L) {
    stop("etaDistCovarSearch: data has no ID column", call. = FALSE)
  }
  .cov <- .vaeCovariateSearch(data, data[[.idCol[1]]])
  if (length(.cov$tvExcl) > 0L) {
    .minfo(paste0("etaDistCovarSearch: not searching time-varying covariate(s): ",
                  paste(.cov$tvExcl, collapse = ", ")))
  }
  .r0 <- .etaDistCovFit(.base, data, est, control, NULL)
  if (!isTRUE(.r0$ok)) {
    stop("etaDistCovarSearch: the base model could not be scored -- ", .r0$why,
         call. = FALSE)
  }
  .objf <- .r0$objf
  .kept <- NULL
  .path <- NULL
  for (.step in seq_len(maxSteps)) {
    .cand <- .etaDistCovCandidates(.base, .cov, .kept)
    if (NROW(.cand) == 0L) break
    .res <- lapply(seq_len(NROW(.cand)), function(.i) {
      .etaDistCovFit(.base, data, est, control,
                     rbind(.kept, .cand[.i, , drop = FALSE]))
    })
    for (.i in seq_along(.res)) {
      .path <- rbind(.path, .etaDistPathRow(.step, .cand[.i, , drop = FALSE],
                                            .res[[.i]], .objf))
    }
    .score <- vapply(.res, function(.x) if (isTRUE(.x$ok)) .x$objf else Inf, numeric(1))
    .best <- which.min(.score)
    if (!is.finite(.score[.best]) || .score[.best] > .objf - dOfvAdd) break
    .objf <- .score[.best]
    .kept <- rbind(.kept, .cand[.best, , drop = FALSE])
  }
  if (!is.null(.kept) && NROW(.kept) > 1L) {
    .bk <- .etaDistCovBackward(.base, data, est, control, .kept, .objf, dOfvDrop)
    .kept <- .bk$kept; .objf <- .bk$objf; .path <- rbind(.path, .bk$path)
  }
  list(ui = .etaDistApplyAll(.base, .kept), relations = .kept,
       objf = .objf, path = .path)
}

#' Backward elimination over the relations forward selection kept
#'
#' Each trial model is the BASE plus the relations that survive, so dropping is
#' the exact inverse of adding (`.etaDistApplyAll`) rather than an edit that
#' would have to undo what the applier wrote.
#' @noRd
.etaDistCovBackward <- function(baseUi, data, est, control, kept, objf, dOfvDrop) {
  .path <- NULL
  repeat {
    if (NROW(kept) <= 1L) break
    .res <- lapply(seq_len(NROW(kept)), function(.i) {
      .etaDistCovFit(baseUi, data, est, control, kept[-.i, , drop = FALSE])
    })
    .score <- vapply(.res, function(.x) if (isTRUE(.x$ok)) .x$objf else Inf, numeric(1))
    .cheapest <- which.min(.score)
    if (!is.finite(.score[.cheapest])) break
    ## a relation is kept only if dropping it costs at least dOfvDrop
    if (.score[.cheapest] - objf >= dOfvDrop) break
    .path <- rbind(.path, .etaDistPathRow(-1L, kept[.cheapest, , drop = FALSE],
                                          .res[[.cheapest]], objf,
                                          note = "dropped in backward elimination"))
    kept <- kept[-.cheapest, , drop = FALSE]
    objf <- .score[.cheapest]
  }
  list(kept = kept, objf = objf, path = .path)
}
