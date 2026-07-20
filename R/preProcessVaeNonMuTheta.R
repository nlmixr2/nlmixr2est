# VAE-only preprocessing: est="vae" only estimates parameters that occupy the
# latent space (structural thetas carrying a random effect, plus omega and
# residual error).  A structural theta with no eta stays frozen at its ini()
# value.  This hook injects a small eta into each such non-mu-referenced theta so
# it enters the latent space and is estimated; the temporary eta is dropped from
# the reported model (collapsed to `theta + mean(eta)`) in the VAE output path.
# Reuses the rxode2 model()/ini() piping to add the eta and rmEta()/.downgradeEtas
# to drop it -- no new UI-surgery machinery.

#' Collapse a character vector to a comma-separated list that (together with
#' `prefix`) fits `budget` characters, replacing the overflow tail with
#' "+k more" so a runInfo note stays single-line and bounded.
#' @param x character vector
#' @param prefix leading text the list is appended to (counted against budget)
#' @param budget maximum total character width
#' @return length-1 character
#' @noRd
.vaeTruncList <- function(x, prefix = "", budget = 74L) {
  .n <- length(x)
  .avail <- budget - nchar(prefix)
  .keep <- .n
  repeat {
    .more <- .n - .keep
    .txt <- paste(x[seq_len(.keep)], collapse = ", ")
    if (.more > 0L) .txt <- paste0(.txt, ", +", .more, " more")
    if (.keep <= 1L || nchar(.txt) <= .avail) break
    .keep <- .keep - 1L
  }
  .txt
}

#' Covariate-coefficient thetas of a VAE model (the parameters that multiply a
#' data covariate in a mu-referenced expression).  Read-only over the SHARED
#' `muRefCovariateDataFrame`/`allCovs` UI fields -- the same covariate identity
#' SAEM (saemMuRefCovariateDataFrame) and FOCEI (muRefClassify groups) consume --
#' so nothing here mutates the shared covariate representation.
#'
#' Two categories are returned:
#'  (a) linear effects rxode2 records in `muRefCovariateDataFrame` (`beta*WT`);
#'  (b) transformed effects it does NOT (`beta*log(WT/70)`) -- detected as a
#'      non-mu, non-error theta whose EVERY referencing model line also mentions
#'      a data covariate (a plain structural theta appears in at least one
#'      covariate-free line, so it is not mis-caught).
#' User-fixed coefficients (`fix=TRUE`) are dropped.
#' @param ui rxode2 ui
#' @return character vector of theta names (possibly empty)
#' @noRd
.vaeCovariateCoefThetas <- function(ui) {
  .idf <- ui$iniDf
  .th <- .idf[!is.na(.idf$ntheta) & is.na(.idf$err) & !isTRUE2(.idf$fix), , drop = FALSE]
  if (nrow(.th) == 0L) return(character(0))
  .thNames <- .th$name
  .mu <- if (is.null(ui$muRefDataFrame)) character(0) else ui$muRefDataFrame$theta
  ## (a) rxode2-recognized linear mu-ref covariate coefficients
  .linear <- if (is.null(ui$muRefCovariateDataFrame)) character(0)
             else as.character(ui$muRefCovariateDataFrame$covariateParameter)
  ## (b) transformed covariate effects: a non-mu theta whose referencing lines
  ## ALL reference a data covariate
  .covData <- ui$allCovs
  .transformed <- character(0)
  if (length(.covData) > 0L) {
    for (.p in setdiff(.thNames, .mu)) {
      .lines <- Filter(function(e) .p %in% all.vars(e), ui$lstExpr)
      if (length(.lines) == 0L) next
      if (all(vapply(.lines, function(e) any(.covData %in% all.vars(e)), logical(1)))) {
        .transformed <- c(.transformed, .p)
      }
    }
  }
  intersect(unique(c(.linear, .transformed)), .thNames)
}

#' Structural population thetas that are NOT mu-referenced (candidates for a
#' VAE-injected eta): a theta that appears in the model, is not a residual-error
#' or covariate-coefficient theta, is not fixed, and has no eta referencing it.
#' @param ui rxode2 ui
#' @return character vector of theta names (possibly empty)
#' @noRd
.vaeNonMuThetas <- function(ui) {
  .idf <- ui$iniDf
  .th <- .idf[!is.na(.idf$ntheta) & is.na(.idf$err) & !isTRUE2(.idf$fix), , drop = FALSE]
  if (nrow(.th) == 0L) return(character(0))
  .mu <- if (is.null(ui$muRefDataFrame)) character(0) else ui$muRefDataFrame$theta
  .cov <- if (is.null(ui$muRefCovariateDataFrame)) character(0) else ui$muRefCovariateDataFrame$theta
  ## covariate coefficients are estimated by the regress M-step (see .vaeDataPrep),
  ## not by nonMuTheta eta/fix injection -- exclude them here
  .covCoef <- .vaeCovariateCoefThetas(ui)
  .cand <- setdiff(.th$name, c(.mu, .cov, .covCoef))
  if (length(.cand) == 0L) return(character(0))
  ## keep only thetas that actually appear in a model expression (so an eta can be
  ## attached to a structural line)
  .modelVars <- unique(unlist(lapply(ui$lstExpr, all.vars)))
  .cand[.cand %in% .modelVars]
}

#' Vectorized isTRUE for a logical/NA column
#' @noRd
isTRUE2 <- function(x) !is.na(x) & x

#' Inject an eta into each of the given non-mu-referenced thetas.
#'
#' For theta `p` with structural line `... p ...`, rewrites the line to
#' `... (p + eta.p) ...` and adds `eta.p ~ nonMuEtaOmega` (fixed when `fix=TRUE`),
#' reusing rxode2 model()/ini() piping.
#' @param ui rxode2 ui
#' @param thetas character theta names to convert
#' @param omega numeric injected variance
#' @param fix logical, hold the injected omega fixed
#' @return list(ui=, etas=) with the modified ui and the injected eta names
#' @noRd
.vaeInjectNonMuEtas <- function(ui, thetas, omega, fix) {
  .ui <- ui
  .etas <- character(0)
  for (.p in thetas) {
    .eta <- .vaeUniqueEtaName(.ui, .p)
    .lines <- .ui$lstExpr
    .idx <- which(vapply(.lines, function(e) .p %in% all.vars(e), logical(1)))
    if (length(.idx) == 0L) next
    ## rewrite EVERY expression that references the theta, not just the first: the
    ## injected eta is a single per-subject random effect, so replacing each `p`
    ## with `p + eta` uses the same eta realization everywhere and consistently
    ## moves the structural parameter into the latent space (a first-only rewrite
    ## would leave the other uses frozen at the population value).
    ## inject FLAT (no wrapping parens) so `exp(theta + eta)` stays the additive
    ## form rxode2 recognizes as a mu-referenced exp() parameter (parens would hide
    ## the exp() back-transform -- see .vaeUpdateModel)
    for (.i in .idx) {
      .newTxt <- gsub(paste0("\\b", .p, "\\b"), paste0(.p, " + ", .eta),
                      deparse1(.lines[[.i]]))
      .ui <- do.call(rxode2::model, list(.ui, str2lang(.newTxt)))
    }
    .ui <- do.call(rxode2::ini, list(.ui, str2lang(paste0(.eta, " ~ ", signif(omega, 12)))))
    .etas <- c(.etas, .eta)
  }
  if (fix && length(.etas) > 0L) {
    ## nonMuTheta="fix" holds the whole parameter fixed: both the injected omega
    ## AND the structural theta (the fixed effect) stay at their ini() value.
    ## rxode2 ini(eta ~ fix(v)) does not flag the omega, so set the diagonal fix
    ## flag directly (read by .vaeDataPrep's omegaFix); also fix the paired theta
    ## rows (read by .vaeDataPrep's zPopFix -> the M-step holds the typical value
    ## at ini and drops it from the iteration print).  Rebuild via the compress
    ## round-trip (as in advi.R) so derived UI fields re-sync.
    .ui <- rxode2::rxUiDecompress(.ui)
    .df <- .ui$iniDf
    .df$fix[!is.na(.df$neta1) & .df$neta1 == .df$neta2 & .df$name %in% .etas] <- TRUE
    .df$fix[!is.na(.df$ntheta) & .df$name %in% thetas] <- TRUE
    assign("iniDf", .df, envir = .ui)
    .ui <- rxode2::rxUiCompress(.ui)
  }
  list(ui = .ui, etas = .etas)
}

#' A fresh eta name for theta `p` (strip a leading structural-theta "t" when that
#' yields an unused name, e.g. tR0 -> eta.R0; else eta.<theta>)
#' @noRd
.vaeUniqueEtaName <- function(ui, p) {
  .used <- ui$iniDf$name
  .base <- if (grepl("^t[A-Z0-9]", p)) sub("^t", "", p) else p
  .cand <- paste0("eta.", .base)
  if (!(.cand %in% .used)) return(.cand)
  .cand <- paste0("eta.", p)
  .i <- 1L
  while (.cand %in% .used) { .cand <- paste0("eta.", p, ".", .i); .i <- .i + 1L }
  .cand
}

#' VAE preprocessing hook: inject etas for non-mu-referenced thetas per
#' vaeControl(nonMuTheta=).
#' @inheritParams nlmixr2
#' @param ui rxode2 ui
#' @return list(ui=) possibly with injected etas
#' @export
#' @author Matthew L. Fidler
.preProcessVaeNonMuTheta <- function(ui, est, data, control) {
  if (!inherits(control, "vaeControl")) return(NULL)
  .mode <- if (is.null(control$nonMuTheta)) "regress" else control$nonMuTheta
  ## reset per-fit record of injected etas (read by the VAE output collapse)
  nlmixr2global$nlmixr2EstEnv$vaeNonMuEtas <- character(0)
  ## covariateSelection=FALSE: model-declared covariate coefficients are estimated
  ## by the regress M-step regardless of nonMuTheta (see .vaeDataPrep) -- note it
  .covCoef <- if (isFALSE(control$covariateSelection)) .vaeCovariateCoefThetas(ui) else character(0)
  if (length(.covCoef) > 0L) {
    .pre <- "estimating covariate coef(s): "
    warning(.pre, .vaeTruncList(.covCoef, prefix = .pre), call. = FALSE)
  }
  if (identical(.mode, "none")) return(NULL)
  .thetas <- .vaeNonMuThetas(ui)
  if (length(.thetas) == 0L) return(NULL)
  if (identical(.mode, "regress")) {
    ## "regress": no eta is injected -- the thetas stay plain fixed effects and are
    ## estimated by a bounded bobyqa regression in the VAE M-step (see
    ## .vaeDataPrep / vaeTrainCpp_).  Only surface a note; leave the UI unchanged.
    .pre <- "regressing non-mu theta(s): "
    warning(.pre, .vaeTruncList(.thetas, prefix = .pre), call. = FALSE)
    return(NULL)
  }
  .omega <- if (is.null(control$nonMuEtaOmega)) 0.01 else control$nonMuEtaOmega
  .inj <- .vaeInjectNonMuEtas(ui, .thetas, .omega, fix = identical(.mode, "fix"))
  nlmixr2global$nlmixr2EstEnv$vaeNonMuEtas <- .inj$etas
  .pre <- paste0("non-mu-referenced theta(s) temp eta [", .mode, "]: ")
  warning(.pre, .vaeTruncList(.thetas, prefix = .pre), call. = FALSE)
  list(ui = .inj$ui)
}
preProcessHooksAdd(".preProcessVaeNonMuTheta", .preProcessVaeNonMuTheta)
