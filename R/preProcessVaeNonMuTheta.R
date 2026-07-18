# VAE-only preprocessing: est="vae" only estimates parameters that occupy the
# latent space (structural thetas carrying a random effect, plus omega and
# residual error).  A structural theta with no eta stays frozen at its ini()
# value.  This hook injects a small eta into each such non-mu-referenced theta so
# it enters the latent space and is estimated; the temporary eta is dropped from
# the reported model (collapsed to `theta + mean(eta)`) in the VAE output path.
# Reuses the rxode2 model()/ini() piping to add the eta and rmEta()/.downgradeEtas
# to drop it -- no new UI-surgery machinery.

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
  .mu <- ui$muRefDataFrame$theta
  .cov <- ui$muRefCovariateDataFrame$theta
  .cand <- setdiff(.th$name, c(.mu, .cov))
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
    .idx <- .idx[1]
    ## inject FLAT (no wrapping parens) so `exp(theta + eta)` stays the additive
    ## form rxode2 recognizes as a mu-referenced exp() parameter (parens would hide
    ## the exp() back-transform -- see .vaeUpdateModel)
    .newTxt <- gsub(paste0("\\b", .p, "\\b"), paste0(.p, " + ", .eta),
                    deparse1(.lines[[.idx]]))
    .ui <- do.call(rxode2::model, list(.ui, str2lang(.newTxt)))
    .ui <- do.call(rxode2::ini, list(.ui, str2lang(paste0(.eta, " ~ ", signif(omega, 12)))))
    .etas <- c(.etas, .eta)
  }
  if (fix && length(.etas) > 0L) {
    ## rxode2 ini(eta ~ fix(v)) does not flag the omega; set the diagonal fix flag
    ## directly (read by .vaeDataPrep's omegaFix detection) and rebuild.
    .ui <- rxode2::rxUiDecompress(.ui)
    .df <- .ui$iniDf
    .df$fix[!is.na(.df$neta1) & .df$neta1 == .df$neta2 & .df$name %in% .etas] <- TRUE
    assign("iniDf", .df, envir = .ui)
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
#' @return list(ui=) possibly with injected etas
#' @export
#' @author Matthew L. Fidler
.preProcessVaeNonMuTheta <- function(ui, est, data, control) {
  if (!inherits(control, "vaeControl")) return(NULL)
  .mode <- if (is.null(control$nonMuTheta)) "eta" else control$nonMuTheta
  ## reset per-fit record of injected etas (read by the VAE output collapse)
  nlmixr2global$nlmixr2EstEnv$vaeNonMuEtas <- character(0)
  if (identical(.mode, "none")) return(NULL)
  .thetas <- .vaeNonMuThetas(ui)
  if (length(.thetas) == 0L) return(NULL)
  .omega <- if (is.null(control$nonMuEtaOmega)) 0.01 else control$nonMuEtaOmega
  .inj <- .vaeInjectNonMuEtas(ui, .thetas, .omega, fix = identical(.mode, "fix"))
  nlmixr2global$nlmixr2EstEnv$vaeNonMuEtas <- .inj$etas
  warning("est=\"vae\": the following non-mu-referenced population parameter(s) were ",
          "temporarily given a small eta (nonMuTheta=\"", .mode, "\") so they can be ",
          "estimated, and are reported as theta+mean(eta) with the eta dropped: ",
          paste(.thetas, collapse = ", "), call. = FALSE)
  list(ui = .inj$ui)
}
preProcessHooksAdd(".preProcessVaeNonMuTheta", .preProcessVaeNonMuTheta)
