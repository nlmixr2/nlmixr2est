#' Normalise `saemControl(etaMat=)`
#'
#' @param etaMat `NULL`, a prior fit, a data frame (optionally with an `ID`
#'   column) or a numeric matrix, one row per subject
#' @return `NULL`, or a numeric matrix of etas whose row names, when present,
#'   are subject IDs
#' @noRd
#' @author Matthew L. Fidler
.saemEtaMatArg <- function(etaMat) {
  if (is.null(etaMat)) {
    return(NULL)
  }
  if (inherits(etaMat, "nlmixr2FitCore")) {
    .m <- etaMat$etaMat
    .id <- etaMat$eta$ID
    if (is.null(.m) || is.null(.id) || length(.id) != nrow(.m)) {
      stop("'etaMat' fit has no etas to start the chains from", call. = FALSE)
    }
    rownames(.m) <- as.character(.id)
    etaMat <- .m
  } else if (is.data.frame(etaMat)) {
    .id <- etaMat$ID
    .m <- as.matrix(.nmDropNonEtaCols(etaMat))
    rownames(.m) <- if (is.null(.id)) NULL else as.character(.id)
    etaMat <- .m
  }
  if (is.matrix(etaMat) && is.integer(etaMat)) {
    storage.mode(etaMat) <- "double"
  }
  checkmate::assertMatrix(
    etaMat,
    mode = "double",
    any.missing = FALSE,
    min.rows = 1,
    min.cols = 1,
    .var.name = "etaMat"
  )
  if (any(!is.finite(etaMat))) {
    stop("'etaMat' must be finite", call. = FALSE)
  }
  .id <- rownames(etaMat)
  if (!is.null(.id) && anyDuplicated(.id)) {
    stop("'etaMat' has duplicated ID(s): ", paste(unique(.id[duplicated(.id)]), collapse = ", "), call. = FALSE)
  }
  etaMat
}

#' Start every SAEM chain at the `saemControl(etaMat=)` etas
#'
#' Each phi1 column becomes `mprior_phi1 + eta` for every chain; a temporary
#' (bounded-transform) eta or an uninformative eta starts at its prior mean.
#'
#' @param cfg `.configsaem()` configuration
#' @param etaMat normalised `etaMat` (`.saemEtaMatArg()`), or `NULL`
#' @param ui rxode2 ui being fitted
#' @param idLvl ID labels of the fitted data, in subject order
#' @param pseudoEtas saem's temporary etas (`.saemPseudoEtaRecord()`)
#' @return `cfg`
#' @noRd
#' @author Matthew L. Fidler
.saemEtaMatCfg <- function(cfg, etaMat, ui, idLvl, pseudoEtas = character(0)) {
  if (is.null(etaMat)) {
    return(cfg)
  }
  if (length(ui$mixProbs) > 0L || isTRUE(cfg$nMix > 1L)) {
    stop("'etaMat' is not supported for mixture models in 'saem'", call. = FALSE)
  }
  if (length(.uiIovEnv$iovVars) > 0L) {
    stop("'etaMat' is not supported for IOV models in 'saem'", call. = FALSE)
  }
  .iniDf <- ui$iniDf
  .eta <- .iniDf[!is.na(.iniDf$neta1) & .iniDf$neta1 == .iniDf$neta2, ]
  .etaNames <- .eta$name
  .need <- setdiff(.etaNames, pseudoEtas)
  .n <- as.integer(cfg$N)
  .cn <- colnames(etaMat)
  if (is.null(.cn)) {
    if (ncol(etaMat) != length(.need)) {
      stop(sprintf("'etaMat' needs %d columns (one per eta), not %d", length(.need), ncol(etaMat)), call. = FALSE)
    }
    colnames(etaMat) <- .need
  } else {
    .bad <- setdiff(.cn, .need)
    if (length(.bad) > 0L) {
      stop("'etaMat' has unknown eta(s): ", paste(.bad, collapse = ", "), call. = FALSE)
    }
    .miss <- setdiff(.need, .cn)
    if (length(.miss) > 0L) {
      stop("'etaMat' has no column for eta(s): ", paste(.miss, collapse = ", "), call. = FALSE)
    }
  }
  .id <- rownames(etaMat)
  if (is.null(.id)) {
    if (nrow(etaMat) != .n) {
      stop(sprintf("'etaMat' needs %d rows (one per subject), not %d", .n, nrow(etaMat)), call. = FALSE)
    }
    .row <- seq_len(.n)
  } else {
    .lvl <- as.character(idLvl)
    .bad <- setdiff(.id, .lvl)
    if (length(.bad) > 0L) {
      stop("'etaMat' has unknown ID(s): ", paste(.bad, collapse = ", "), call. = FALSE)
    }
    .row <- match(.lvl, .id)
    if (anyNA(.row)) {
      stop("'etaMat' has no row for ID(s): ", paste(.lvl[is.na(.row)], collapse = ", "), call. = FALSE)
    }
  }
  .trans <- ui$saemOmegaTrans
  .i1 <- cfg$i1 + 1L
  .phi1 <- cfg$mprior_phi1
  for (.k in seq_along(.etaNames)) {
    .j <- .trans[.k]
    if (is.na(.j) || !(.etaNames[.k] %in% .need)) {
      next
    }
    .phi1[, .j] <- .phi1[, .j] + etaMat[.row, .etaNames[.k]] * cfg$ue[seq_len(.n), .i1[.j]]
  }
  .phiM <- cfg$phiM
  .phiM[, .i1] <- .phi1[rep(seq_len(.n), cfg$nmc), , drop = FALSE]
  .i0 <- cfg$i0 + 1L
  .mean <- unname(rowsum(.phiM, rep(seq_len(.n), cfg$nmc), reorder = TRUE)) / cfg$nmc
  cfg$phiM <- .phiM
  cfg$statphi11 <- .mean[, .i1, drop = FALSE]
  cfg$statphi01 <- .mean[, .i0, drop = FALSE]
  cfg$statphi12 <- crossprod(.phiM[, .i1, drop = FALSE])
  cfg$statphi02 <- crossprod(.phiM[, .i0, drop = FALSE])
  cfg
}
