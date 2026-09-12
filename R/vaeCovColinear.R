# Colinearity clusters for the VAE covariate M-step.  A cluster NAMES a set of
# interchangeable covariates; it never constrains what may be selected.  Two
# consumers use it: the hysteresis pass, which refuses to displace last
# iteration's incumbent for a cluster mate that does not clearly beat it, and the
# near-tie diagnostic, which reports the mates that came close.
#
# covMat is constant across a fit, so this is computed once in R.

#' Default `abs(cor)` at which two covariates cluster.
#'
#' Shared by [vaeControl()] and [vaeCovariates()] so the two cannot drift.
#' Measured justification is in `tools/vaeColinearPreflight.R`.
#' @noRd
.vaeColinearCut <- 0.9

#' Colinearity cluster id per covariate column.
#'
#' Clusters are built at the GROUP level and broadcast back to columns, so the
#' result is always a COARSENING of `group`: every group lies entirely inside one
#' cluster.  Only CROSS-group column pairs are ever compared, which is what keeps
#' the two shapes of one covariate out of a cluster -- they share a group, and
#' the existing mutual-exclusion machinery already arbitrates them.  The arms of
#' a hockey block are excluded by the same rule, since a block never spans groups.
#'
#' Groups are joined by single linkage (connected components).  At a threshold
#' this high, chaining is rare, and because a cluster is a label rather than a
#' constraint, an over-generous component costs a few extra swap scores -- never
#' an unreachable model.  Components are order-independent, which a greedy
#' pairwise assignment would not be.
#'
#' @param covMat covariate design, one column per shape column
#' @param group `covGroup`, one entry per column
#' @param cut `abs(cor)` at which two groups join
#' @return integer cluster id per column, numbered by first appearance
#' @noRd
.vaeCovCluster <- function(covMat, group, cut = .vaeColinearCut) {
  .n <- if (is.null(covMat)) 0L else ncol(covMat)
  if (.n == 0L) return(integer(0))
  .g <- as.integer(group)
  if (length(.g) != .n) {
    stop("'group' must have one entry per covariate column", call. = FALSE)
  }
  ## first appearance order, matching how covGroup/covBlock are numbered
  .idOf <- function(x) match(x, unique(x))
  ## Nothing to join: hand back the groups themselves, which is the identity
  ## coarsening and makes .vaeClusterBinds() FALSE.
  if (!.vaeCanCluster(covMat, .g, cut)) return(.idOf(.g))
  ## a constant column correlates with nothing; excluding it here also keeps
  ## stats::cor from emitting a zero-variance warning at the user
  .sd <- apply(covMat, 2, stats::sd)
  .ok <- which(is.finite(.sd) & .sd > 0)
  if (length(.ok) < 2L) return(.idOf(.g))
  .r <- abs(stats::cor(covMat[, .ok, drop = FALSE]))
  .r[!is.finite(.r)] <- 0
  .gu <- unique(.g)
  .idOf(.vaeGroupLink(.r, .g[.ok], .gu, cut)[match(.g, .gu)])
}

#' Is there anything for the correlation pass to join?
#'
#' Fewer than two columns, fewer than three rows, a single group or a
#' non-finite cut/design all mean the identity coarsening is the answer.
#'
#' @inheritParams .vaeCovCluster
#' @param g integer `group`
#' @return single logical
#' @noRd
.vaeCanCluster <- function(covMat, g, cut) {
  ncol(covMat) >= 2L && nrow(covMat) >= 3L && length(unique(g)) >= 2L &&
    is.finite(cut) && all(is.finite(covMat))
}

#' Single-linkage components over the covariate GROUPS.
#'
#' Only cross-group pairs are edges, so a group never links to itself.  A union
#' relabels the whole component rather than one entry, which is what makes the
#' result independent of the order the edges are visited.
#'
#' @param r `abs(cor)` among the usable columns
#' @param gok group id of each usable column
#' @param gu unique group ids, first-appearance order
#' @param cut `abs(cor)` at which two groups join
#' @return component id, one per entry of `gu`
#' @noRd
.vaeGroupLink <- function(r, gok, gu, cut) {
  .comp <- seq_along(gu)
  for (.a in seq_along(gok)) {
    for (.b in seq_len(.a - 1L)) {
      if (gok[.a] == gok[.b] || r[.a, .b] < cut) next
      .ia <- .comp[match(gok[.a], gu)]
      .ib <- .comp[match(gok[.b], gu)]
      if (.ia == .ib) next
      .comp[.comp == max(.ia, .ib)] <- min(.ia, .ib)
    }
  }
  .comp
}

#' Does any cluster actually merge two covariate groups?
#'
#' The gate for sending the cluster vector to C++.  It is deliberately NOT
#' `anyDuplicated(cluster) > 0L`, the idiom `covGroup`/`covBlock` use: a cluster
#' id repeats on every multi-shape covariate even when no cross-covariate
#' colinearity exists at all, so that predicate would ship the vector on ordinary
#' designs and silently switch the new mechanisms on everywhere.
#' @param cluster from [.vaeCovCluster()]
#' @param group `covGroup`
#' @return single logical
#' @noRd
.vaeClusterBinds <- function(cluster, group) {
  if (is.null(cluster) || is.null(group)) return(FALSE)
  if (!length(cluster) || length(cluster) != length(group)) return(FALSE)
  length(unique(cluster)) < length(unique(group))
}

#' Run-time note when correlated latent dims were found under a diagonal omega.
#'
#' With a diagonal omega the covariate objective is separable across latent
#' dimensions, so the cross-parameter refinement provably cannot improve anything
#' and is skipped.  Declaring the correlation is what enables it.
#' @param nPair correlated pairs found
#' @param omOff whether the model declares off-diagonal omega, as REPORTED by the
#'   fit rather than re-derived: the C++ flag comes from the omega selection
#'   structure, so a declared block whose ini covariance is exactly 0 would make
#'   an R-side check disagree with the gate that actually ran
#' @param anySel whether any covariate was selected at all
#' @return character(0) or a single message
#' @noRd
.vaePhiDiagMsg <- function(nPair, omOff, anySel) {
  if (!length(nPair) || is.na(nPair) || nPair <= 0L) return(character(0))
  if (isTRUE(omOff) || !isTRUE(anySel)) return(character(0))
  "correlated etas found; declare an omega block to refine jointly"
}
