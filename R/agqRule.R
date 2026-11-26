# Quadrature rules for the agq family.  Every rule is returned in the
# Gauss-Hermite layout inner.cpp consumes: nodes x for the exp(-x'x) kernel
# (eta = center + sqrt(2)*L*x) and per-row weights whose log-sum is the node
# weight.  A rule for the standard normal measure, nodes z and weights p, maps
# onto that layout exactly as x = z/sqrt(2) with p in the first weight column.

# Refuse (or warn about) a grid before building it
.agqNodeGuard <- function(nodes, what, warnNodes = 1e4, maxNodes = 1e6) {
  if (nodes > maxNodes) {
    stop(
      sprintf(
        "quadrature grid of %s = %s nodes exceeds the limit of %s; ",
        what,
        format(nodes, big.mark = ",", scientific = FALSE),
        format(maxNodes, big.mark = ",", scientific = FALSE)
      ),
      "use a smaller grid, nAGQ=1 (laplace) or est=\"focei\"",
      call. = FALSE
    )
  }
  if (nodes > warnNodes) {
    warning(
      sprintf(
        "large quadrature grid: %s = %s nodes",
        what,
        format(nodes, big.mark = ",", scientific = FALSE)
      ),
      call. = FALSE
    )
  }
  invisible(nodes)
}

#' Gauss-Legendre nodes and weights on \[-1, 1\] (Golub-Welsch)
#'
#' @param n number of nodes
#' @return list with `x` (nodes, increasing) and `w` (weights)
#' @noRd
.gaussLegendre <- function(n) {
  if (n == 1L) {
    return(list(x = 0, w = 2))
  }
  .k <- seq_len(n - 1L)
  .b <- .k / sqrt(4 * .k^2 - 1)
  .j <- matrix(0, n, n)
  .j[cbind(.k, .k + 1L)] <- .b
  .j[cbind(.k + 1L, .k)] <- .b
  .e <- eigen(.j, symmetric = TRUE)
  .o <- order(.e$values)
  .x <- .e$values[.o]
  if (n %% 2L == 1L) {
    .x[(n + 1L) / 2L] <- 0 # exact zero so the mode node can be cached
  }
  list(x = .x, w = 2 * .e$vectors[1, .o]^2)
}

# Convert a standard-normal rule (z, p) to the inner.cpp layout
.agqFromNormalRule <- function(z, p, neta, nAGQ, first = FALSE) {
  z <- matrix(z, ncol = neta)
  .w <- matrix(1, nrow(z), neta)
  .w[, 1] <- p
  list(
    x = z / sqrt(2),
    w = .w,
    n = nrow(z),
    neta = neta,
    nAGQ = nAGQ,
    first = first
  )
}

#' Adaptive Gauss-Legendre grid
#'
#' Product Gauss-Legendre rule over `[-nsd, nsd]` standard deviations of the
#' conditional (Laplace) distribution of each eta, weights normalised to 1.
#'
#' @inheritParams .agq
#' @param nsd half-width of the integration box in standard deviations
#' @return same list as [.agq()]
#' @noRd
.agqLegendre <- function(neta = 2, nAGQ = 3, nsd = 3, warnNodes = 1e4, maxNodes = 1e6) {
  .agqNodeGuard(
    as.double(nAGQ)^as.double(neta),
    sprintf("nAGQ^neta = %d^%d", as.integer(nAGQ), as.integer(neta)),
    warnNodes,
    maxNodes
  )
  .gl <- .gaussLegendre(as.integer(nAGQ))
  .z <- nsd * .gl$x
  # box rule for f(z)*phi(z), normalised to sum to 1 so a constant ratio to the
  # Laplace Gaussian (a model linear in eta) integrates exactly
  .p <- .gl$w * stats::dnorm(.z)
  .p <- .p / sum(.p)
  .first <- FALSE
  if (nAGQ %% 2L == 1L) {
    .zero <- (nAGQ + 1L) / 2L
    .z <- c(.z[.zero], .z[-.zero])
    .p <- c(.p[.zero], .p[-.zero])
    .first <- TRUE
  }
  .zg <- as.matrix(do.call("expand.grid", lapply(seq_len(neta), function(i) .z)))
  .pg <- as.matrix(do.call("expand.grid", lapply(seq_len(neta), function(i) .p)))
  .agqFromNormalRule(.zg, apply(.pg, 1, prod), neta, nAGQ, first = .first)
}

#' Unit directions of the Stieltjes grid
#'
#' `ns` directions per quadrant of every coordinate plane, at the midpoints of
#' `ns` equal angles; the set's second moment is proportional to the identity.
#'
#' @param neta dimension
#' @param ns directions per quadrant
#' @return matrix, one unit direction per row
#' @noRd
.agqStieltjesDirections <- function(neta, ns = 1L) {
  if (neta == 1L) {
    return(matrix(c(1, -1), ncol = 1L))
  }
  .th <- (seq_len(ns) - 0.5) * pi / (2 * ns)
  .pairs <- utils::combn(neta, 2L)
  .s <- rbind(c(1, 1), c(1, -1), c(-1, 1), c(-1, -1))
  .rows <- lapply(seq_len(ncol(.pairs)), function(k) {
    .i <- .pairs[1, k]
    .j <- .pairs[2, k]
    do.call(
      rbind,
      lapply(.th, function(t) {
        .u <- matrix(0, 4L, neta)
        .u[, .i] <- .s[, 1] * cos(t)
        .u[, .j] <- .s[, 2] * sin(t)
        .u
      })
    )
  })
  do.call(rbind, .rows)
}

#' Stieltjes (radial) grid
#'
#' Partitions the standard normal measure into ellipsoidal shells by the
#' probability of the squared radius: a central region `[0, r0]`, `nr` equal
#' shells over `[r0, r1]`, and a tail `[r1, 1]` unless `r1 > 0.9999`.  Each
#' shell carries its probability, split evenly over the directions of
#' `.agqStieltjesDirections()`, at the radius matching the shell's second
#' moment, so the rule integrates polynomials of degree 3 exactly.
#'
#' @param neta dimension
#' @param grid `c(nr, ns, r0, r1)`, as in NONMEM's `GRID=`
#' @inheritParams .agq
#' @return same list as [.agq()]
#' @noRd
.agqStieltjes <- function(neta = 2, grid = c(3, 1, 0.6, 0.9), warnNodes = 1e4, maxNodes = 1e6) {
  .nr <- as.integer(grid[1])
  .ns <- as.integer(grid[2])
  .r0 <- grid[3]
  .r1 <- grid[4]
  .p <- c(0, .r0 + (.r1 - .r0) * (0:.nr) / .nr)
  if (.r1 <= 0.9999) {
    .p <- c(.p, 1)
  }
  .nd <- if (neta == 1L) 2 else 2 * .ns * neta * (neta - 1)
  .agqNodeGuard(
    (length(.p) - 1) * .nd,
    sprintf("%d shells x %d directions", length(.p) - 1L, as.integer(.nd)),
    warnNodes,
    maxNodes
  )
  .u <- .agqStieltjesDirections(neta, .ns)
  .q <- stats::qchisq(.p, neta)
  .pr <- diff(.p)
  # E[R^2 | shell] for R^2 ~ chisq(neta)
  .r2 <- neta * diff(stats::pchisq(.q, neta + 2)) / .pr
  .z <- do.call(rbind, lapply(sqrt(.r2), function(r) r * .u))
  .w <- rep(.pr / nrow(.u), each = nrow(.u))
  .agqFromNormalRule(.z, .w / sum(.w), neta, nrow(.z))
}

#' Quadrature grid for the agq family from a control
#'
#' @param neta number of etas
#' @param control a `foceiControl()`-family control
#' @inheritParams .agq
#' @return same list as [.agq()]
#' @noRd
.agqRuleGrid <- function(neta, control, warnNodes = 1e4, maxNodes = 1e6) {
  .rule <- control$agqRule
  if (is.null(.rule)) {
    .rule <- "hermite"
  }
  switch(
    .rule,
    hermite = .agq(neta, control$nAGQ, warnNodes = warnNodes, maxNodes = maxNodes),
    legendre = .agqLegendre(neta, control$nAGQ, control$agqNsd, warnNodes, maxNodes),
    stieltjes = .agqStieltjes(neta, control$agqGrid, warnNodes, maxNodes)
  )
}

# Objective label for the quadrature rule
.agqRuleLabel <- function(control) {
  switch(
    if (is.null(control$agqRule)) "hermite" else control$agqRule,
    hermite = "",
    legendre = paste0("AGQ-GL", control$nAGQ),
    stieltjes = "AGQ-Stieltjes"
  )
}
