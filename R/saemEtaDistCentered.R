# Centered (interweaved) M-step for declared ETA distributions.
#
# BHHH moves a declared family's thetas with the unit-normal latents held
# fixed, so every decoded eta moves with them; when the data pin the
# individual parameters that step is tiny and the family's location barely
# converges.  After it, this step holds the DECODED etas fixed instead: it
# maximizes their family density (plus the Gaussian copula of their latents)
# over the family thetas and rxCor, then re-encodes the latents so the decoded
# etas, and therefore every prediction, are unchanged (Yu and Meng 2011).

.saemEtaDistRFamilies <- c("dgamma", "dbeta", "dlnorm", "dweibull", "dexp",
                           "dnorm", "dlogis", "dcauchy", "dunif", "dchisq")

#' Build the centered M-step for SAEM's declared ETA distributions
#'
#' @param ui expanded rxode2 ui (`rxEtaDistExpand(param = "cdf")`)
#' @param cfg `.configsaem()` configuration
#' @param families named character: the declared family of each eta
#' @return `NULL` when the model has no family this step handles, otherwise a
#'   function of `(phiM, phi0, gain)` returning `NULL` or a list with the new
#'   `theta` (phi0 order), the re-encoded latent `z` and their 0-based `zcols`
#' @noRd
.saemEtaDistCenteredFn <- function(ui, cfg, families) {
  if (length(families) == 0L || !all(families %in% .saemEtaDistRFamilies)) return(NULL)
  .asg <- Filter(function(x) {
    is.call(x) && (identical(x[[1]], quote(`<-`)) || identical(x[[1]], quote(`=`))) &&
      is.name(x[[2]])
  }, as.list(ui$lstExpr))
  .lhs <- vapply(.asg, function(x) as.character(x[[2]]), character(1))
  # the parameter prelude: every generated rx* line (latents, copula, family
  # arguments, and any bounded-transform back-mapping such as rxCor)
  .pre <- .asg[startsWith(.lhs, "rx")]
  .preLhs <- .lhs[startsWith(.lhs, "rx")]
  .etas <- names(families)
  .pn <- ui$saemParamsToEstimate
  .zcol <- match(paste0("rxz.", .etas), .pn)
  .ncol <- paste0("rxN.", .etas)
  if (anyNA(.zcol) || !all(.ncol %in% .preLhs)) return(NULL)
  # a copula that also carries an ORDINARY eta (its latent rxz.* has no declared
  # family here) is outside this step: decline, and the fit keeps the other updates
  .zUsed <- unique(unlist(lapply(.pre, function(x) grep("^rxz\\.", all.vars(x[[3]]), value = TRUE))))
  if (!all(.zUsed %in% paste0("rxz.", .etas))) return(NULL)
  .ini <- ui$iniDf
  .th <- !is.na(.ini$ntheta)
  .thetaAll <- stats::setNames(.ini$est[.th], .ini$name[.th])
  .names0 <- .pn[cfg$i0 + 1L]
  .paramLines <- .pre[!startsWith(.preLhs, "rxN.")]
  .vars <- setdiff(unique(unlist(lapply(.paramLines, function(x) all.vars(x[[3]])))), .preLhs)
  # arguments that depend on covariates or etas are not handled here; the
  # current values come from phi0 (rxCor is a phi0 parameter, not an ini theta)
  if (!all(.vars %in% union(names(.thetaAll), .names0))) return(NULL)
  .free <- setdiff(intersect(.vars, .names0), .ini$name[.th & .ini$fix])
  if (length(.free) == 0L) return(NULL)
  .args <- lapply(.etas, function(e) {
    .a <- .preLhs[startsWith(.preLhs, paste0("rxEdA.", e, "."))]
    stats::setNames(.a, substring(.a, nchar(paste0("rxEdA.", e, ".")) + 1L))
  })
  .dfun <- lapply(families, function(f) get(f, envir = asNamespace("stats")))
  .pfun <- lapply(families, function(f) get(sub("^d", "p", f), envir = asNamespace("stats")))
  .qfun <- lapply(families, function(f) get(sub("^d", "q", f), envir = asNamespace("stats")))
  .k <- length(.etas)
  .eval <- function(theta, z) {
    .env <- list2env(as.list(.thetaAll), parent = asNamespace("rxode2"))
    for (.n in names(theta)) assign(.n, theta[[.n]], envir = .env)
    for (.j in seq_len(.k)) assign(paste0("rxz.", .etas[.j]), z[, .j], envir = .env)
    for (.x in .pre) assign(as.character(.x[[2]]), eval(.x[[3]], .env), envir = .env)
    .env
  }
  # the latents enter linearly: N = L z
  .chol <- function(theta) {
    .env <- .eval(theta, diag(.k))
    L <- t(vapply(.ncol, function(n) rep_len(get(n, envir = .env), .k), numeric(.k)))
    list(L = matrix(L, .k, .k), env = .env)
  }
  .famArgs <- function(env, j) lapply(.args[[j]], function(n) get(n, envir = env)[1L])
  .encode <- function(theta, eta) {
    .c <- .chol(theta)
    N <- vapply(seq_len(.k), function(j) {
      .p <- do.call(.pfun[[j]], c(list(eta[, j]), .famArgs(.c$env, j)))
      stats::qnorm(pmin(pmax(.p, 1e-15), 1 - 1e-15))
    }, numeric(nrow(eta)))
    N <- matrix(N, nrow(eta), .k)
    list(N = N, z = t(solve(.c$L, t(N))), L = .c$L, env = .c$env)
  }
  # per-subject complete-data log-likelihood of the decoded etas: family
  # densities plus the Gaussian copula log-density of their latent normals
  .rowLik <- function(theta, eta) {
    .e <- .encode(theta, eta)
    .v <- -0.5 * rowSums(.e$z^2) + 0.5 * rowSums(.e$N^2) - log(abs(det(.e$L)))
    for (j in seq_len(.k)) {
      .v <- .v + do.call(.dfun[[j]], c(list(eta[, j]), .famArgs(.e$env, j), list(log = TRUE)))
    }
    .v[!is.finite(.v)] <- -Inf
    .v
  }
  # warm-started BHHH ascent on the per-subject contributions: the optimum
  # moves little between SAEM iterations, so a few steps reach it where a
  # general optimizer with numeric gradients spends ~100 full evaluations
  .bhhh <- function(theta, eta, maxit = 20L) {
    .p <- length(.free)
    .l0 <- .rowLik(theta, eta)
    if (!all(is.finite(.l0))) return(NULL)
    .s0 <- sum(.l0)
    for (it in seq_len(maxit)) {
      .h <- .Machine$double.eps^(1 / 3) * pmax(1, abs(theta[.free]))
      G <- vapply(seq_len(.p), function(q) {
        .tp <- .tm <- theta
        .tp[.free[q]] <- .tp[.free[q]] + .h[q]
        .tm[.free[q]] <- .tm[.free[q]] - .h[q]
        (.rowLik(.tp, eta) - .rowLik(.tm, eta)) / (2 * .h[q])
      }, numeric(length(.l0)))
      G <- matrix(G, length(.l0), .p)
      if (!all(is.finite(G))) return(NULL)
      .g <- colSums(G)
      .d <- try(solve(crossprod(G), .g), silent = TRUE)
      if (inherits(.d, "try-error") || !all(is.finite(.d))) return(NULL)
      .a <- 1
      repeat {
        .t <- theta
        .t[.free] <- .t[.free] + .a * .d
        .l <- .rowLik(.t, eta)
        .s <- sum(.l)
        if (is.finite(.s) && .s >= .s0) break
        .a <- .a / 2
        if (.a < 1e-4) return(theta)
      }
      .dec <- .s - .s0
      theta <- .t
      .l0 <- .l
      .s0 <- .s
      if (.dec < 1e-9 * length(.l0) || max(abs(.a * .d)) < 1e-8) break
    }
    theta
  }
  function(phiM, phi0, gain) {
    theta <- stats::setNames(as.numeric(phi0), .names0)[unique(c(.free, intersect(.vars, .names0)))]
    z <- phiM[, .zcol, drop = FALSE]
    .c <- .chol(theta)
    N <- z %*% t(.c$L)
    eta <- vapply(seq_len(.k), function(j) {
      do.call(.qfun[[j]], c(list(stats::pnorm(N[, j])), .famArgs(.c$env, j)))
    }, numeric(nrow(z)))
    eta <- matrix(eta, nrow(z), .k)
    if (!all(is.finite(eta))) return(NULL)
    .opt <- .bhhh(theta, eta)
    if (is.null(.opt)) return(NULL)
    .new <- theta
    .new[.free] <- theta[.free] + gain * (.opt[.free] - theta[.free])
    .z <- .encode(.new, eta)$z
    if (!all(is.finite(.z))) return(NULL)
    .out <- stats::setNames(as.numeric(phi0), .names0)
    .out[.free] <- .new[.free]
    list(theta = unname(.out), z = .z, zcols = as.integer(.zcol - 1L))
  }
}

#' Declared family of each eta, from the declaration stash
#'
#' @param ui expanded rxode2 ui
#' @return named character (eta name -> family function name, e.g.
#'   `"dgamma"`); empty when nothing is declared
#' @noRd
.saemEtaDistFamilies <- function(ui) {
  .st <- .etaDistDeclGet(ui)
  if (is.null(.st) || length(.st$etaDist) == 0L) return(character(0))
  .fam <- vapply(.st$etaDist, function(.t) {
    .cl <- try(str2lang(.t), silent = TRUE)
    if (inherits(.cl, "try-error") || !is.call(.cl)) NA_character_ else as.character(.cl[[1]])
  }, character(1), USE.NAMES = FALSE)
  if (anyNA(.fam)) return(character(0))
  stats::setNames(.fam, .st$name)
}
