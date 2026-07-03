# Analytic FOCEI observed-information covariance R-matrix for nlmixr2est.
#
# Replaces the finite-difference covariance Hessian with the exact analytic
# observed information: a data term (2nd-order sensitivities) plus a
# log-determinant term (3rd-order sensitivities, via Shi(2021) finite differences
# of the 2nd-order sensitivities).  Any per-subject solve failure returns NULL and
# the caller falls back to the finite-difference covariance.

#' Are the rxode2 sensitivity primitives available for the analytic covariance
#'
#' Gate for the analytic-covariance scope: needs rxExpandSens2_ from rxode2 and
#' the symengine package.
#' @return logical, TRUE when the primitives are present
#' @author Hidde van de Beek
#' @noRd
.hasRxExpandSens2 <- function() {
  exists("rxExpandSens2_", envir = asNamespace("rxode2"), inherits = FALSE) &&
    requireNamespace("symengine", quietly = TRUE)
}

#' Omega variance-covariance derivatives for the analytic Omega block
#'
#' First and second derivatives of Omega^{-1} and log|Omega| with respect to the
#' free variance-covariance elements, via rxode2's C++ rxOmegaVarCovDeriv.
#' @param Om the estimated Omega matrix
#' @param pairs matrix of free lower-triangle element indices, each row c(a, b)
#'   with a >= b
#' @return list with nom, dOi (list of dOmega^{-1}/dw), d2Oi (list-of-lists) and
#'   d2LD (second derivative of log|Omega|)
#' @author Hidde van de Beek
#' @noRd
.omegaVarCovDeriv <- function(Om, pairs) {
  .nom <- nrow(pairs)
  .d <- rxode2::rxOmegaVarCovDeriv(Om, order = 2L)
  .key <- function(.p) paste(.p[, 1], .p[, 2], sep = "-")
  .idx <- match(.key(pairs), .key(.d$elements))
  list(nom = .nom,
       dOi = .d$dOmegaInv[.idx],
       d2Oi = lapply(.idx, function(.a) .d$d2OmegaInv[[.a]][.idx]),
       d2LD = .d$d2LogDet[.idx, .idx, drop = FALSE])
}

#' Residual error machinery for the theta + sigma + Omega R-matrix
#'
#' Builds the symbolic f-derivatives of rho and p plus the sigma partials, in
#' (f, y, sa, sp) with sa the additive SD and sp the proportional SD.  Supports
#' add / prop / combined1 / combined2; returns NULL for lnorm / propT / propF.
#' @param ui the rxode2 UI object
#' @return list of symbolic derivative expressions and an evaluator, or NULL for
#'   out-of-scope error models (lnorm / propT / propF)
#' @author Hidde van de Beek
#' @noRd
.foceiAnalyticErrFull <- function(ui) {
  # Single endpoint only: a multi-endpoint model pools error rows into one R against
  # one rx_pred_, the wrong likelihood -- bow out to the finite-difference cov.
  if (!is.null(ui$predDf) && nrow(ui$predDf) != 1L) {
    return(NULL)
  }
  .ini <- ui$iniDf
  .er <- .ini[!is.na(.ini$err), , drop = FALSE]
  # only additive / proportional error is supported; any other rider (lnorm, propT/
  # propF, or a DV/variance transform such as boxCox / yeoJohnson / pow) is out of
  # scope -> finite differences.
  if (!all(.er$err %in% c("add", "prop"))) {
    return(NULL)
  }
  # model-declared addProp wins; the control applies only when the model says "default"
  .addPr <- as.character(ui$predDf$addProp)
  if (length(.addPr) != 1L || is.na(.addPr) || .addPr == "default") {
    .addPr <- tryCatch(rxode2::rxGetControl(ui, "addProp", "combined2"),
                       error = function(e) "combined2")
  }
  .addN <- .er$name[.er$err == "add"]
  .propN <- .er$name[.er$err == "prop"]
  .hasA <- length(.addN) == 1L
  .hasP <- length(.propN) == 1L
  if (!.hasA && !.hasP) {
    return(NULL)
  }
  # combined1 is (sa + sp f)^2, combined2 is sa^2 + sp^2 f^2.
  .Rstr <- if (.hasA && .hasP) {
    if (identical(.addPr, "combined1")) {
      "(sa+sp*f)^2"
    } else {
      "sa^2+sp^2*f^2"
    }
  } else if (.hasP) {
    "sp^2*f^2"
  } else {
    "sa^2"
  }
  .Rq <- parse(text = .Rstr)[[1]]
  .rhoE <- bquote(0.5 * ((y - f)^2 / .(.Rq) + log(.(.Rq))))
  .pE <- bquote(1 / .(.Rq) + 0.5 * (.(D(.Rq, "f")) / .(.Rq))^2)
  .DD <- function(.e, ...) {
    for (.v in c(...)) {
      .e <- D(.e, .v)
    }
    .e
  }
  .sgVar <- c(if (.hasA) "sa", if (.hasP) "sp")
  .sgName <- c(if (.hasA) .addN, if (.hasP) .propN)
  .val <- setNames(c(if (.hasA) .er$est[.er$name == .addN],
                     if (.hasP) .er$est[.er$name == .propN]),
                   .sgVar)
  .sc <- list(r1 = .DD(.rhoE, "f"),
              r2 = .DD(.rhoE, "f", "f"),
              r3 = .DD(.rhoE, "f", "f", "f"),
              p = .pE,
              p1 = .DD(.pE, "f"),
              p2 = .DD(.pE, "f", "f"))
  .per <- list()
  for (.s in .sgVar) {
    .per[[.s]] <- list(rf = .DD(.rhoE, "f", .s),
                       rff = .DD(.rhoE, "f", "f", .s),
                       ps = .DD(.pE, .s),
                       pf = .DD(.pE, "f", .s))
  }
  .pair <- list()
  for (.i in seq_along(.sgVar)) {
    for (.j in .i:length(.sgVar)) {
      .a <- .sgVar[.i]
      .b <- .sgVar[.j]
      .pair[[paste0(.a, .b)]] <- list(rss = .DD(.rhoE, .a, .b),
                                      rfss = .DD(.rhoE, "f", .a, .b),
                                      pss = .DD(.pE, .a, .b))
    }
  }
  list(sgVar = .sgVar,
       sgName = .sgName,
       sc = .sc,
       per = .per,
       pair = .pair,
       ev = function(.e, f, y) eval(.e, c(list(f = f, y = y), as.list(.val))))
}

#' Augmented sensitivity model over an arbitrary direction set
#'
#' Builds an rxode2 model carrying prediction sensitivities over a mixed set of
#' ETA_i_ and THETA_j_ directions (a mu-referenced theta reuses its eta's
#' direction; any other structural theta gets its own THETA_j_ direction).  State
#' sensitivities use rxode2's .rxSens; the higher-order prediction chain f1/f2 is
#' built here.  Returns NULL if the 1st-order sensitivities are unavailable.
#' @param ui the rxode2 UI object
#' @param dirs character vector of direction names (ETA_i_ / THETA_j_)
#' @return list(augMod, dirs, ndir, st, P2), or NULL on any build failure
#' @author Hidde van de Beek
#' @noRd
.foceiAnalyticAugModelDirs <- function(ui, dirs) {
  tryCatch({
    .s <- ui$loadPruneSens
    .st <- rxode2::rxStateOde(.s)
    rxode2::.rxJacobian(.s, c(.st, dirs))
    # 1st-order sensitivities are already expanded for the gradient (free); if
    # unavailable the model is not differentiable, so bail before the 2nd-order build.
    .s1 <- rxode2::.rxSens(.s, dirs)
    if (length(.s1) == 0L) {
      return(NULL)
    }
    .s2 <- rxode2::.rxSens(.s, dirs, dirs)          # 2nd order: the expensive expansion
    .pred <- get("rx_pred_", .s)
    .Dn <- function(.e, .v) symengine::D(.e, symengine::S(.v))
    .sn1 <- function(.j, ...) {
      symengine::S(paste0("rx__sens_", .j, "_BY_", paste(c(...), collapse = "_BY_"), "__"))
    }
    .toRx <- function(.l) rxode2::rxFromSE(.l)
    .g1 <- function(.ex, .p) {
      .e <- .Dn(.ex, .p)
      for (.j in .st) {
        .e <- .e + .Dn(.ex, .j) * .sn1(.j, .p)
      }
      .e
    }
    .g2 <- function(.ex, .p, .q) {
      .gq <- .g1(.ex, .q)
      .e <- .Dn(.gq, .p)
      for (.k in .st) {
        .e <- .e + .Dn(.gq, .k) * .sn1(.k, .p)
      }
      for (.j in .st) {
        .e <- .e + .Dn(.ex, .j) * .sn1(.j, .p, .q)
      }
      .e
    }
    .P2 <- expand.grid(i = dirs, j = dirs, stringsAsFactors = FALSE)
    .fL1 <- vapply(dirs, function(.p) {
      paste0("f1_", .p, "=", .toRx(.g1(.pred, .p)))
    }, character(1))
    .fL2 <- vapply(seq_len(nrow(.P2)), function(.r) {
      paste0("f2_", .P2$i[.r], "_", .P2$j[.r], "=",
             .toRx(.g2(.pred, .P2$i[.r], .P2$j[.r])))
    }, character(1))
    .baseOde <- vapply(.st, function(.x) {
      paste0("d/dt(", .x, ")=", .toRx(get(paste0("rx__d_dt_", .x, "__"), .s)))
    }, character(1))
    .modTxt <- paste(c(.baseOde, .s1, .s2, paste0("predf=", .toRx(.pred)),
                       .fL1, .fL2),
                     collapse = "\n")
    .modTxt <- gsub("ETA\\[([0-9]+)\\]", "ETA_\\1_", .modTxt)
    .modTxt <- gsub("THETA\\[([0-9]+)\\]", "THETA_\\1_", .modTxt)
    list(augMod = rxode2::rxode2(.modTxt),
         dirs = dirs,
         ndir = length(dirs),
         st = .st,
         P2 = .P2)
  }, error = function(e) NULL)
}

#' Solve the direction-set augmented model (order 2) with a finite-difference 3rd order
#'
#' Solves the 2nd-order model and recovers the 3rd-order tensor Ath by Shi (2021)
#' adaptive central differences of the exact 2nd-order sensitivities.  The result
#' is symmetrized and returned as list(f, a, A, Ath) for the R-matrix assembly.
#' @param aug the augmented model list (order 2) from .foceiAnalyticAugModelDirs
#' @param params named parameter vector (thetas and eta directions)
#' @param ev the subject's event data frame
#' @param times observation times to keep
#' @param ef expected function error for the Shi adaptive step
#' @return list(f, a, A, Ath), or NULL on any failure
#' @author Hidde van de Beek
#' @noRd
.foceiAnalyticSolveDirFD3 <- function(aug, params, ev, times, ef = 7e-7, solveOpts = NULL) {
  .dirs <- aug$dirs
  .nd <- length(.dirs)
  .solveA <- function(.p) {
    .args <- c(list(aug$augMod, params = .p, ev, returnType = "data.frame",
                    atol = 1e-10, rtol = 1e-10),         # tight tol for the sensitivity states
               solveOpts)                                # + the fit's covsInterpolation / method
    .d <- tryCatch(
      withCallingHandlers(
        as.data.frame(do.call(rxode2::rxSolve, .args)),
        warning = function(w) invokeRestart("muffleWarning")),  # muffle benign warnings, don't abort the tier
      error = function(e) NULL)
    if (is.null(.d)) {
      return(NULL)
    }
    .d <- .d[.d$time %in% times, , drop = FALSE]
    if (nrow(.d) != length(times) || !all(c("predf", paste0("f1_", .dirs)) %in% names(.d))) {
      return(NULL)   # extra solve rows (EVID=2, dose at an obs time) misalign f vs y -> bow out to FD
    }
    .no <- nrow(.d)
    .a <- matrix(vapply(.dirs, function(.q) .d[[paste0("f1_", .q)]], numeric(.no)),
                 .no, .nd)
    .A <- array(0, c(.no, .nd, .nd))
    for (.r in seq_len(nrow(aug$P2))) {
      .A[, match(aug$P2$i[.r], .dirs), match(aug$P2$j[.r], .dirs)] <-
        .d[[paste0("f2_", aug$P2$i[.r], "_", aug$P2$j[.r])]]
    }
    list(f = .d$predf, a = .a, A = .A)
  }
  .E0 <- .solveA(params)
  if (is.null(.E0)) {
    return(NULL)
  }
  .no <- length(.E0$f)
  .f0 <- as.vector(.E0$A)
  .Aflat <- function(.p) {
    .e <- .solveA(.p)
    if (is.null(.e) || !all(is.finite(.e$A))) {
      return(NULL)
    }
    as.vector(.e$A)
  }
  .Ath <- array(0, c(.no, .nd, .nd, .nd))
  for (.d in seq_len(.nd)) {
    .idx <- match(.dirs[.d], names(params))
    if (is.na(.idx)) {
      return(NULL)
    }
    .ffun <- function(.tt) .Aflat(setNames(.tt, names(params)))
    .shi <- .shi21CentralR(.ffun, params, .idx, .f0, ef = ef)
    if (is.null(.shi$gr) || !all(is.finite(.shi$gr))) {
      return(NULL)
    }
    .Ath[, , , .d] <- array(.shi$gr, c(.no, .nd, .nd))
  }
  .S <- array(0, c(.no, .nd, .nd, .nd))
  .pr <- list(c(1, 2, 3), c(1, 3, 2), c(2, 1, 3), c(2, 3, 1), c(3, 1, 2), c(3, 2, 1))
  for (.l in seq_len(.nd)) {
    for (.m in seq_len(.nd)) {
      for (.n in seq_len(.nd)) {
        .ix <- c(.l, .m, .n)
        .S[, .l, .m, .n] <-
          Reduce(`+`, lapply(.pr, function(.p) .Ath[, .ix[.p[1]], .ix[.p[2]], .ix[.p[3]]])) / 6
      }
    }
  }
  .out <- list(f = .E0$f, a = .E0$a, A = .E0$A, Ath = .S)
  if (!all(is.finite(.out$f)) || !all(is.finite(.out$a)) ||
      !all(is.finite(.out$A)) || !all(is.finite(.out$Ath))) {
    return(NULL)
  }
  .out
}

#' Harmonic-mean gate for the Shi (2021) adaptive step
#'
#' R port of nlmixr2's shi21Central (src/shi21.cpp), which is not exported to R and
#' so cannot difference the R closure used here.  Collapses the per-component
#' third-difference ratios to one value by their harmonic mean, with a correction
#' for locally linear (zero third difference) components.
#' @param allv vector of per-component third-difference ratios
#' @return the corrected harmonic mean
#' @author Hidde van de Beek
#' @noRd
.shiHarmonicMean <- function(allv) {
  if (length(allv) == 1L) {
    return(allv[1L])
  }
  .zero <- allv == 0
  .n <- sum(!.zero)
  .nzero <- sum(.zero)
  if (.n == 0L) {
    return(0)
  }
  .correction <- (.n - .nzero) / .n
  if (.correction <= 0) {
    .correction <- 1
  }
  (.n / sum(1 / allv[!.zero])) * .correction
}

#' One Shi (2021) central-difference probe at step `h`
#'
#' Paper Algorithm 3.1: the harmonic-mean third-difference ratio plus the +/-h
#' evaluations, or a non-finite marker naming the evaluation that failed.
#' @param h the step size
#' @param ffun the vector-valued function being differenced
#' @param ef expected function error
#' @param t the base evaluation point
#' @param idx the coordinate being perturbed
#' @return list with the ratio r (NA on failure), a `bad` marker, and the +/-h
#'   evaluations
#' @author Hidde van de Beek
#' @noRd
.shiRC <- function(h, ffun, ef, t, idx) {
  .pert <- function(.delta) {
    .tt <- t
    .tt[idx] <- .tt[idx] + .delta
    ffun(.tt)
  }
  .fp1 <- .pert(h)
  if (is.null(.fp1)) {
    return(list(r = NA_real_, bad = "p1"))
  }
  .fm1 <- .pert(-h)
  if (is.null(.fm1)) {
    return(list(r = NA_real_, bad = "m1", fp1 = .fp1))
  }
  .fp3 <- .pert(3 * h)
  if (is.null(.fp3)) {
    return(list(r = NA_real_, bad = "p3", fp1 = .fp1, fm1 = .fm1))
  }
  .fm3 <- .pert(-3 * h)
  if (is.null(.fm3)) {
    return(list(r = NA_real_, bad = "m3", fp1 = .fp1, fm1 = .fm1))
  }
  .allv <- abs(.fp3 - 3 * .fp1 + 3 * .fm1 - .fm3) / (8 * ef)
  list(r = .shiHarmonicMean(.allv), fp1 = .fp1, fm1 = .fm1)
}

#' Shi (2021) adaptive central-difference step and gradient
#'
#' R port of nlmixr2 shi21Central for a vector-valued function `ffun` along
#' coordinate `idx`.  The step is driven down/up by the third-difference ratio
#' until it sits in [rl, ru]; returns the accepted step and the central-difference
#' gradient at it.
#' @param ffun the vector-valued function being differenced
#' @param t the base evaluation point
#' @param idx the coordinate being perturbed
#' @param f0 the base function value
#' @param ef expected function error
#' @param rl,ru acceptance bounds on the third-difference ratio
#' @param nu step scale factor
#' @param maxiter maximum step iterations
#' @return list(h, gr) with the accepted step and the gradient
#' @author Hidde van de Beek
#' @noRd
.shi21CentralR <- function(ffun, t, idx, f0, ef = 7e-7, rl = 1.5, ru = 4.5,
                           nu = 3.0, maxiter = 15L) {
  .h <- (3 * ef)^(1 / 3)
  .l <- 0
  .u <- Inf
  .hlast <- .h
  .gr <- NULL
  for (.iter in seq_len(maxiter)) {
    .shi <- .shiRC(.h, ffun, ef, t, idx)
    if (is.na(.shi$r)) {
      if (.shi$bad == "p1") {
        .h <- .h * 0.5 / 3
        next
      }
      if (.shi$bad == "m1") {
        if (is.null(.gr)) {
          .gr <- (.shi$fp1 - f0) / .h
        }
        .h <- .h * 0.5 / 3
        next
      }
      if (is.null(.gr)) {
        .gr <- (.shi$fp1 - .shi$fm1) / (2 * .h)   # use the .h fp1/fm1 were evaluated at, before shrinking
        .hlast <- .h
      }
      .h <- .h * 2 / 3
      next
    }
    .gr <- (.shi$fp1 - .shi$fm1) / (2 * .h)
    .hlast <- .h
    if (.shi$r < rl) {
      .l <- .h
    } else if (.shi$r > ru) {
      .u <- .h
    } else {
      break
    }
    if (!is.finite(.u)) {
      .h <- nu * .h
    } else if (.l == 0) {
      .h <- .h / nu
    } else {
      .h <- (.l + .u) / 2
    }
  }
  list(h = .hlast, gr = .gr)
}

#' Per-subject observed-information R-matrix over all population parameters
#'
#' Pure assembly (no ODE solve) from already-evaluated sensitivities (E$a/A/Ath,
#' the residual partials in `ef`, the Omega derivatives `omd`).  Each entry is the
#' data term plus the log-determinant term, dispatched over the parameter types
#' (structural theta, sigma, Omega).  The inner Hessian H is over the etas; each
#' structural theta differentiates in its own direction-slot (dirTh maps theta
#' index -> direction), so a fully mu-referenced model has ndir == neta.
#' @param E per-observation sensitivities from the solve (f, a, A, Ath) plus y
#' @param ehat the subject's EBEs
#' @param Om the estimated Omega matrix
#' @param ef the residual error machinery from .foceiAnalyticErrFull
#' @param neta number of etas
#' @param nth number of structural thetas
#' @param nsg number of residual sigma parameters
#' @param sgVar residual sigma variable names
#' @param omd the Omega derivatives from .omegaVarCovDeriv
#' @param covFull FALSE for the structural-theta block only, TRUE for theta +
#'   sigma + Omega
#' @param ndir number of sensitivity directions
#' @param dirTh integer map from structural theta index to direction index
#' @return the per-subject observed-information R-matrix (npE x npE)
#' @author Hidde van de Beek
#' @noRd
.foceiAnalyticSubjectR <- function(E, ehat, Om, ef, neta, nth, nsg, sgVar, omd,
                                   covFull = FALSE, ndir = neta,
                                   dirTh = seq_len(nth)) {
  tr <- function(M) sum(diag(M))
  a <- E$a
  A <- E$A
  Ath <- E$Ath
  f <- E$f
  y <- E$y
  Oi <- solve(Om)
  evf <- function(e) ef$ev(e, f, y)
  rd <- list(r1 = evf(ef$sc$r1), r2 = evf(ef$sc$r2), r3 = evf(ef$sc$r3))
  pf <- list(p = evf(ef$sc$p), p1 = evf(ef$sc$p1), p2 = evf(ef$sc$p2))
  np <- nth + nsg + omd$nom
  ei <- seq_len(neta)
  di <- seq_len(ndir)
  ae <- a[, ei, drop = FALSE]                          # eta-cols for sigma/Omega
  .dirOf <- function(p) dirTh[p]
  H <- Oi
  for (l in ei) {
    for (m in ei) {
      H[l, m] <- H[l, m] + sum(rd$r2 * a[, l] * a[, m] + rd$r1 * A[, l, m])
    }
  }
  N <- matrix(0, neta, ndir)
  for (l in ei) {
    for (d in di) {
      N[l, d] <- sum(rd$r2 * a[, l] * a[, d] + rd$r1 * A[, l, d])
    }
  }
  HiM <- solve(H)
  Ht <- Oi
  for (l in ei) {
    for (m in ei) {
      Ht[l, m] <- Ht[l, m] + sum(pf$p * a[, l] * a[, m])
    }
  }
  Hti <- solve(Ht)
  ouAA <- function(v) {
    M <- matrix(0, neta, neta)
    for (l in ei) {
      for (m in ei) {
        M[l, m] <- sum(v * a[, l] * a[, m])
      }
    }
    M
  }
  dHtD <- lapply(di, function(s) {
    D <- matrix(0, neta, neta)
    for (l in ei) {
      for (m in ei) {
        D[l, m] <- sum(pf$p1 * a[, s] * a[, l] * a[, m] +
                         pf$p * A[, l, s] * a[, m] + pf$p * a[, l] * A[, m, s])
      }
    }
    D
  })
  d2HtDD <- lapply(di, function(s) {
    lapply(di, function(t) {
      D <- matrix(0, neta, neta)
      for (l in ei) {
        for (m in ei) {
          D[l, m] <- sum(pf$p2 * a[, s] * a[, t] * a[, l] * a[, m] +
            pf$p1 * (A[, s, t] * a[, l] * a[, m] + a[, s] * A[, l, t] * a[, m] +
                       a[, s] * a[, l] * A[, m, t] + a[, t] * A[, l, s] * a[, m] +
                       a[, t] * a[, l] * A[, m, s]) +
            pf$p * (Ath[, l, s, t] * a[, m] + A[, l, s] * A[, m, t] +
                      A[, l, t] * A[, m, s] + a[, l] * Ath[, m, s, t]))
        }
      }
      D
    })
  })
  Cen <- vapply(ei, function(l) 0.5 * tr(Hti %*% dHtD[[l]]), numeric(1))
  Cee <- matrix(0, neta, neta)
  for (s in ei) {
    for (t in ei) {
      Cee[s, t] <- 0.5 * (tr(Hti %*% d2HtDD[[s]][[t]]) -
                            tr(Hti %*% dHtD[[s]] %*% Hti %*% dHtD[[t]]))
    }
  }
  Tn <- array(0, c(neta, ndir, ndir))
  for (l in ei) {
    for (s in di) {
      for (t in di) {
        Tn[l, s, t] <- sum(rd$r3 * a[, l] * a[, s] * a[, t] +
                             rd$r2 * (A[, l, s] * a[, t] + A[, l, t] * a[, s] +
                                        A[, s, t] * a[, l]) +
                             rd$r1 * Ath[, l, s, t])
      }
    }
  }
  Ndd <- function(da, db) sum(rd$r2 * a[, da] * a[, db] + rd$r1 * A[, da, db])
  typ <- function(p) {
    if (p <= nth) {
      "th"
    } else if (p <= nth + nsg) {
      "sg"
    } else {
      "om"
    }
  }
  sgi <- function(p) sgVar[p - nth]
  omc <- function(p) p - nth - nsg
  PVper <- function(p) lapply(ef$per[[sgi(p)]], evf)
  PVpair <- function(aa, bb) {
    .s1 <- sgi(aa)
    .s2 <- sgi(bb)
    .k <- if (paste0(.s1, .s2) %in% names(ef$pair)) paste0(.s1, .s2) else paste0(.s2, .s1)
    lapply(ef$pair[[.k]], evf)
  }
  # Omega enters only the prior and Ht's +Omega^-1 term, so every Omega-block
  # quantity is an E-basis contraction from `omd`.
  Mcol <- function(p) {
    .pt <- typ(p)
    if (.pt == "th") {
      return(N[, .dirOf(p)])
    }
    if (.pt == "sg") {
      return(as.numeric(crossprod(ae, PVper(p)$rf)))
    }
    as.numeric(omd$dOi[[omc(p)]] %*% ehat)
  }
  dHt_p <- function(p) {
    .pt <- typ(p)
    if (.pt == "th") {
      return(dHtD[[.dirOf(p)]])
    }
    if (.pt == "sg") {
      return(ouAA(PVper(p)$ps))
    }
    omd$dOi[[omc(p)]]
  }
  d2HtEtaP <- function(p, l) {
    .pt <- typ(p)
    if (.pt == "th") {
      return(d2HtDD[[.dirOf(p)]][[l]])
    }
    if (.pt == "om") {
      return(matrix(0, neta, neta))
    }
    .pv <- PVper(p)
    D <- matrix(0, neta, neta)
    for (s in ei) {
      for (m in ei) {
        D[s, m] <- sum(.pv$pf * a[, l] * a[, s] * a[, m] +
                         .pv$ps * (A[, s, l] * a[, m] + a[, s] * A[, m, l]))
      }
    }
    D
  }
  d2Ht_pp <- function(aa, bb) {
    .ta <- typ(aa)
    .tb <- typ(bb)
    if (.ta == "th" && .tb == "th") {
      return(d2HtDD[[.dirOf(aa)]][[.dirOf(bb)]])
    }
    if (.ta == "om" && .tb == "om") {
      return(omd$d2Oi[[omc(aa)]][[omc(bb)]])
    }
    if (.ta == "om" || .tb == "om") {
      return(matrix(0, neta, neta))
    }
    if (.ta == "sg" && .tb == "sg") {
      return(ouAA(PVpair(aa, bb)$pss))
    }
    .thp <- if (.ta == "th") aa else bb
    .sg <- if (.ta == "th") bb else aa
    d2HtEtaP(.sg, .dirOf(.thp))
  }
  Smat <- function(p) {
    .pt <- typ(p)
    if (.pt == "th") {
      M <- matrix(0, neta, ndir)
      for (l in ei) {
        for (s in di) {
          M[l, s] <- Tn[l, .dirOf(p), s]
        }
      }
      return(M)
    }
    if (.pt == "om") {
      return(omd$dOi[[omc(p)]])
    }
    .pv <- PVper(p)
    M <- matrix(0, neta, ndir)
    for (l in ei) {
      for (s in di) {
        M[l, s] <- sum(.pv$rff * a[, s] * a[, l] + .pv$rf * A[, l, s])
      }
    }
    M
  }
  Svec <- function(aa, bb) {
    .ta <- typ(aa)
    .tb <- typ(bb)
    if (.ta == "th" && .tb == "th") {
      return(Tn[, .dirOf(aa), .dirOf(bb)])
    }
    if (.ta == "om" && .tb == "om") {
      return(as.numeric(omd$d2Oi[[omc(aa)]][[omc(bb)]] %*% ehat))
    }
    if (.ta == "om" || .tb == "om") {
      return(rep(0, neta))
    }
    if (.ta == "sg" && .tb == "sg") {
      return(as.numeric(crossprod(ae, PVpair(aa, bb)$rfss)))
    }
    .thp <- if (.ta == "th") aa else bb
    .sg <- if (.ta == "th") bb else aa
    Smat(.sg)[, .dirOf(.thp)]
  }
  d2Phi <- function(aa, bb) {
    .ta <- typ(aa)
    .tb <- typ(bb)
    if (.ta == "th" && .tb == "th") {
      return(Ndd(.dirOf(aa), .dirOf(bb)))
    }
    if (.ta == "om" && .tb == "om") {
      return(0.5 * as.numeric(t(ehat) %*% omd$d2Oi[[omc(aa)]][[omc(bb)]] %*% ehat) +
               0.5 * omd$d2LD[omc(aa), omc(bb)])
    }
    if (.ta == "om" || .tb == "om") {
      return(0)
    }
    if (.ta == "sg" && .tb == "sg") {
      return(sum(PVpair(aa, bb)$rss))
    }
    .thp <- if (.ta == "th") aa else bb
    .sg <- if (.ta == "th") bb else aa
    as.numeric(crossprod(a[, .dirOf(.thp)], PVper(.sg)$rf))
  }
  # covFull=FALSE (theta-only): the theta-theta block still needs the 3rd-order
  # sensitivities, so only the sigma/Omega columns of the assembly are skipped.
  npE <- if (covFull) np else nth
  etaP <- matrix(vapply(1:npE, function(p) as.numeric(-HiM %*% Mcol(p)), numeric(neta)),
                 nrow = neta)                           # neta x npE (neta == 1 safe)
  eta2 <- function(aa, bb) {
    .b <- Svec(aa, bb) + Smat(aa)[, ei, drop = FALSE] %*% etaP[, bb] +
      Smat(bb)[, ei, drop = FALSE] %*% etaP[, aa]
    for (l in ei) {
      .b[l] <- .b[l] + as.numeric(t(etaP[, aa]) %*% Tn[l, ei, ei] %*% etaP[, bb])
    }
    as.numeric(-HiM %*% .b)
  }
  Cpe <- function(p, l) {
    0.5 * (tr(Hti %*% d2HtEtaP(p, l)) - tr(Hti %*% dHt_p(p) %*% Hti %*% dHtD[[l]]))
  }
  Cpp <- function(aa, bb) {
    0.5 * (tr(Hti %*% d2Ht_pp(aa, bb)) - tr(Hti %*% dHt_p(aa) %*% Hti %*% dHt_p(bb)))
  }
  R <- matrix(0, npE, npE)
  for (aa in 1:npE) {
    for (bb in aa:npE) {                                # R symmetric: fill upper, mirror
      .dat <- d2Phi(aa, bb) - as.numeric(t(Mcol(aa)) %*% HiM %*% Mcol(bb))
      .ld <- Cpp(aa, bb) +
        sum(vapply(ei, function(l) Cpe(aa, l), numeric(1)) * etaP[, bb]) +
        sum(vapply(ei, function(mm) Cpe(bb, mm), numeric(1)) * etaP[, aa]) +
        as.numeric(t(etaP[, aa]) %*% Cee %*% etaP[, bb]) +
        sum(Cen * eta2(aa, bb))
      R[aa, bb] <- R[bb, aa] <- .dat + .ld
    }
  }
  R
}

#' Full analytic FOCEI covariance for a fitted nlmixr2 object
#'
#' Computes the exact observed-information R-matrix over the structural fixed
#' effects, residual sigma, and Omega variances and covariances, and returns the
#' covariance solve(R).  Returns NULL (caller falls back to finite differences)
#' when the model is out of scope or any per-subject augmented solve fails.
#' @param fit a fitted nlmixr2 focei object
#' @param covFull FALSE (default) returns the structural-theta block only; TRUE
#'   returns the full theta + sigma + Omega covariance
#' @return list(cov, se, R, params), with `params` covering theta, sigma, the
#'   `om.<theta>` Omega variances and any `cov.<theta>.<theta>` covariances, or
#'   NULL
#' @author Hidde van de Beek
#' @noRd
.foceiCovAnalytic <- function(fit, covFull = FALSE) {
  .ui <- fit$finalUi
  if (!.hasRxExpandSens2()) {
    return(NULL)                                        # need rxExpandSens2_ + symengine
  }
  # reuse the fit's covariate-interpolation and integration method for the augmented
  # solves; only these scalar options are forwarded, since rxControl's per-compartment
  # tolerance vectors are sized to the original model, not the augmented one.
  .rxc <- tryCatch(rxode2::rxGetControl(.ui, "rxControl", NULL), error = function(e) NULL)
  .solveOpts <- list()
  if (is.list(.rxc)) {
    for (.n in c("covsInterpolation", "method")) {
      if (!is.null(.rxc[[.n]])) .solveOpts[[.n]] <- .rxc[[.n]]
    }
  }
  if (length(.solveOpts) == 0L) .solveOpts <- NULL
  if (isTRUE(as.logical(rxode2::rxGetControl(.ui, "fo", FALSE)))) {
    return(NULL)                                        # FO/FOI: first-order marginal, not conditional
  }
  # Analytic scope is FOCEI only.  FOCE (no interaction) has a non-smooth objective,
  # so it is left to its own finite-difference cov path.
  .interaction <- !identical(as.integer(rxode2::rxGetControl(.ui, "interaction", 1L)), 0L)
  if (!.interaction) {
    return(NULL)                                        # FOCE -> finite-difference fallback
  }
  # Censored (BLOQ, M2/M3/M4) data: the engine scores every DV as an exact Gaussian
  # observation, disagreeing with the FD censored likelihood -- bow out to FD.
  .cens <- fit$dataSav$CENS
  if (!is.null(.cens) && any(.cens != 0, na.rm = TRUE)) {
    return(NULL)
  }
  .ef <- .foceiAnalyticErrFull(.ui)
  if (is.null(.ef)) {
    return(NULL)                                        # unsupported error model -> fallback
  }
  .ini <- .ui$iniDf
  .map <- .foceiEtaThetaMap(.ui)                        # theta <-> eta pairing
  .etaNames <- .map$etaNames
  .neta <- length(.etaNames)
  if (.neta == 0L) {
    return(NULL)
  }
  .thetaForEta <- .map$thetaForEta
  # A non-mu-ref eta (thetaForEta NA) is still a valid random effect: it keeps its
  # own ETA_i_ direction and Omega variance, named by the eta rather than a theta.
  if (any(.iniIsFixed(.ini, .thetaForEta))) {
    return(NULL)                                        # fixed structural theta breaks eta indexing
  }
  .keep <- !.iniIsFixed(.ini, .ef$sgName)
  .ef$sgVar <- .ef$sgVar[.keep]
  .ef$sgName <- .ef$sgName[.keep]                       # drop fixed sigma
  .Om <- fit$omega
  .blocks <- .omegaBlocks(.Om)                          # free Omega lower-triangle
  .pairs <- do.call(rbind, lapply(.blocks, function(.b) {
    do.call(rbind, lapply(seq_along(.b), function(.a) cbind(.b[.a], .b[seq_len(.a)])))
  }))
  .pairs <- .pairs[!.omegaFixed(.ini, .pairs), , drop = FALSE]   # drop fixed Omega elements
  .omd <- .omegaVarCovDeriv(.Om, .pairs)
  .thRows <- .ini[!is.na(.ini$ntheta), , drop = FALSE]
  .thRows <- .thRows[order(.thRows$ntheta), , drop = FALSE]
  # Uniform direction assembly: every estimated structural theta is one direction.
  # A mu-ref theta reuses its eta's direction; any other structural theta gets its
  # own THETA_j_ direction carrying its true sensitivity (nothing is scaled).
  .thStructRows <- .thRows[!.iniIsFixed(.ini, .thRows$name) & !(.thRows$name %in% .ef$sgName), , drop = FALSE]
  .thStruct <- .thStructRows$name
  .nth <- length(.thStruct)
  if (.nth == 0L) {
    return(NULL)
  }
  .etaDirs <- paste0("ETA_", seq_len(.neta), "_")
  .dirTh <- integer(.nth)
  .nonMuTheta <- character(0)
  for (.p in seq_len(.nth)) {
    .k <- match(.thStruct[.p], .thetaForEta)            # eta this theta mu-references (NA if none)
    if (!is.na(.k)) {
      .dirTh[.p] <- .k                                  # reuse eta k's direction (free)
    } else {
      .nonMuTheta <- c(.nonMuTheta, paste0("THETA_", .thStructRows$ntheta[.p], "_"))
      .dirTh[.p] <- .neta + length(.nonMuTheta)         # its own true-sensitivity direction
    }
  }
  .dirs <- c(.etaDirs, .nonMuTheta)
  .ndir <- length(.dirs)
  .th <- setNames(.thRows$est, paste0("THETA_", seq_len(nrow(.thRows)), "_"))
  .ebes <- as.matrix(fit$eta[, .etaNames, drop = FALSE])
  .nsg <- length(.ef$sgVar)
  .np <- .nth + .nsg + .omd$nom
  # The direction-set augmented model carries the 2nd-order prediction chain; the
  # 3rd-order term is recovered by finite-differencing it (.foceiAnalyticSolveDirFD3).
  .am <- .foceiAnalyticAugModelDirs(.ui, .dirs)
  if (is.null(.am) || .am$ndir != .ndir) {
    return(NULL)
  }
  .ds <- fit$dataSav
  # dataSav$ID is etTrans's integer 1..N renumbering; fit$eta$ID is a factor whose
  # labels are the original IDs, so join on the integer factor code, not the label.
  .dsId <- if (is.factor(.ds$ID)) as.integer(.ds$ID) else .ds$ID
  .ids <- as.integer(fit$eta$ID)
  .npOut <- if (covFull) .np else .nth                  # covFull=FALSE: structural thetas only
  .R <- matrix(0, .npOut, .npOut)
  for (.i in seq_along(.ids)) {
    .s <- .ds[.dsId == .ids[.i], , drop = FALSE]        # solve over the subject's actual events
    .obs <- .s[.s$EVID == 0, , drop = FALSE]            # (CMT/EVID/II/SS/ADDL/covariates preserved)
    .p <- c(.th, setNames(.ebes[.i, ], .etaDirs))
    .E <- .foceiAnalyticSolveDirFD3(.am, .p, .s, .obs$TIME, solveOpts = .solveOpts)
    if (is.null(.E)) {
      return(NULL)                                      # solve failure -> finite-difference fallback
    }
    .E$y <- .obs$DV
    .Ri <- tryCatch(.foceiAnalyticSubjectR(.E, .ebes[.i, ], .Om, .ef, .neta, .nth,
                                           .nsg, .ef$sgVar, .omd, covFull = covFull,
                                           ndir = .ndir, dirTh = .dirTh),
                    error = function(e) NULL)
    if (is.null(.Ri) || !all(is.finite(.Ri))) {
      return(NULL)
    }
    .R <- .R + .Ri
  }
  .cov <- tryCatch(solve(.R), error = function(e) NULL)
  if (is.null(.cov)) {
    return(NULL)
  }
  .onm <- ifelse(is.na(.thetaForEta), .etaNames, .thetaForEta)   # non-mu-ref eta named by the eta
  .omNm <- .foceiOmegaCovNames(.pairs, .onm)
  .nm <- if (covFull) c(.thStruct, .ef$sgName, .omNm) else .thStruct
  dimnames(.R) <- dimnames(.cov) <- list(.nm, .nm)
  list(cov = .cov,
       se = setNames(suppressWarnings(sqrt(diag(.cov))), .nm),   # NaN flags non-PD
       R = .R,
       params = .nm)
}
