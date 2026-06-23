# Finite-difference covariance fallback for nlmixr2 FOCEI/FOCE, in NON-CHOLESKY
# (variance-covariance) Omega parameterization.
#
# nlmixr2's covMethod="r" only yields the structural-theta block (no Omega /
# residual SEs), and its Omega is optimized as a Cholesky factor.  This computes
# the FULL covariance (theta + Omega + residual) as the inverse FD Hessian of the
# FOCEI/FOCE objective over the natural variance-covariance Omega elements, the
# residual sigma, and the structural thetas -- re-optimizing the EBEs at each
# perturbation.  It is the robust fallback for the analytic path and the
# mechanism for Omega/residual SEs on the interpretable variance scale.

#' Lightweight 1st-order sensitivity model (predf and f1_ETA_k = df/deta),
#' cheap to re-solve for the FD covariance.
#' @noRd
.foceiSens1Model <- function(ui) {
  tryCatch({
    s <- ui$foceiEtaS
    st <- rxode2::rxStateOde(s); neta <- s$..maxEta; etav <- paste0("ETA_", seq_len(neta), "_")
    pred <- get("rx_pred_", s)
    Dn <- function(e, v) symengine::D(e, symengine::S(v))
    sn1 <- function(j, p) symengine::S(paste0("rx__sens_", j, "_BY_", p, "__"))
    toRx <- function(.l) rxode2::rxFromSE(.l)
    g1 <- function(p) { e <- Dn(pred, p); for (j in st) e <- e + Dn(pred, j) * sn1(j, p); e }
    fL1 <- vapply(etav, function(p) { .l <- g1(p); paste0("f1_", p, "=", toRx(.l)) }, character(1))
    baseOde <- vapply(st, function(x) { .l <- get(paste0("rx__d_dt_", x, "__"), s)
      paste0("d/dt(", x, ")=", toRx(.l)) }, character(1))
    modTxt <- paste(c(baseOde, s$..sens, paste0("predf=", toRx(pred)), fL1), collapse = "\n")
    modTxt <- gsub("ETA\\[([0-9]+)\\]", "ETA_\\1_", modTxt)
    modTxt <- gsub("THETA\\[([0-9]+)\\]", "THETA_\\1_", modTxt)
    list(mod = rxode2::rxode2(modTxt), etav = etav, neta = neta)
  }, error = function(e) NULL)
}

#' Residual-variance R(f) and dR/df closures from add/prop error values.
#' @noRd
.foceiRfun <- function(errType) {
  function(ev) {
    add <- ev[["add"]]; prop <- ev[["prop"]]
    if (errType == "add")  list(R = function(f) rep(add^2, length(f)),        d = function(f) rep(0, length(f)))
    else if (errType == "prop") list(R = function(f) (prop * f)^2,            d = function(f) 2 * prop^2 * f)
    else list(R = function(f) add^2 + (prop * f)^2,                           d = function(f) 2 * prop^2 * f)  # combined2
  }
}

#' Finite-difference covariance (full theta + Omega + residual, non-Cholesky)
#'
#' Computes the full FOCEI/FOCE covariance as the inverse finite-difference
#' Hessian of the objective over the structural thetas, the residual sigma, and
#' the random-effect `Omega` elements in the natural variance-covariance
#' (non-Cholesky) parameterization, re-optimizing the EBEs at each perturbation.
#' This produces the `Omega` and residual SEs that `covMethod = "r"` does not,
#' and is the robust fallback for the analytic path ([foceiCovAnalytic()]).
#'
#' @param fit a fitted nlmixr2 focei object.
#' @param h relative finite-difference step.
#' @return list with `cov`, `se`, `params`, and the Hessian `R`, or `NULL`.
#' @author Matthew Fidler, Hidde van de Beek
#' @export
foceiCovFD <- function(fit, h = 1e-4) {
  ui <- fit$finalUi
  m1 <- .foceiSens1Model(ui); if (is.null(m1)) return(NULL)
  neta <- m1$neta; etav <- m1$etav; mod <- m1$mod
  inter <- as.integer(rxode2::rxGetControl(ui, "interaction", 1L))
  ini <- ui$iniDf; muRef <- ui$muRefDataFrame
  etaRows <- ini[!is.na(ini$neta1) & ini$neta1 == ini$neta2, , drop = FALSE]
  etaRows <- etaRows[order(etaRows$neta1), , drop = FALSE]
  thetaForEta <- muRef$theta[match(etaRows$name, muRef$eta)]
  if (anyNA(thetaForEta)) return(NULL)
  thRows <- ini[!is.na(ini$ntheta), , drop = FALSE]; thRows <- thRows[order(thRows$ntheta), , drop = FALSE]
  errRows <- ini[!is.na(ini$err), , drop = FALSE]
  hasAdd <- any(errRows$err %in% c("add", "lnorm")); hasProp <- any(errRows$err %in% c("prop", "propT", "propF"))
  errType <- if (hasAdd && hasProp) "combined" else if (hasProp) "prop" else "add"
  errName <- errRows$name; addName <- errRows$name[errRows$err %in% c("add","lnorm")]; propName <- errRows$name[errRows$err %in% c("prop","propT","propF")]
  Rbuild <- .foceiRfun(errType)
  Om0 <- fit$omega; isDiag <- all(abs(Om0[lower.tri(Om0)]) < 1e-10)
  omIdx <- if (isDiag) cbind(1:neta, 1:neta) else which(lower.tri(Om0, diag = TRUE), arr.ind = TRUE)
  omNames <- if (isDiag) paste0("om.", thetaForEta) else apply(omIdx, 1, function(z) paste0("om", z[1], z[2]))
  ebes0 <- as.matrix(fit$eta[, etaRows$name, drop = FALSE]); ids <- fit$eta$ID; ds <- fit$dataSav
  subs <- lapply(ids, function(id) { s <- ds[ds$ID == id, ]; dose <- s[s$EVID != 0, ]; obs <- s[s$EVID == 0, ]
    ev <- rxode2::et(); for (k in seq_len(nrow(dose))) ev <- ev |> rxode2::et(amt = dose$AMT[k], cmt = "depot", time = dose$TIME[k])
    list(ev = ev |> rxode2::et(obs$TIME), times = obs$TIME, y = obs$DV) })

  # full theta value vector (named THETA_k_) with structural + error overrides
  thBase <- setNames(thRows$est, paste0("THETA_", seq_len(nrow(thRows)), "_"))
  unpack <- function(p) {
    nT <- length(thetaForEta); nE <- length(errName)
    th <- thBase
    th[paste0("THETA_", match(thetaForEta, thRows$name), "_")] <- p[seq_len(nT)]
    th[paste0("THETA_", match(errName, thRows$name), "_")] <- p[nT + seq_len(nE)]
    ev <- list(add = if (length(addName)) p[nT + match(addName, errName)] else 0,
               prop = if (length(propName)) p[nT + match(propName, errName)] else 0)
    ov <- p[nT + nE + seq_len(nrow(omIdx))]
    Om <- matrix(0, neta, neta); for (r in seq_len(nrow(omIdx))) { Om[omIdx[r,1],omIdx[r,2]] <- ov[r]; Om[omIdx[r,2],omIdx[r,1]] <- ov[r] }
    list(th = th, ev = ev, Om = Om)
  }
  par0 <- c(thRows$est[match(thetaForEta, thRows$name)], thRows$est[match(errName, thRows$name)],
            if (isDiag) diag(Om0) else Om0[omIdx])
  np <- length(par0)

  Fi <- function(z, sb, eta) {
    Oi <- solve(z$Om); rr <- Rbuild(z$ev)
    sol <- function(e) { d <- as.data.frame(rxode2::rxSolve(mod, params = c(z$th, setNames(e, etav)), sb$ev,
        returnType = "data.frame", atol = 1e-12, rtol = 1e-12))
      d <- d[d$time %in% sb$times, , drop = FALSE]; list(f = d$predf, a = as.matrix(d[, paste0("f1_", etav), drop = FALSE])) }
    for (it in 1:40) { E <- sol(eta); R <- rr$R(E$f); dR <- rr$d(E$f); eps <- sb$y - E$f
      gradCoef <- -eps / R + inter * 0.5 * dR * (R - eps^2) / R^2
      g <- as.numeric(crossprod(E$a, gradCoef)) + as.numeric(Oi %*% eta)
      pw <- 1 / R + inter * 0.5 * (dR / R)^2
      H <- Oi; for (l in 1:neta) for (mm in 1:neta) H[l, mm] <- H[l, mm] + sum(pw * E$a[, l] * E$a[, mm])
      st <- solve(H, g); eta <- eta - st; if (max(abs(st)) < 1e-11) break }
    E <- sol(eta); R <- rr$R(E$f); dR <- rr$d(E$f); eps <- sb$y - E$f
    pw <- 1 / R + inter * 0.5 * (dR / R)^2
    Ht <- Oi; for (l in 1:neta) for (mm in 1:neta) Ht[l, mm] <- Ht[l, mm] + sum(pw * E$a[, l] * E$a[, mm])
    sum(0.5 * (eps^2 / R + log(R))) + 0.5 * as.numeric(t(eta) %*% Oi %*% eta) + 0.5 * log(det(z$Om)) + 0.5 * log(det(Ht))
  }
  Ftot <- function(p) { z <- unpack(p); s <- 0; for (i in seq_along(subs)) s <- s + Fi(z, subs[[i]], ebes0[i, ]); s }
  M <- matrix(0, np, np)
  for (i in 1:np) for (j in i:np) { hi <- h * (abs(par0[i]) + h); hj <- h * (abs(par0[j]) + h)
    pp <- par0; pp[i] <- pp[i]+hi; pp[j] <- pp[j]+hj; pm <- par0; pm[i] <- pm[i]+hi; pm[j] <- pm[j]-hj
    mp <- par0; mp[i] <- mp[i]-hi; mp[j] <- mp[j]+hj; mm <- par0; mm[i] <- mm[i]-hi; mm[j] <- mm[j]-hj
    M[i,j] <- M[j,i] <- (Ftot(pp)-Ftot(pm)-Ftot(mp)+Ftot(mm))/(4*hi*hj) }
  cov <- tryCatch(solve(M), error = function(e) NULL); if (is.null(cov)) return(NULL)
  nm <- c(thetaForEta, errName, omNames); dimnames(cov) <- list(nm, nm)
  list(cov = cov, se = setNames(sqrt(abs(diag(cov))), nm), params = nm, R = M)
}
