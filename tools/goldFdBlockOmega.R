## Phase 8a: gold brute-force FD for a BLOCK Omega under FOCE.
## The existing gold tests perturb only diag(Om) -- no off-diagonal element has ever
## been validated.  Compare the analytic R against an independent FD Hessian.
suppressMessages(devtools::load_all("/home/matt-fidler/src/nlmixr2est", quiet=TRUE))
blk <- function() {
  ini({ tka <- log(1.5); tcl <- log(2.7); tv <- log(31.5)
        eta.ka ~ 0.6; eta.cl + eta.v ~ c(0.3, 0.03, 0.1); add.sd <- 0.7 })
  model({ ka <- exp(tka+eta.ka); cl <- exp(tcl+eta.cl); v <- exp(tv+eta.v)
          d/dt(depot) <- -ka*depot; d/dt(center) <- ka*depot - cl/v*center
          cp <- center/v; cp ~ add(add.sd) })
}
d0 <- nlmixr2data::theo_sd
dat <- do.call(rbind, lapply(1:4, function(k){ .x<-d0; .x$ID <- .x$ID+(k-1)*100; .x}))
fitF <- suppressMessages(suppressWarnings(nlmixr(blk, dat, "focei",
          foceiControl(print=0L, covMethod="", interaction=FALSE, sigdig=6))))
rF <- nlmixr2est:::foceiCovAnalytic(fitF)
fitI <- suppressMessages(suppressWarnings(nlmixr(blk, dat, "focei",
          foceiControl(print=0L, covMethod="", sigdig=6))))
rI <- nlmixr2est:::foceiCovAnalytic(fitI)
cat("analytic R params:", paste(colnames(rF$R), collapse=" "), "\n")

ui <- fitF$finalUi; neta <- 3L; etav <- paste0("ETA_", 1:neta, "_")
am <- nlmixr2est:::.foceiAnalyticAugModelDirs(ui, etav)
thNames <- names(fitF$theta)
thBase <- setNames(as.numeric(fitF$theta[thNames]), paste0("THETA_", seq_along(thNames), "_"))
iTh <- match(c("tka","tcl","tv","add.sd"), thNames)
Om0 <- fitF$omega
byId <- split(fitF$dataSav, as.character(fitF$dataSav$ID))
ids <- fitF$eta$ID; idCode <- as.integer(ids)
eta0m <- as.matrix(fitF$eta[, c("eta.ka","eta.cl","eta.v")])
subj <- lapply(seq_along(ids), function(i) {
  s <- byId[[as.character(idCode[i])]]; obs <- s[s$EVID==0, , drop=FALSE]
  list(s=s, times=obs$TIME, y=obs$DV, eta0=eta0m[i,]) })
.fa <- function(th, eta, s, times)
  nlmixr2est:::.foceiAnalyticSolveFA(am, c(th, setNames(eta, etav)), s, times, tol=1e-12)

mkOm <- function(p) { Om <- matrix(0,3,3); Om[1,1]<-p[1]; Om[2,2]<-p[2]
                      Om[3,2]<-Om[2,3]<-p[3]; Om[3,3]<-p[4]; Om }
objFOCE <- function(psi) {
  th <- thBase; th[iTh] <- psi[1:4]; sa <- psi[4]
  Om <- mkOm(psi[5:8]); Oi <- tryCatch(solve(Om), error=function(e) NULL)
  if (is.null(Oi)) return(NA_real_)
  dOm <- det(Om); if (!is.finite(dOm) || dOm <= 0) return(NA_real_)
  ldOm <- log(dOm); tot <- 0
  for (sj in subj) {
    y <- sj$y; s <- sj$s; times <- sj$times
    E0 <- .fa(th, rep(0,neta), s, times); if (is.null(E0)) return(NA_real_)
    R0 <- rep(sa^2, length(y))
    eta <- sj$eta0
    for (it in 1:100) {
      E <- .fa(th, eta, s, times); if (is.null(E)) return(NA_real_)
      q0 <- -(y - E$f)/R0
      S <- as.numeric(Oi %*% eta); for (l in 1:neta) S[l] <- S[l] + sum(q0*E$a[,l])
      if (max(abs(S)) < 1e-12) break
      Hf <- Oi; for (l in 1:neta) for (m in 1:neta)
        Hf[l,m] <- Hf[l,m] + sum((1/R0)*E$a[,l]*E$a[,m] + q0*E$A[,l,m])
      eta <- eta - solve(Hf, S)
    }
    E <- .fa(th, eta, s, times); f <- E$f
    Phi <- 0.5*sum((y-f)^2/R0 + log(R0)) + 0.5*as.numeric(t(eta)%*%Oi%*%eta) + 0.5*ldOm
    Ht <- Oi; for (l in 1:neta) for (m in 1:neta)
      Ht[l,m] <- Ht[l,m] + sum((1/R0)*E$a[,l]*E$a[,m])
    tot <- tot + Phi + 0.5*log(det(Ht))
  }
  tot
}
psi0 <- c(as.numeric(fitF$theta[c("tka","tcl","tv","add.sd")]),
          Om0[1,1], Om0[2,2], Om0[3,2], Om0[3,3])
np <- length(psi0); h <- pmax(abs(psi0), 0.5)*5e-5
H <- matrix(0,np,np); f0 <- objFOCE(psi0)
cat("gold f0 =", sprintf("%.8f", f0), "\n")
for (i in 1:np) { ei <- numeric(np); ei[i] <- h[i]
  H[i,i] <- (objFOCE(psi0+2*ei) - 2*f0 + objFOCE(psi0-2*ei))/(4*h[i]^2) }
for (i in 1:(np-1)) for (j in (i+1):np) {
  ei <- numeric(np); ei[i]<-h[i]; ej <- numeric(np); ej[j]<-h[j]
  H[i,j] <- H[j,i] <- (objFOCE(psi0+ei+ej)-objFOCE(psi0+ei-ej)-
                       objFOCE(psi0-ei+ej)+objFOCE(psi0-ei-ej))/(4*h[i]*h[j]) }
pn <- c("tka","tcl","tv","add.sd","om.eta.ka","om.eta.cl","cov.eta.v.eta.cl","om.eta.v")
dimnames(H) <- list(pn,pn)
Ran <- rF$R[pn,pn]
big <- abs(H) > 0.01*max(abs(H))
rel <- abs(Ran-H)/pmax(abs(H),1e-12)
cat("\nmax rel err on significant entries:", sprintf("%.3e", max(rel[big])), "\n")
w <- which(big & rel > 3e-3, arr.ind=TRUE)
if (nrow(w)) { cat("entries where analytic disagrees with GOLD FD (>0.3%):\n")
  for (i in seq_len(nrow(w))) cat(sprintf("  [%-17s, %-17s] analytic=%12.5g gold=%12.5g rel=%.3f\n",
    pn[w[i,1]], pn[w[i,2]], Ran[w[i,1],w[i,2]], H[w[i,1],w[i,2]], rel[w[i,1],w[i,2]]))
} else cat("analytic R MATCHES the gold FD on all significant entries\n")
## same gold H is the reference for FOCEI too (additive => identical objectives)
RanI <- rI$R[pn,pn]
relI <- abs(RanI-H)/pmax(abs(H),1e-12)
cat("\n=== FOCEI vs the SAME gold FD ===\n")
cat("max rel err on significant entries:", sprintf("%.3e", max(relI[big])), "\n")
wI <- which(big & relI > 3e-3, arr.ind=TRUE)
if (nrow(wI)) { cat("FOCEI entries disagreeing with gold (>0.3%):\n")
  for (i in seq_len(nrow(wI))) cat(sprintf("  [%-17s, %-17s] focei=%12.5g gold=%12.5g rel=%.3f\n",
    pn[wI[i,1]], pn[wI[i,2]], RanI[wI[i,1],wI[i,2]], H[wI[i,1],wI[i,2]], relI[wI[i,1],wI[i,2]]))
} else cat("FOCEI matches gold on all significant entries\n")
cat("\ngold eigen:", paste(sprintf("%.3e", eigen(H, symmetric=TRUE, only.values=TRUE)$values), collapse=" "), "\n")
cat("anlyt eigen:", paste(sprintf("%.3e", eigen(Ran, symmetric=TRUE, only.values=TRUE)$values), collapse=" "), "\n")
