## FOCEI vs FOCE analytic R at a PINNED common point (removes the convergence confound).
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
ref <- suppressMessages(suppressWarnings(nlmixr(blk, dat, "focei",
         foceiControl(print=0L, covMethod="", interaction=FALSE, sigdig=6))))
th <- ref$theta; em <- as.matrix(ref$eta[,-1]); om <- ref$omega
mdl <- blk |> ini(th)
for (nm in rownames(om)) mdl <- do.call(ini, list(mdl, stats::setNames(list(om[nm,nm]), nm)))
pin <- function(inter) suppressMessages(suppressWarnings(nlmixr(mdl, dat, "focei",
  foceiControl(print=0L, covMethod="", interaction=inter, sigdig=6,
               etaMat=em, maxInnerIterations=0L, maxOuterIterations=0L))))
fI <- pin(TRUE); fF <- pin(FALSE)
cat("pinned thetas identical:", isTRUE(all.equal(unname(fI$theta), unname(fF$theta))), "\n")
cat("pinned etas identical  :", isTRUE(all.equal(as.matrix(fI$eta[,-1]), as.matrix(fF$eta[,-1]))), "\n")
rI <- nlmixr2est:::foceiCovAnalytic(fI); rF <- nlmixr2est:::foceiCovAnalytic(fF)
pn <- c("tka","tcl","tv","add.sd","om.eta.ka","om.eta.cl","cov.eta.v.eta.cl","om.eta.v")
eig <- function(m) paste(sprintf("%.3e", eigen(m, symmetric=TRUE, only.values=TRUE)$values), collapse=" ")
cat("\nFOCEI eigen(R)  :", eig(rI$R[pn,pn]), "\n")
cat("FOCE  eigen(R)  :", eig(rF$R[pn,pn]), "\n")
cat("FOCEI eigen(cov):", eig(rI$cov), "\n")
cat("FOCE  eigen(cov):", eig(rF$cov), "\n")
A <- rI$R[pn,pn]; B <- rF$R[pn,pn]
big <- abs(B) > 0.01*max(abs(B)); rel <- abs(A-B)/pmax(abs(B),1e-12)
cat("\nmax rel |FOCEI R - FOCE R| at the SAME point:", sprintf("%.3e", max(rel[big])), "\n")
w <- which(big & rel > 1e-3, arr.ind=TRUE)
if (nrow(w)) { cat("disagreeing entries (>0.1%) -- FOCEI assembly defect:\n")
  for (i in seq_len(nrow(w))) cat(sprintf("  [%-17s, %-17s] focei=%12.5g foce=%12.5g rel=%.3f\n",
    pn[w[i,1]], pn[w[i,2]], A[w[i,1],w[i,2]], B[w[i,1],w[i,2]], rel[w[i,1],w[i,2]]))
} else cat("FOCEI R == FOCE R at the same point -> earlier 3.8x was the convergence confound\n")
