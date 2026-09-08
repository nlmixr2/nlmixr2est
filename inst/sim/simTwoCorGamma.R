## Simulate a TWO-CORRELATION declared-distribution dataset.
##
## Bauer's four datasets all have ONE copula correlation, so the multi-
## correlation path -- the one that must never collapse to the single-parameter
## uniroot -- has never been exercised.  Two INDEPENDENT correlated pairs give
## exactly two correlations, which is the smallest case that distinguishes
## "optimize one scalar" from "optimize a system".
.libPaths(c("~/src/.rx-lib", "~/src/.diag-lib", "~/src/.etadist-lib4", .libPaths()))
suppressMessages(library(rxode2))
set.seed(20260908)

N   <- 300                      # subjects, matching Bauer's
tim <- c(0.25,0.5,1,2,4,6,8,12,16,24)

## truth -- gamma BSV on all four structural parameters, in two correlated
## pairs.  Means chosen near Bauer's so the fits are comparable in scale.
CLm <- 5.0;  CLrv <- 0.25       # pair 1
V1m <- 4.7;  V1rv <- 0.25
Qm  <- 2.5;  Qrv  <- 0.20       # pair 2
V2m <- 60 ;  V2rv <- 0.20
rho1 <- 0.6                     # eta.cl  ~ eta.v1
rho2 <- -0.4                    # eta.q   ~ eta.v2   (opposite sign on purpose)

## latent standard normals, correlated WITHIN each pair only
z <- matrix(rnorm(N*4), N, 4)
L1 <- chol(matrix(c(1,rho1,rho1,1),2,2))
L2 <- chol(matrix(c(1,rho2,rho2,1),2,2))
w  <- cbind(z[,1:2] %*% L1, z[,3:4] %*% L2)

## Bauer's construction: eta = Q(Phi(z); args), args from mean/relative variance
qg <- function(u, m, rv) qgamma(u, shape = 1/rv, rate = 1/(rv*m))
u  <- pnorm(w)
eta <- cbind(qg(u[,1], CLm, CLrv), qg(u[,2], V1m, V1rv),
             qg(u[,3], Qm , Qrv ), qg(u[,4], V2m, V2rv))
colnames(eta) <- c("cl","v","q","v2")

mod <- rxode2({
  d/dt(central)  = -(cl/v)*central - (q/v)*central + (q/v2)*periph
  d/dt(periph)   =  (q/v)*central - (q/v2)*periph
  cp = central/v
})
ev <- et(amt = 100, cmt = "central")
ev <- et(ev, tim)
sim <- rxSolve(mod, ev, params = as.data.frame(eta), returnType = "data.frame")
names(sim) <- tolower(names(sim))
cat("sim cols:", paste(names(sim), collapse=","), " rows:", nrow(sim), "\n")
sim$DV <- sim$cp * (1 + rnorm(nrow(sim), 0, 0.15))    # 15% proportional
d <- data.frame(ID = sim$sim.id, TIME = sim$time, AMT = 0, RATE = 0,
                EVID = 0, MDV = 0, DV = sim$DV, IPRED = sim$cp)
dose <- data.frame(ID = unique(d$ID), TIME = 0, AMT = 100, RATE = 0,
                   EVID = 1, MDV = 1, DV = 0, IPRED = 0)
out <- rbind(dose, d); out <- out[order(out$ID, out$TIME, -out$EVID), ]
f <- "~/src/gamma_indpar/gamma2cor_clv1.dat"
writeLines(sprintf(
  "; two-correlation declared-distribution set: CL/V1 rho=%.2f, Q/V2 rho=%.2f; CL=%.1f V1=%.1f Q=%.1f V2=%.0f; rv=%.2f/%.2f",
  rho1, rho2, CLm, V1m, Qm, V2m, CLrv, Qrv), f)
suppressWarnings(write.table(out, f, append = TRUE, row.names = FALSE,
                             quote = FALSE, sep = " "))
cat("wrote", f, "\n  subjects:", length(unique(out$ID)),
    " rows:", nrow(out), "\n")
cat("  empirical cor(eta.cl,eta.v1) =", round(cor(eta[,1], eta[,2]), 3),
    "  cor(eta.q,eta.v2) =", round(cor(eta[,3], eta[,4]), 3), "\n")
