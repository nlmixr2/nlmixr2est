## Phase 4.2 T3: a declared gamma whose RATE carries a known covariate effect.
##
## Writes t3cov.rds.  The pairing is the point: ONE declaration carries the
## covariate (so the family M-step must stand down for it) and the other does
## not (so the M-step must still run for that one).  A single declaration
## cannot test that.
##
## Output directory is $CLAUDE_JOB_DIR/tmp when set, else the working directory.
.outDir <- Sys.getenv("CLAUDE_JOB_DIR")
.outDir <- if (nzchar(.outDir)) file.path(.outDir, "tmp") else getwd()
## Phase 4.2 / T3: simulate a declared gamma whose RATE carries a known
## covariate effect, so a fit can be asked to recover it.
##
## The point of the dataset is engagement, not power: one declaration carries
## the covariate (so the family M-step must stand down for it) and the other
## does not (so the M-step must still run for that one).  That pairing is what
## the per-declaration classification claims to handle and what a single
## declaration could not test.
suppressMessages(library(rxode2))
set.seed(20260910)

nSub <- 120L
## Relative variance 0.09 (about 30% CV), which is Bauer's g1 arm -- a spread
## the 0.25-24 h schedule can actually observe.
##
## An earlier version used rv = 0.30 (55% CV) on BOTH cl and v.  That makes the
## elimination rate cl/v span 0.03-57 /h, i.e. half-lives from 23 h to 45 s, and
## no single sampling schedule covers it: fast subjects decay past the solver's
## absolute tolerance (central goes NEGATIVE, measured -2.7e-10 at 0.5 h) and
## slow ones never leave the peak.  Both repairs tried on that design failed --
## an LLOQ at 1e-3 of max removed 40% of records, all of them the low late ones
## that identify clearance, and biased lclm 1.63 -> 2.63; keeping every positive
## record instead drove it to 9.44, because under prop() a 1e-12 observation
## against a 1e-5 prediction is a relative residual of 1e7 and those records
## dominate.  The design was the problem, not the filter.
.lclm <- 1.63; .lv1m <- 1.55; .lclrv <- -2.4; .lv1rv <- -2.4
.bWT  <- 0.75                      # TRUTH: allometric-style exponent on WT/70
.rho  <- 0.50
WT <- round(stats::rnorm(nSub, 70, 12), 1)

## Bauer's construction, written out directly so the simulation does not depend
## on the machinery under test: z ~ N(0,1), the copula pairs them, and the
## gamma quantile decodes each.
z1 <- stats::rnorm(nSub); z2 <- stats::rnorm(nSub)
w2 <- .rho * z1 + sqrt(1 - .rho^2) * z2
u  <- function(z) pmin(pmax(stats::pnorm(z), 1e-15), 1 - 1e-15)
shCL <- 1/exp(.lclrv); rtCL <- 1/(exp(.lclrv) * exp(.lclm + .bWT * log(WT/70)))
shV1 <- 1/exp(.lv1rv); rtV1 <- 1/(exp(.lv1rv) * exp(.lv1m))
CL <- stats::qgamma(u(z1), shape = shCL, rate = rtCL)
V1 <- stats::qgamma(u(w2), shape = shV1, rate = rtV1)

m <- rxode2({ d/dt(central) <- -(cl/v)*central; cp <- central/v })
tim <- c(0.25, 0.5, 1, 2, 4, 8, 12, 24)
ev <- do.call(rbind, lapply(seq_len(nSub), function(i) {
  data.frame(id = i, time = c(0, tim), amt = c(100, rep(NA_real_, length(tim))),
             evid = c(1L, rep(0L, length(tim))), cmt = 1L,
             cl = CL[i], v = V1[i])
}))
s <- rxSolve(m, ev, returnType = "data.frame")
s <- s[!is.na(s$cp) & s$time > 0, ]
## Assay limit.  Applied to the TRUE concentration rather than the measured
## one, so this is a limit of quantification and not a selection on the noise
## draw (which would bias the retained values upward at the limit).
##
## This matters more than it looks.  Keeping every positive record makes the
## declared M-step run AWAY from truth: started exactly at truth it drifted
## lclm 1.63 -> 2.12 in 5 iterations and -> 5.55 in 20, with the step firing
## every iteration.  With the limit below it sits at 1.641 after 5 and 1.641
## after 20 -- truth is a fixed point.  The M-step fits the family to the EBEs,
## and a subject whose only records are ~1e-12 has an essentially undetermined
## EBE; a plain log-normal fit of the same data tolerates them (it recovers
## CL 4.93 / V 4.77) which is why this took a controlled experiment to find.
.LLOQ <- 0.01
s <- s[s$cp > .LLOQ, ]
s$DV <- s$cp * (1 + stats::rnorm(nrow(s), 0, 0.10))     # 10% proportional
stopifnot(all(s$DV > 0))

d <- rbind(
  data.frame(ID = seq_len(nSub), TIME = 0, DV = NA_real_, AMT = 100,
             EVID = 1L, CMT = 1L, WT = WT),
  data.frame(ID = s$id, TIME = s$time, DV = s$DV, AMT = NA_real_,
             EVID = 0L, CMT = 1L, WT = WT[s$id]))
d <- d[order(d$ID, d$TIME, -d$EVID), ]
rownames(d) <- NULL
attr(d, "truth") <- c(lclm = .lclm, lv1m = .lv1m, lclrv = .lclrv,
                      lv1rv = .lv1rv, bWT = .bWT, rho = .rho)
saveRDS(d, file.path(.outDir, "t3cov.rds"))
cat(sprintf("T3 simulated: %d subjects, %d rows; CL %.2f-%.2f, V1 %.2f-%.2f\n",
            nSub, nrow(d), min(CL), max(CL), min(V1), max(V1)))
cat(sprintf("  truth bWT=%.2f rho=%.2f; WT %.1f-%.1f\n", .bWT, .rho, min(WT), max(WT)))
