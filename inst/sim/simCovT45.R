## Phase 4.2 T4 and T5.  Writes t4bal.rds, t4imb.rds, t5zero.rds.
##
## Output directory is $CLAUDE_JOB_DIR/tmp when set, else the working directory.
.outDir <- Sys.getenv("CLAUDE_JOB_DIR")
.outDir <- if (nzchar(.outDir)) file.path(.outDir, "tmp") else getwd()
## Phase 4.2 T4 and T5.
##
## T4  TIME-VARYING covariate.  Two datasets from the SAME subjects and the same
##     underlying etas: one balanced, one with a 4:1 observation-count imbalance.
##     bWT must agree between them.  That is what PROVES the record weighting
##     rather than merely choosing it -- a per-record weighting would let the
##     heavily-sampled subjects dominate and the two answers would part.
##
## T5  DEGENERATE: the covariate effect is truly zero.  Must return ~0, and must
##     not DRIFT -- the estimate at 2n iterations must be no further from 0 than
##     at n.  A slow drift away from zero is the failure a single fit hides.
suppressMessages(library(rxode2))

.lclm <- 1.63; .lv1m <- 1.55; .lclrv <- -1.2; .lv1rv <- -1.2; .rho <- 0.50
u <- function(z) pmin(pmax(stats::pnorm(z), 1e-15), 1 - 1e-15)
tim <- c(0.25, 0.5, 1, 2, 4, 8, 12, 24)

build <- function(bWT, seed, nSub = 120L, tag) {
  set.seed(seed)
  z1 <- stats::rnorm(nSub); z2 <- stats::rnorm(nSub)
  w2 <- .rho * z1 + sqrt(1 - .rho^2) * z2
  ## WT varies WITHIN subject: a baseline plus a per-visit drift, so the
  ## declaration's rate is a different number at every record.
  wtBase <- stats::rnorm(nSub, 70, 12)
  wtRec <- lapply(seq_len(nSub), function(i)
    round(wtBase[i] + cumsum(stats::rnorm(length(tim), 0, 1.5)), 1))
  shCL <- 1/exp(.lclrv); shV1 <- 1/exp(.lv1rv)
  V1 <- stats::qgamma(u(w2), shape = shV1, rate = 1/(exp(.lv1rv)*exp(.lv1m)))
  m <- rxode2({ d/dt(central) <- -(cl/v)*central; cp <- central/v })
  rows <- lapply(seq_len(nSub), function(i) {
    ## CL is a per-RECORD quantity here: same latent z1[i], covariate per visit
    cl_i <- stats::qgamma(u(z1[i]), shape = shCL,
                          rate = 1/(exp(.lclrv)*exp(.lclm + bWT*log(wtRec[[i]]/70))))
    ev <- data.frame(id = i, time = c(0, tim), amt = c(100, rep(NA_real_, length(tim))),
                     evid = c(1L, rep(0L, length(tim))), cmt = 1L,
                     cl = c(cl_i[1], cl_i), v = V1[i])
    s <- rxSolve(m, ev, returnType = "data.frame")
    s <- s[!is.na(s$cp) & s$time > 0, ]
    .w <- wtRec[[i]][seq_len(nrow(s))]
    data.frame(ID = i, TIME = s$time, CP = s$cp,
               AMT = NA_real_, EVID = 0L, CMT = 1L, WT = .w)
  })
  obs <- do.call(rbind, rows)
  ## LLOQ on the TRUE concentration, before the noise -- see simT3.R.  Applied
  ## across the whole dataset so it is one assay limit, not a per-subject one.
  obs <- obs[obs$CP >= 1e-3 * max(obs$CP), ]
  obs$DV <- obs$CP * (1 + stats::rnorm(nrow(obs), 0, 0.10))
  stopifnot(all(obs$DV > 0))
  obs <- obs[, c("ID", "TIME", "DV", "AMT", "EVID", "CMT", "WT")]
  dose <- data.frame(ID = seq_len(nSub), TIME = 0, DV = NA_real_, AMT = 100,
                     EVID = 1L, CMT = 1L, WT = round(wtBase, 1))
  d <- rbind(dose, obs); d <- d[order(d$ID, d$TIME, -d$EVID), ]; rownames(d) <- NULL
  attr(d, "truth") <- c(lclm = .lclm, lv1m = .lv1m, lclrv = .lclrv,
                        lv1rv = .lv1rv, bWT = bWT, rho = .rho)
  saveRDS(d, file.path(.outDir, paste0(tag, ".rds")))
  cat(sprintf("%-12s %d subj, %d rows, bWT=%.2f, WT %.1f-%.1f\n",
              tag, nSub, nrow(d), bWT, min(d$WT), max(d$WT)))
  d
}

## T4 balanced, then the SAME data thinned 4:1 -- same subjects, same etas, so
## any disagreement in bWT is the weighting and nothing else.
d4 <- build(0.75, 20260911L, tag = "t4bal")
set.seed(99)
heavy <- sample(unique(d4$ID), length(unique(d4$ID)) %/% 2)
## Thin by POSITION, after the LLOQ.  Selecting on time values instead left the
## ratio at the mercy of which records the assay limit had already removed --
## it came out 7:1 with one subject on a single sample, which confounds the
## weighting question with a near-uninformative subject.
.lightKeep <- 2L
keep <- unlist(lapply(split(seq_len(nrow(d4)), d4$ID), function(ix) {
  .obs <- ix[d4$EVID[ix] == 0L]
  .dose <- ix[d4$EVID[ix] == 1L]
  if (d4$ID[ix[1]] %in% heavy) c(.dose, .obs)
  else c(.dose, utils::head(.obs, .lightKeep))
}), use.names = FALSE)
keep <- seq_len(nrow(d4)) %in% keep
d4i <- d4[keep, ]; rownames(d4i) <- NULL
attr(d4i, "truth") <- attr(d4, "truth")
saveRDS(d4i, file.path(.outDir, "t4imb.rds"))
.n <- table(d4i$ID[d4i$EVID == 0L])
cat(sprintf("t4imb       %d rows; obs/subject %d-%d (ratio %.1f:1)\n",
            nrow(d4i), min(.n), max(.n), max(.n)/min(.n)))

## T5 degenerate
invisible(build(0.0, 20260912L, tag = "t5zero"))
