# Shared simulated arm for the declared-distribution covariate tests.
#
# A helper rather than a definition in one test file: testthat sources each test
# file in its own frame, so a generator defined in test-etaDistCovariateFit.R is
# invisible to test-etaDistCovSearch.R -- which is exactly how the search test
# first failed ("could not find function \".edT5Data\"").  Both files need the
# same arm, and two copies would drift.
#
#
# `timeVarying` picks which of the two arms is generated, and they behave very
# differently: with the covariate constant within a subject saem recovers the
# coefficient, and with it varying within a subject saem freezes the whole
# declaration.  Both arms carry essentially the SAME information about the
# coefficient -- measured on 120 subjects, log(WT/70) has a between-subject sd
# of 0.185 (constant arm) against 0.174 (time-varying arm), while the
# time-varying arm's WITHIN-subject sd is only 0.022 -- so the difference
# between the two is structural, not statistical.
.edT5Data <- function(nSub = 60L, bWT = 0.0, seed = 20260912L,
                      timeVarying = TRUE) {
  set.seed(seed)
  .lclm <- 1.63; .lv1m <- 1.55; .lclrv <- -2.4; .lv1rv <- -2.4; .rho <- 0.5
  .tim <- c(0.25, 0.5, 1, 2, 4, 8, 12, 24)
  .u <- function(z) pmin(pmax(stats::pnorm(z), 1e-15), 1 - 1e-15)
  .z1 <- stats::rnorm(nSub); .z2 <- stats::rnorm(nSub)
  .w2 <- .rho * .z1 + sqrt(1 - .rho^2) * .z2
  .wtBase <- stats::rnorm(nSub, 70, 12)
  ## the random walk is drawn either way so the two arms consume the SAME
  ## stream and differ ONLY in whether the covariate moves within a subject
  .wtWalk <- lapply(seq_len(nSub), function(i)
    round(.wtBase[i] + cumsum(stats::rnorm(length(.tim), 0, 1.5)), 1))
  .wtRec <- if (timeVarying) .wtWalk else {
    lapply(seq_len(nSub), function(i) rep(round(.wtBase[i], 1), length(.tim)))
  }
  .shCL <- 1 / exp(.lclrv); .shV1 <- 1 / exp(.lv1rv)
  .V1 <- stats::qgamma(.u(.w2), shape = .shV1,
                       rate = 1 / (exp(.lv1rv) * exp(.lv1m)))
  .m <- rxode2::rxode2({
    d/dt(central) <- -(cl / v) * central
    cp <- central / v
  })
  .rows <- lapply(seq_len(nSub), function(i) {
    .cl <- stats::qgamma(.u(.z1[i]), shape = .shCL,
                         rate = 1 / (exp(.lclrv) *
                                     exp(.lclm + bWT * log(.wtRec[[i]] / 70))))
    .ev <- data.frame(id = i, time = c(0, .tim),
                      amt = c(100, rep(NA_real_, length(.tim))),
                      evid = c(1L, rep(0L, length(.tim))), cmt = 1L,
                      cl = c(.cl[1], .cl), v = .V1[i])
    .s <- rxode2::rxSolve(.m, .ev, returnType = "data.frame")
    .s <- .s[!is.na(.s$cp) & .s$time > 0, ]
    data.frame(ID = i, TIME = .s$time, CP = .s$cp, AMT = NA_real_,
               EVID = 0L, CMT = 1L,
               WT = .wtRec[[i]][seq_len(nrow(.s))])
  })
  .obs <- do.call(rbind, .rows)
  # assay limit on the TRUE concentration: without it a 1e-12 observation
  # against a 1e-5 prediction is a relative residual of 1e7 under prop(), and
  # those records dominate the residual parameter
  .obs <- .obs[.obs$CP > 0.01, ]
  .obs$DV <- .obs$CP * (1 + stats::rnorm(nrow(.obs), 0, 0.10))
  .obs <- .obs[.obs$DV > 0, c("ID", "TIME", "DV", "AMT", "EVID", "CMT", "WT")]
  .dose <- data.frame(ID = seq_len(nSub), TIME = 0, DV = NA_real_, AMT = 100,
                      EVID = 1L, CMT = 1L, WT = round(.wtBase, 1))
  .d <- rbind(.dose, .obs)
  .d[order(.d$ID, .d$TIME, -.d$EVID), ]
}
