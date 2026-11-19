# An ODE model with a subject-level eta (eta.v) and an occasion eta
# (iov.cl ~ 0.1 | occ).  An ODE model, not linCmt(): under LOCF, rxode2's
# linCmt() takes a time-varying covariate's value at the end of each
# interval, so augPred()'s added grid points change its solution where an
# occasion changes.
.iovOdeMod <- function() {
  ini({
    tka <- 0.45
    tcl <- 1
    tv <- 3.45
    add.sd <- 0.7
    eta.v ~ 0.1
    iov.cl ~ 0.1 | occ
  })
  model({
    ka <- exp(tka)
    cl <- exp(tcl + iov.cl)
    v <- exp(tv + eta.v)
    d / dt(depot) <- -ka * depot
    d / dt(center) <- ka * depot - cl / v * center
    cp <- center / v
    cp ~ add(add.sd)
  })
}

# .iovOdeMod() with a second occasion variable, occ2, for v
.iovOdeTwoOccMod <- function() {
  ini({
    tka <- 0.45
    tcl <- 1
    tv <- 3.45
    add.sd <- 0.7
    eta.v ~ 0.1
    iov.cl ~ 0.1 | occ
    iov.v ~ 0.1 | occ2
  })
  model({
    ka <- exp(tka)
    cl <- exp(tcl + iov.cl)
    v <- exp(tv + eta.v + iov.v)
    d / dt(depot) <- -ka * depot
    d / dt(center) <- ka * depot - cl / v * center
    cp <- center / v
    cp ~ add(add.sd)
  })
}

# theo_sd with two occasions; "subject" puts the boundary at 3 h for even
# subjects and 8 h for odd ones, "same" at 6 h for all.  occ2 changes at 4 h.
.iovOdeData <- function(occ = c("same", "subject")) {
  occ <- match.arg(occ)
  .d <- nlmixr2data::theo_sd
  .cut <- if (occ == "same") 6 else ifelse(.d$ID %% 2 == 0, 3, 8)
  .d$occ <- ifelse(.d$TIME < .cut, 1L, 2L)
  .d$occ2 <- ifelse(.d$TIME < 4, 1L, 2L)
  .d
}

# .iovOdeData() with its subjects relabelled or reordered: "shifted" labels
# them 101-112, "character" S01-S12, "descending" puts them in descending
# order, and "dose-only" keeps only subject 3's dose, so the fit drops it and
# the later subjects' codes shift by one
.iovOdeIdData <- function(layout = c("shifted", "character", "descending", "dose-only")) {
  layout <- match.arg(layout)
  .d <- .iovOdeData()
  if (layout == "shifted") {
    .d$ID <- .d$ID + 100L
  } else if (layout == "character") {
    .d$ID <- sprintf("S%02d", .d$ID)
  } else if (layout == "descending") {
    .d <- .d[order(-.d$ID, seq_len(nrow(.d))), ]
  } else {
    .d <- .d[!(.d$ID == 3 & .d$EVID == 0), ]
  }
  rownames(.d) <- NULL
  .d
}

# theo_sd with a second dose at each subject's observation nearest 12 h, at
# the same time, so the occasion changes between two tied records: the
# observation is occasion 1 and the dose occasion 2.  `order` is the data's
# row order of the tied pair.
.iovOdeTieData <- function(order = c("obs-dose", "dose-obs")) {
  order <- match.arg(order)
  .d <- nlmixr2data::theo_sd
  .d$occ <- 1L
  .rows <- list()
  for (.id in unique(.d$ID)) {
    .s <- .d[.d$ID == .id, ]
    .obs <- which(.s$EVID == 0)
    .t <- .s$TIME[.obs[which.min(abs(.s$TIME[.obs] - 12))]]
    .s$occ[.s$TIME > .t] <- 2L
    .dose <- .s[1, ]
    .dose$TIME <- .t
    .dose$occ <- 2L
    .at <- which(.s$TIME == .t)
    .before <- if (order == "obs-dose") .at else .at - 1L
    .rows[[length(.rows) + 1L]] <- rbind(.s[seq_len(.before), ], .dose, .s[-seq_len(.before), ])
  }
  .d <- do.call(rbind, .rows)
  rownames(.d) <- NULL
  .d
}

# theo_sd with its dose repeated once 12 h later (addl = 1, ii = 12) and the
# occasion changing at 12 h, so the dose rxode2 adds starts occasion 2
.iovOdeAddlData <- function() {
  .d <- nlmixr2data::theo_sd
  .d$ADDL <- ifelse(.d$EVID != 0, 1L, 0L)
  .d$II <- ifelse(.d$EVID != 0, 12, 0)
  .d$occ <- ifelse(.d$TIME < 12, 1L, 2L)
  .d
}

# A quiet fit with no covariance step.  The solver tolerances are tight, so
# a solve of the fit differs from its table by far less than an occasion
# eta does; saem runs few iterations, which leaves its occasion variance
# large enough to see.
.iovOdeFit <- function(
  est,
  occ = "same",
  data = .iovOdeData(occ),
  model = .iovOdeMod,
  interpolation = "locf"
) {
  .rx <- rxode2::rxControl(atol = 1e-10, rtol = 1e-10, covsInterpolation = interpolation)
  .ctl <- if (est == "saem") {
    saemControl(print = 0, covMethod = "", nBurn = 20, nEm = 20, rxControl = .rx)
  } else {
    foceiControl(print = 0, covMethod = "", rxControl = .rx)
  }
  suppressWarnings(suppressMessages(nlmixr2(model, data, est, control = .ctl)))
}
