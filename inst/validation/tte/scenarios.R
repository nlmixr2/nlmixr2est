# Time-to-event validation scenarios (nlmixr2/nlmixr2est#1461).
#
# Each scenario holds the nlmixr2 model (with the TRUE values in ini()), the
# starting values used by every tool, the simulation design and the columns of
# its data set.  The data are simulated by simulate.R; the NONMEM control
# streams and Monolix projects in nonmem/ and monolix/ fit the same models.
#
# The random-effect variances start at 0.3 (0.2 for the PK ones): FOCEi bounds
# a variance to 5 times its initial value (foceiControl(diagOmegaBoundUpper=)),
# and NONMEM and Monolix have no such bound.

# Weibull hazard in calendar time, with the scale on the log scale
.tteWeibull <- function() {
  ini({
    lsc <- log(10)
    lsh <- log(1.5)
    eta.sc ~ 0.6
  })
  model({
    sc <- exp(lsc + eta.sc)
    sh <- exp(lsh)
    haz <- sh / sc * (time / sc)^(sh - 1)
    ev ~ tte(haz, maxEventNumber = 1)
  })
}

.tteScenarioList <- function() {
  list(
    exponential = list(
      description = "constant hazard, single event, right-censored at 20",
      model = function() {
        ini({
          lh <- log(0.05)
        })
        model({
          haz <- exp(lh)
          ev ~ tte(haz, maxEventNumber = 1)
        })
      },
      start = c(lh = log(0.1)),
      n = 200L,
      seed = 14611L,
      design = function(n) .tteWindow(n, 0, 20)
    ),
    weibull = list(
      description = "Weibull hazard, random effect on the scale, single event",
      model = .tteWeibull,
      start = c(lsc = log(8), lsh = log(1.2), eta.sc = 0.3),
      n = 400L,
      seed = 14612L,
      design = function(n) .tteWindow(n, 0, 30)
    ),
    rtte = list(
      description = "repeated events, constant hazard with a subject frailty",
      model = function() {
        ini({
          lh <- log(0.15)
          eta.h ~ 0.4
        })
        model({
          haz <- exp(lh + eta.h)
          ev ~ tte(haz)
        })
      },
      start = c(lh = log(0.1), eta.h = 0.3),
      n = 150L,
      seed = 14613L,
      design = function(n) .tteWindow(n, 0, 24)
    ),
    gapTime = list(
      description = "repeated events, Weibull hazard in the time since the last event",
      model = function() {
        ini({
          lte <- log(5)
          lsh <- log(2)
          eta.te ~ 0.2
        })
        model({
          te <- exp(lte + eta.te)
          sh <- exp(lsh)
          haz <- sh / te * ((time - tteTime(ev)) / te)^(sh - 1)
          ev ~ tte(haz)
        })
      },
      start = c(lte = log(4), lsh = log(1.5), eta.te = 0.3),
      n = 150L,
      seed = 14614L,
      design = function(n) .tteWindow(n, 0, 40),
      # NONMEM reads the hazard clock from a data item: the last event before
      # each record (NONMEM data items in $DES come from the record ending the
      # interval)
      addColumns = function(d) {
        d$TLAST <- .tteLastEventBefore(d)
        d
      }
    ),
    pkFixed = list(
      description = "hazard inhibited by a one-compartment PK model with known individual PK",
      model = function() {
        ini({
          llam <- log(0.1)
          lic50 <- log(2)
          eta.lam ~ 0.3
        })
        model({
          d / dt(central) <- -CLI / VI * central
          cp <- central / VI
          haz <- exp(llam + eta.lam) * (1 - cp / (exp(lic50) + cp))
          ev ~ tte(haz)
        })
      },
      start = c(llam = log(0.05), lic50 = log(1), eta.lam = 0.3),
      n = 150L,
      seed = 14615L,
      design = function(n) {
        .d <- .tteWindow(n, 0, 48)
        .pk <- data.frame(
          id = seq_len(n),
          CLI = exp(stats::rnorm(n, 0, sqrt(0.1))),
          VI = 10 * exp(stats::rnorm(n, 0, sqrt(0.1)))
        )
        .tteAddDose(merge(.d, .pk, by = "id"), 100)
      },
      covariates = c("CLI", "VI")
    ),
    pkJoint = list(
      description = "PK endpoint estimated jointly with the PK-inhibited hazard",
      model = function() {
        ini({
          tcl <- log(1)
          tv <- log(10)
          prop.sd <- 0.1
          llam <- log(0.1)
          lic50 <- log(2)
          eta.cl ~ 0.1
          eta.v ~ 0.1
          eta.lam ~ 0.3
        })
        model({
          cl <- exp(tcl + eta.cl)
          v <- exp(tv + eta.v)
          d / dt(central) <- -cl / v * central
          cp <- central / v
          haz <- exp(llam + eta.lam) * (1 - cp / (exp(lic50) + cp))
          cp ~ prop(prop.sd)
          ev ~ tte(haz)
        })
      },
      start = c(
        tcl = log(1.5),
        tv = log(8),
        prop.sd = 0.2,
        llam = log(0.05),
        lic50 = log(1),
        eta.cl = 0.2,
        eta.v = 0.2,
        eta.lam = 0.3
      ),
      n = 150L,
      seed = 14616L,
      design = function(n) {
        .d <- rbind(
          .tteWindow(n, 0, 48),
          data.frame(
            id = rep(seq_len(n), each = 6),
            time = c(0.5, 2, 4, 8, 12, 24),
            evid = 0L,
            amt = 0,
            cmt = "cp"
          )
        )
        .tteAddDose(.d, 100)
      }
    ),
    tvCov = list(
      description = "repeated events with a covariate that switches on during follow-up",
      model = function() {
        ini({
          lh <- log(0.1)
          bx <- -0.7
          eta.h ~ 0.3
        })
        model({
          haz <- exp(lh + eta.h + bx * X)
          ev ~ tte(haz)
        })
      },
      # nested model for the objective function difference
      reduced = function(model) {
        rxode2::ini(model, bx = fix(0))
      },
      start = c(lh = log(0.2), bx = 0, eta.h = 0.3),
      n = 150L,
      seed = 14617L,
      design = function(n) {
        .d <- .tteWindow(n, 0, 30)
        .d$X <- 0
        .switch <- data.frame(
          id = seq_len(n),
          time = round(stats::runif(n, 5, 15), 4),
          evid = 2L,
          amt = 0,
          cmt = "ev",
          X = 1,
          TSW = 0
        )
        .d$TSW <- .switch$time[.d$id]
        .switch$TSW <- .switch$time
        .d <- rbind(.d, .switch)
        .d$X <- as.numeric(.d$time >= .d$TSW)
        .d[order(.d$id, .d$time), ]
      },
      covariates = c("X", "TSW"),
      # the simulated event records carry the covariate in force at their time
      fixSim = function(d) {
        d$X <- as.numeric(d$time >= d$TSW)
        d
      },
      # NONMEM data items in $DES come from the record ending the interval, so
      # each record carries the value that held just before it
      addColumns = function(d) {
        d$XNM <- .tteValueBefore(d, d$X)
        d
      }
    ),
    interval = list(
      description = "interval-censored single event (inspections every 4 time units)",
      model = function() {
        ini({
          lsc <- log(15)
          lsh <- log(1.5)
          eta.sc ~ 0.6
        })
        model({
          sc <- exp(lsc + eta.sc)
          sh <- exp(lsh)
          haz <- sh / sc * (time / sc)^(sh - 1)
          ev ~ tte(haz, eventType = "intervalCensored", intervalLength = 4, maxEventNumber = 1)
        })
      },
      start = c(lsc = log(10), lsh = log(1.2), eta.sc = 0.3),
      n = 400L,
      seed = 14618L,
      design = function(n) .tteWindow(n, 0, 40)
    ),
    lateEntry = list(
      description = "Weibull single event with late entry (observation starts after time 0)",
      model = function() {
        ini({
          lsc <- log(15)
          lsh <- log(1.5)
          eta.sc ~ 0.6
        })
        model({
          sc <- exp(lsc + eta.sc)
          sh <- exp(lsh)
          haz <- sh / sc * (time / sc)^(sh - 1)
          ev ~ tte(haz, maxEventNumber = 1)
        })
      },
      start = c(lsc = log(10), lsh = log(1.2), eta.sc = 0.3),
      n = 400L,
      seed = 14619L,
      design = function(n) {
        .start <- stats::runif(n, 0, 10)
        data.frame(
          id = rep(seq_len(n), each = 2),
          time = as.vector(rbind(.start, 40)),
          evid = 0L,
          amt = 0,
          cmt = "ev"
        )
      }
    ),
    dropout = list(
      description = "single event with dropout as a competing process",
      model = function() {
        ini({
          lh1 <- log(0.05)
          lh2 <- log(0.03)
          eta.h1 ~ 0.6
        })
        model({
          h1 <- exp(lh1 + eta.h1)
          h2 <- exp(lh2)
          ev ~ tte(h1, maxEventNumber = 1)
          drop ~ tte(h2, maxEventNumber = 1, dropout = TRUE)
        })
      },
      start = c(lh1 = log(0.1), lh2 = log(0.1), eta.h1 = 0.3),
      n = 400L,
      seed = 14620L,
      design = function(n) {
        .d <- rbind(.tteWindow(n, 0, 30), .tteWindow(n, 0, 30, "drop"))
        .d[order(.d$id, .d$time), ]
      },
      # a dropout ends the observation of `ev` too, but the simulation leaves
      # no record for it there: right-censor `ev` at the dropout
      fixSim = .tteCensorAtDropout
    )
  )
}

#' Time-to-event validation scenarios
#'
#' @param name scenario name(s); all scenarios when NULL
#' @return a named list of scenarios
#' @noRd
tteValidationScenarios <- function(name = NULL) {
  .s <- .tteScenarioList()
  if (is.null(name)) {
    return(.s)
  }
  .bad <- setdiff(name, names(.s))
  if (length(.bad) > 0) {
    stop("unknown tte validation scenario: ", paste(.bad, collapse = ", "), call. = FALSE)
  }
  .s[name]
}

# observation window [from, to] of one endpoint for n subjects
.tteWindow <- function(n, from, to, cmt = "ev") {
  data.frame(
    id = rep(seq_len(n), each = 2),
    time = rep(c(from, to), n),
    evid = 0L,
    amt = 0,
    cmt = cmt
  )
}

# a bolus into `central` at time 0 for every subject
.tteAddDose <- function(d, amt) {
  .dose <- d[!duplicated(d$id), ]
  .dose$time <- 0
  .dose$evid <- 1L
  .dose$amt <- amt
  .dose$cmt <- "central"
  .d <- rbind(.dose, d)
  .d[order(.d$id, .d$time, -.d$evid), ]
}

# per record, the value of `x` on the subject's previous record (its own value
# on the first record)
.tteValueBefore <- function(d, x) {
  .first <- !duplicated(d$ID)
  .prev <- c(x[1], x[-length(x)])
  ifelse(.first, x, .prev)
}

# per record, the time of the last tte event strictly before it (or the
# subject's first record time)
.tteLastEventBefore <- function(d) {
  .out <- numeric(nrow(d))
  for (.id in unique(d$ID)) {
    .w <- which(d$ID == .id)
    .last <- d$TIME[.w[1]]
    for (.i in .w) {
      .out[.i] <- .last
      if (d$EVID[.i] == 0 && d$DV[.i] == 1) .last <- d$TIME[.i]
    }
  }
  .out
}

# right-censor the other tte endpoints at a dropout event
.tteCensorAtDropout <- function(d, dropName = "drop") {
  .drop <- d[d$cmt == dropName & d$dv == 1, ]
  .add <- lapply(seq_len(nrow(.drop)), function(i) {
    .other <- d[d$id == .drop$id[i] & d$cmt != dropName & d$evid == 0, ]
    .ends <- split(.other, .other$cmt)
    do.call(
      rbind,
      lapply(.ends, function(e) {
        .last <- e[nrow(e), ]
        # an endpoint that already ended (its single event happened) needs no row
        if (.last$dv == 1 || .last$time >= .drop$time[i]) {
          return(NULL)
        }
        .last$time <- .drop$time[i]
        .last$dv <- 0
        .last
      })
    )
  })
  .d <- rbind(d, do.call(rbind, .add))
  .d[order(.d$id, .d$time, .d$cmt != dropName), ]
}

# ---- data set: simulation and loading --------------------------------------

# endpoint names of a scenario's model, in DVID order
.tteEndpoints <- function(scenario) {
  rxode2::rxode2(scenario$model)$predDf$cond
}

# NONMEM compartment of each endpoint / dosing compartment: the rxode2 state
# order, with each tte endpoint observed in its cumulative-hazard state and a
# PK endpoint in the compartment its concentration comes from
.tteNmCmt <- function(scenario, cmt) {
  .state <- rxode2::rxModelVars(rxode2::rxode2(scenario$model))$state
  vapply(
    cmt,
    function(x) {
      .w <- match(c(x, paste0("rx_tteH_", x)), .state)
      .w <- .w[!is.na(.w)]
      if (length(.w) == 1L) {
        return(.w)
      }
      .nm <- scenario$nmCmt[[x]]
      if (is.null(.nm)) {
        .nm <- 1L
      }
      as.integer(.nm)
    },
    integer(1),
    USE.NAMES = FALSE
  )
}

#' Simulate one scenario's data set
#'
#' @param scenario a scenario from `tteValidationScenarios()`
#' @return the NONMEM-style data set written to `data/<name>.csv`
#' @noRd
tteValidationSimulate <- function(scenario) {
  .ui <- rxode2::rxode2(scenario$model)
  .cov <- scenario$covariates
  .design <- NULL
  .sim <- rxode2::rxWithSeed(scenario$seed, {
    .design <- scenario$design(scenario$n)
    as.data.frame(rxode2::rxSolve(.ui, .design, addDosing = TRUE, returnType = "data.frame"))
  })
  .ep <- .tteEndpoints(scenario)
  .state <- rxode2::rxModelVars(.ui)$state
  .isDose <- .sim$evid != 0
  .simCmt <- if (is.null(.sim$CMT)) rep(.ui$predDf$cmt[1], nrow(.sim)) else .sim$CMT
  .cmt <- ifelse(.isDose, .state[.simCmt], .ep[match(.simCmt, .ui$predDf$cmt)])
  .d <- data.frame(
    id = .sim$id,
    time = round(.sim$time, 4),
    evid = ifelse(.isDose, 1L, 0L),
    amt = ifelse(.isDose, .sim$amt, 0),
    cmt = .cmt,
    dv = ifelse(.isDose, 0, signif(.sim$sim, 6))
  )
  # the solve returns neither the covariate-only (evid = 2) records nor the
  # design columns: take both from the design
  .other <- .design[.design$evid == 2L, , drop = FALSE]
  if (nrow(.other) > 0) {
    .other$time <- round(.other$time, 4)
    .other$dv <- 0
    .d <- merge(.d, .other[, c("id", "time", "evid", "amt", "cmt", "dv")], all = TRUE, sort = FALSE)
  }
  .first <- .design[!duplicated(.design$id), , drop = FALSE]
  for (.c in .cov) {
    .d[[.c]] <- .first[[.c]][match(.d$id, .first$id)]
  }
  .d <- .d[order(.d$id, .d$time, .d$evid != 1L), ]
  if (is.function(scenario$fixSim)) {
    .d <- scenario$fixSim(.d)
  }
  .out <- data.frame(
    ID = .d$id,
    TIME = .d$time,
    EVID = .d$evid,
    AMT = .d$amt,
    CMT = .tteNmCmt(scenario, .d$cmt),
    DV = .d$dv,
    MDV = as.integer(.d$evid != 0),
    DVID = ifelse(.d$evid == 0, match(.d$cmt, .ep), 0L)
  )
  for (.c in .cov) {
    .out[[.c]] <- signif(.d[[.c]], 6)
  }
  if (is.function(scenario$addColumns)) {
    .out <- scenario$addColumns(.out)
  }
  .out
}

#' Read a scenario's data set for nlmixr2
#'
#' @param name scenario name
#' @param dir the validation directory (holding `data/`)
#' @return a data frame with `cmt` set to the endpoint (or dosing) name
#' @noRd
tteValidationData <- function(name, dir) {
  .s <- tteValidationScenarios(name)[[1]]
  .d <- utils::read.csv(file.path(dir, "data", paste0(name, ".csv")))
  .ep <- .tteEndpoints(.s)
  .doseCmt <- if (is.null(.s$doseCmt)) "central" else .s$doseCmt
  .out <- data.frame(
    id = .d$ID,
    time = .d$TIME,
    evid = .d$EVID,
    amt = .d$AMT,
    cmt = ifelse(.d$EVID == 1, .doseCmt, .ep[pmax(.d$DVID, 1L)]),
    dv = .d$DV
  )
  for (.c in .s$covariates) {
    .out[[.c]] <- .d[[.c]]
  }
  .out
}

#' A scenario's model at its starting values
#'
#' @param scenario a scenario from `tteValidationScenarios()`
#' @return rxode2 ui
#' @noRd
tteValidationStart <- function(scenario) {
  .ui <- rxode2::rxode2(scenario$model)
  .st <- as.list(scenario$start)
  do.call(rxode2::ini, c(list(.ui), .st))
}
