# SAEM refactor evaluation harness
#
# Simulate-from-truth + refit, then score parameter recovery (bias / RMSE),
# convergence, objective, and runtime.  Reused across the SAEM-unification
# phases to (a) assert the shared-inner path is numerically equivalent to the
# classic path and (b) gate any default flip on the "strictly better + no
# regressions" bar.
#
# Not part of the package test suite; source() this file after
# devtools::load_all(".").  Example:
#
#   devtools::load_all(".")
#   source("tools/saemRefactorEval.R")
#   res <- seCompareConfigs("oneCmt",
#            configs = list(classic = saemControl(sharedInner = "classic"),
#                           shared  = saemControl(sharedInner = "shared")),
#            nsub = 40, seeds = 1:3)
#   print(res$summary)

## ---------------------------------------------------------------------------
## Curated model registry
## ---------------------------------------------------------------------------
## Each entry is a list with:
##   ui       function() returning a nlmixr2 model (ini/model blocks)
##   truth    named numeric of the population values used to simulate; names
##            match c(fixef(fit), omega-diagonal eta names)
##   dose     dose amount for the default simulator
##   times    observation times for the default simulator
##   simData  (optional) function(entry, nsub, seed) overriding the default
##            rxSolve-based simulator (used for the mixture model)

.seModels <- function() {
  list(
    ## 1-cmt oral, proportional error, one mu-referenced eta.  Baseline
    ## equivalence + recovery model.
    oneCmt = list(
      ui = function() {
        function() {
          ini({
            tka <- log(1.2)
            tcl <- log(0.25)
            tv  <- log(5)
            eta.cl ~ 0.09
            prop.sd <- 0.15
          })
          model({
            ka <- exp(tka)
            cl <- exp(tcl + eta.cl)
            v  <- exp(tv)
            d/dt(depot)   <- -ka * depot
            d/dt(central) <-  ka * depot - cl / v * central
            cp <- central / v
            cp ~ prop(prop.sd)
          })
        }
      },
      truth = c(tka = log(1.2), tcl = log(0.25), tv = log(5),
                prop.sd = 0.15, eta.cl = 0.09),
      dose = 100,
      times = c(0.25, 0.5, 1, 2, 4, 6, 8, 12, 24)
    ),

    ## 1-cmt with a non-mu-referenced structural theta (v enters through a
    ## Hill-like nonlinearity so it is not a clean theta+eta).  Exercises the
    ## nonMuTheta path.
    nonMu = list(
      ui = function() {
        function() {
          ini({
            tka  <- log(1.2)
            tcl  <- log(0.25)
            tv   <- log(5)
            pow  <- 0.75
            eta.cl ~ 0.09
            prop.sd <- 0.15
          })
          model({
            ka <- exp(tka)
            cl <- exp(tcl + eta.cl)
            v  <- exp(tv) * (1 + pow)
            d/dt(depot)   <- -ka * depot
            d/dt(central) <-  ka * depot - cl / v * central
            cp <- central / v
            cp ~ prop(prop.sd)
          })
        }
      },
      truth = c(tka = log(1.2), tcl = log(0.25), tv = log(5),
                pow = 0.75, prop.sd = 0.15, eta.cl = 0.09),
      dose = 100,
      times = c(0.25, 0.5, 1, 2, 4, 6, 8, 12, 24)
    ),

    ## 2-component mixture on clearance.  Exercises mixProbMethod (incl.
    ## "regress").  Custom simulator draws membership from the true prob.
    mix2 = list(
      ui = function() {
        function() {
          ini({
            tka  <- log(1.2)
            tcl1 <- log(0.15)
            tcl2 <- log(0.45)
            tv   <- log(5)
            p1   <- 0.6
            eta.cl ~ 0.05
            prop.sd <- 0.12
          })
          model({
            ka <- exp(tka)
            cl <- mix(exp(tcl1 + eta.cl), p1, exp(tcl2 + eta.cl))
            v  <- exp(tv)
            d/dt(depot)   <- -ka * depot
            d/dt(central) <-  ka * depot - cl / v * central
            cp <- central / v
            cp ~ prop(prop.sd)
          })
        }
      },
      truth = c(tka = log(1.2), tcl1 = log(0.15), tcl2 = log(0.45),
                tv = log(5), p1 = 0.6, prop.sd = 0.12, eta.cl = 0.05),
      dose = 100,
      times = c(0.25, 0.5, 1, 2, 4, 6, 8, 12, 24),
      simData = function(entry, nsub, seed) {
        set.seed(seed)
        tr <- entry$truth
        comp <- ifelse(stats::runif(nsub) < tr[["p1"]], 1L, 2L)
        tcl <- ifelse(comp == 1L, tr[["tcl1"]], tr[["tcl2"]])
        eta <- stats::rnorm(nsub, 0, sqrt(tr[["eta.cl"]]))
        ka <- exp(tr[["tka"]]); v <- exp(tr[["tv"]]); cl <- exp(tcl + eta)
        do.call(rbind, lapply(seq_len(nsub), function(i) {
          tt <- entry$times
          ke <- cl[i] / v
          ## 1-cmt oral closed form (ka != ke)
          f <- entry$dose / v * ka[1] / (ka[1] - ke) *
            (exp(-ke * tt) - exp(-ka[1] * tt))
          dv <- f * (1 + stats::rnorm(length(tt), 0, tr[["prop.sd"]]))
          rbind(
            data.frame(ID = i, TIME = 0, DV = NA_real_, AMT = entry$dose, EVID = 1),
            data.frame(ID = i, TIME = tt, DV = dv, AMT = 0, EVID = 0)
          )
        }))
      }
    )
  )
}

## ---------------------------------------------------------------------------
## Default simulator (rxSolve on the true UI; sim = error-included DV)
## ---------------------------------------------------------------------------
.seSimDefault <- function(entry, nsub, seed) {
  ui <- entry$ui()
  et <- rxode2::et(amt = entry$dose, time = 0, cmt = "depot")
  et <- rxode2::et(et, entry$times)
  et <- rxode2::et(et, id = seq_len(nsub))
  rxode2::rxSetSeed(seed)
  set.seed(seed)
  s <- as.data.frame(rxode2::rxSolve(ui, et, addDosing = FALSE))
  obs <- data.frame(ID = s$id, TIME = s$time, DV = s$sim, AMT = 0, EVID = 0)
  dose <- data.frame(ID = seq_len(nsub), TIME = 0, DV = NA_real_,
                     AMT = entry$dose, EVID = 1)
  d <- rbind(dose, obs)
  d[order(d$ID, d$TIME, -d$EVID), ]
}

seSimData <- function(entry, nsub, seed) {
  if (is.function(entry$simData)) return(entry$simData(entry, nsub, seed))
  .seSimDefault(entry, nsub, seed)
}

## ---------------------------------------------------------------------------
## Fit + score
## ---------------------------------------------------------------------------
.seEstVec <- function(fit) {
  est <- tryCatch(nlmixr2est::fixef(fit), error = function(e) stats::coef(fit))
  om <- tryCatch({
    o <- fit$omega
    if (is.null(o)) numeric(0) else stats::setNames(diag(o), rownames(o))
  }, error = function(e) numeric(0))
  c(est, om[!(names(om) %in% names(est))])
}

seFitScore <- function(entry, data, control, seed) {
  ui <- entry$ui()
  t0 <- proc.time()[["elapsed"]]
  fit <- tryCatch(
    suppressWarnings(nlmixr2est::nlmixr2(ui, data, est = "saem", control = control)),
    error = function(e) e)
  elapsed <- proc.time()[["elapsed"]] - t0
  if (inherits(fit, "error")) {
    return(list(converged = FALSE, error = conditionMessage(fit),
                time = elapsed, seed = seed))
  }
  est <- .seEstVec(fit)
  truth <- entry$truth
  nm <- intersect(names(truth), names(est))
  objf <- tryCatch(as.numeric(fit$objf), error = function(e) NA_real_)
  bias <- est[nm] - truth[nm]
  list(converged = all(is.finite(est[nm])),
       est = est[nm], truth = truth[nm], bias = bias,
       objf = objf, time = elapsed, seed = seed)
}

## ---------------------------------------------------------------------------
## Compare configs across seeds
## ---------------------------------------------------------------------------
## Returns a list with $perParam (bias/RMSE per parameter per config) and
## $summary (convergence rate, mean objf, mean runtime per config).
seCompareConfigs <- function(entryName, configs, nsub = 40, seeds = 1:3,
                             models = .seModels()) {
  entry <- models[[entryName]]
  if (is.null(entry)) stop("unknown model '", entryName, "'")
  runs <- list()
  for (cn in names(configs)) {
    for (sd in seeds) {
      data <- seSimData(entry, nsub, sd)
      sc <- seFitScore(entry, data, configs[[cn]], sd)
      sc$config <- cn
      runs[[length(runs) + 1L]] <- sc
    }
  }
  ok <- Filter(function(x) isTRUE(x$converged), runs)
  ## per-parameter bias/RMSE
  perParam <- do.call(rbind, lapply(names(configs), function(cn) {
    cok <- Filter(function(x) x$config == cn, ok)
    if (length(cok) == 0L) return(NULL)
    pn <- names(entry$truth)
    do.call(rbind, lapply(pn, function(p) {
      b <- vapply(cok, function(x) unname(x$bias[p]), numeric(1))
      b <- b[is.finite(b)]
      if (length(b) == 0L) return(NULL)
      data.frame(config = cn, param = p, truth = unname(entry$truth[p]),
                 meanBias = mean(b), rmse = sqrt(mean(b^2)), n = length(b))
    }))
  }))
  summary <- do.call(rbind, lapply(names(configs), function(cn) {
    cruns <- Filter(function(x) x$config == cn, runs)
    cok <- Filter(function(x) isTRUE(x$converged), cruns)
    data.frame(config = cn,
               convRate = length(cok) / length(cruns),
               meanObjf = mean(vapply(cok, function(x) x$objf, numeric(1))),
               meanTime = mean(vapply(cruns, function(x) x$time, numeric(1))))
  }))
  list(perParam = perParam, summary = summary, runs = runs)
}

## Convenience: is config B strictly better than A (lower RMSE on every
## parameter and no worse convergence)?  Encodes the default-flip bar.
seStrictlyBetter <- function(cmp, a, b, tol = 1e-8) {
  pp <- cmp$perParam
  ra <- pp[pp$config == a, c("param", "rmse")]
  rb <- pp[pp$config == b, c("param", "rmse")]
  m <- merge(ra, rb, by = "param", suffixes = c(".a", ".b"))
  sa <- cmp$summary[cmp$summary$config == a, "convRate"]
  sb <- cmp$summary[cmp$summary$config == b, "convRate"]
  all(m$rmse.b <= m$rmse.a + tol) && sb + tol >= sa &&
    (any(m$rmse.b < m$rmse.a - tol) || sb > sa + tol)
}
