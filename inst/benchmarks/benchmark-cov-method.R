# Covariance-step cost of covMethod = "analytic" versus the finite-difference "r,s"
# (covFull = TRUE) as the model grows: ODE states, etas, thetas and residual
# parameters (#1217).  Each case is simulated from its own initial estimates and
# fitted with maxOuterIterations = 0, so both methods run at the same point and
# only the covariance step differs.
#
#   OMP_NUM_THREADS=1 taskset -c 2 Rscript inst/benchmarks/benchmark-cov-method.R [out.csv] [grid]
#
# `grid` is "axes" (one axis at a time from a base, the default), "large", "ode"
# (transit chains with combined error) or "all".
# library(nlmixr2est) loads the INSTALLED package: install the tree first.
suppressPackageStartupMessages({
  library(nlmixr2est)
  library(rxode2)
})
setRxThreads(1L)

.args <- commandArgs(trailingOnly = TRUE)
.out <- if (length(.args) >= 1L) .args[1] else "cov-method-benchmark.csv"
.grid <- if (length(.args) >= 2L) .args[2] else "axes"

# A depot -> nTransit transit states -> central (+ nPeriph peripherals) linear model, or
# Michaelis-Menten elimination (mm).  Structural parameters in order: ka, cl, v, ktr,
# then q/vp per peripheral (and km for mm); the first nEta of them get an eta.  nCov
# covariate effects on cl add thetas only.  nSigma 1 is proportional, 2 combined.
.covBenchModel <- function(nTransit = 0L, nPeriph = 0L, nEta = 3L, nCov = 0L, nSigma = 1L, mm = FALSE) {
  .par <- c("ka", "cl", "v")
  .val <- c(ka = log(1), cl = log(4), v = log(40))
  if (nTransit > 0L) {
    .par <- c(.par, "ktr")
    .val <- c(.val, ktr = log(3))
  }
  for (.p in seq_len(nPeriph)) {
    .q <- paste0("q", .p)
    .vp <- paste0("vp", .p)
    .par <- c(.par, .q, .vp)
    .val[[.q]] <- log(2 / .p)
    .val[[.vp]] <- log(60 * .p)
  }
  if (mm) {
    .par <- c(.par, "km")
    .val[["km"]] <- log(2)
  }
  stopifnot(nEta <= length(.par))
  .ini <- c(
    sprintf("t%s <- %s", .par, format(.val[.par], digits = 6)),
    sprintf("cov%d <- 0.1", seq_len(nCov)),
    sprintf("eta.%s ~ 0.1", .par[seq_len(nEta)]),
    if (nSigma >= 2L) "add.sd <- 0.1",
    "prop.sd <- 0.1"
  )
  .lp <- function(p) {
    .e <- if (p %in% .par[seq_len(nEta)]) sprintf(" + eta.%s", p) else ""
    .c <- if (p == "cl" && nCov > 0L) {
      paste0(" + ", paste0(sprintf("cov%d * x%d", seq_len(nCov), seq_len(nCov)), collapse = " + "))
    } else {
      ""
    }
    sprintf("%s <- exp(t%s%s%s)", p, p, .e, .c)
  }
  .into <- if (nTransit > 0L) "tr1" else "center"
  .model <- c(
    vapply(.par, .lp, character(1)),
    sprintf("d/dt(depot) <- -ka * depot"),
    if (nTransit > 0L) {
      vapply(
        seq_len(nTransit),
        function(i) {
          .from <- if (i == 1L) "ka * depot" else sprintf("ktr * tr%d", i - 1L)
          sprintf("d/dt(tr%d) <- %s - ktr * tr%d", i, .from, i)
        },
        character(1)
      )
    },
    {
      .in <- if (nTransit > 0L) sprintf("ktr * tr%d", nTransit) else "ka * depot"
      .el <- if (mm) "cl * center / v / (km + center / v)" else "cl / v * center"
      .pe <- vapply(
        seq_len(nPeriph),
        function(p) {
          sprintf(" - q%d / v * center + q%d / vp%d * periph%d", p, p, p, p)
        },
        character(1)
      )
      sprintf("d/dt(center) <- %s - %s%s", .in, .el, paste(.pe, collapse = ""))
    },
    vapply(
      seq_len(nPeriph),
      function(p) {
        sprintf("d/dt(periph%d) <- q%d / v * center - q%d / vp%d * periph%d", p, p, p, p, p)
      },
      character(1)
    ),
    "cp <- center / v",
    if (nSigma >= 2L) "cp ~ add(add.sd) + prop(prop.sd)" else "cp ~ prop(prop.sd)"
  )
  .f <- function() {}
  body(.f) <- str2lang(paste0(
    "{ini({",
    paste(.ini, collapse = "\n"),
    "})\nmodel({",
    paste(.model, collapse = "\n"),
    "})}"
  ))
  .f
}

.covBenchData <- function(ui, nCov, nSub = 40L, seed = 1217L) {
  set.seed(seed)
  .t <- c(0.25, 0.5, 1, 2, 3, 4, 6, 8, 12, 16, 24, 36)
  .ev <- et(amt = 100, cmt = "depot") |> et(.t) |> et(id = seq_len(nSub))
  .ev <- as.data.frame(.ev)
  for (.k in seq_len(nCov)) {
    .x <- stats::setNames(stats::rnorm(nSub, 0, 0.3), seq_len(nSub))
    .ev[[paste0("x", .k)]] <- .x[as.character(.ev$id)]
  }
  .s <- rxSolve(ui, .ev, returnType = "data.frame", addDosing = TRUE, seed = seed)
  .d <- .ev
  .d$dv <- NA_real_
  .obs <- .d$evid == 0
  .d$dv[.obs] <- .s$sim[.s$evid == 0]
  .d$dv[.obs & .d$dv <= 0] <- 1e-3
  .std <- names(.d) %in% c("id", "time", "amt", "evid", "cmt", "dv")
  names(.d)[.std] <- toupper(names(.d)[.std])
  .d
}

.covBenchOne <- function(case, mod, dat, covMethod) {
  .ctl <- foceiControl(print = 0L, maxOuterIterations = 0L, covMethod = covMethod)
  .w <- character(0)
  .t0 <- proc.time()[["elapsed"]]
  .fit <- withCallingHandlers(
    tryCatch(
      suppressMessages(nlmixr2(mod, dat, "focei", control = .ctl)),
      error = function(e) e
    ),
    warning = function(w) {
      .w <<- c(.w, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  .wall <- proc.time()[["elapsed"]] - .t0
  if (inherits(.fit, "error")) {
    return(data.frame(
      case,
      covMethod = covMethod,
      wall = .wall,
      covTime = NA_real_,
      estTime = NA_real_,
      installed = paste("error:", conditionMessage(.fit)),
      evals = NA_integer_,
      se = I(list(NULL))
    ))
  }
  .tm <- .fit$time
  .ev <- .fit$env$covEvals
  .cov <- .fit$cov
  data.frame(
    case,
    covMethod = covMethod,
    wall = .wall,
    covTime = if ("covariance" %in% names(.tm)) .tm$covariance else NA_real_,
    estTime = sum(unlist(.tm[setdiff(names(.tm), c("covariance", "table"))])),
    installed = if (length(.fit$covMethod)) as.character(.fit$covMethod)[1] else NA_character_,
    evals = if (is.null(.ev)) NA_integer_ else sum(.ev),
    se = I(list(if (is.matrix(.cov)) sqrt(diag(.cov)) else NULL))
  )
}

.covBenchCases <- function(grid) {
  .base <- list(nTransit = 0L, nPeriph = 0L, nEta = 3L, nCov = 0L, nSigma = 1L, mm = FALSE)
  .vary <- function(...) {
    .l <- list(...)
    .n <- names(.l)[1]
    lapply(.l[[1]], function(v) {
      .c <- .base
      .c[[.n]] <- v
      for (.o in names(.l)[-1]) {
        .c[[.o]] <- .l[[.o]]
      }
      .c
    })
  }
  .axes <- c(
    list(.base),
    .vary(nTransit = c(1L, 3L, 6L, 10L, 16L)),
    .vary(nEta = c(1L, 2L, 4L, 6L, 8L), nTransit = 1L, nPeriph = 2L),
    .vary(nCov = c(2L, 4L, 8L, 16L)),
    .vary(nSigma = 2L),
    .vary(mm = TRUE, nEta = 4L)
  )
  .large <- list(
    list(nTransit = 6L, nPeriph = 2L, nEta = 8L, nCov = 8L, nSigma = 2L, mm = FALSE),
    list(nTransit = 16L, nPeriph = 2L, nEta = 8L, nCov = 16L, nSigma = 2L, mm = FALSE),
    list(nTransit = 10L, nPeriph = 2L, nEta = 9L, nCov = 4L, nSigma = 2L, mm = TRUE),
    list(nTransit = 30L, nPeriph = 1L, nEta = 4L, nCov = 0L, nSigma = 1L, mm = FALSE)
  )
  # combined error: a proportional-only model declines the analytic step when a long
  # chain predicts near zero at the first samples
  .ode <- c(
    .vary(nTransit = c(3L, 6L, 10L, 16L, 30L), nSigma = 2L),
    .vary(nTransit = c(10L, 30L), nSigma = 2L, nEta = 4L, nCov = 8L)
  )
  switch(grid, axes = .axes, large = .large, ode = .ode, all = c(.axes, .large, .ode), stop("unknown grid"))
}

.res <- list()
for (.c in .covBenchCases(.grid)) {
  .mod <- do.call(.covBenchModel, .c)
  .ui <- nlmixr2(.mod)
  .dat <- .covBenchData(.ui, .c$nCov)
  .case <- data.frame(
    nOde = 2L + .c$nTransit + .c$nPeriph,
    nTransit = .c$nTransit,
    nPeriph = .c$nPeriph,
    nEta = .c$nEta,
    nTheta = sum(!is.na(.ui$iniDf$ntheta)),
    nSigma = .c$nSigma,
    mm = .c$mm,
    nObs = sum(.dat$EVID == 0)
  )
  # one warm-up compiles both models before either is timed
  invisible(.covBenchOne(.case, .mod, .dat, ""))
  # the analytic augmented model compiles on its first use: run it cold, then warm
  for (.i in seq_along(.ms <- c("r,s", "analytic", "analytic"))) {
    .m <- .ms[.i]
    .r <- .covBenchOne(.case, .mod, .dat, .m)
    .r$run <- sum(.ms[seq_len(.i)] == .m)
    .res[[length(.res) + 1L]] <- .r
    message(sprintf(
      "nOde=%d nEta=%d nTheta=%d nSigma=%d mm=%s %-8s run %d cov=%.2fs est=%.2fs %s",
      .case$nOde,
      .case$nEta,
      .case$nTheta,
      .case$nSigma,
      .case$mm,
      .m,
      .r$run,
      .r$covTime,
      .r$estTime,
      .r$installed
    ))
  }
  .all <- do.call(rbind, .res)
  .all$se <- vapply(.all$se, function(x) if (is.numeric(x)) paste(signif(x, 6), collapse = ";") else "", character(1))
  utils::write.csv(.all, .out, row.names = FALSE)
}
