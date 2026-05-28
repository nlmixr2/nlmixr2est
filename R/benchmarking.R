.nlmixr2BenchmarkStageMap <- c(
  preprocess = "preprocess",
  setup = "setup",
  optimize = "estimate",
  saem = "estimate",
  covariance = "covariance",
  table = "table",
  compress = "compress",
  CWRES = "CWRES",
  NPDE = "NPDE",
  other = "other"
)

.nlmixr2BenchmarkOneCompartment <- function() {
  ini({
    tka <- 0.45
    tcl <- 1
    tv <- 3.45
    eta.ka ~ 0.6
    eta.cl ~ 0.3
    eta.v ~ 0.1
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    d/dt(depot) = -ka * depot
    d/dt(center) = ka * depot - cl / v * center
    cp = center / v
    cp ~ add(add.sd)
  })
}

.nlmixr2BenchmarkOneCompartmentLag <- function() {
  ini({
    ltlag <- log(0.2)
    lka <- log(1.15)
    lcl <- log(0.135)
    lv <- log(8)
    prop.err <- 0.15
    add.err <- 0.6
    eta.tlag ~ 0.5
    eta.ka ~ 0.5
    eta.cl ~ 0.1
    eta.v ~ 0.1
  })
  model({
    tlag <- exp(ltlag + eta.tlag)
    ka <- exp(lka + eta.ka)
    cl <- exp(lcl + eta.cl)
    v <- exp(lv + eta.v)
    d/dt(gut) <- -ka * gut
    d/dt(central) <- ka * gut - (cl / v) * central
    lag(gut) <- tlag
    cp <- central / v
    cp ~ prop(prop.err) + add(add.err)
  })
}

.nlmixr2BenchmarkDefaultCases <- function(calcTables = TRUE) {
  list(
    one_compartment_saem = list(
      est = "saem",
      model = .nlmixr2BenchmarkOneCompartment,
      data = nlmixr2data::theo_sd,
      control = function() {
        saemControl(print = 0L, calcTables = calcTables)
      }
    ),
    one_compartment_focei = list(
      est = "focei",
      model = .nlmixr2BenchmarkOneCompartment,
      data = nlmixr2data::theo_sd,
      control = function() {
        foceiControl(print = 0L, calcTables = calcTables)
      }
    ),
    lag_focei = list(
      est = "focei",
      model = .nlmixr2BenchmarkOneCompartmentLag,
      data = subset(nlmixr2data::warfarin, as.character(dvid) == "cp"),
      control = function() {
        foceiControl(print = 0L, calcTables = calcTables)
      }
    )
  )
}

.nlmixr2BenchmarkRxThreads <- function() {
  .ret <- try(rxode2::getRxThreads(), silent = TRUE)
  if (inherits(.ret, "try-error")) {
    return(NA_integer_)
  }
  as.integer(.ret)
}

.nlmixr2BenchmarkSetRxThreads <- function(threads) {
  if (is.null(threads) || length(threads) != 1L || is.na(threads)) {
    return(invisible())
  }
  try(rxode2::setRxThreads(as.integer(threads)), silent = TRUE)
  invisible()
}

.nlmixr2BenchmarkNormalizeTime <- function(timeDf) {
  checkmate::assertDataFrame(timeDf, min.rows = 1, max.rows = 1, .var.name = "timeDf")
  .raw <- names(timeDf)
  .elapsed <- vapply(.raw, function(.stage) {
    as.numeric(timeDf[[.stage]][1])
  }, numeric(1), USE.NAMES = FALSE)
  .stage <- ifelse(.raw %in% names(.nlmixr2BenchmarkStageMap),
                   unname(.nlmixr2BenchmarkStageMap[.raw]),
                   .raw)
  .ret <- data.frame(
    raw_stage = .raw,
    stage = .stage,
    elapsed = .elapsed,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  .ret <- stats::aggregate(elapsed ~ stage, data = .ret, FUN = sum)
  .rawMap <- stats::aggregate(raw_stage ~ stage,
                              data = data.frame(stage = .stage, raw_stage = .raw,
                                                stringsAsFactors = FALSE),
                              FUN = function(x) paste(x, collapse = ","))
  .ret <- merge(.ret, .rawMap, by = "stage", sort = FALSE)
  .ret <- .ret[match(unique(.stage), .ret$stage), c("stage", "raw_stage", "elapsed"), drop = FALSE]
  row.names(.ret) <- NULL
  .ret
}

.nlmixr2BenchmarkExtractFit <- function(fit, case, threads = .nlmixr2BenchmarkRxThreads()) {
  .time <- fit$time
  .stages <- .nlmixr2BenchmarkNormalizeTime(.time)
  .origData <- fit$origData
  .fitEst <- try(trimws(fit$est), silent = TRUE)
  if (inherits(.fitEst, "try-error") || length(.fitEst) != 1L) .fitEst <- NA_character_
  .method <- try(trimws(fit$method), silent = TRUE)
  if (inherits(.method, "try-error") || length(.method) != 1L) .method <- NA_character_
  .nTheta <- try(length(fit$fixef), silent = TRUE)
  if (inherits(.nTheta, "try-error")) .nTheta <- NA_integer_
  .omega <- try(fit$omega, silent = TRUE)
  if (inherits(.omega, "try-error")) {
    .nEta <- NA_integer_
  } else if (is.list(.omega) && !is.null(.omega$id)) {
    .nEta <- nrow(.omega$id)
  } else if (is.matrix(.omega)) {
    .nEta <- nrow(.omega)
  } else {
    .nEta <- NA_integer_
  }
  .subjects <- if (inherits(.origData, "data.frame") && "ID" %in% names(.origData)) {
    length(unique(.origData$ID))
  } else if (inherits(.origData, "data.frame")) {
    1L
  } else {
    NA_integer_
  }
  .nObs <- if (inherits(.origData, "data.frame") && "EVID" %in% names(.origData)) {
    sum(is.na(.origData$EVID) | .origData$EVID == 0)
  } else if (inherits(.origData, "data.frame") && "DV" %in% names(.origData)) {
    sum(!is.na(.origData$DV))
  } else {
    NA_integer_
  }
  .totalElapsed <- sum(.stages$elapsed)
  data.frame(
    case = rep_len(case, nrow(.stages)),
    estimator = rep_len(.fitEst, nrow(.stages)),
    method = rep_len(.method, nrow(.stages)),
    subjects = rep_len(.subjects, nrow(.stages)),
    observations = rep_len(.nObs, nrow(.stages)),
    ntheta = rep_len(as.integer(.nTheta), nrow(.stages)),
    neta = rep_len(as.integer(.nEta), nrow(.stages)),
    rx_threads = rep_len(as.integer(threads), nrow(.stages)),
    omp_threads = rep_len(Sys.getenv("OMP_NUM_THREADS", unset = NA_character_), nrow(.stages)),
    total_elapsed = rep_len(.totalElapsed, nrow(.stages)),
    .stages,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

.nlmixr2BenchmarkRun <- function(cases = .nlmixr2BenchmarkDefaultCases(),
                                 threads = NULL,
                                 file = NULL,
                                 quiet = TRUE) {
  checkmate::assertList(cases, min.len = 1, names = "named")
  .oldThreads <- .nlmixr2BenchmarkRxThreads()
  on.exit(.nlmixr2BenchmarkSetRxThreads(.oldThreads), add = TRUE)
  .nlmixr2BenchmarkSetRxThreads(threads)
  .threads <- .nlmixr2BenchmarkRxThreads()
  .ret <- lapply(names(cases), function(.caseName) {
    .case <- cases[[.caseName]]
    .control <- if (is.function(.case$control)) .case$control() else .case$control
    .fit <- if (isTRUE(quiet)) {
      suppressWarnings(suppressMessages(
        nlmixr2(.case$model, .case$data, est = .case$est, control = .control)
      ))
    } else {
      nlmixr2(.case$model, .case$data, est = .case$est, control = .control)
    }
    .nlmixr2BenchmarkExtractFit(.fit, .caseName, threads = .threads)
  })
  .ret <- do.call(rbind, .ret)
  row.names(.ret) <- NULL
  if (!is.null(file)) {
    utils::write.csv(.ret, file = file, row.names = FALSE)
  }
  .ret
}
