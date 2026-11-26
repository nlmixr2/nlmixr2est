# Shared checks of the tte validation scenarios (inst/validation/tte).

.tteValidationEnv <- new.env()

# the validation directory with its scripts sourced into one environment
.tteValidation <- function() {
  if (is.null(.tteValidationEnv$dir)) {
    .dir <- system.file("validation", "tte", package = "nlmixr2est")
    if (.dir == "") {
      return(NULL)
    }
    .env <- new.env()
    suppressMessages({
      for (.f in c("scenarios.R", "reference.R", "exact.R")) {
        sys.source(file.path(.dir, .f), envir = .env)
      }
    })
    .tteValidationEnv$dir <- .dir
    .tteValidationEnv$env <- .env
    .tteValidationEnv$exact <- list()
  }
  .tteValidationEnv
}

# exact marginal ML of a scenario file, computed once per session
.tteValidationExactFit <- function(file) {
  .v <- .tteValidation()
  if (is.null(.v$exact[[file]])) {
    .name <- sub("-reduced$", "", file)
    .fixed <- if (grepl("-reduced$", file)) c(bx = 0) else NULL
    .v$exact[[file]] <- suppressMessages(.v$env$tteValidationExact(.name, .v$dir, fixed = .fixed))
  }
  .v$exact[[file]]
}

# fit a scenario file with one estimation method
.tteValidationFit <- function(file, est, control) {
  .v <- .tteValidation()
  .name <- sub("-reduced$", "", file)
  .s <- .v$env$tteValidationScenarios(.name)[[1]]
  .ui <- suppressMessages(.v$env$tteValidationStart(.s))
  if (grepl("-reduced$", file)) {
    .ui <- .s$reduced(.ui)
  }
  .d <- suppressMessages(.v$env$tteValidationData(.name, .v$dir))
  suppressMessages(suppressWarnings(.nlmixr(.ui, .d, est = est, control = control)))
}

# estimates of a fit by nlmixr2 parameter name (eta variances)
.tteValidationEst <- function(fit) {
  .ini <- fit$iniDf
  .ini <- .ini[is.na(.ini$neta1) | .ini$neta1 == .ini$neta2, ]
  stats::setNames(.ini$est, .ini$name)
}

# compare a fit with a reference (est, se) within `tol` reference standard
# errors; fixed parameters and parameters without a standard error are skipped
.tteValidationExpectClose <- function(fit, est, se, tol, info) {
  .f <- .tteValidationEst(fit)
  .p <- intersect(names(est), names(.f))
  .p <- .p[!is.na(se[.p]) & !fit$iniDf$fix[match(.p, fit$iniDf$name)]]
  for (.n in .p) {
    expect_lt(abs(.f[[.n]] - est[[.n]]), tol * se[[.n]], label = paste(info, .n, "difference / se"))
  }
}

# check one scenario fitted by `est` against the exact likelihood and the
# NONMEM / Monolix references stored with it
.tteValidationCheck <- function(file, est, control, tol, llTol = NA, fit = NULL) {
  if (is.null(fit)) {
    fit <- .tteValidationFit(file, est, control)
  }
  .info <- paste0(file, " (", est, ")")
  .ex <- .tteValidationExactFit(file)
  if (!is.null(.ex)) {
    .tteValidationExpectClose(fit, .ex$est, .ex$se, tol, paste(.info, "vs exact"))
    if (!is.na(llTol)) {
      expect_lt(abs(-2 * as.numeric(logLik(fit)) - .ex$m2ll), llTol, label = paste(.info, "-2LL"))
    }
  }
  .v <- .tteValidation()
  for (.tool in c("nonmem", "monolix")) {
    .ref <- .v$env$tteValidationReference(.v$dir, file, .tool)
    if (is.null(.ref)) {
      next
    }
    .ref <- .ref[.ref$parameter != "OFV", ]
    .tteValidationExpectClose(
      fit,
      stats::setNames(.ref$estimate, .ref$parameter),
      stats::setNames(.ref$se, .ref$parameter),
      tol,
      paste(.info, "vs", .tool)
    )
  }
  invisible(fit)
}

# the objective function difference between the full and reduced tvCov
# models, against the exact one and each stored reference
.tteValidationCheckDofv <- function(full, reduced, est, tol) {
  .d <- as.numeric(reduced$objf - full$objf)
  .ex <- .tteValidationExactFit("tvCov-reduced")$m2ll - .tteValidationExactFit("tvCov")$m2ll
  expect_lt(abs(.d - .ex), tol, label = paste("tvCov dOFV (", est, ") vs exact"))
  .v <- .tteValidation()
  for (.tool in c("nonmem", "monolix")) {
    .f <- .v$env$tteValidationReference(.v$dir, "tvCov", .tool)
    .r <- .v$env$tteValidationReference(.v$dir, "tvCov-reduced", .tool)
    if (is.null(.f) || is.null(.r)) {
      next
    }
    .ref <- .r$estimate[.r$parameter == "OFV"] - .f$estimate[.f$parameter == "OFV"]
    expect_lt(abs(.d - .ref), tol, label = paste("tvCov dOFV (", est, ") vs", .tool))
  }
}

# thetas of a fit within `tol` of its standard errors from the simulation truth
.tteValidationExpectTruth <- function(fit, file, tol = 3.5) {
  .v <- .tteValidation()
  .s <- .v$env$tteValidationScenarios(sub("-reduced$", "", file))[[1]]
  .ini <- suppressMessages(rxode2::rxode2(.s$model))$iniDf
  .se <- fit$parFixedDf$SE
  names(.se) <- rownames(fit$parFixedDf)
  .f <- .tteValidationEst(fit)
  for (.n in .ini$name[!is.na(.ini$ntheta)]) {
    if (is.na(.se[.n])) {
      next
    }
    expect_lt(abs(.f[[.n]] - .ini$est[.ini$name == .n]), tol * .se[[.n]], label = paste(file, .n, "vs truth / se"))
  }
}
