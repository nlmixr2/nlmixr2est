# Reading NONMEM and Monolix results into reference files, and reading those
# back.  A reference file reference/<tool>/<scenario>.csv has one row per
# parameter (nlmixr2 names; random-effect variances, not SDs) and an `OFV`
# row holding -2 log-likelihood, with the tool version and run date.

# a scenario's model without its tte() endpoints (tteTime() read as 0): the
# parameter names do not depend on them, and an rxode2 without tte() support
# can still parse it
.tteRefUi <- function(model) {
  if (!is.function(model)) {
    return(rxode2::rxode2(model))
  }
  .sub <- function(x) {
    if (!is.call(x)) {
      return(x)
    }
    if (identical(x[[1]], quote(tteTime))) {
      return(0)
    }
    as.call(lapply(as.list(x), .sub))
  }
  .isTte <- function(x) {
    is.call(x) && identical(x[[1]], quote(`~`)) && is.call(x[[3]]) && identical(x[[3]][[1]], quote(tte))
  }
  .b <- body(model)
  for (.i in seq_along(.b)[-1]) {
    if (is.call(.b[[.i]]) && identical(.b[[.i]][[1]], quote(model))) {
      .m <- as.list(.b[[.i]][[2]])
      .m <- .m[!vapply(.m, .isTte, logical(1))]
      .b[[.i]][[2]] <- .sub(as.call(.m))
    }
  }
  body(model) <- .b
  rxode2::rxode2(model)
}

# theta names (in ini() order) and eta names of a scenario's model
.tteParNames <- function(scenario) {
  .ini <- .tteRefUi(scenario$model)$iniDf
  list(
    theta = .ini$name[!is.na(.ini$ntheta)][order(.ini$ntheta[!is.na(.ini$ntheta)])],
    eta = .ini$name[is.na(.ini$ntheta) & .ini$neta1 == .ini$neta2]
  )
}

# the model of a scenario file: `tvCov-reduced` is the reduced tvCov model
.tteRefScenario <- function(file) {
  .name <- sub("-reduced$", "", file)
  .s <- tteValidationScenarios(.name)[[1]]
  if (grepl("-reduced$", file)) {
    .s$model <- .s$reduced(.tteRefUi(.s$model))
  }
  .s
}

#' Read a NONMEM run of a scenario
#'
#' @param ext path of the run's .ext file
#' @param lst path of the run's .lst file (for the version)
#' @param scenario the scenario the run fitted
#' @return reference data frame
#' @noRd
tteNonmemRead <- function(ext, lst, scenario) {
  .l <- readLines(ext)
  .tab <- grep("^TABLE NO", .l)
  .l <- .l[(max(.tab) + 1):length(.l)]
  .e <- utils::read.table(text = .l, header = TRUE, check.names = FALSE)
  .est <- .e[.e$ITERATION == -1000000000, ]
  .se <- .e[.e$ITERATION == -1000000001, ]
  .nm <- .tteParNames(scenario)
  .col <- c(
    paste0("THETA", seq_along(.nm$theta)),
    sprintf("OMEGA(%d,%d)", seq_along(.nm$eta), seq_along(.nm$eta))
  )
  .val <- function(row, col) if (nrow(row) == 0) NA_real_ else as.numeric(row[[col]])
  .ver <- grep("VERSION", grep("\\(NONMEM\\)", readLines(lst), value = TRUE), value = TRUE)
  .ver <- if (length(.ver) == 0) NA_character_ else sub(".*VERSION\\s+", "", .ver[1])
  data.frame(
    tool = "nonmem",
    version = trimws(.ver),
    date = format(file.mtime(ext), "%Y-%m-%d"),
    parameter = c(.nm$theta, .nm$eta, "OFV"),
    estimate = c(vapply(.col, function(c) .val(.est, c), numeric(1)), .val(.est, "OBJ")),
    se = c(
      vapply(
        .col,
        function(c) {
          .v <- .val(.se, c)
          # NONMEM writes 1e10 for a fixed parameter's standard error
          if (!is.na(.v) && .v >= 1e10) NA_real_ else .v
        },
        numeric(1)
      ),
      NA_real_
    ),
    row.names = NULL
  )
}

# a lixoftConnectors result as a named numeric vector: older versions return
# named vectors, newer ones data frames with a parameter (or criteria) column
# and the values in `col` (by default the first numeric column)
.tteMlxVector <- function(x, col = NULL) {
  if (is.data.frame(x)) {
    .key <- intersect(c("parameter", "criteria"), names(x))[1]
    if (is.null(col) || !col %in% names(x)) {
      col <- names(x)[vapply(x, is.numeric, logical(1))][1]
    }
    return(stats::setNames(as.numeric(x[[col]]), x[[.key]]))
  }
  unlist(x)
}

#' Reference rows of a Monolix fit of a scenario
#'
#' @param pop estimated population parameters, named by Monolix parameter
#' @param se their standard errors, named the same way
#' @param ofv -2 log-likelihood
#' @param scenario the scenario the run fitted
#' @param version Monolix version
#' @param date run date
#' @return reference data frame
#' @noRd
tteMonolixReference <- function(
  pop,
  se,
  ofv,
  scenario,
  version = NA_character_,
  date = format(Sys.Date(), "%Y-%m-%d")
) {
  .nm <- .tteParNames(scenario)
  .mu <- .tteRefUi(scenario$model)$muRefDataFrame
  .get <- function(x, v) if (x %in% names(v)) as.numeric(v[[x]]) else NA_real_
  # residual-error parameters keep their nlmixr2 name; Monolix names them b
  .mlx <- vapply(
    .nm$theta,
    function(x) {
      if (!is.null(scenario$monolixNames[[x]])) scenario$monolixNames[[x]] else paste0(x, "_pop")
    },
    character(1)
  )
  .est <- vapply(.mlx, .get, numeric(1), v = pop)
  .se <- vapply(.mlx, .get, numeric(1), v = se)
  # omega_<theta> is the SD of the eta on <theta>: report the variance
  .om <- paste0("omega_", .mu$theta[match(.nm$eta, .mu$eta)])
  .sd <- vapply(.om, .get, numeric(1), v = pop)
  .sdSe <- vapply(.om, .get, numeric(1), v = se)
  data.frame(
    tool = "monolix",
    version = version,
    date = date,
    parameter = c(.nm$theta, .nm$eta, "OFV"),
    estimate = c(.est, .sd^2, as.numeric(ofv)),
    se = c(.se, 2 * .sd * .sdSe, NA_real_),
    row.names = NULL
  )
}

#' Read the Monolix fit of a scenario from the lixoftConnectors session
#'
#' The project must have been loaded and run in this session.
#'
#' @param scenario the scenario the run fitted
#' @param version Monolix version
#' @return reference data frame
#' @noRd
tteMonolixRead <- function(scenario, version = NA_character_) {
  .pop <- .tteMlxVector(lixoftConnectors::getEstimatedPopulationParameters(), "value")
  .se <- .tteMlxVector(lixoftConnectors::getEstimatedStandardErrors()$stochasticApproximation, "se")
  .ll <- .tteMlxVector(lixoftConnectors::getEstimatedLogLikelihood()$importanceSampling)
  tteMonolixReference(.pop, .se, .ll[["OFV"]], scenario, version)
}

#' Read a stored reference fit
#'
#' @param dir the validation directory
#' @param file scenario file name (`tvCov-reduced` for the nested model)
#' @param tool "nonmem" or "monolix"
#' @return the reference data frame, or NULL when there is none
#' @noRd
tteValidationReference <- function(dir, file, tool) {
  .f <- file.path(dir, "reference", tool, paste0(file, ".csv"))
  if (!file.exists(.f)) {
    return(NULL)
  }
  utils::read.csv(.f, stringsAsFactors = FALSE)
}
