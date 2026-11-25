# Reading NONMEM and Monolix results into reference files, and reading those
# back.  A reference file reference/<tool>/<scenario>.csv has one row per
# parameter (nlmixr2 names; random-effect variances, not SDs) and an `OFV`
# row holding -2 log-likelihood, with the tool version and run date.

# theta names (in ini() order) and eta names of a scenario's model
.tteParNames <- function(scenario) {
  .ini <- rxode2::rxode2(scenario$model)$iniDf
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
    .s$model <- .s$reduced(rxode2::rxode2(.s$model))
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

#' Read a Monolix run of a scenario
#'
#' @param dir the run's export directory
#' @param scenario the scenario the run fitted
#' @param version Monolix version
#' @return reference data frame
#' @noRd
tteMonolixRead <- function(dir, scenario, version = NA_character_) {
  .p <- utils::read.csv(file.path(dir, "populationParameters.txt"))
  .ll <- utils::read.csv(file.path(dir, "LogLikelihood", "logLikelihood.txt"))
  .nm <- .tteParNames(scenario)
  .mu <- rxode2::rxode2(scenario$model)$muRefDataFrame
  .seCol <- grep("^se_", names(.p), value = TRUE)[1]
  .get <- function(x, col) {
    .w <- match(x, .p$parameter)
    if (is.na(.w) || is.na(col)) NA_real_ else as.numeric(.p[[col]][.w])
  }
  # residual-error parameters keep their nlmixr2 name; Monolix names them b
  .mlx <- vapply(
    .nm$theta,
    function(x) {
      if (!is.null(scenario$monolixNames[[x]])) scenario$monolixNames[[x]] else paste0(x, "_pop")
    },
    character(1)
  )
  .est <- vapply(.mlx, .get, numeric(1), col = "value")
  .se <- vapply(.mlx, .get, numeric(1), col = .seCol)
  # omega_<theta> is the SD of the eta on <theta>: report the variance
  .om <- paste0("omega_", .mu$theta[match(.nm$eta, .mu$eta)])
  .sd <- vapply(.om, .get, numeric(1), col = "value")
  .sdSe <- vapply(.om, .get, numeric(1), col = .seCol)
  .ofv <- .ll[[2]][match("OFV", .ll[[1]])]
  data.frame(
    tool = "monolix",
    version = version,
    date = format(file.mtime(file.path(dir, "populationParameters.txt")), "%Y-%m-%d"),
    parameter = c(.nm$theta, .nm$eta, "OFV"),
    estimate = c(.est, .sd^2, .ofv),
    se = c(.se, 2 * .sd * .sdSe, NA_real_),
    row.names = NULL
  )
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
