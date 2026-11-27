# Regenerate the tte validation data sets (data/<scenario>.csv).
#
#   Rscript simulate.R [validation-dir]
#
# The data sets are simulated with rxode2's tte() from the true values in
# each scenario's model, under a fixed seed per scenario.
.args <- commandArgs(trailingOnly = TRUE)
.dir <- if (length(.args) > 0) {
  .args[1]
} else {
  .f <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
  if (length(.f) == 1) dirname(normalizePath(.f)) else getwd()
}
source(file.path(.dir, "scenarios.R"))
dir.create(file.path(.dir, "data"), showWarnings = FALSE)
for (.name in names(tteValidationScenarios())) {
  .d <- tteValidationSimulate(tteValidationScenarios(.name)[[1]])
  utils::write.csv(.d, file.path(.dir, "data", paste0(.name, ".csv")), row.names = FALSE, quote = FALSE)
  .obs <- .d[.d$EVID == 0, ]
  message(sprintf(
    "%-12s %4d subjects, %5d records, %4d events",
    .name,
    length(unique(.d$ID)),
    nrow(.d),
    sum(.obs$DV == 1 & .obs$DVID == max(.obs$DVID))
  ))
}
