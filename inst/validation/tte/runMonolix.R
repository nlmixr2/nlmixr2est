# Fit the tte validation scenarios with Monolix and store the reference fits.
#
#   Rscript runMonolix.R [validation-dir] [scenario ...]
#
# Needs a licensed Monolix and its lixoftConnectors R package (the Monolix
# install directory can be given in $MONOLIX_PATH).  Each project is run in a
# copy of monolix/ and data/ in runs/monolix/ (or under $TTE_RUN_DIR), which
# keeps Monolix's output, and its estimates, standard errors and
# -2LL are read through lixoftConnectors into reference/monolix/<name>.csv.
.args <- commandArgs(trailingOnly = TRUE)
.dir <- if (length(.args) > 0 && dir.exists(.args[1])) {
  normalizePath(.args[1])
} else {
  .f <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
  if (length(.f) == 1) dirname(normalizePath(.f)) else getwd()
}
source(file.path(.dir, "reference.R"))
source(file.path(.dir, "scenarios.R"))
.which <- .args[!dir.exists(.args)]
.path <- Sys.getenv("MONOLIX_PATH")
if (.path == "") {
  lixoftConnectors::initializeLixoftConnectors(software = "monolix")
} else {
  lixoftConnectors::initializeLixoftConnectors(software = "monolix", path = .path)
}
.version <- lixoftConnectors::getLixoftConnectorsState()$version
.prj <- sub("[.]mlxtran$", "", list.files(file.path(.dir, "monolix"), "[.]mlxtran$"))
if (length(.which) > 0) {
  .prj <- .prj[sub("-reduced$", "", .prj) %in% .which]
}
.out <- file.path(.dir, "reference", "monolix")
dir.create(.out, recursive = TRUE, showWarnings = FALSE)
.run <- file.path(Sys.getenv("TTE_RUN_DIR", file.path(.dir, "runs")), "monolix")
unlink(.run, recursive = TRUE)
dir.create(.run, recursive = TRUE)
file.copy(file.path(.dir, c("monolix", "data")), .run, recursive = TRUE)
for (.name in .prj) {
  message("Monolix: ", .name, " in ", .run)
  lixoftConnectors::loadProject(file.path(.run, "monolix", paste0(.name, ".mlxtran")))
  .ok <- lixoftConnectors::runScenario()
  .ref <- if (isFALSE(.ok)) NULL else try(tteMonolixRead(.tteRefScenario(.name), .version))
  if (!is.data.frame(.ref) || anyNA(.ref$estimate)) {
    warning("Monolix run failed: ", .name, call. = FALSE)
    next
  }
  utils::write.csv(.ref, file.path(.out, paste0(.name, ".csv")), row.names = FALSE)
}
