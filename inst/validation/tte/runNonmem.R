# Fit the tte validation scenarios with NONMEM and store the reference fits.
#
#   Rscript runNonmem.R [validation-dir] [scenario ...]
#
# Needs a licensed NONMEM: the command is $NMFE (default "nmfe75"), called as
# `<nmfe> <name>.ctl <name>.lst` in a scratch directory holding the control
# stream and its data set.  Results go to reference/nonmem/<name>.csv.
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
.nmfe <- Sys.getenv("NMFE", "nmfe75")
.ctl <- sub("[.]ctl$", "", list.files(file.path(.dir, "nonmem"), "[.]ctl$"))
if (length(.which) > 0) {
  .ctl <- .ctl[sub("-reduced$", "", .ctl) %in% .which]
}
.out <- file.path(.dir, "reference", "nonmem")
dir.create(.out, recursive = TRUE, showWarnings = FALSE)
for (.name in .ctl) {
  .run <- tempfile(paste0("nm-", .name, "-"))
  dir.create(.run)
  .data <- paste0(sub("-reduced$", "", .name), ".csv")
  file.copy(file.path(.dir, "nonmem", paste0(.name, ".ctl")), .run)
  file.copy(file.path(.dir, "data", .data), .run)
  message("NONMEM: ", .name, " in ", .run)
  .old <- setwd(.run)
  .status <- system2(.nmfe, c(paste0(.name, ".ctl"), paste0(.name, ".lst")))
  setwd(.old)
  .ext <- file.path(.run, paste0(.name, ".ext"))
  if (.status != 0 || !file.exists(.ext)) {
    warning("NONMEM run failed: ", .name, call. = FALSE)
    next
  }
  .ref <- try(tteNonmemRead(.ext, file.path(.run, paste0(.name, ".lst")), .tteRefScenario(.name)))
  if (!is.data.frame(.ref)) {
    warning("could not read the NONMEM results: ", .name, call. = FALSE)
    next
  }
  utils::write.csv(.ref, file.path(.out, paste0(.name, ".csv")), row.names = FALSE)
}
