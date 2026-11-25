# Check the NONMEM control streams and Monolix projects without a licence:
# nonmem2rx parses each control stream, and monolix2rx parses each project and
# its model file (it stops at the event endpoint, which it does not translate).
#
#   Rscript checkSyntax.R [validation-dir]
.args <- commandArgs(trailingOnly = TRUE)
.dir <- if (length(.args) > 0) {
  .args[1]
} else {
  .f <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
  if (length(.f) == 1) dirname(normalizePath(.f)) else getwd()
}
.bad <- character(0)
for (.f in list.files(file.path(.dir, "nonmem"), "[.]ctl$", full.names = TRUE)) {
  .x <- try(
    suppressMessages(nonmem2rx::nonmem2rx(.f, determineError = FALSE, load = FALSE, validate = FALSE)),
    silent = TRUE
  )
  if (inherits(.x, "try-error")) .bad <- c(.bad, paste(basename(.f), conditionMessage(attr(.x, "condition"))))
}
.old <- setwd(file.path(.dir, "monolix"))
for (.f in list.files(".", "[.]mlxtran$")) {
  .x <- try(suppressMessages(suppressWarnings(monolix2rx::monolix2rx(.f, update = FALSE))), silent = TRUE)
  .msg <- if (inherits(.x, "try-error")) conditionMessage(attr(.x, "condition")) else ""
  if (!grepl("'event' endpoint not supported", .msg, fixed = TRUE)) .bad <- c(.bad, paste(.f, .msg))
}
setwd(.old)
if (length(.bad) > 0) {
  stop(paste(.bad, collapse = "\n"), call. = FALSE)
}
message("all control streams and projects parse")
