## This is only for nlmixr
for (f in c("src/RcppExports.cpp", "inst/include/nlmixr2est_RcppExports.h")) {
  l <- readLines(f)
  w <- which(regexpr("^[#]include <RcppArmadillo.h>", l) != -1)
  if (length(w) == 1) {
    l <- l[-w]
    message("Excluding RcppArmadillo from ", f)
    writeLines(l, f)
  }
}

unlink("R/nlmixr2_md5.R")

cpp <- list.files("src", pattern = ".(c|h|cpp|f)$")
include <- list.files("inst/include")
Rfiles <- list.files("R/", pattern = ".R")
md5 <- digest::digest(lapply(c(paste0("src/", cpp),
                               paste0("inst/include/", include),
                               paste0("R/", Rfiles)), digest::digest, file = TRUE))

md5file <- file("R/nlmixr2_md5.R", "wb")
writeLines(sprintf("nlmixr2.md5 <- \"%s\"\n", md5), md5file)
close(md5file)

.in <- suppressWarnings(readLines("src/Makevars.in"))

.in <- gsub("@ARMA@", file.path(find.package("RcppArmadillo"),"include"), .in)
.in <- gsub("@BH@", file.path(find.package("BH"),"include"), .in)
.in <- gsub("@RCPP@", file.path(find.package("Rcpp"),"include"), .in)
.in <- gsub("@RXP@", file.path(find.package("rxode2"),"include"), .in)
.in <- gsub("@EG@", file.path(find.package("RcppEigen"),"include"), .in)

## COMPATIBILITY PROBE -- remove with the layers it gates (see the follow-up issue).
##
## rxode2 5.1.6 added setIndCmt() and the event-sensitivity shape save/restore entry
## points to the linked function-pointer table.  nlmixr2est uses them when they are
## there and falls back to what it did before when they are not, so one source tree
## builds against both 5.1.5 and 5.1.6 -- which is what lets the two packages be
## submitted to CRAN at the same time.
##
## Probe the HEADER rather than packageVersion(): the header is what the compile
## actually sees, and a source install of an in-between rxode2 can carry either.
.rxh <- file.path(find.package("rxode2"), "include", "rxode2ptr.h")
.rxsrc <- if (file.exists(.rxh)) paste(readLines(.rxh, warn = FALSE), collapse = "\n") else ""
.hasCmtW <- grepl("setIndCmt_t", .rxsrc, fixed = TRUE)
.hasShape <- grepl("rxode2EventSensShapeSave_t", .rxsrc, fixed = TRUE)
.compat <- character(0)
if (.hasCmtW)  .compat <- c(.compat, "-DNLMIXR2EST_HAS_SETINDCMT")
if (.hasShape) .compat <- c(.compat, "-DNLMIXR2EST_HAS_ESSHAPE")
message("nlmixr2est: rxode2 API probe -- setIndCmt=", .hasCmtW, " esShape=", .hasShape)
.in <- gsub("@RXCOMPAT@", paste(.compat, collapse = " "), .in)


if (.Platform$OS.type == "windows") {
  .makevars <- file("src/Makevars.win", "wb")
  .i <- "I"
} else {
  .makevars <- file("src/Makevars", "wb")
  if (any(grepl("Pop!_OS", utils::osVersion, fixed=TRUE)) ||
        any(grepl("Ubuntu", utils::osVersion, fixed=TRUE))) {
    .i <- "isystem"
  } else {
    .i <- "I"
  }
}

writeLines(gsub("@ISYSTEM@", .i, .in),
           .makevars)
close(.makevars)
