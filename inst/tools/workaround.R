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

if (.Platform$OS.type == "windows" && !file.exists("src/Makevars.win")) {
  writeLines(gsub("@ISYSTEM@", "I",
                  gsub("@CXX14STD@", "CXX14STD = -std=c++1y", .in)),
             "src/Makevars.win")
} else {
  writeLines(gsub("@ISYSTEM@", "isystem",
                  gsub("@CXX14STD@", "CXX14STD = -std=gnu++14", .in)),
             "src/Makevars")
}
