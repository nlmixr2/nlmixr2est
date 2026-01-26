## nocov start
.genAgq <- function(n=64) {
  .c <- as.call(c(list(quote(`switch`), quote(`n`)),
                  lapply(seq_len(n), function(n) {
                    .v <- fastGHQuad::gaussHermiteData(n)
                    .v$w <- .v$w/sqrt(pi)
                    str2lang(deparse1(.v))})))
  .f <- function(n) {}
  body(.f) <- as.call(c(list(quote(`{`)),
                        str2lang(paste0("if (is.na(n)) return(", n, ")")),
                        .c))
  .f
}

.nlmixr2estbuild <- function() {
  # This function is used to build some code in the nlmixr2est package
  # it is called with devtools::document()
  message("Pregenerate gaussHermiteData for nlmixr2est")
  .tmp <- deparse(.genAgq())
  .tmp[1]<- paste0(".nlmixr2estAgq <- ", .tmp[1])
  writeLines(c("## created by .nlmixr2estbuild() in build.R edit there",
               .tmp),
             devtools::package_file("R/agqGen.R"))
  message("done")
  message("Clear testing model cache")
  .fitCacheDir <- file.path(testthat::test_path(), "fixtures")
  if (dir.exists(.fitCacheDir)) {
    unlink(list.files(.fitCacheDir, pattern=".rds$", full.names = TRUE))
  }
  message("done")
  invisible("")
}

## nocov end
