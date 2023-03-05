.resetCacheIfNeeded <- function() {
  .wd <- rxode2::rxTempDir()
  if (.wd != "") {
    .md5File <- file.path(.wd, "nlmixr2est.md5")
    if (file.exists(.md5File)) {
      .md5 <- readLines(.md5File)
      if (.md5 != nlmixr2.md5) {
        packageStartupMessage("detected new version of nlmixr2est, cleaning rxode2 cache")
        rxode2::rxClean()
      }
    } else {
      writeLines(nlmixr2.md5, .md5File)
    }
  }
}

.onLoad <- function(libname, pkgname) {
  backports::import(pkgname)
  if (requireNamespace("generics", quietly = TRUE)) {
    rxode2::.s3register("generics::tidy", "nlmixr2FitCore")
    rxode2::.s3register("generics::tidy", "nlmixr2FitCoreSilent")
    rxode2::.s3register("generics::glance", "nlmixr2FitCore")
    rxode2::.s3register("generics::glance", "nlmixr2FitCoreSilent")
    rxode2::.s3register("generics::augment", "nlmixr2FitCore")
    rxode2::.s3register("generics::augment", "nlmixr2FitCoreSilent")
  }
  rxode2::.s3register("rxode2::rxModelVarsS3", "nlmixr2FitCore")
  rxode2::.s3register("rxode2::rxModelVarsS3", "nlmixr2FitCoreSilent")
  rxode2::.s3register("rxode2::getBaseSimModel", "nlmixr2FitCoreSilent")
  rxode2::.s3register("rxode2::getBaseSimModel", "nlmixr2FitCore")
  rxode2::.s3register("rxode2::getBaseSimModel", "nlmixr2FitData")
  .resetCacheIfNeeded()
}

compiled.rxode2.md5 <- rxode2::rxMd5()

.onAttach <- function(libname, pkgname) {
  ## nocov start
  ## Setup rxode2.prefer.tbl
  if (compiled.rxode2.md5 != rxode2::rxMd5()) {
    stop("nlmixr2 compiled against different version of rxode2, cannot run nlmixr2\ntry `install.packages(\"nlmixr2\", type = \"source\")` to recompile", call.=FALSE)
  }
  requireNamespace("rxode2parse", quietly=TRUE)  
  requireNamespace("rxode2random", quietly=TRUE)
  ## nlmixr2SetupMemoize()
  ## options(keep.source = TRUE)
  ## nocov end
}
