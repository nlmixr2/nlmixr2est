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

.iniLotriPtr <- function() {
  .Call(`_nlmixr2est_iniLotriPtr`,
        lotri::.lotriPointers(),
        PACKAGE = "nlmixr2est")
}

# This will be saved when compiled
rxode2.api <- names(rxode2::.rxode2ptrs())

.iniRxode2Ptr <- function() {
  .ptr <- rxode2::.rxode2ptrs()
  .nptr <- names(.ptr)
  if (length(rxode2.api) > length(.nptr)) {
    stop("nlmixr2est requires a newer version of rxode2 api, cannot run nlmixr2est\ntry `install.packages(\"rxode2\")` to get a newer version of rxode2", call.=FALSE)
  } else {
    .nptr <- .nptr[seq_along(rxode2.api)]
    if (!identical(rxode2.api, .nptr)) {
      .bad <- TRUE
      stop("nlmixr2est needs a different version of rxode2 api, cannot run nlmixr2est\ntry `install.packages(\"rxode2\")` to get a newer version of rxode2, or update both packages", call.=FALSE)
    }
  }
  .Call(`_nlmixr2est_iniRxodePtrs`, .ptr,
        PACKAGE = "nlmixr2est")
}

.iniN1qn1ptr <- function() {
  .Call(`_nlmixr2est_iniN1qn1cPtrs`,
        n1qn1::.n1qn1ptr(),
        PACKAGE = "nlmixr2est")
}

.iniLbfgsb3c <- function() {
  .Call(`_nlmixr2est_iniLbfgsb3ptr`,
        lbfgsb3c::.lbfgsb3cPtr(),
        PACKAGE = "nlmixr2est")
}

.iniPtrs <- function() {
  .iniLotriPtr()
  .iniRxode2Ptr()
  .iniN1qn1ptr()
  .iniLbfgsb3c()
}

.iniS3 <- function() {
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

  rxode2::.s3register("rxode2::rxUiDeparse", "foceiControl")
  rxode2::.s3register("rxode2::rxUiDeparse", "saemControl")
  rxode2::.s3register("rxode2::rxUiDeparse", "bobyqaControl")
  rxode2::.s3register("rxode2::rxUiDeparse", "lbfgsb3cControl")
  rxode2::.s3register("rxode2::rxUiDeparse", "n1qn1Control")
  rxode2::.s3register("rxode2::rxUiDeparse", "newuoaControl")
  rxode2::.s3register("rxode2::rxUiDeparse", "nlmeControl")
  rxode2::.s3register("rxode2::rxUiDeparse", "nlminbControl")
  rxode2::.s3register("rxode2::rxUiDeparse", "nlmControl")
  rxode2::.s3register("rxode2::rxUiDeparse", "nlsControl")
  rxode2::.s3register("rxode2::rxUiDeparse", "optimControl")
  rxode2::.s3register("rxode2::rxUiDeparse", "uobyqaControl")
  rxode2::.s3register("rxode2::rxUiDeparse", "tableControl")
  .resetCacheIfNeeded()
}

.onLoad <- function(libname, pkgname) {
  .nlmixr2globalReset(TRUE)
  backports::import(pkgname)
  .iniPtrs()
  .iniS3()

}

compiled.rxode2.md5 <- rxode2::rxMd5()

.onAttach <- function(libname, pkgname) {
  ## nocov start
  ## Setup rxode2.prefer.tbl
  .nlmixr2globalReset(TRUE)
  .iniPtrs()
  .iniS3()
  ## nlmixr2SetupMemoize()
  ## options(keep.source = TRUE)
  ## nocov end
}
