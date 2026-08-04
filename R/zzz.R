.resetCacheIfNeeded <- function() {
  .wd <- rxode2::rxTempDir()
  if (.wd != "") {
    .md5File <- file.path(.wd, "nlmixr2est.md5")
    if (file.exists(.md5File)) {
      .md5 <- readLines(.md5File)
      if (.md5 != nlmixr2.md5) {
        ## Deliberately NOT rxClean() here.  It wipes rxTempDir(), and deleting a
        ## generated model's compiled artifact out from under a live model object makes
        ## rxode2's deferred-compile thunk rebuild it -- emitting different code than
        ## the build it replaces, which silently corrupts that model (see
        ## nlmixr2/rxode2#1171 and .nlmixr2estRxode2()).  Every model nlmixr2est
        ## generates now carries a role-tagged, content-hashed artifact name and builds
        ## outside rxTempDir(), so an upgrade cannot make this package reuse a stale
        ## artifact and there is nothing here to clean.
        writeLines(nlmixr2.md5, .md5File)
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

## DO NOT re-add an rxode2 API version/name check here.
##
## There used to be one: the NAMES of rxode2's function-pointer table were captured
## into `rxode2.api` when nlmixr2est was compiled, and .onLoad refused to start
## ("nlmixr2est needs a different version of rxode2 api") whenever the live table's
## names did not match that snapshot.
##
## It is the wrong check, and it was the approach that caused the ABI problems the
## current design was requested to solve.  The table is read POSITIONALLY and is
## ADDITIVE -- new entry points are appended, so an older consumer reads the lower
## indices and ignores the rest, which is what makes an rxode2 upgrade safe without
## rebuilding every downstream package.  Names are documentation, not the contract:
## rxode2 5.1.6 corrected two drifted slot LABELS (slot 8 was mislabeled
## `getSolvingOptionsInd` when it is `getSolvingOptions`, slot 9 named an export that
## no longer exists) and the guard then refused a combination that was completely
## sound, taking down every fresh-session test with it.
##
## The reverse case is not a reason to bring it back either: a downstream package
## built against MORE slots than the installed rxode2 provides is caught where it
## matters, by rxode2's own build-time check that no slot is left unset and by the
## DESCRIPTION version requirement -- not by string-comparing a snapshot at load.
.iniRxode2Ptr <- function() {
  .Call(`_nlmixr2est_iniRxodePtrs`, rxode2::.rxode2ptrs(),
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
#' Initialize S3 methods
#'
#' @return nothing, called for side effects
#' @export
#' @keywords internal
#' @author Matthew L. Fidler
#' @examples
#'
#' # Lets tools like nlmixr2shiny import S3 registration without
#' # depending on nlmixr2est
#'
#' .iniS3() # run to register S3 methods
#'
.iniS3 <- function() {
  if (requireNamespace("generics", quietly = TRUE)) {
    rxode2::.s3register("generics::tidy", "nlmixr2FitCore")
    rxode2::.s3register("generics::tidy", "nlmixr2FitCoreSilent")
    rxode2::.s3register("generics::glance", "nlmixr2FitCore")
    rxode2::.s3register("generics::glance", "nlmixr2FitCoreSilent")
    rxode2::.s3register("generics::augment", "nlmixr2FitCore")
    rxode2::.s3register("generics::augment", "nlmixr2FitCoreSilent")
  }
  if (exists("dim.rxEt", envir = asNamespace("nlmixr2est"), inherits = FALSE)) {
    rxode2::.s3register("base::dim", "rxEt")
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
  rxode2::.s3register("rxode2::rxUiDeparse", "agqControl")
  rxode2::.s3register("rxode2::rxUiDeparse", "laplaceControl")
  rxode2::.s3register("rxode2::rxUiGet", "foceiOuter")
  rxode2::.s3register("rxode2::rxUiGet", "impmapThetaSens")
  .resetCacheIfNeeded()
  .Call(`_rxode2version4`, as.integer(utils::packageVersion("rxode2") >= "4.0.0"))
}

.onLoad <- function(libname, pkgname) {
  .nlmixr2globalReset(TRUE)
  backports::import(pkgname)
  .iniPtrs()
  .iniS3()
}

## Stamp type/description attrs onto the built-in nlmixr2Est.* methods at
## namespace-build time (bindings not yet locked, so no unlockBinding needed).
.nlmixr2EstTypeApply(environment())

compiled.rxode2.md5 <- rxode2::rxMd5()

.onAttach <- function(libname, pkgname) {
  ## nocov start
  ## Setup rxode2.prefer.tbl
  .nlmixr2globalReset(TRUE)
  .iniPtrs()
  .iniS3()
  ## nocov end
}
