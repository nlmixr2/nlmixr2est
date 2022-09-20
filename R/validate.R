##' Validate nlmixr2
##'
##' This allows easy validation/qualification of nlmixr2 by running the
##' testing suite on your system.
##' @param type of test to be run
##' @param skipOnCran when `TRUE` skip the test on CRAN.
##' @author Matthew L. Fidler
##' @return Nothing, called for its side effects
##' @export
nlmixr2Validate <- function(type = NULL, skipOnCran=TRUE) {
  if (is(substitute(type), "{")) {
    if (isTRUE(skipOnCran)) {
      if (!identical(Sys.getenv("NOT_CRAN"), "true") ||
            !identical(Sys.getenv("nmTest"), "")) {
        return(invisible())
      }
    }
    gc()
    rxode2::rxUnloadAll()
    return(force(type))
  }

  pt <- proc.time()
  .filter <- NULL
  if (is.null(type)) type <- FALSE
  if (is.character(type)) {
    .filter <- type
    type <- TRUE
  }
  if (type == TRUE) {
    .oldCran <- Sys.getenv("NOT_CRAN")
    .oldNmTest <- Sys.getenv("nmTest")
    Sys.setenv("NOT_CRAN" = "true") # nolint
    Sys.setenv("nmTest" = "") # nolint
    on.exit(Sys.setenv("NOT_CRAN" = .oldCran, "nmTest"=.oldNmTest)) # nolint
  } else if (type == FALSE) {
    .oldCran <- Sys.getenv("NOT_CRAN")
    .oldNmTest <- Sys.getenv("nmTest")
    Sys.setenv("NOT_CRAN" = "false") # nolint
    Sys.setenv("rxTest" = "false") # nolint
    on.exit(Sys.setenv("NOT_CRAN" = .oldCran, "rxTest"=.oldNmTest)) # nolint
  }
  rxode2::.rxWithOptions(list(testthat.progress.max_fails=10000000000), {
    path <- file.path(system.file("tests", package = "nlmixr2est"), "testthat")
    rxode2::.rxWithWd(path, {
      try(testthat::test_dir(path, filter = .filter))
      message("================================================================================")
      print(proc.time() - pt)
      message("================================================================================")
    })
  })
}

##' @rdname nlmixr2Validate
##' @export
nmTest <- nlmixr2Validate
