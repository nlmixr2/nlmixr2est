#' Get abberviate package information for nlmixr2 fit object
#'
#' @param pkg package to add
#' @return list of package information
#' @noRd
#' @author Matthew L. Fidler
.pkgInfo <- function(pkg) {
  .pkg <- suppressWarnings(utils::packageDescription(pkg))
  if (length(.pkg) == 1L && is.na(.pkg)) {
    return(list(
      Package = pkg,
      Version = NA_character_,
      dev = NA,
      installed = FALSE,
      install = NA_character_))
  }
  .ret <- list(
    Package = pkg,
    Version = .pkg$Version)
  if (!is.null(.pkg$GithubUsername)) {
    .ret <- c(.ret,
              list(
                dev=TRUE,
                installed = TRUE,
                install=deparse1(bquote(remotes::install_github(
                  .(paste0(.pkg$GithubUsername,"/", .pkg$GithubRepo)),
                  ref = .(.pkg$GithubSHA1))))))

  } else {
    .ret <- c(.ret,
              list(
                dev=FALSE,
                installed = TRUE,
                install=deparse1(bquote(remotes::install_version(.(pkg), version=.(.pkg$Version))))))
  }
  class(.ret) <- "nlmixr2estPkgInfo"
  .ret
}

.sessionInfoEnv <- new.env(parent=emptyenv())
.sessionInfoEnv$pkg <- c("dparser",
                         "lotri",
                         "PreciseSums",
                         "rxode2ll",
                         "rxode2",
                         "lbfgsb3c",
                         "n1qn1",
                         "nlmixr2est",
                         "nlmixr2extra",
                         "nlmixr2lib",
                         "nlmixr2",
                         "nonemem2rx",
                         "monolix2rx",
                         "babelmixr2",
                         "PopED",
                         "PKNCA",
                         "lotri",
                         "nlmixr2data",
                         "nlmixr2est",
                         "nlmixr2extra",
                         "nlmixr2plot",
                         "rxode2",
                         "ggPMX",
                         "shinyMixR",
                         "xpose.nlmixr2")

#' Adds a package to the nlmixr2's $sessioninfo inside the fit
#'
#'
#' @param pkg character vector of the package to add
#' @return nothing, called for side effects
#' @export
#' @author Matthew L. Fidler
#' @keywords internal
#' @examples
#' .addPkgNlmixr2("nlmixr2") # already present
.addPkgNlmixr2 <- function(pkg) {
  .sessionInfoEnv$pkg <- unique(c(.sessionInfoEnv$pkg, pkg))
}
#' Create the sessionInfo for nlmixr2 fit
#'
#' @return a nlmixr2 abbreviated session information  object
#' @noRd
#' @author Matthew L. Fidler
.sessionInfo <- function() {
  # Get the minimum information for the nlmixr2 related packages
  .extra <- character(0)

  .ret <- setNames(lapply(.sessionInfoEnv$pkg, .pkgInfo), .sessionInfoEnv$pkg)

  .os <- suppressWarnings(utils::sessionInfo("base")$running)
  if (is.null(.os))
    return(NA_character_)
  .os <- gsub("Service Pack", "SP", .os)
  if (is.null(.os)) {
    .os <- NA_character_
  }
  if (!is.na(.os)) {
    .extra <- paste0(.extra, paste0("# OS: ", .os, "\n"))
  }
  .la <- tryCatch(base::La_library(), error = function(err) NA_character_)
  if (!is.na(.la)) {
    .extra <- paste0(.extra, paste0("# LAPACK: ", .la, "\n"))
  }
  .lv <- tryCatch(base::La_version(), error = function(err) NA_character_)
  if (!is.na(.lv)) {
    .extra <- paste0(.extra, paste0("# LAPACK Version: ", .lv, "\n"))
  }
  .r <- R.version.string
  if (!is.na(.r)) {
    .extra <- paste0(.extra, paste0("# R Version: ", .r, "\n"))
  }
  class(.ret) <- "nlmixr2estSessionInfo"
  attr(.ret, "extra") <- .extra
  .ret
}

#' @export
print.nlmixr2estSessionInfo <- function(x, ...) {
  cat("## ==============================\n")
  cat("## nlmixr2est Session Information\n")
  cat("## ==============================\n")
  cat(paste(attr(x, "extra"), collapse="\n"), sep="\n")
  for (pkg in names(x)) {
    print.nlmixr2estPkgInfo(x[[pkg]])
  }
  invisible()
}

#' @export
print.nlmixr2estPkgInfo <- function(x, ...) {
  if (x$installed) {
    cat("\n# Install ")
    if (x$dev) {
      cat("Development version of '", x$Package, "' from GitHub (shows ver ", x$Version, ")\n", sep="")
      cat(x$install, "\n")
    } else {
      cat("Package version ", x$Version, " of '", x$Package, "'\n", sep="")
      cat(x$install, "\n")
    }
  } else {
    cat("\n# Package '", x$Package, "' is not installed, but known to enhance nlmixr2/babelmixr2\n", sep="")
  }
  invisible(x)
}
