.h2 <- function(x) {
  cli::cli_text(crayon::bold(paste0(cli::symbol$line, cli::symbol$line, " ", x, " ", cli::symbol$line, cli::symbol$line)))
}

.getR <- function(x, sd = FALSE) {
  if (is.null(x)) {
    return(x)
  }
  .rs <- x
  .lt <- lower.tri(.rs)
  .dn1 <- dimnames(x)[[2]]
  .nms <- apply(
    which(.lt, arr.ind = TRUE), 1,
    function(x) {
      sprintf(
        "cor%s%s", getOption("broom.mixed.sep1", "__"),
        paste(.dn1[x], collapse = getOption("broom.mixed.sep2", "."))
      )
    }
  )
  .lt <- structure(.rs[.lt], .Names = .nms)
  .lt <- .lt[.lt != 0]
  if (sd) {
    .d <- dim(x)
    if (.d[1] > 0) {
      .lt <- c(setNames(diag(x), paste0("sd", getOption("broom.mixed.sep1", "__"), .dn1)), .lt)
    }
  }
  # Omit NA values that may be present with one eta fixed as zero and other etas
  # with estimated correlations.
  .lt[!is.na(.lt)]
}

.getCorPrint <- function(x) {
  .op <- options()
  on.exit(options(.op))
  options(
    "broom.mixed.sep1" = ":",
    "broom.mixed.sep2" = ","
  )
  .strong <- getOption("nlmixr2.strong.corr", 0.7)
  .moderate <- getOption("nlmixr2.moderate.corr", 0.3)
  .lt <- .getR(x)
  .digs <- 3
  .lts <- sapply(.lt, function(x) {
    .x <- abs(x)
    .ret <- "<"
    if (.x > .strong) {
      .ret <- ">" ## Strong
    } else if (.x > .moderate) {
      .ret <- "=" ## Moderate
    }
    return(.ret)
  })
  .nms <- names(.lt)
  .lt <- sprintf("%s%s", formatC(signif(.lt, digits = .digs), digits = .digs, format = "fg", flag = "#"), .lts)
  names(.lt) <- .nms
  .lt <- gsub(rex::rex("\""), "", paste0("    ", .captureOutput(print(.lt))))
  if (crayon::has_color()) {
    .lt <- gsub(rex::rex(capture(.regNum), ">"), "  \033[1m\033[31m\\1 \033[39m\033[22m", .lt, perl = TRUE)
    .lt <- gsub(rex::rex(capture(.regNum), "="), "  \033[1m\033[32m\\1 \033[39m\033[22m", .lt, perl = TRUE)
    .lt <- gsub(rex::rex(capture(.regNum), "<"), "  \\1  ", .lt, perl = TRUE)
  } else {
    .lt <- gsub(rex::rex(capture(.regNum), or(">", "=", "<")), "  \\1 ", .lt, perl = TRUE)
  }
  cat(paste(.lt, collapse = "\n"), "\n\n")
}

##' Print data-frame for Rstudio notebooks
##'
##' @param x Data frame
##' @param name Name of data  frame
##' @param bound What the nlmixr2 object is bound to
##' @return If this is printing to a notebook object
##' @author Matthew Fidler
##' @noRd
.pagedPrint <- function(x, name, bound) {
  if (is.null(x)) {
    return(invisible())
  }
  access <- as.character(substitute(x))
  if (length(access) == 3) {
    if (access[1] == "$") {
      access <- paste0("$", access[3])
    }
  } else {
    access <- ""
  }
  .df <- x
  if (inherits(.df, "matrix")) {
    .df <- as.data.frame(.df)
  } else if (inherits(.df, "character")) {
    .df <- data.frame(Information = .df)
  }
  .cls <- c(
    paste0(bound, access, ": ", gsub(" +", " ", name)),
    "paged_df", "data.frame")
  class(.df) <- .cls
  return(length(utils::capture.output(print(.df))) == 0)
}

##' @export
print.nlmixr2Class <- function(x, ...) {
  tmp <- x
  attr(tmp, ".foceiEnv") <- NULL
  class(tmp) <- NULL
  print(tmp)
  return(invisible(x))
}

##' @export
print.nlmixr2FitCoreSilent <- function(x, ...) {
  return(invisible(x))
}

##' @export
print.nlmixr2LstSilent <- function(x, ...) {
  return(invisible(x))
}


##' @export
print.nlmixr2Gill83 <- function(x, ...) {
  cat(sprintf(
    "Gill83 Derivative/Forward Difference\n  (rtol=%s; K=%s, step=%s, ftol=%s)\n\n",
    x$gillRtol, x$gillK, x$gillStep, x$gillFtol
  ))
  NextMethod(x)
}

##' @export
print.nlmixr2FitCore <- function(x, ...) {
  .parent <- parent.frame(2)
  .bound <- do.call("c", lapply(ls(.parent), function(.cur) {
    if (identical(.parent[[.cur]], x)) {
      return(.cur)
    }
    return(NULL)
  }))
  if (length(.bound) == 0) {
    .bound <- ""
  } else if (length(.bound) >= 2) {
    .bound <- .bound[order(sapply(.bound, nchar))]
    if (.bound[1] == "x") {
      .bound <- .bound[-1]
    }
    .bound <- .bound[1]
  } else if (.bound == "x") {
    .parent <- globalenv()
    .bound2 <- do.call("c", lapply(ls(.parent), function(.cur) {
      if (identical(.parent[[.cur]], x)) {
        return(.cur)
      }
      return(NULL)
    }))
    if (length(.bound2) > 0) {
      .bound <- .bound2[order(sapply(.bound2, nchar))]
      .bound <- .bound[1]
    }
  }
  .nb <- TRUE
  if (!is.na(get("objective", x$env))) {
    .nb <- .pagedPrint(x$objDf, "Objective", .bound)
  }
  if (.nb) .nb <- .pagedPrint(x$time, "Time (sec)", .bound)
  if (.nb) {
    .pagedPrint(x$parFixedDf, "Pop. Pars", .bound)
    .pagedPrint(x$omega, "BSV Cov", .bound)
    .pagedPrint(x$omegaR, "BSV Corr", .bound)
    .pagedPrint(x$shrink, "Dist. Stats", .bound)
    .pagedPrint(x$notes, "Fit notes", .bound)
    .pagedPrint(x, "Fit Data", .bound)
  } else {
    .width <- getOption("width")
    .parent <- parent.frame(2)

    cat(cli::cli_format_method({
      .h2(paste0(
        crayon::bold$blue("nlmix"),
        crayon::bold$red(paste0("r", ifelse(use.utf(), "\u00B2", "2"))), " ",
        crayon::bold(ifelse(any(x$ui$predDf$distribution != "norm"), "log-likelihood ", "")),
        crayon::bold$yellow(x$method),
        x$extra, x$posthoc
      ))
    }), sep = "\n")
    cat("\n")
    if (length(.bound) == 0) {
      .bound <- ""
    } else if (length(.bound) >= 2) {
      .bound <- .bound[order(sapply(.bound, nchar))]
      if (.bound[1] == "x") {
        .bound <- .bound[-1]
      }
      .bound <- .bound[1]
    } else if (.bound == "x") {
      .parent <- globalenv()
      .bound2 <- do.call("c", lapply(ls(.parent), function(.cur) {
        if (identical(.parent[[.cur]], x)) {
          return(.cur)
        }
        return(NULL)
      }))
      if (length(.bound2) > 0) {
        .bound <- .bound2[order(sapply(.bound2, nchar))]
        .bound <- .bound[1]
      }
    }
    if (is.na(get("objective", x$env))) {
      cat(sprintf(
        " Gaussian/Laplacian Likelihoods: AIC(%s) or %s etc.",
        crayon::yellow(.bound),
        paste0(crayon::yellow(.bound), crayon::bold$blue("$objf"))
      ), "\n")
      cat(sprintf(
        " FOCEi CWRES & Likelihoods: addCwres(%s)",
        crayon::yellow(.bound)
      ), "\n")
    } else {
      print(x$objDf)
    }
    cat("\n")
    .fmt3("Time", .bound, "time", "sec ")
    cat("\n")
    print(x$time)
    cat("\n")
    .boundChar <- nchar(.bound)
    .tmp <- x$omega
    .noEta <- (dim(.tmp)[1] == 0)

    .populationParameters <- ifelse(.noEta, "Parameters", "Population Parameters")
    if (2 * .boundChar + 54 < .width) {
      .fmt3(.populationParameters, .bound, c("parFixed", "parFixedDf"))
    } else if (.boundChar + 54 < .width) {
      .fmt3(.populationParameters, .bound, c("parFixed", "parFixedDf"),
        on = c(TRUE, FALSE)
      )
    } else {
      .fmt3(.populationParameters, .bound, c("parFixed", "parFixedDf"),
        on = c(FALSE, FALSE)
      )
    }
    cat("\n")
    .file <- raw(0L)
    .pf <- .captureOutput(print(x$parFixed))
    if (crayon::has_color()) {
      .pf <- gsub(rex::rex(capture(.regNum), "%>"), "\033[1;31m\\1%\033[0m ", .pf, perl = TRUE)
      .pf <- gsub(rex::rex(capture(.regNum), "%="), "\033[1;32m\\1%\033[0m ", .pf, perl = TRUE)
      .pf <- gsub(rex::rex(capture(.regNum), "="), "\033[1;32m\\1\033[0m ", .pf, perl = TRUE)
      .pf <- gsub(rex::rex(capture(.regNum), "%<"), "\\1% ", .pf, perl = TRUE)
      .tmp <- c(row.names(x$parFixed), names(x$parFixed))
      .tmp <- .tmp[order(-sapply(.tmp, nchar))]
      .pf <- gsub(rex::rex(boundary, capture(or(.tmp)), boundary), "\033[1m\\1\033[0m", .pf, perl = TRUE)
      .pf <- gsub(rex::rex(capture(or(.tmp))), "\033[1m\\1\033[0m", .pf, perl = TRUE)
      .pf <- gsub(rex::rex("FIXED"), "\033[1;32mFIXED\033[0m", .pf, perl = TRUE)
      .pf <- gsub(rex::rex("fix(", capture(.regNum), ")"), "\033[1;32mfix(\\1)\033[0m", .pf, perl = TRUE)
    } else {
      .pf <- gsub(rex::rex(capture(.regNum), "%", or(">", "=", "<")), "\\1% ", .pf, perl = TRUE)
      .pf <- gsub(rex::rex(capture(.regNum), "="), "\\1 ", .pf, perl = TRUE)
    }
    cat(paste(.pf, collapse = "\n"), "\n")
    ## Correlations
    .covMethod <- x$covMethod
    if (!checkmate::testCharacter(.covMethod, len=1)) .covMethod <- ""
    if (.covMethod != "") {
      cat(paste0(
        "  Covariance Type (", crayon::yellow(.bound), crayon::bold$blue("$covMethod"), "): ",
        crayon::bold(x$covMethod), "\n"
      ))
    }
    if (exists("covList", x$env)) {
      cat("    other calculated covs (", crayon::bold$blue("setCov()"), "): ",
        paste(crayon::bold(names(x$env$covList)), collapse = ", "),
        "\n",
        sep = ""
      )
    }
    if (exists("cor", x$env)) {
      .tmp <- .getR(x$cor)
      if (any(abs(.tmp) >= getOption("nlmixr2.strong.corr", 0.7))) {
        cat(paste0("  Some strong fixed parameter correlations exist (", crayon::yellow(.bound), crayon::bold$blue("$cor"), ") :\n"))
        .getCorPrint(x$cor)
      } else {
        cat(paste0("  Fixed parameter correlations in ", crayon::yellow(.bound), crayon::bold$blue("$cor"), "\n"))
      }
    }
    .tmp <- x$omega
    if (!is.null(.tmp)) {
      diag(.tmp) <- 0
      if (!.noEta) {
        if (all(.tmp == 0)) {
          cat("  No correlations in between subject variability (BSV) matrix\n")
        } else {
          cat("  Correlations in between subject variability (BSV) matrix:\n")
          .getCorPrint(x$omegaR)
        }
        if (.boundChar * 2 + 70 < .width && !.noEta) {
          cat(paste0("  Full BSV covariance (", crayon::yellow(.bound), crayon::bold$blue("$omega"), ") or correlation (", crayon::yellow(.bound), crayon::bold$blue("$omegaR"), "; diagonals=SDs)"), "\n")
        } else if (!.noEta) {
          if (.boundChar + 43 < .width) {
            cat(paste0("  Full BSV covariance (", crayon::yellow(.bound), crayon::bold$blue("$omega"), ")"), "\n")
            cat(paste0("    or correlation (", crayon::yellow(.bound), crayon::bold$blue("$omegaR"), "; diagonals=SDs)", "\n"))
          } else {
            cat(paste0("  Full BSV covariance (", crayon::bold$blue("$omega"), ")\n"))
            cat(paste0("    or correlation (", crayon::bold$blue("$omegaR"), "; diagonals=SDs)\n"))
          }
        }
      }
      if (.boundChar + 74 < .width && !.noEta) {
        cat(paste0(
          "  Distribution stats (mean/skewness/kurtosis/p-value) available in ",
          crayon::yellow(.bound), crayon::bold$blue("$shrink")
        ), "\n")
      } else if (!.noEta) {
        cat(paste0(
          "  Distribution stats (mean/skewness/kurtosis/p-value) available in ",
          crayon::bold$blue("$shrink")
        ), "\n")
      }
    }

    if (length(x$runInfo) > 0) {
      cat(paste0("  Information about run found (", crayon::yellow(.bound),
                 crayon::bold$blue("$runInfo"), "):\n"))
      lapply(x$runInfo, function(msg) {
        cat("  ", cli::cli_format_method({
          cli::cli_li(msg)
        }), "\n")
      })
    }
    cat(paste0(
      "  Censoring (", crayon::yellow(.bound), crayon::bold$blue("$censInformation"), "): ",
      as.character(x$censInformation), "\n"
    ))
    .msg <- x$message
    if (length(.msg) >= 1L) {
      if (length(.msg) == 1L && x$message == "") {
        .msg <- NULL
      } else {
        .msg <- suppressWarnings(try(gsub("^ *", "", .msg), silent=TRUE))
        if (inherits(.msg, "try-error")) {
          .msg <- "    $message cannot be displayed, examine manually"
        } else {
          .msg <- paste(paste0("    ", .msg), collapse="\n")
        }
      }
    }
    if (length(.msg) > 0L) {
      cat(paste0("  Minimization message (", crayon::yellow(.bound), crayon::bold$blue("$message"), "): "), "\n")
      cat(.msg, "\n")
      if (x$message[1] == "false convergence (8)") {
        cat("  In an ODE system, false convergence may mean \"useless\" evaluations were performed.\n")
        cat("  See https://tinyurl.com/yyrrwkce\n")
        cat("  It could also mean the convergence is poor, check results before accepting fit\n")
        cat("  You may also try a good derivative free optimization:\n")
        cat("    nlmixr2(...,control=list(outerOpt=\"bobyqa\"))\n")
      }
    }
    if (rxode2::rxIs(x, "nlmixr2FitData")) {
      .dfName <- "data.frame"
      if (rxode2::rxIs(x, "tbl")) .dfName <- "tibble"
      if (rxode2::rxIs(x, "data.table")) .dfName <- "data.table"
      cat("\n")
      cat(cli::cli_format_method({
        .h2(paste0(
          crayon::bold("Fit Data"),
          " (object",
          ifelse(.bound == "", "", " "),
          crayon::yellow(.bound),
          " is a modified ",
          crayon::blue(.dfName), "):"
        ))
      }), sep = "\n")
      if (rxode2::rxIs(x, "tbl") || rxode2::rxIs(x, "data.table")) {
        .oldOpts <- options("tibble.print_max", "tibble.print_min")
        on.exit(options(
          tibble.print_max = .oldOpts$tibble.print_max,
          tibble.print_min = .oldOpts$tibble.print_min
        ))
        options(tibble.print_max = 3, tibble.print_min = 3)
        NextMethod()
        options(
          tibble.print_max = .oldOpts$tibble.print_max,
          tibble.print_min = .oldOpts$tibble.print_min
        )
      } else {
        print(head(x))
      }
    }
  }

  return(invisible(x))
}


.notesFit <- function(x) {
  .parent <- globalenv()
  .bound2 <- do.call("c", lapply(ls(.parent), function(.cur) {
    if (identical(.parent[[.cur]], x)) {
      return(.cur)
    }
    return(NULL)
  }))
  if (length(.bound2) > 0) {
    .bound <- .bound2[order(sapply(.bound2, nchar))]
    .bound <- .bound[1]
  } else {
    .bound <- ""
  }
  .c <- NULL
  if (x$covMethod != "") {
    .c <- c(.c, paste0(
      "  Covariance Type (", .bound, "$covMethod): ",
      x$covMethod
    ))
  }
  if (is.na(get("objective", x$env))) {
    .c <- c(
      .c,
      "Missing Objective function; Can add by:",
      sprintf(
        " Gaussian/Laplacian Likelihoods: AIC(%s) or %s etc.",
        .bound, .bound, "$objf"
      ),
      sprintf(
        " FOCEi CWRES & Likelihoods: addCwres(%s)",
        .bound
      )
    )
  }
  .c <- c(.c, paste0(
    "  Censoring: ",
    as.character(x$censInformation)))
  if (length(x$runInfo) > 0) {
    .c <- c(.c, paste0("  Information about run found in (", .bound, "$runInfo):"),
            paste0("    ", x$runInfo))
  }
  .msg <- x$message
  if (length(.msg) >= 1L) {
    if (length(.msg) == 1L && x$message == "") {
      .msg <- NULL
    } else {
      .msg <- gsub("^ *", "", .msg)
      .msg <- paste(paste0("    ", .msg), collapse="\n")
    }
  }
  if (length(.msg) > 0L) {
    .c <- c(
      .c, paste0("  Minimization message (", .bound, "$message): "),
      paste0("    ", x$message)
    )
    if (.msg == "false convergence (8)") {
      .c <- c(
        .c, "  In an ODE system, false convergence may mean \"useless\" evaluations were performed.",
        "  See https://tinyurl.com/yyrrwkce",
        "  It could also mean the convergence is poor, check results before accepting fit",
        "  You may also try a good derivative free optimization:",
        "    nlmixr2(...,control=list(outerOpt=\"bobyqa\"))"
      )
    }
  }
  .c
}


.fmt3 <- function(name, bound, access, extra = "",
                  on = c(TRUE, TRUE)) {
  if (length(access) == 1) {
    on <- on[1]
  }
  cat(cli::cli_format_method({
    .h2(paste0(
      crayon::bold(name), " (", extra,
      paste(crayon::bold$blue(paste0(ifelse(on, crayon::yellow(bound), ""), "$", access)), collapse = " or "), "):"
    ))
   }), sep = "\n")
}
