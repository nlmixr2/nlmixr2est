# Adding a covariate to a declared distribution (plan phase 3.4).
#
# rxode2's expansion hoists each family argument onto its own role anchor,
# `rxEdA.<eta>.<role>`, so a covariate on a declaration is a term on ONE named
# line rather than a substitution buried inside a quantile call.  This adds that
# term at the source -- the declaration in `ini({})` -- so both spellings and
# every downstream consumer see it, and the anchor picks it up on expansion.
#
# The term is multiplicative on the argument: `arg * exp(beta * f(cov))`.  That
# is the one form that means something for every role, and it is why `rate` is a
# separate role from `scale`: scaling a gamma's MEAN by exp(+b*x) scales its
# RATE by exp(-b*x), so the same covariate carries the opposite sign depending
# on which the family parameterizes.  The caller (or the search) reads the sign
# off the role; this function does not try to be clever about it.

#' Coefficient name for a covariate on a declaration's role
#'
#' `beta.<eta>.<role>.<cov>` plus the shape tag when the shape is not the
#' default, mirroring the `beta.<par>.<COV>.<tag>` convention the vae covariate
#' search already uses so the two read the same in `ini()`.
#'
#' @param eta declared random effect
#' @param role the family argument's role
#' @param cov covariate column
#' @param shape covariate shape
#' @return length-one character
#' @noRd
.etaDistBetaName <- function(eta, role, cov, shape) {
  .tag <- if (identical(shape, "power")) "" else paste0(".", .vaeShapeBeta(shape))
  paste0("beta.", eta, ".", role, ".", cov, .tag)
}

#' Add a covariate to one role of one declared distribution
#'
#' @param ui rxode2 ui with at least one `dist()` declaration
#' @param eta the declared random effect's name
#' @param role which family argument to put the covariate on, by lotri role
#'   (`shape`/`rate` for a gamma, `location`/`scale` for a normal)
#' @param cov covariate column name
#' @param shape covariate shape, as `.vaeShapeExpr()` spells them
#' @param center centering value the shape needs (`power`, `lin`, `center`,
#'   the hockey arms)
#' @param est starting value for the coefficient.  NOT zero by default, and
#'   deliberately: a coefficient started at exactly 0 has no magnitude for the
#'   outer search to scale by, and measured on a known effect of 0.75 it comes
#'   back as 0.0017 from a start of 0 against 0.7353 from a start of 0.1.
#'   `foceiControl(zeroThetaRetry=)` recovers from it, but starting somewhere
#'   sane costs nothing.
#' @return the ui, with the declaration rewritten and the coefficient added
#' @noRd
.etaDistAddCovariate <- function(ui, eta, role, cov, shape = "power",
                                 center = NA_real_, est = 0.1) {
  .ui <- rxode2::rxUiDecompress(ui)
  .iniDf <- .ui$iniDf
  .w <- which(.iniDf$name == eta & !is.na(.iniDf$etaDist))
  if (length(.w) != 1L) {
    stop("'", eta, "' is not a declared random effect in this model",
         call. = FALSE)
  }
  .call <- str2lang(.iniDf$etaDist[.w])
  .fam <- as.character(.call[[1]])
  ## Feature-detected, never version-detected: the lotri that carries this
  ## reports the SAME 1.0.5 as the CRAN one that does not, so a DESCRIPTION
  ## requirement cannot express it and a version test would pass while the call
  ## still failed.  Without the guard this died on lotri's own namespace error,
  ## which names neither what was wanted nor what to install.
  .tab <- try(lotri::lotriEtaDists(), silent = TRUE)
  if (inherits(.tab, "try-error") || is.null(.tab$name)) {
    stop("adding a covariate to a dist() declaration needs a 'lotri' that ",
         "describes declared distributions, and the installed one does not ",
         "provide 'lotriEtaDists()'\n",
         "  install the development 'lotri':\n",
         "    remotes::install_github(\"nlmixr2/lotri\")",
         call. = FALSE)
  }
  .f <- which(.tab$name == .fam)
  if (length(.f) != 1L || !any(names(.tab) == "roles") ||
        !nzchar(.tab$roles[.f])) {
    stop("the installed 'lotri' has no role table for '", .fam, "'",
         call. = FALSE)
  }
  .roles <- strsplit(.tab$roles[.f], ",", fixed = TRUE)[[1]]
  .i <- which(.roles == role)
  if (length(.i) != 1L) {
    stop("'", .fam, "' has no role '", role, "'; it has: ",
         paste(.roles, collapse = ", "), call. = FALSE)
  }
  if (role %in% .etaDistRoleNoCovariate) {
    stop("role '", role, "' is a support endpoint and cannot carry a ",
         "covariate: a subject-varying bound makes the density ",
         "discontinuous in theta", call. = FALSE)
  }
  .args <- as.list(.call)[-1]
  if (.i > length(.args)) {
    stop("'", eta, "' does not supply the '", role, "' argument of '", .fam,
         "'", call. = FALSE)
  }
  .beta <- .etaDistBetaName(eta, role, cov, shape)
  if (.beta %in% .iniDf$name) {
    stop("'", .beta, "' is already in the model: ", cov,
         " is already on the '", role, "' role of dist(", eta, ")",
         call. = FALSE)
  }
  .expr <- .vaeShapeExpr(shape, cov, center = center)
  .args[[.i]] <- str2lang(paste0("(", deparse1(.args[[.i]]), ") * exp(",
                                 .beta, " * ", .expr, ")"))
  .newDecl <- deparse1(as.call(c(.call[[1]], .args)))

  ## Rebuild rather than mutate iniDf.  The coefficient is only USED inside the
  ## declaration, and a declaration does not become a model line until
  ## rxEtaDistExpand() runs -- so a theta added straight to iniDf trips the ui's
  ## "in the ini block but not in the model block" check.  The declaration and
  ## its coefficient have to be written at the same time, which is the same
  ## reason rxEtaDistMuRef() rebuilds.
  .iniTxt <- deparse(.ui$iniFun)
  ## The ui handed in may already carry the expansion, whose generated lines
  ## would then be re-emitted alongside a second expansion of the same
  ## declaration -- two sets of anchors and two decoders for one eta.  Recover
  ## the user's own model block by dropping every line the expansion writes:
  ## the copula/latent intermediates, the role anchors, and the assignment to a
  ## declared eta itself.
  .gen <- c(paste0("^rx[NTLSUuc][.]"), "^rxEdA[.]")
  .decl <- .iniDf$name[!is.na(.iniDf$etaDist)]
  .keep <- vapply(.ui$lstExpr, function(.l) {
    if (!(is.call(.l) && length(.l) > 2L && identical(.l[[1]], quote(`<-`)) &&
            is.name(.l[[2]]))) {
      return(TRUE)
    }
    .lhs <- as.character(.l[[2]])
    !(any(vapply(.gen, function(.p) grepl(.p, .lhs), logical(1))) ||
        .lhs %in% .decl)
  }, logical(1), USE.NAMES = FALSE)
  .modTxt <- vapply(.ui$lstExpr[.keep], function(.l) paste0("  ", deparse1(.l)),
                    character(1), USE.NAMES = FALSE)
  .declRe <- paste0("dist\\(", gsub("([.\\[\\]])", "\\\\\\1", eta), "\\)")
  .inIni <- grepl(.declRe, .iniTxt)
  .inMod <- grepl(.declRe, .modTxt)
  .repl <- paste0("  dist(", eta, ") ~ ", .newDecl)
  if (any(.inIni)) {
    .iniTxt[which(.inIni)[1]] <- .repl
  } else if (any(.inMod)) {
    ## the model({}) spelling declares on the parameter, not the eta
    .modTxt[which(.inMod)[1]] <- .repl
  } else {
    stop("cannot find the 'dist(", eta, ")' line to add '", cov, "' to",
         call. = FALSE)
  }
  .iniTxt <- c(.iniTxt[-length(.iniTxt)],            # drop the closing "})"
               paste0("  ", .beta, " <- ", est),
               "})")
  .txt <- paste0("function() {\n",
                 paste(.iniTxt, collapse = "\n"), "\n",
                 "model({\n", paste(.modTxt, collapse = "\n"), "\n})\n}")
  .fun <- try(eval(parse(text = .txt)), silent = TRUE)
  if (inherits(.fun, "try-error")) {
    message(.txt)
    stop("could not add '", cov, "' to the '", role, "' role of dist(", eta,
         "); the model this tried to build is echoed above", call. = FALSE)
  }
  rxode2::rxUiDecompress(rxode2::rxode2(.fun))
}
