# The per-model context and the final gradient composition of the linCmt()
# carry emission (foceiLinCmtCarry.R holds the per-row lines).  The context
# is parsed once per pair row from rx_pred_'s structural call
# (foceiLinCmtCarryPredCall.R); the composition turns the carried
# compartment sensitivity into d(rx_pred_)/d(eta) -- directly for a bare
# prediction, through the outer chain rule around rx_lcConc_ for a wrapped
# ll() likelihood (#1004).

#' The control's covariate interpolation
#'
#' @param ui rxode2 ui
#' @return "locf", "nocb", "midpoint" or "linear" (rxode2's default, locf,
#'   when the control sets none)
#' @noRd
.rxFoceiLinCmtCarryGlobalInterp <- function(ui) {
  .i <- rxode2::rxGetControl(ui, "rxControl", NULL)$covsInterpolation
  if (is.null(.i) || length(.i) != 1L) {
    return("locf")
  }
  if (is.numeric(.i)) {
    # rxControl's integer codes
    return(c("linear", "locf", "nocb", "midpoint")[as.integer(.i) + 1L])
  }
  as.character(.i)
}

#' The covariate interpolation of a carry build, for each covariate
#'
#' A covariate's own declaration in the model (`nocb(wt)`, in `ui$mv0$interp`)
#' wins, as it does in rxode2; the others take the control's.
#'
#' @param ui rxode2 ui
#' @param covs the model's covariates
#' @return named character, "locf", "nocb", "midpoint" or "linear" for each of
#'   `covs`
#' @noRd
.rxFoceiLinCmtCarryInterp <- function(ui, covs) {
  .ret <- stats::setNames(rep(.rxFoceiLinCmtCarryGlobalInterp(ui), length(covs)), covs)
  .decl <- ui$mv0$interp
  if (is.factor(.decl) && length(covs) > 0L) {
    .d <- as.character(.decl[match(covs, names(.decl))])
    .w <- !is.na(.d) & .d != "default"
    .ret[.w] <- .d[.w]
  }
  .ret
}

#' The covariates the carry reads through hidden locf and nocb copies
#'
#' Under midpoint rxode2 advances an interval with the mean of the two
#' records that carry the covariate on either side of it.  A record with no
#' value of its own (an infusion's end, an addl dose, a missing value) lies
#' between them, so neither the record nor its neighbour stands for them.
#' The carry reads them instead as the covariate's locf view at the
#' interval's start and its nocb view at its end, from two hidden copies of
#' the column (`rx_lcL_<cov>`, `rx_lcN_<cov>`) that the generated models
#' declare `locf()` and `nocb()`.
#'
#' @param ui rxode2 ui with its control
#' @return the covariates under midpoint, by the control or by their own
#'   declaration, of a `linCmt()` model the carry is on for; `character(0)`
#'   otherwise
#' @noRd
.rxFoceiLinCmtCarryMidCovs <- function(ui) {
  if (!.foceiUsesLinCmt(ui) || !.rxFoceiLinCmtCarryBuildEnabled(ui)) {
    return(character(0))
  }
  .covs <- ui$allCovs
  .covs[.rxFoceiLinCmtCarryInterp(ui, .covs) == "midpoint"]
}

#' Data with the hidden copies the carry reads under midpoint
#'
#' @param data data frame or event table with the covariate columns
#' @param covs `.rxFoceiLinCmtCarryMidCovs()`
#' @return `data` with `rx_lcL_<cov>` and `rx_lcN_<cov>`, copies of each of
#'   `covs`
#' @noRd
.rxFoceiLinCmtCarryAliasData <- function(data, covs) {
  for (.c in covs) {
    data[[paste0("rx_lcL_", .c)]] <- data[[.c]]
    data[[paste0("rx_lcN_", .c)]] <- data[[.c]]
  }
  data
}

#' The interval value of each covariate the carry reads
#'
#' rxode2 advances the linear compartments over each interval between
#' records with that interval's covariates, each read at the interval's
#' middle by its interpolation: under locf the start record's, under nocb
#' the end record's, under midpoint the mean of the two records that carry
#' it on either side.  The carry's generated lines run at the record that
#' ends the interval, so a covariate under locf or midpoint is read in the
#' advance's parameters and in g through its interval value:
#' - locf: its value at the previous record, lagged from a model copy;
#' - midpoint: the mean of its locf view at the previous record and its nocb
#'   view at this one, from the hidden copies of
#'   `.rxFoceiLinCmtCarryAliasData()`;
#' - linear: the mean of its values at the two records, its value at the
#'   middle of an interval between two records that carry it.
#'
#' The lag is taken from a model copy (`copies`, real lhs ahead of
#' rx_pred_, as lag() needs), not from the column: a record with no value of
#' its own holds NA, which lag0() turns into 0, where the model evaluation
#' fills it as rxode2 fills any NA covariate.  A subject's first record has
#' no previous one and an empty interval; its own value stands in (newind is
#' 1 only there).  lag0() reads the last lhs computed, which at a subject's
#' first record is NA only in a fresh solve: FOCEi solves a subject again for
#' each eta, and there it holds the previous pass's last record, so the first
#' record ignores it.
#'
#' @param covs the model's covariates
#' @param interp `.rxFoceiLinCmtCarryInterp()`, named by covariate
#' @return NULL when every covariate's record value is its interval's (nocb,
#'   or no covariates), else list(copies = the copy and declaration lines, to
#'   go ahead of rx_pred_, lines = the lines defining each interval value, map
#'   = covariate -> the symbol holding its interval value)
#' @noRd
.rxFoceiLinCmtCarryCovInterval <- function(covs, interp) {
  covs <- covs[interp[covs] %in% c("locf", "midpoint", "linear")]
  if (length(covs) == 0L) {
    return(NULL)
  }
  .copies <- character(0)
  .lines <- character(0)
  for (.c in covs) {
    .copy <- paste0("rx_lcCov_", .c, "_")
    .lag <- paste0("rx_lcLag_", .c, "_")
    .int <- paste0("rx_lcCovI_", .c, "_")
    .mid <- interp[[.c]] == "midpoint"
    # under midpoint the lag reads the locf copy of the column
    .src <- if (.mid) paste0("rx_lcL_", .c) else .c
    .copies <- c(
      .copies,
      if (.mid) c(paste0("locf(rx_lcL_", .c, ")"), paste0("nocb(rx_lcN_", .c, ")")),
      paste0(.copy, "=", .src)
    )
    .lines <- c(
      .lines,
      paste0(.lag, "~lag0(", .copy, ")*(newind!=1)+", .src, "*(newind==1)"),
      switch(
        interp[[.c]],
        locf = paste0(.int, "~", .lag),
        midpoint = paste0(.int, "~0.5*(", .lag, "+rx_lcN_", .c, ")"),
        linear = paste0(.int, "~0.5*(", .c, "+", .lag, ")")
      )
    )
  }
  list(
    copies = .copies,
    lines = .lines,
    map = stats::setNames(paste0("rx_lcCovI_", covs, "_"), covs)
  )
}

#' Model text with each covariate read through its interval value
#'
#' @param txt model text (rxFromSE() output)
#' @param covInterval `.rxFoceiLinCmtCarryCovInterval()`
#' @return `txt` with every covariate symbol replaced
#' @noRd
.rxFoceiLinCmtCarryCovSub <- function(txt, covInterval) {
  if (is.null(covInterval)) {
    return(txt)
  }
  .map <- covInterval$map
  for (.c in names(.map)) {
    .pat <- paste0("(?<![A-Za-z0-9._])", gsub("([.])", "\\\\\\1", .c), "(?![A-Za-z0-9._])")
    txt <- gsub(.pat, .map[[.c]], txt, perl = TRUE)
  }
  txt
}

#' Shared per-model context for the carry emission (parsed once per pair
#' row; cheap)
#' @noRd
.rxFoceiLinCmtCarryCtx <- function(pairs, s) {
  .pc <- .rxFoceiCarryPredCall(s)
  .a <- .pc$args
  # rxFromSE() deparses a non-character argument's expression
  # (substitute()-based), so hand it the repr string
  .txt <- vapply(
    seq_along(.a),
    function(i) {
      .r <- paste(.a[[i]])
      rxode2::rxFromSE(.r)
    },
    character(1)
  )
  .callRepr <- paste(.pc$call)
  .ncmt <- as.integer(as.numeric(.txt[4]))
  .oral0 <- as.integer(as.numeric(.txt[5]))
  .nP <- nrow(pairs)
  .anyJump <- any(!is.na(pairs$fD) | !is.na(pairs$lagD))
  .anyLag <- any(!is.na(pairs$lagD))
  .m <- .ncmt + .oral0
  list(
    ncmt = .ncmt,
    oral0 = .oral0,
    m = .m,
    central = .oral0,
    trans = .txt[8],
    pfx = paste(.txt[1:5], collapse = ","),
    thetas = .txt[9:15],
    # the advance's parameters over the interval ending at the row
    covInterval = attr(pairs, "carryCovInterval"),
    thetasInterval = .rxFoceiLinCmtCarryCovSub(.txt[9:15], attr(pairs, "carryCovInterval")),
    zero = rep("0", 7),
    nP = .nP,
    anyJump = .anyJump,
    anyLag = .anyLag,
    aCol = 2L * .nP,
    lCol = 2L * .nP + 1L,
    slotExpr = lapply(1:7, function(k) .a[[k + 8L]]),
    rows = seq_len(.m) - 1L,
    # an ll() endpoint embeds the concentration call in a larger
    # expression (#1004): the carry differentiates the concentration,
    # read back as rx_lcConc_, and symengine supplies the outer chain rule
    bare = .pc$bare,
    predSym = .pc$predSym,
    conc = if (.pc$bare) "rx_pred_" else "rx_lcConc_",
    concLine = if (.pc$bare) {
      character(0)
    } else {
      paste0("rx_lcConc_~", rxode2::rxFromSE(.callRepr))
    }
  )
}

#' The pair's final `dfe=` composition: s_central / Vc plus, when the
#' slot enters the observation scaling, the row-local direct term
#' -(dVc/dslot) g rx_pred_ / Vc (Vc from rxode2's own micro-constant
#' builder for this trans); a wrapped prediction multiplies the outer
#' derivative in and adds the direct eta term of the likelihood
#' @noRd
.rxFoceiLinCmtCarryFinal <- function(cx, pairs, w, dfe) {
  .p <- w - 1L
  .sC <- paste0("rx_lcCarryS", .p, "r", cx$central, "_")
  .vc <- .rxFoceiCarryVc(
    cx$ncmt,
    cx$oral0,
    as.numeric(cx$trans),
    pairs$slot[w],
    cx$slotExpr
  )
  .vcRepr <- paste(.vc$vc)
  .vcTxt <- paste0("(", rxode2::rxFromSE(.vcRepr), ")")
  .conc <- paste0(.sC, "/", .vcTxt)
  if (!is.null(.vc$dVc)) {
    .dRepr <- paste(.vc$dVc)
    .dTxt <- rxode2::rxFromSE(.dRepr)
    # lead with the symbol: rxode2's expression optimizer folds "-(1)*x"
    # to "-x", and "-" + "-x" would read as a C decrement
    # the observation's own scaling: g at the record, not over the interval
    .conc <- paste0(.conc, "-rx_lcCarryGr", .p, "_*(", .dTxt, ")*", cx$conc, "/", .vcTxt)
  }
  if (cx$bare) {
    return(paste0(dfe, "=", .conc))
  }
  # ll() endpoint: d(pred)/d(eta) = d(pred)/d(conc) * d(conc)/d(eta) plus
  # whatever eta dependence the expression has with the concentration held
  # fixed (symengine sees rx_lcConc_ as a free symbol for that)
  .outer <- paste(symengine::D(cx$predSym, symengine::S("rx_lcConc_")))
  .direct <- symengine::D(cx$predSym, symengine::S(pairs$eta[w]))
  .fin <- paste0(dfe, "=(", rxode2::rxFromSE(.outer), ")*(", .conc, ")")
  if (!.rxFoceiCarryIsZero(.direct)) {
    .directRepr <- paste(.direct)
    .fin <- paste0(.fin, "+(", rxode2::rxFromSE(.directRepr), ")")
  }
  .fin
}
