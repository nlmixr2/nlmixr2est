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
#' Under midpoint and linear rxode2 advances an interval with the covariate
#' at the interval's middle, from the nearest records that carry it on
#' either side: their mean under midpoint, the line between them under
#' linear.  A record with no value of its own (an infusion's end, an addl
#' dose, a missing value) lies between them, so neither the record nor its
#' neighbour stands for them.  The carry reads them instead as the
#' covariate's locf view at the interval's start and its nocb view at its
#' end, from two hidden copies of the column (`rx_lcL_<cov>`,
#' `rx_lcN_<cov>`) that the generated models declare `locf()` and `nocb()`;
#' under linear, two more hidden columns (`rx_lcTL_<cov>`, `rx_lcTN_<cov>`)
#' give those records' times the same way.
#'
#' @param ui rxode2 ui with its control
#' @return the interpolation ("midpoint" or "linear") of each covariate under
#'   one of them, by the control or by its own declaration, named by
#'   covariate, of a `linCmt()` model the carry is on for; `character(0)`
#'   otherwise
#' @noRd
.rxFoceiLinCmtCarryAliasCovs <- function(ui) {
  if (!.foceiUsesLinCmt(ui) || !.rxFoceiLinCmtCarryBuildEnabled(ui)) {
    return(character(0))
  }
  .covs <- ui$allCovs
  .interp <- .rxFoceiLinCmtCarryInterp(ui, .covs)
  .interp[.interp %in% c("midpoint", "linear")]
}

#' Data with the hidden copies the carry reads under midpoint and linear
#'
#' @param data data frame or event table with the covariate columns and a
#'   time column
#' @param covs `.rxFoceiLinCmtCarryAliasCovs()`
#' @return `data` with `rx_lcL_<cov>` and `rx_lcN_<cov>`, copies of each of
#'   `covs`, and for each linear one `rx_lcTL_<cov>` and `rx_lcTN_<cov>`,
#'   the record's time where the covariate has a value and NA where it has
#'   none
#' @noRd
.rxFoceiLinCmtCarryAliasData <- function(data, covs) {
  .time <- NULL
  for (.c in names(covs)) {
    data[[paste0("rx_lcL_", .c)]] <- data[[.c]]
    data[[paste0("rx_lcN_", .c)]] <- data[[.c]]
    if (covs[[.c]] == "linear") {
      if (is.null(.time)) {
        .time <- data[[which(tolower(names(data)) == "time")[1]]]
      }
      .t <- ifelse(is.na(data[[.c]]), NA_real_, .time)
      data[[paste0("rx_lcTL_", .c)]] <- .t
      data[[paste0("rx_lcTN_", .c)]] <- .t
    }
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
#' - linear: the line between its locf view at the previous record and its
#'   nocb view at this one, at the interval's middle, with those records'
#'   times from the hidden time copies.
#'
#' The lag is taken from a model copy (`copies`, real lhs ahead of
#' rx_pred_, as lag() needs), not from the column: a record with no value of
#' its own holds NA, which lag0() turns into 0, where the model evaluation
#' fills it as rxode2 fills any NA covariate.  A subject's first record ends
#' an empty interval, so the carry there does not depend on the lag, which is
#' 0 in a fresh solve and the previous pass's last record when FOCEi solves
#' the subject again for another eta.
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
  if (any(interp[covs] == "linear")) {
    # the previous record's time, for the interval's middle
    .copies <- "rx_lcT_=t"
    .lines <- "rx_lcTPrev_~lag0(rx_lcT_)"
  }
  for (.c in covs) {
    .copy <- paste0("rx_lcCov_", .c, "_")
    .lag <- paste0("rx_lcLag_", .c, "_")
    .int <- paste0("rx_lcCovI_", .c, "_")
    .alias <- interp[[.c]] %in% c("midpoint", "linear")
    .lin <- interp[[.c]] == "linear"
    # under midpoint and linear the lag reads the locf copy of the column
    .src <- if (.alias) paste0("rx_lcL_", .c) else .c
    .copies <- c(
      .copies,
      if (.alias) c(paste0("locf(rx_lcL_", .c, ")"), paste0("nocb(rx_lcN_", .c, ")")),
      if (.lin) {
        c(
          paste0("locf(rx_lcTL_", .c, ")"),
          paste0("nocb(rx_lcTN_", .c, ")"),
          paste0("rx_lcTim_", .c, "_=rx_lcTL_", .c)
        )
      },
      paste0(.copy, "=", .src)
    )
    .tLag <- paste0("rx_lcTLag_", .c, "_")
    .tN <- paste0("rx_lcTN_", .c)
    .lines <- c(
      .lines,
      paste0(.lag, "~lag0(", .copy, ")"),
      if (.lin) paste0(.tLag, "~lag0(rx_lcTim_", .c, "_)"),
      switch(
        interp[[.c]],
        locf = paste0(.int, "~", .lag),
        midpoint = paste0(.int, "~0.5*(", .lag, "+rx_lcN_", .c, ")"),
        # equal times: the interval lies past the covariate's last value (or
        # before its first), where rxode2 holds that value
        linear = paste0(
          .int,
          "~",
          .lag,
          "+(rx_lcN_",
          .c,
          "-",
          .lag,
          ")*(0.5*(rx_lcTPrev_+t)-",
          .tLag,
          ")/(",
          .tN,
          "-",
          .tLag,
          "+(",
          .tN,
          "==",
          .tLag,
          "))"
        )
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
