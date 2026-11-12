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
#' Around a lagged dose (`.rxFoceiLinCmtCarryLagKnot()`) every locf, nocb
#' and midpoint covariate is read through both copies, and through its next
#' data knot (`.rxFoceiLinCmtCarryNextKnot()`).
#'
#' @param ui rxode2 ui with its control
#' @return the covariates under midpoint, by the control or by their own
#'   declaration, of a `linCmt()` model the carry is on for, with an rxode2
#'   that advances an interval with its own covariates
#'   (`.rxFoceiLinCmtIntervalCovs()`); `character(0)` otherwise.  Attribute
#'   `knot` is `TRUE` around a lagged dose.
#' @noRd
.rxFoceiLinCmtCarryMidCovs <- function(ui) {
  if (!.foceiUsesLinCmt(ui) || !.rxFoceiLinCmtCarryBuildEnabled(ui) || !.rxFoceiLinCmtIntervalCovs()) {
    return(character(0))
  }
  .covs <- ui$allCovs
  .interp <- .rxFoceiLinCmtCarryInterp(ui, .covs)
  if (.rxFoceiLinCmtCarryLagKnot(ui)) {
    return(structure(.covs[.interp %in% c("locf", "nocb", "midpoint")], knot = TRUE))
  }
  .covs[.interp == "midpoint"]
}

#' Data with the hidden copies the carry reads under midpoint
#'
#' @param data data frame or event table with the covariate columns
#' @param covs `.rxFoceiLinCmtCarryMidCovs()`
#' @return `data` with `rx_lcL_<cov>` and `rx_lcN_<cov>`, copies of each of
#'   `covs`, and, when `covs` has `knot`, each one's next knot
#'   (`.rxFoceiLinCmtCarryNextKnot()`)
#' @noRd
.rxFoceiLinCmtCarryAliasData <- function(data, covs) {
  for (.c in covs) {
    data[[paste0("rx_lcL_", .c)]] <- data[[.c]]
    data[[paste0("rx_lcN_", .c)]] <- data[[.c]]
  }
  if (isTRUE(attr(covs, "knot"))) {
    data <- .rxFoceiLinCmtCarryNextKnot(data, covs)
  }
  data
}

#' Each row's next data knot of each covariate
#'
#' A knot is a later time at which the subject has a non-missing value.  Each
#' row gets that time (`rx_lcTn_<cov>`, a huge time past the last), the value
#' a locf read there takes (`rx_lcNl_<cov>`, the time's last row) and the one
#' a nocb read takes (`rx_lcNn_<cov>`, its first row).
#'
#' @param data data sorted by subject and time
#' @param covs covariates
#' @return `data` with the three columns per covariate
#' @noRd
.rxFoceiLinCmtCarryNextKnot <- function(data, covs) {
  .t <- data[[intersect(c("TIME", "time"), names(data))[1]]]
  .id <- data[[intersect(c("ID", "id"), names(data))[1]]]
  .rows <- split(seq_along(.t), factor(.id, levels = unique(.id)))
  for (.c in covs) {
    .v <- data[[.c]]
    .tn <- rep(1e300, length(.t))
    .nl <- .v
    .nn <- .v
    for (.r in .rows) {
      .ok <- .r[!is.na(.v[.r])]
      if (length(.ok) == 0L) {
        next
      }
      .kt <- unique(.t[.ok])
      .last <- .ok[!duplicated(.t[.ok], fromLast = TRUE)]
      .first <- .ok[!duplicated(.t[.ok])]
      # index of the first knot strictly after each row's time
      .k <- findInterval(.t[.r], .kt) + 1L
      .has <- .k <= length(.kt)
      .tn[.r[.has]] <- .kt[.k[.has]]
      .nl[.r[.has]] <- .v[.last[.k[.has]]]
      .nn[.r[.has]] <- .v[.first[.k[.has]]]
    }
    data[[paste0("rx_lcTn_", .c)]] <- .tn
    data[[paste0("rx_lcNl_", .c)]] <- .nl
    data[[paste0("rx_lcNn_", .c)]] <- .nn
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
#' With `knot` (`.rxFoceiLinCmtCarryLagKnot()`) a lagged dose's data time is
#' a knot but not a record, so it can lie inside an interval: the one that
#' ends at its arrival, or one that ends at a record before it.  The
#' previous record's next knot `F` (`.rxFoceiLinCmtCarryNextKnot()`) splits
#' an interval it lies inside: with `m` the interval's middle, `wP` the
#' covariate's locf view at the previous record, `wN` its nocb view here and
#' `wF` its locf (nocb) value at `F`, locf reads `m >= F ? wF : wP`, nocb
#' `m < F ? wF : wN` and midpoint their mean.  It assumes at most one knot
#' inside an interval (two lagged doses still to arrive at one record).
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
#' @param knot `.rxFoceiLinCmtCarryLagKnot()`
#' @return NULL when every covariate's record value is its interval's (nocb,
#'   no covariates, or an rxode2 that advances every interval with the end
#'   record's covariates, `.rxFoceiLinCmtIntervalCovs()`), else list(copies =
#'   the copy and declaration lines, to go ahead of rx_pred_, lines = the lines
#'   defining each interval value, map = covariate -> the symbol holding its
#'   interval value)
#' @noRd
.rxFoceiLinCmtCarryCovInterval <- function(covs, interp, knot = FALSE) {
  if (!.rxFoceiLinCmtIntervalCovs()) {
    return(NULL)
  }
  covs <- covs[interp[covs] %in% c(if (knot) "nocb", "locf", "midpoint", "linear")]
  if (length(covs) == 0L) {
    return(NULL)
  }
  .copies <- character(0)
  .lines <- character(0)
  if (knot) {
    .copies <- "rx_lcT_=t"
    .lines <- c("rx_lcTp_~lag0(rx_lcT_)", "rx_lcMid_~0.5*(rx_lcTp_+t)")
  }
  for (.c in covs) {
    .copy <- paste0("rx_lcCov_", .c, "_")
    .lag <- paste0("rx_lcLag_", .c, "_")
    .int <- paste0("rx_lcCovI_", .c, "_")
    .two <- knot || interp[[.c]] == "midpoint"
    # with the hidden copies the lag reads the locf copy of the column
    .src <- if (.two) paste0("rx_lcL_", .c) else .c
    .copies <- c(
      .copies,
      if (.two) c(paste0("locf(rx_lcL_", .c, ")"), paste0("nocb(rx_lcN_", .c, ")")),
      paste0(.copy, "=", .src)
    )
    .lines <- c(.lines, paste0(.lag, "~lag0(", .copy, ")"))
    if (!knot) {
      .lines <- c(
        .lines,
        switch(
          interp[[.c]],
          locf = paste0(.int, "~", .lag),
          midpoint = paste0(.int, "~0.5*(", .lag, "+rx_lcN_", .c, ")"),
          linear = paste0(.int, "~0.5*(", .c, "+", .lag, ")")
        )
      )
      next
    }
    # the previous record's next knot, and the locf/nocb values there
    .k <- paste0(c("rx_lcTn_", "rx_lcNl_", "rx_lcNn_"), .c)
    .kc <- paste0(.k, "C_")
    .kp <- paste0(.k, "P_")
    .copies <- c(.copies, paste0("locf(", .k, ")"), paste0(.kc, "=", .k))
    .lines <- c(.lines, paste0(.kp, "~lag0(", .kc, ")"))
    .in <- paste0(.kp[1], "<t")
    .lo <- paste0("ifelse(", .in, "&&rx_lcMid_>=", .kp[1], ",", .kp[2], ",", .lag, ")")
    .no <- paste0("ifelse(", .in, "&&rx_lcMid_<", .kp[1], ",", .kp[3], ",rx_lcN_", .c, ")")
    .lines <- c(
      .lines,
      paste0(
        .int,
        "~",
        switch(interp[[.c]], locf = .lo, nocb = .no, midpoint = paste0("0.5*(", .lo, "+", .no, ")"))
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

#' The carry's tracker columns, one per slot
#'
#' @param slot the pairs' slots (1-7, `NA` for an event-only pair)
#' @return list(col = each pair's 0-based tracker column after the `nP`
#'   carry columns, `NA` without a slot; last = whether the pair is the
#'   last of its slot, which updates the tracker once every pair of the slot
#'   has read it; n = the number of tracker columns)
#' @noRd
.rxFoceiLinCmtCarryTrackers <- function(slot) {
  .slots <- unique(slot[!is.na(slot)])
  list(
    col = length(slot) + match(slot, .slots) - 1L,
    last = !is.na(slot) & !duplicated(slot, fromLast = TRUE),
    n = length(.slots)
  )
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
  # one tracker column per slot: it holds M_i J_{i-1}(:, k), which belongs
  # to the slot, so the pairs of one slot share it
  .trk <- .rxFoceiLinCmtCarryTrackers(pairs$slot)
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
    trk = .trk$col,
    trkLast = .trk$last,
    anyJump = .anyJump,
    anyMultCov = any(pairs$shape %in% "multCov"),
    anyLag = .anyLag,
    aCol = .nP + .trk$n,
    lCol = .nP + .trk$n + 1L,
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
