#' Outer-gradient columns of a time-varying `linCmt()` model through the carry
#'
#' A `linCmt()` slot that reads a covariate changing within a subject changes the
#' kernel's parameter vector between records; the parameter-space reads then expand
#' around the wrong point and read `NA` (nlmixr2/rxode2lincmt#19), which sent every
#' subject to the per-subject finite difference.  rxode2lincmt's direction carry
#' (`linCmtB()` which1 = -30 .. -34) follows parameters that change between records,
#' so it is driven here over the outer gradient's model directions: per row each
#' slot's first and second direction derivatives go in (-30/-31), the carry advances
#' (-32) and the concentration's direction derivatives come back (-33/-34).
#' @param ui rxode2 ui
#' @param s pruned symengine env of the augmented build
#' @param dirs model directions (integrated, not covariate-scaled)
#' @return NULL when out of scope (no time-varying slot), a decline reason (character),
#'   or list(lines, call, conc, c1, c2) with `c1[dir]` the first-order symbol name and
#'   `c2(p, q)` the second-order one
#' @noRd
.foceiOuterLinCmtCarry <- function(ui, s, dirs) {
  if (!.foceiUsesLinCmt(ui) || !.foceiLinCmtHasSecondOrder()) {
    return(NULL)
  }
  .a <- .rxFoceiCarryPredArgs(ui, s)
  if (is.null(.a)) {
    return(NULL)
  }
  .txt <- vapply(seq_along(.a), function(i) .foceiEcSe(.a[[i]]), character(1))
  .slotVars <- unique(unlist(lapply(9:15, function(k) all.vars(str2lang(.txt[k])))))
  .const <- tryCatch(rxode2::rxGetControl(ui, "foceiConstCovs", NULL), error = function(e) NULL)
  .tv <- intersect(.slotVars, c(setdiff(ui$allCovs, .const), "t", "time"))
  if (length(.tv) == 0L) {
    return(NULL)
  }
  if (any(grepl("^rx_(f|lag|alag|rate|dur)_.+_$", ls(envir = s, all.names = TRUE)))) {
    return("a modeled dose with a time-varying linCmt() parameter")
  }
  .E <- 16L # rxode2lincmt's RX_LINCMT_ETA_MAX
  if (length(dirs) > .E) {
    return("more than 16 directions for the linCmt() carry")
  }
  .g <- .foceiEcGradLines(.a, .txt, dirs, .E)
  if (length(.g) == 0L) {
    return(NULL)
  }
  .n <- length(dirs)
  .c1 <- setNames(sprintf("rx_ecC%d_", seq_len(.n)), dirs)
  .lines <- c(
    .g,
    paste0("rx_ecAdv_=", .foceiEcCall(.txt, "-32", 0L)),
    paste0("rx_lcConc_~", .foceiEcCall(.txt, "-1", -1L)),
    sprintf("%s~%s", .c1, vapply(seq_len(.n), function(.i) .foceiEcCall(.txt, "-33", .i - 1L), character(1)))
  )
  for (.j in seq_len(.n)) {
    for (.i in seq_len(.j)) {
      .lines <- c(.lines, sprintf("rx_ecC%d_%d_~%s", .i, .j, .foceiEcCall(.txt, "-34", (.i - 1L) * .E + (.j - 1L))))
    }
  }
  .pc <- .rxFoceiCarryPredCall(s)
  list(
    lines = .lines,
    call = .pc$call,
    conc = "rx_lcConc_",
    c1 = .c1,
    c2 = function(p, q) {
      .i <- sort(match(c(p, q), dirs))
      sprintf("rx_ecC%d_%d_", .i[1], .i[2])
    }
  )
}
