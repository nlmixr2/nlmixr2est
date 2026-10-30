# What the loaded rxode2lincmt does with the records that modify a linCmt()
# compartment outright: a reset (EVID 3, and the reset an EVID 4 begins
# with), a replace (EVID 5) and a multiply (EVID 6).  The data-aware carry
# gate in .foceiFamilyReturn() reads these.  Each answer comes from a small
# solve, made once per session and only when the data has such a record;
# they are behavioural probes because no version or exported symbol tells
# the builds apart.

.linCmtEventEnv <- new.env(parent = emptyenv())

#' A cached probe answer: TRUE only when `probe()` ran and returned TRUE
#' @param name the cache key
#' @param probe a function of no arguments
#' @noRd
.rxFoceiLinCmtEventProbe <- function(name, probe) {
  if (is.null(.linCmtEventEnv[[name]])) {
    .linCmtEventEnv[[name]] <- isTRUE(try(probe(), silent = TRUE))
  }
  .linCmtEventEnv[[name]]
}

#' A linCmtB() carry model: per record, the -5 advance, a -7 add of 0.5 and
#' the -6 read of the column
#' @noRd
.rxFoceiLinCmtEventCarryModel <- function() {
  .a <- "rx__PTR__, t, 1, 1, 0, %d, %d, 1, 1, 10, 0.5, 0, 0, 0, 0"
  suppressWarnings(suppressMessages(rxode2::rxode2(paste0(
    "cp=linCmtB(", sprintf(.a, -1L, -1L), ")\n",
    "ad=linCmtB(", sprintf(.a, -5L, 0L), ")\n",
    "a0=linCmtB(", sprintf(.a, -7L, 0L), ")\n",
    "r0=linCmtB(", sprintf(.a, -6L, 0L), ")"
  ))))
}

#' The carry column of .rxFoceiLinCmtEventCarryModel() at each observation,
#' with the constant-theta fast path pinned on (and restored)
#' @param evid,amt the record at 2 h, after a dose at 0 and an observation
#'   at 1, before one at 3
#' @noRd
.rxFoceiLinCmtEventCarryRead <- function(evid, amt) {
  .fast <- rxode2::linCmtCarrySetFast(TRUE)
  on.exit(rxode2::linCmtCarrySetFast(.fast), add = TRUE)
  .ev <- data.frame(
    id = 1L,
    time = 0:3,
    evid = c(1L, 0L, evid, 0L),
    amt = c(100, 0, amt, 0),
    cmt = 1L
  )
  .s <- suppressWarnings(rxode2::rxSolve(.rxFoceiLinCmtEventCarryModel(), .ev, returnType = "data.frame"))
  setNames(.s$r0, .s$time)
}

#' Does the carry restart at a reset (rxode2lincmt PR 4)?  The column reads
#' the same after the reset as after the dose only if it does.
#' @noRd
.rxFoceiLinCmtCarryResetProbe <- function() {
  .r <- .rxFoceiLinCmtEventCarryRead(3L, 0)
  isTRUE(all.equal(.r[["3"]], .r[["1"]]))
}

#' Does the carry apply a replace or multiply to its columns?  After a
#' replace with 0 the column restarts from 0, so the next row reads 0.5; a
#' carry that ignores the event reads 2.
#' @noRd
.rxFoceiLinCmtCarryEventProbe <- function() {
  .r <- .rxFoceiLinCmtEventCarryRead(5L, 0)
  isTRUE(abs(.r[["3"]] - 0.5) < 1e-12)
}

#' Does linCmt()'s own Jacobian follow a replace?  After replacing the
#' central amount with 50 at 24 h, d(A)/d(cl) at 30 h is the replaced
#' amount's alone, -50 * (6 / v) * exp(-6 * cl / v).
#' @noRd
.rxFoceiLinCmtJacEventProbe <- function() {
  .a <- "rx__PTR__, t, 1, 1, 0, %d, %d, 1, cl, 20, 0, 0, 0, 0, 0"
  .m <- suppressWarnings(suppressMessages(rxode2::rxode2(paste0(
    "cp=linCmtB(", sprintf(.a, -1L, -1L), ")\n",
    "jr=linCmtB(", sprintf(.a, 0L, 0L), ")"
  ))))
  .ev <- data.frame(
    id = 1L,
    time = c(0, 12, 24, 30),
    evid = c(1L, 0L, 5L, 0L),
    amt = c(100, 0, 50, 0),
    cmt = 1L
  )
  .s <- suppressWarnings(rxode2::rxSolve(.m, .ev, params = c(cl = 2), returnType = "data.frame"))
  .want <- -50 * (6 / 20) * exp(-6 * 2 / 20)
  isTRUE(abs(.s$jr[.s$time == 30] - .want) < 1e-8 * abs(.want))
}

#' The three probes, cached
#' @noRd
.rxFoceiLinCmtCarryResetOk <- function() {
  .rxFoceiLinCmtEventProbe("reset", .rxFoceiLinCmtCarryResetProbe)
}

#' @noRd
.rxFoceiLinCmtCarryEventOk <- function() {
  .rxFoceiLinCmtEventProbe("carryEvent", .rxFoceiLinCmtCarryEventProbe)
}

#' @noRd
.rxFoceiLinCmtJacEventOk <- function() {
  .rxFoceiLinCmtEventProbe("jacEvent", .rxFoceiLinCmtJacEventProbe)
}

#' Which records in the data modify a compartment outright
#'
#' EVID 3 and 4 reset, 5 replaces and 6 multiplies; rxode2's encoded EVIDs
#' (100 and up) carry the type in their fourth digit (4 replace, 5
#' multiply).
#' @param evid the data's EVID column (NULL when it has none)
#' @return list(reset =, repMult =) of logicals
#' @noRd
.rxFoceiLinCmtEventKinds <- function(evid) {
  .e <- evid[!is.na(evid)]
  .whI <- floor(.e[.e >= 100] / 1e4) %% 10
  list(
    reset = any(.e %in% c(3, 4)),
    repMult = any(.e %in% c(5, 6)) || any(.whI %in% c(4, 5))
  )
}

#' Why the carry cannot run on this data, or NULL
#'
#' Steady-state doses are not in the carry's recurrence.  A reset needs an
#' rxode2lincmt whose carry restarts there, and a replace or multiply one
#' whose carry and whose own Jacobian both follow the event.
#' @param rd the data, column names upper case
#' @noRd
.rxFoceiLinCmtCarryEventWhy <- function(rd) {
  if ("SS" %in% names(rd) && any(rd[["SS"]] > 0, na.rm = TRUE)) {
    return("ss rows")
  }
  .ev <- .rxFoceiLinCmtEventKinds(rd[["EVID"]])
  if (.ev$reset && !.rxFoceiLinCmtCarryResetOk()) {
    return("reset rows")
  }
  if (.ev$repMult && !(.rxFoceiLinCmtCarryEventOk() && .rxFoceiLinCmtJacEventOk())) {
    return("replace/multiply rows")
  }
  NULL
}

#' Are a linCmt() fit's eta gradients approximate on this data?
#'
#' TRUE for replace or multiply records when the loaded rxode2lincmt's own
#' Jacobian does not follow them: then the standard gradient is as wrong as
#' the carry.
#' @param ui the model
#' @param rd the data, column names upper case
#' @noRd
.rxFoceiLinCmtJacApprox <- function(ui, rd) {
  .flg <- rxode2::rxModelVars(ui)$flags
  if (!("linCmtFlg" %in% names(.flg)) || .flg[["linCmtFlg"]] <= 0) {
    return(FALSE)
  }
  .rxFoceiLinCmtEventKinds(rd[["EVID"]])$repMult && !.rxFoceiLinCmtJacEventOk()
}
