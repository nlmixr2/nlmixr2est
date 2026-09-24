# Does the installed rxode2 parameterize a block-internal omega zero itself
# (rxode2#1391), or refuse it so nlmixr2est has to fill it first (#1128)?
.rxSymInvCholHoldsBlockZeros <- function() {
  .m <- matrix(c(1, 0.1, 0.1, 0.1, 1, 0, 0.1, 0, 1), 3, 3)
  !inherits(try(rxode2::rxSymInvCholCreate(mat = .m, diag.xform = "sqrt"), silent = TRUE), "try-error")
}
