# TRUE when rxSymInvCholCreate() cannot parameterize `m` by its own nonzero
# pattern: older rxode2 refuses it, rxode2 >= #1365 makes a zero a free cell.
.rxInvPatternMismatch <- function(m, same = NULL) {
  .r <- try(rxode2::rxSymInvCholCreate(mat = m, diag.xform = "sqrt", same = same), silent = TRUE)
  inherits(.r, "try-error") || length(.r$theta) != sum(m[upper.tri(m, diag = TRUE)] != 0)
}
