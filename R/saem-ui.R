.setSaemExtra <- function(.env, type) {
  .uif <- .env$uif
  .txt <- paste0("(", crayon::italic(ifelse(is.null(.uif$nmodel$lin.solved), ifelse(.uif$predSys, "PRED", "ODE"), "Solved")), "); ")
  if (tolower(type) == "focei") {
    .txt <- paste0(.txt, crayon::blurred$italic("OBJF by FOCEi approximation"))
  } else if (tolower(type) == "foce") {
    .txt <- paste0(.txt, crayon::blurred$italic("OBJF by FOCE approximation"))
  } else if (tolower(type) == "fo") {
    .txt <- paste0(.txt, crayon::blurred$italic("OBJF by FO approximation"))
  } else if (type == "") {
    .txt <- paste0(.txt, crayon::blurred$italic("OBJF not calculated"))
  } else {
    .reg <- rex::rex(start, "laplace", capture(.regNum), end)
    .regG <- rex::rex(start, "gauss", capture(.regNum), "_", capture(.regNum), end)
    if (regexpr(.reg, type, perl = TRUE) != -1) {
      .nnode <- 1
      .nsd <- as.numeric(sub(.reg, "\\1", type, perl = TRUE))
    } else if (regexpr(.regG, type, perl = TRUE) != -1) {
      .nnode <- as.numeric(sub(.regG, "\\1", type, perl = TRUE))
      .nsd <- as.numeric(sub(.regG, "\\2", type, perl = TRUE))
    } else {
      stop("unknown error")
    }
    .txt <- paste0(.txt, crayon::blurred$italic(sprintf("OBJF by %s", paste0(ifelse(.nnode == 1, "Lapalcian (n.sd=", sprintf("Gaussian Quadrature (n.nodes=%s, n.sd=", .nnode)), .nsd, ")"))))
  }
  .env$extra <- .txt
  return(invisible(.txt))
}
