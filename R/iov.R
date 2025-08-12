.uiApplyIov <- function(env) {
  .ui <- env$ui
  .iniDf <- .ui$iniDf
  .lvls <- .iniDf$condition[which(!is.na(.iniDf$condition) &
                                    .iniDf$condition != "id" &
                                     is.na(.iniDf$err))]
  if (length(.lvls) > 0) {
    .n <- .iniDf[which(.iniDf$condition %in% .lvls), "name"]
    .ui <- eval(str2lang(paste0("rxode2::rxRename(.ui, ",
                                paste(paste0("rx.", .n, "=", .n),
                                      collapse=", "), ")")))
    .ui <- rxode2::rxUiDecompress(.ui)
    # For the new iniDf, we will take out all the level variables and
    # then renumber the etas
    .thetas <- .iniDf[is.na(.iniDf$neta1),, drop=FALSE]
    .etas <- .iniDf[is.na(.iniDf$ntheta),, drop=FALSE]
    if (length(.thetas$name) > 0) {
      .maxtheta <- max(.thetas$ntheta, na.rm = TRUE)
      .theta1 <- .thetas[1,]
      .theta1$ntheta <- .maxtheta
      .theta1$lower <- 0
    } else {
      .maxtheta <- 0L
      .theta1 <- .etas[1,]
      .theta1$ntheta <- 0L
      .theta1$neta1 <- NA_real_
      .theta1$neta2 <- NA_real_
      .theta1$lower <- 0
    }
    .eta1 <- .etas[1, ]
    .eta1$fix <- TRUE
    .eta1$neta1 <- .eta1$neta2 <- 0
    .eta1$est <- 1

    .etas <- .etas[which(!(.etas$condition %in% .lvls)), , drop=FALSE]
    if (length(.etas$name) > 0) {
      .etas$neta1 <- factor(.etas$neta1, levels = sort(unique(.etas$neta1)))
      .etas$neta2 <- factor(.etas$neta2, levels = sort(unique(.etas$neta2)))
      .etas$neta1 <- as.integer(.etas$neta1)
      .etas$neta2 <- as.integer(.etas$neta2)
      .maxeta <- max(.etas$neta1, na.rm = TRUE)
    } else {
      .maxeta <- 0L
      .theta1 <- .etas[1,]
    }

    .data <- env$data
    .lvls <- setNames(lapply(.lvls, function(l) {
      sort(unique(.data[[l]]))
    }), .lvls)
    .env <- new.env(parent = emptyenv())
    .env$thetas <- .thetas
    .env$etas <- .etas
    .env$maxtheta <- .maxtheta
    .env$maxeta <- .maxeta
    # Now we have enough information to create the IOV variables
    # changed to etas on id
    .lines <- lapply(names(.lvls),
                     function(l1) {
                       .var <- .iniDf$name[which(.iniDf$condition == l1)]
                       .lst <- lapply(.var, function(v) {
                         # Add theta to dataset; represents variance of iov
                         .curTheta <- .theta1
                         .curTheta$est <- .iniDf[which(.iniDf$name == v &
                                                         is.na(.iniDf$ntheta)), "est"]
                         .curTheta$name <- v
                         .env$maxtheta <- .curTheta$ntheta <- .env$maxtheta + 1L
                         .env$thetas <- rbind(.env$thetas, .curTheta)
                         for (n in .lvls[[l1]]) {
                           .curEta <- .eta1
                           .curEta$name <- paste0("rx.", v, ".", n)
                           .env$maxeta <- .curEta$neta1 <-
                             .curEta$neta2 <- .env$maxeta + 1L
                           .env$etas <- rbind(.env$etas, .curEta)
                         }
                         str2lang(paste0("rx.", v, " <- sqrt(", v, ")*(",
                                         paste(paste0("rx.", v, ".", .lvls[[l1]],
                                                      "*(", l1,
                                                      " == ", .lvls[[l1]], ")"),
                                               collapse="+"),
                                         ")"))
                       })
                       .lst
                     })
    .lines <- do.call(`c`, c(.lines, list(.ui$lstExpr)))
    # Now the lines can be added to the model
    assign("iniDf", rbind(.env$thetas,.env$etas), envir = .ui)
    assign("lstExpr", .lines, envir = .ui)
    env$iov <- env$ui
    env$ui <- .ui$fun()
  }
}
