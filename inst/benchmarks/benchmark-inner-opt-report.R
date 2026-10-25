# Report for benchmark-inner-opt.R: the #1153 comparison tables, with refArm as
# the reference.  Gap = OFV minus the best any arm reached for that case and
# start; better/worse = an OFV difference of more than 0.01 against refArm.

.innerOptReport <- function(outDir, refArm = "trust") {
  .files <- list.files(outDir, "^arm-.*\\.csv$", full.names = TRUE)
  .d <- do.call(rbind, lapply(.files, utils::read.csv, stringsAsFactors = FALSE))
  .d$key <- paste(.d$case, .d$start)
  .arms <- unique(.d$arm)
  .arms <- c(refArm, setdiff(.arms, refArm))
  # keep rows every arm finished with a finite objective
  .ok <- .d[.d$ok & is.finite(.d$objf), ]
  .keys <- Reduce(intersect, lapply(.arms, function(a) .ok$key[.ok$arm == a]))
  .ok <- .ok[.ok$key %in% .keys, ]
  .w <- function(col) {
    .m <- stats::reshape(.ok[, c("key", "arm", col)], idvar = "key", timevar = "arm", direction = "wide")
    names(.m) <- sub(paste0("^", col, "\\."), "", names(.m))
    .m[order(.m$key), c("key", .arms)]
  }
  .of <- .w("objf")
  .tm <- .w("time")
  .meta <- .ok[.ok$arm == refArm, ]
  .meta <- .meta[match(.of$key, .meta$key), ]
  .real <- !(.meta$est %in% "posthoc") & !is.na(.meta$maxOuter) & .meta$maxOuter > 0
  .best <- apply(.of[, .arms, drop = FALSE], 1, min)
  .gm <- function(x) exp(mean(log(pmax(x, 1e-3))))
  .fmt <- function(x, d = 3) formatC(x, digits = d, format = "f")
  .out <- character()
  .add <- function(...) .out <<- c(.out, paste0(...))
  .tab <- function(rows, hdr) {
    .add("| ", paste(hdr, collapse = " | "), " |")
    .add("|", paste(rep("---", length(hdr)), collapse = "|"), "|")
    for (r in rows) {
      .add("| ", paste(r, collapse = " | "), " |")
    }
    .add("")
  }
  .summ <- function(sel, label) {
    .add("**", label, " (", sum(sel), ")**")
    .add("")
    .rows <- lapply(.arms, function(a) {
      .gap <- .of[[a]][sel] - .best[sel]
      .dl <- .of[[a]][sel] - .of[[refArm]][sel]
      c(
        a,
        paste0(.fmt(100 * mean(.gap <= 0.01), 1), "%"),
        sum(.gap > 1),
        if (a == refArm) "" else paste0(sum(.dl < -0.01), " / ", sum(.dl > 0.01)),
        paste0(.fmt(.gm(.tm[[a]][sel] / .tm[[refArm]][sel]), 2), "x"),
        paste0(round(sum(.tm[[a]][sel])), " s")
      )
    })
    .tab(
      .rows,
      c(
        "arm",
        "within 0.01 of best",
        "> 1 OFV off",
        paste0("better / worse vs ", refArm),
        "time (geometric mean)",
        "total time"
      )
    )
  }
  .add(
    "Cases: ",
    length(unique(.meta$case)),
    " from ",
    length(unique(.meta$source)),
    " sources; ",
    nrow(.of),
    " paired rows (",
    sum(.real),
    " real fits, ",
    sum(!.real),
    " fixed evaluations).  Rows that errored in any arm: ",
    length(unique(.d$key)) - nrow(.of),
    "."
  )
  .add("")
  .summ(.real, "Real fits")
  .big <- .real & .tm[[refArm]] > 5
  .summ(.big, paste0("Real fits taking more than 5 s under ", refArm))
  .summ(!.real, "Fixed-parameter evaluations")
  # by type
  .feat <- c(
    errAdd = "additive error",
    errProp = "proportional",
    errAddProp = "add+prop",
    linCmt = "linCmt",
    ode = "ODE",
    blockOmega = "block omega",
    multiEndpoint = "multiple endpoints",
    cens = "censoring",
    covariates = "covariates",
    fast = "`fast = TRUE`",
    eta3 = "3+ ETAs"
  )
  .add("**By FOCEi type (real fits): better/worse vs ", refArm, ", > 1 OFV off, time**")
  .add("")
  .rows <- list()
  for (.f in names(.feat)) {
    .sel <- .real & as.logical(.meta[[.f]]) %in% TRUE
    if (!any(.sel)) {
      next
    }
    .r <- c(.feat[[.f]], sum(.sel))
    for (a in setdiff(.arms, refArm)) {
      .dl <- .of[[a]][.sel] - .of[[refArm]][.sel]
      .r <- c(
        .r,
        paste0(
          sum(.dl < -0.01),
          " / ",
          sum(.dl > 0.01),
          "; ",
          sum(.of[[a]][.sel] - .best[.sel] > 1),
          "; ",
          .fmt(.gm(.tm[[a]][.sel] / .tm[[refArm]][.sel]), 2),
          "x"
        )
      )
    }
    .rows[[length(.rows) + 1]] <- .r
  }
  for (.e in sort(unique(.meta$est[.real]))) {
    .sel <- .real & .meta$est == .e
    .r <- c(.e, sum(.sel))
    for (a in setdiff(.arms, refArm)) {
      .dl <- .of[[a]][.sel] - .of[[refArm]][.sel]
      .r <- c(
        .r,
        paste0(
          sum(.dl < -0.01),
          " / ",
          sum(.dl > 0.01),
          "; ",
          sum(.of[[a]][.sel] - .best[.sel] > 1),
          "; ",
          .fmt(.gm(.tm[[a]][.sel] / .tm[[refArm]][.sel]), 2),
          "x"
        )
      )
    }
    .rows[[length(.rows) + 1]] <- .r
  }
  .tab(.rows, c("group", "n", setdiff(.arms, refArm)))
  # size of the changes
  .add("**Size of the changes against ", refArm, " (real fits)**")
  .add("")
  .rows <- list()
  .stat <- list(
    "fits better / worse" = function(dl) paste0(sum(dl < -0.01), " / ", sum(dl > 0.01)),
    "median improvement" = function(dl) .fmt(stats::median(-dl[dl < -0.01])),
    "mean improvement" = function(dl) .fmt(mean(-dl[dl < -0.01]), 2),
    "largest improvement" = function(dl) .fmt(max(c(0, -dl)), 1),
    "median regression" = function(dl) .fmt(stats::median(dl[dl > 0.01])),
    "mean regression" = function(dl) .fmt(mean(dl[dl > 0.01]), 2),
    "largest regression" = function(dl) .fmt(max(c(0, dl)), 1),
    "net OFV over all fits" = function(dl) .fmt(sum(dl), 1),
    "net, each fit capped at +/-10" = function(dl) .fmt(sum(pmin(pmax(dl, -10), 10)), 1)
  )
  for (.s in names(.stat)) {
    .r <- .s
    for (a in setdiff(.arms, refArm)) {
      .r <- c(.r, .stat[[.s]](.of[[a]][.real] - .of[[refArm]][.real]))
    }
    .rows[[length(.rows) + 1]] <- .r
  }
  .tab(.rows, c("", setdiff(.arms, refArm)))
  # largest movements
  for (a in setdiff(.arms, refArm)) {
    .dl <- .of[[a]][.real] - .of[[refArm]][.real]
    .k <- .of$key[.real]
    .o <- order(.dl)
    .imp <- head(.o[.dl[.o] < -0.1], 5)
    .reg <- head(rev(.o)[.dl[rev(.o)] > 0.1], 5)
    .lbl <- function(i) {
      sprintf("`%s` %s -> %s (%+.1f)", .k[i], .fmt(.of[[refArm]][.real][i]), .fmt(.of[[a]][.real][i]), .dl[i])
    }
    .add(
      "**",
      a,
      "** largest improvements: ",
      if (length(.imp)) paste(vapply(.imp, .lbl, ""), collapse = "; ") else "none"
    )
    .add("")
    .add(
      "**",
      a,
      "** largest regressions: ",
      if (length(.reg)) paste(vapply(.reg, .lbl, ""), collapse = "; ") else "none"
    )
    .add("")
  }
  writeLines(.out, file.path(outDir, "report.md"))
  cat(.out, sep = "\n")
  invisible(.out)
}
