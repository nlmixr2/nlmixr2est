## vaeCovariates(): exported view of the VAE automatic covariate search
## (fast, data-only -- no model build or training)

nmTest({
  test_that("vaeCovariates discovers subject-level covariates with fit rules", {
    d <- data.frame(
      id = rep(1:4, each = 3),
      time = rep(0:2, 4),
      dv = 1:12,
      wt = rep(c(70, 80, 60, 75), each = 3),
      sex = rep(c(0, 1, 0, 1), each = 3))
    res <- vaeCovariates(d)
    expect_s3_class(res, "data.frame")
    ## a continuous covariate contributes one column per shape family, and the
    ## hockey family contributes one per ARM
    expect_equal(res$covariate,
                 c("WT_power", "WT_lin", "WT_hockeyLow", "WT_hockeyHi", "SEX"))
    expect_equal(res$raw, c("WT", "WT", "WT", "WT", "SEX"))
    expect_equal(res$shape, c("power", "lin", "hockeyLow", "hockeyHi", "cat"))
    expect_equal(res$type, c(rep("continuous", 4), "categorical"))
    ## WT is centered at its median; SEX is a 0/1 indicator so it is left RAW
    ## (center 0) -- its coefficient is the level-1 shift and the structural
    ## theta stays the reference (SEX=0) value
    expect_equal(res$center, c(rep(stats::median(c(70, 80, 60, 75)), 4), 0))
    ## every WT shape competes for one slot; SEX is its own group
    expect_equal(res$group, c(1L, 1L, 1L, 1L, 2L))
    ## but the two arms are ONE block, so they enter together or not at all
    expect_equal(res$block, c(1L, 2L, 3L, 3L, 4L))
  })

  ## ---- hockey stick ---------------------------------------------------------

  ## UPPER-CASE names: .vaeCovariateSearch is called below directly, and it
  ## expects the normalized frame its caller builds -- with lower-case names it
  ## finds no ID column and silently discovers no covariates at all
  hkData <- function(wt = c(50, 60, 70, 80, 90, 100, 110, 120, 130, 140)) {
    data.frame(ID = rep(seq_along(wt), each = 2), TIME = rep(0:1, length(wt)),
               DV = seq_len(2 * length(wt)), WT = rep(wt, each = 2))
  }

  test_that("hockey emits one column per arm, blocked together", {
    res <- vaeCovariates(hkData(), shapes = c("lin", "hockey"))
    expect_equal(res$covariate, c("WT_lin", "WT_hockeyLow", "WT_hockeyHi"))
    expect_equal(res$shape, c("lin", "hockeyLow", "hockeyHi"))
    ## all three compete for one slot (one GROUP), but the two arms are one
    ## all-or-none BLOCK -- that pairing is what stops the search selecting
    ## `lin + one arm`, which spans the same space for the same penalty
    expect_equal(res$group, c(1L, 1L, 1L))
    expect_equal(res$block, c(1L, 2L, 2L))
    ## the knot is the covariate's centering value, nothing new
    expect_equal(res$center, rep(stats::median(hkData()$WT), 3))
  })

  test_that("the arms partition the subjects and sum to the lin column", {
    d <- hkData()
    s <- .vaeCovariateSearch(d, unique(d$ID), .vaeResolveShapes(c("lin", "hockey"))$rules)
    m <- s$covMat
    expect_equal(unname(m[, "WT_hockeyLow"] + m[, "WT_hockeyHi"]),
                 unname(m[, "WT_lin"]))
    ## disjoint: no subject contributes to both arms
    expect_true(all(m[, "WT_hockeyLow"] == 0 | m[, "WT_hockeyHi"] == 0))
    ## the design column is exactly what the written text evaluates to
    WT <- unique(d$WT)
    for (arm in c("hockeyLow", "hockeyHi")) {
      .j <- match(paste0("WT_", arm), s$covNames)
      expect_equal(eval(str2lang(s$covExpr[.j])), unname(m[, .j]), info = arm)
    }
  })

  test_that("hockey is searched by default and can be turned off", {
    expect_true(any(grepl("hockey", vaeCovariates(hkData())$covariate)))
    expect_true("hockey" %in% .vaeDefaultShapes)
    expect_true("hockey" %in% .vaeContShapes)
    expect_silent(.vaeAssertContShapes(c("lin", "hockey")))
    ## naming shapes without it is how a user opts out
    expect_false(any(grepl("hockey",
                           vaeCovariates(hkData(), shapes = c("power", "lin"))$covariate)))
    ## the literal defaults in the exported signatures must not drift from
    ## .vaeDefaultShapes -- they are three copies of one decision
    expect_equal(eval(formals(vaeControl)$shapes), .vaeDefaultShapes)
    expect_equal(eval(formals(vaeCovariates)$shapes), .vaeDefaultShapes)
  })

  test_that("every non-hockey column is its own block", {
    ## block ids must reduce to the historic per-column search when nothing asks
    ## for a multi-column shape, or the search changes for everyone
    d <- data.frame(id = rep(1:6, each = 2), time = rep(0:1, 6), dv = 1:12,
                    wt = rep(c(70, 80, 60, 75, 90, 65), each = 2),
                    sex = rep(c(0, 1, 0, 1, 1, 0), each = 2))
    res <- vaeCovariates(d, shapes = c("power", "lin", "log", "identity", "center"))
    expect_false(any(grepl("hockey", res$covariate)))
    expect_equal(res$block, seq_len(nrow(res)))
  })

  test_that("a knot with (almost) nothing on one side drops hockey", {
    ## the arm would be a column of zeros, making the least-squares M-step
    ## singular.  The median splits the subjects in half, so this needs a
    ## covCenter= override to reach.
    expect_warning(res <- vaeCovariates(hkData(), shapes = c("lin", "hockey"),
                                        covCenter = c(WT = 200)),
                   "hockey skipped")
    expect_equal(res$covariate, "WT_lin")
    ## a knot below every subject is the same failure on the other arm
    expect_warning(vaeCovariates(hkData(), shapes = c("lin", "hockey"),
                                 covCenter = c(WT = 10)), "hockey skipped")
    ## asking for hockey ALONE falls back to the linear family rather than
    ## dropping the covariate from the search entirely
    expect_warning(res <- vaeCovariates(hkData(), shapes = "hockey",
                                        covCenter = c(WT = 200)),
                   "hockey skipped")
    expect_equal(res$covariate, "WT_lin")
    expect_equal(res$shape, "lin")
    ## a usable knot warns about nothing
    expect_silent(vaeCovariates(hkData(), shapes = c("lin", "hockey")))
  })

  test_that("an empty arm is dropped even when catCutoff admits it", {
    ## catCutoff = 0 is documented as "test every level", and against a bare
    ## proportion `0 >= 0` is TRUE -- so a knot outside the data range would
    ## emit an all-zero arm and make the least-squares M-step singular
    expect_warning(res <- vaeCovariates(hkData(), shapes = "hockey",
                                        covCenter = c(WT = 200),
                                        catCutoff = 0),
                   "hockey skipped")
    expect_equal(res$covariate, "WT_lin")
    ## catCutoff still sets the threshold above that floor: 10 subjects, 2 below
    ## a knot of 70 -- kept at 20%, dropped at 30%
    expect_silent(vaeCovariates(hkData(), shapes = c("lin", "hockey"),
                                covCenter = c(WT = 70), catCutoff = 0.2))
    expect_warning(vaeCovariates(hkData(), shapes = c("lin", "hockey"),
                                 covCenter = c(WT = 70), catCutoff = 0.3),
                   "hockey skipped")
  })

  test_that("a list-format shapes= may ask for hockey", {
    res <- vaeCovariates(hkData(), shapes = list(WT = c("lin", "hockey")))
    expect_equal(res$shape, c("lin", "hockeyLow", "hockeyHi"))
    res <- vaeCovariates(hkData(),
                         shapes = list(list(var = "ka", covar = "wt",
                                            shapes = "hockey")))
    expect_equal(res$shape, c("hockeyLow", "hockeyHi"))
  })

  test_that("the hockey runInfo note fits on one line", {
    ## $runInfo renders one bullet per warning; CLAUDE.md caps these at 75 chars
    d <- hkData()
    .cov <- .vaeCovariateSearch(d, unique(d$ID),
                                .vaeResolveShapes(c("lin", "hockey"))$rules,
                                covCenter = c(WT = 200))
    expect_equal(.cov$hockeyDrop, "WT")
    .pre <- "<5% of subjects one side of knot, hockey skipped: "
    expect_lte(nchar(paste0(.pre, .vaeTruncList(.cov$hockeyDrop, prefix = .pre))),
               75L)
  })

  test_that("restricting shapes= restricts the candidate columns", {
    d <- data.frame(
      id = rep(1:4, each = 3), time = rep(0:2, 4), dv = 1:12,
      wt = rep(c(70, 80, 60, 75), each = 3))
    res <- vaeCovariates(d, shapes = "power", covCenterType = "mean")
    ## a single shape family reproduces the historic one-column-per-covariate
    ## design, mean-centered
    expect_equal(res$covariate, "WT_power")
    expect_equal(res$center, mean(c(70, 80, 60, 75)))
    ## a linear-family shape names its column after the shape asked for
    expect_equal(vaeCovariates(d, shapes = "identity")$covariate, "WT_identity")
  })

  test_that("covCenter overrides the centering statistic per covariate", {
    d <- data.frame(
      id = rep(1:4, each = 3), time = rep(0:2, 4), dv = 1:12,
      wt = rep(c(70, 80, 60, 75), each = 3))
    expect_equal(unique(vaeCovariates(d, covCenter = c(wt = 70))$center), 70)
    ## an unmatched name leaves the statistic in place
    expect_equal(unique(vaeCovariates(d, covCenter = c(age = 40))$center),
                 stats::median(c(70, 80, 60, 75)))
  })

  test_that("a non-0/1 two-level covariate becomes an indicator", {
    d <- data.frame(
      id = rep(1:4, each = 3), time = rep(0:2, 4), dv = 1:12,
      grp = rep(c(1, 2, 1, 2), each = 3))       # 2 levels, but not 0/1
    res <- vaeCovariates(d)
    expect_equal(res$type, "categorical")
    ## reference is the modal level (a tie resolves to the first alphabetically,
    ## here "1"), so the indicator is for level 2 and is never centered
    expect_equal(res$covariate, "GRP_2")
    expect_equal(res$level, "2")
    expect_equal(res$center, 0)
  })

  test_that("a factor covariate expands to indicators off the modal level", {
    d <- data.frame(
      id = rep(1:10, each = 2), time = rep(0:1, 10), dv = 1:20,
      race = rep(c("W", "W", "W", "W", "W", "B", "B", "B", "A", "A"), each = 2),
      stringsAsFactors = FALSE)
    res <- vaeCovariates(d)
    ## W holds 5/10 subjects -> reference; B and A each clear catCutoff
    expect_equal(res$covariate, c("RACE_B", "RACE_A"))
    expect_equal(res$level, c("B", "A"))
    ## every level is its own exclusion group -- several levels may enter together
    expect_equal(length(unique(res$group)), 2L)
  })

  test_that("levels below catCutoff are lumped with the reference", {
    d <- data.frame(
      id = rep(1:10, each = 2), time = rep(0:1, 10), dv = 1:20,
      race = rep(c("W", "W", "W", "W", "W", "W", "B", "B", "B", "A"), each = 2),
      stringsAsFactors = FALSE)
    ## A holds 1/10 = 0.10 of subjects, below a 0.2 cutoff
    res <- vaeCovariates(d, catCutoff = 0.2)
    expect_equal(res$covariate, "RACE_B")
    ## with no cutoff both non-reference levels are tested
    expect_equal(vaeCovariates(d, catCutoff = 0)$covariate,
                 c("RACE_B", "RACE_A"))
  })

  test_that("log shapes are skipped for a non-positive covariate", {
    d <- data.frame(
      id = rep(1:4, each = 3), time = rep(0:2, 4), dv = 1:12,
      score = rep(c(-1, 0, 2, 3), each = 3))
    res <- vaeCovariates(d)
    ## log(score/ctr) is undefined, so only the linear family survives
    expect_false(any(res$shape %in% c("power", "log")))
    expect_true(all(res$type == "continuous"))
    ## and asking ONLY for a log shape still leaves the covariate searchable
    expect_equal(nrow(vaeCovariates(d, shapes = "power")), 1L)
  })

  test_that("vaeCovariates matches what .vaeDataPrep explores (theo_sd)", {
    res <- vaeCovariates(nlmixr2data::theo_sd)
    expect_equal(res$raw, rep("WT", 4L))
    expect_equal(res$type, rep("continuous", 4L))
    ## the exported view must agree with the pool the fit actually searches --
    ## asserting only on vaeCovariates() would miss any divergence
    mod <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; add.err <- 0.7; eta.cl ~ 0.1 })
      model({ ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v; cp ~ add(add.err) })
    }
    prep <- suppressWarnings(.vaeDataPrep(rxode2::assertRxUi(mod),
                                          nlmixr2data::theo_sd, vaeControl()))
    expect_equal(res$covariate, prep$covNames)
    expect_equal(res$raw, prep$covRaw)
    expect_equal(res$shape, prep$covShape)
    expect_equal(res$group, prep$covGroup)
    expect_equal(res$center, prep$covPop)
  })

  test_that("vaeCovariates warns on time-varying columns and drops them", {
    d <- data.frame(
      id = rep(1:2, each = 3),
      time = rep(0:2, 2),
      dv = 1:6,
      crcl = c(90, 85, 80, 100, 95, 90),   # varies within subject
      wt = rep(c(70, 80), each = 3))
    expect_warning(res <- vaeCovariates(d), "time-varying")
    expect_equal(unique(res$raw), "WT")
    ## warn=FALSE drops them silently
    expect_silent(res2 <- vaeCovariates(d, warn = FALSE))
    expect_equal(res2, res)
  })

  test_that("vaeCovariates returns zero rows when nothing qualifies", {
    d <- data.frame(id = rep(1:2, each = 2), time = rep(0:1, 2), dv = 1:4)
    res <- vaeCovariates(d)
    expect_equal(nrow(res), 0L)
    expect_equal(names(res), c("covariate", "raw", "shape", "level", "group",
                               "block", "type", "center"))
  })

  test_that("vaeCovariates requires an ID column", {
    expect_error(vaeCovariates(data.frame(time = 0:1, dv = 1:2)), "ID")
  })
})

## Regressions for the seventh review pass (discovery and encoding).
nmTest({
  .mk <- function(extra) {
    cbind(data.frame(id = rep(1:6, each = 3), time = rep(0:2, 6), dv = 1:18),
          extra)
  }

  test_that("a covariate missing for any subject is excluded, not crashed on", {
    ## an NA reached `all(v > 0)` and aborted the whole fit with
    ## "missing value where TRUE/FALSE needed"
    d <- .mk(data.frame(wt = rep(c(70, 80, NA, 90, 75, 65), each = 3)))
    expect_warning(res <- vaeCovariates(d), "missing values")
    expect_equal(nrow(res), 0L)
    ## a complete covariate alongside an incomplete one is still searched
    d2 <- .mk(data.frame(wt = rep(c(70, 80, NA, 90, 75, 65), each = 3),
                         age = rep(c(30, 40, 50, 60, 35, 45), each = 3)))
    res2 <- suppressWarnings(vaeCovariates(d2))
    expect_equal(unique(res2$raw), "AGE")
  })

  test_that("a non-finite covariate is excluded", {
    ## Inf survived the numeric checks and put Inf into the design column, which
    ## silently breaks the least-squares M-step
    d <- .mk(data.frame(wt = rep(c(70, 80, Inf, 90, 75, 65), each = 3)))
    expect_warning(res <- vaeCovariates(d), "missing values")
    expect_equal(nrow(res), 0L)
  })

  test_that("a two-level covariate with a missing value is excluded", {
    d <- .mk(data.frame(grp = rep(c(1, 2, NA, 1, 2, 1), each = 3)))
    expect_warning(res <- vaeCovariates(d), "missing values")
    expect_equal(nrow(res), 0L)
  })

  test_that("a covariate carrying no information yields no column", {
    ## constant -> a design column of zeros would make the OLS singular
    expect_equal(nrow(vaeCovariates(.mk(data.frame(wt = rep(70, 18))))), 0L)
    ## a single-level factor has no non-reference level
    expect_equal(nrow(vaeCovariates(.mk(data.frame(
      rc = rep("A", 18), stringsAsFactors = FALSE)))), 0L)
  })

  test_that("a shape suffix cannot be confused with a same-named data column", {
    ## continuous WT emits WT_lin; a 0/1 data column named WT_LIN must stay
    ## distinct so their coefficients do not collide
    d <- .mk(data.frame(wt = rep(c(70, 80, 60, 90, 75, 65), each = 3),
                        wt_lin = rep(c(0, 1, 0, 1, 0, 1), each = 3)))
    res <- vaeCovariates(d)
    expect_equal(anyDuplicated(res$covariate), 0L)
    expect_true(all(c("WT_lin", "WT_LIN") %in% res$covariate))
  })
})
