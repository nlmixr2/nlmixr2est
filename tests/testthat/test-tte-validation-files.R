nmTest({
  skip_if_not(.nlmixr2estRxode2HasTte(), "rxode2 without tte() endpoints")
  .dir <- system.file("validation", "tte", package = "nlmixr2est")
  skip_if(.dir == "", "tte validation scenarios not installed")
  .env <- new.env()
  suppressMessages({
    sys.source(file.path(.dir, "scenarios.R"), envir = .env)
    sys.source(file.path(.dir, "reference.R"), envir = .env)
  })
  .files <- function(name) {
    .s <- .env$tteValidationScenarios(name)[[1]]
    .out <- name
    if (is.function(.s$reduced)) {
      .out <- c(.out, paste0(name, "-reduced"))
    }
    .out
  }
  # starting values of a scenario file, in NONMEM order
  .start <- function(file) {
    .s <- suppressMessages(.env$.tteRefScenario(file))
    .nm <- suppressMessages(.env$.tteParNames(.s))
    .ini <- suppressMessages(rxode2::rxode2(.s$model))$iniDf
    .st <- .s$start
    .fix <- .ini$name[.ini$fix]
    .st[.fix] <- .ini$est[match(.fix, .ini$name)]
    list(theta = .st[.nm$theta], eta = .st[.nm$eta], fix = .fix, s = .s)
  }

  test_that("NONMEM control streams start from the scenario starting values", {
    for (.name in names(.env$tteValidationScenarios())) {
      for (.file in .files(.name)) {
        .l <- readLines(file.path(.dir, "nonmem", paste0(.file, ".ctl")))
        .st <- .start(.file)
        .num <- function(rec) {
          .x <- sub(";.*$", "", grep(paste0("^\\$", rec), .l, value = TRUE))
          .fixed <- grepl("FIX", .x)
          .x <- gsub("\\$[A-Z]+|FIX|[()]", "", .x)
          # (lower, init) -> init
          .x <- vapply(strsplit(.x, ","), function(v) as.numeric(trimws(v[length(v)])), numeric(1))
          attr(.x, "fixed") <- .fixed
          .x
        }
        .th <- .num("THETA")
        expect_equal(as.vector(.th), unname(.st$theta), tolerance = 1e-6, info = .file)
        expect_equal(attr(.th, "fixed"), names(.st$theta) %in% .st$fix, info = .file)
        .om <- .num("OMEGA")
        if (length(.st$eta) == 0) {
          expect_equal(as.vector(.om), 0, info = .file)
        } else {
          expect_equal(as.vector(.om), unname(.st$eta), tolerance = 1e-6, info = .file)
        }
        expect_true(any(grepl(paste0("^\\$DATA ", sub("-reduced$", "", .file), ".csv"), .l)), info = .file)
      }
    }
  })

  test_that("Monolix projects start from the scenario starting values", {
    for (.name in names(.env$tteValidationScenarios())) {
      .f <- file.path(.dir, "monolix", paste0(.name, ".mlxtran"))
      if (isFALSE(.env$tteValidationScenarios(.name)[[1]]$monolix)) {
        expect_false(file.exists(.f), info = .name)
        next
      }
      for (.file in .files(.name)) {
        .l <- readLines(file.path(.dir, "monolix", paste0(.file, ".mlxtran")))
        .st <- .start(.file)
        .mu <- suppressMessages(rxode2::rxode2(.st$s$model))$muRefDataFrame
        .p <- grep("^[A-Za-z0-9_]+ = \\{value=", .l, value = TRUE)
        .val <- stats::setNames(
          as.numeric(sub(".*value=([^,]+),.*", "\\1", .p)),
          sub(" = .*", "", .p)
        )
        .method <- stats::setNames(sub(".*method=([A-Z]+)\\}.*", "\\1", .p), names(.val))
        for (.t in names(.st$theta)) {
          .m <- .st$s$monolixNames[[.t]]
          if (is.null(.m)) {
            .m <- paste0(.t, "_pop")
          }
          expect_equal(.val[[.m]], unname(.st$theta[[.t]]), tolerance = 1e-6, info = paste(.file, .t))
          expect_equal(.method[[.m]] == "FIXED", .t %in% .st$fix, info = paste(.file, .t))
        }
        for (.e in names(.st$eta)) {
          .m <- paste0("omega_", .mu$theta[.mu$eta == .e])
          expect_equal(.val[[.m]]^2, unname(.st$eta[[.e]]), tolerance = 1e-6, info = paste(.file, .e))
        }
        expect_true(any(grepl(paste0("../data/", sub("-reduced$", "", .file), ".csv"), .l, fixed = TRUE)), info = .file)
      }
    }
  })

  test_that("NONMEM and Monolix results read into nlmixr2 parameter names", {
    .s <- suppressMessages(.env$.tteRefScenario("tvCov-reduced"))
    .tmp <- withr::local_tempdir()
    .ext <- file.path(.tmp, "run.ext")
    writeLines(
      c(
        "TABLE NO.     1: Laplacian Conditional Estimation: Goal Function=MINIMUM VALUE OF OBJECTIVE FUNCTION",
        " ITERATION    THETA1       THETA2       OMEGA(1,1)   OBJ",
        "          0 -1.60944E+00  0.00000E+00  1.00000E-01    2500.0",
        " -1000000000 -2.27000E+00  0.00000E+00  2.90000E-01    2443.9",
        " -1000000001  9.50000E-02  1.00000E+10  8.40000E-02    0"
      ),
      .ext
    )
    .lst <- file.path(.tmp, "run.lst")
    writeLines("1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.5.1", .lst)
    .r <- suppressMessages(.env$tteNonmemRead(.ext, .lst, .s))
    expect_equal(.r$parameter, c("lh", "bx", "eta.h", "OFV"))
    expect_equal(.r$estimate, c(-2.27, 0, 0.29, 2443.9))
    expect_equal(.r$se, c(0.095, NA, 0.084, NA))
    expect_equal(unique(.r$version), "7.5.1")

    .s <- suppressMessages(.env$tteValidationScenarios("pkJoint")[[1]])
    # lixoftConnectors returns data frames in newer versions, named vectors in older
    .se <- data.frame(
      parameter = c("tcl_pop", "tv_pop", "b", "llam_pop", "lic50_pop", "omega_tcl", "omega_tv", "omega_llam"),
      se = c(0.03, 0.03, 0.003, 0.08, 0.2, 0.02, 0.03, 0.05),
      rse = 1
    )
    .se <- .env$.tteMlxVector(.se, "se")
    expect_equal(.env$.tteMlxVector(as.list(.se)), .se)
    .pop <- c(
      tcl_pop = 0.01,
      tv_pop = 2.3,
      b = 0.1,
      llam_pop = -2.2,
      lic50_pop = 0.6,
      omega_tcl = 0.3,
      omega_tv = 0.4,
      omega_llam = 0.5
    )
    .ll <- .env$.tteMlxVector(data.frame(criteria = c("OFV", "AIC"), importanceSampling = c(4000.5, 4016.5)))
    expect_equal(.ll, c(OFV = 4000.5, AIC = 4016.5))
    .r <- suppressMessages(.env$tteMonolixReference(.pop, .se, .ll[["OFV"]], .s, "2024R1"))
    expect_equal(.r$parameter, c("tcl", "tv", "prop.sd", "llam", "lic50", "eta.cl", "eta.v", "eta.lam", "OFV"))
    expect_equal(.r$estimate, c(0.01, 2.3, 0.1, -2.2, 0.6, 0.09, 0.16, 0.25, 4000.5))
    # the SE of a variance from the SE of its SD
    expect_equal(.r$se[6:8], 2 * c(0.3, 0.4, 0.5) * c(0.02, 0.03, 0.05))

    # no random effects
    .s <- suppressMessages(.env$tteValidationScenarios("exponential")[[1]])
    .r <- suppressMessages(.env$tteMonolixReference(c(lh_pop = -2.97), c(lh_pop = 0.09), 1001.5, .s, "2024R1"))
    expect_equal(.r$parameter, c("lh", "OFV"))
    expect_equal(.r$estimate, c(-2.97, 1001.5))
  })
})
