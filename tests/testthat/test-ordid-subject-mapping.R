nmTest({
  # rx->ordId maps a POSITION in a parallel region to a SUBJECT ID; the region
  # then has to use that id for everything -- the solve AND the per-subject
  # arrays it writes.  Using the position for one and the id for the other
  # silently swaps subjects, which is the bug class already found on the rxode2
  # side, where `solveid` carried both meanings across ~100 ind_*() drivers
  # (nlmixr2est issue 1020, guarded there by test-ind-solve-subject-id.R).
  #
  # nmForEachSubject() (src/nmParallel.h) is now the only way to write one of
  # these regions, and it hands the body the ID -- the position is never in
  # scope, so the mixing bug cannot be expressed.  What is left to assert is
  # that nobody hand-rolls a region around it, and that the one mapping rule
  # keeps its guard.

  .src <- file.path("..", "..", "src")
  .cpp <- function() list.files(.src, "\\.(cpp|h)$", full.names = TRUE)

  test_that("no per-subject parallel region maps ordId by hand", {
    skip_if(!dir.exists(.src), "source tree not available (installed package)")
    .bad <- character(0)
    for (.f in .cpp()) {
      if (basename(.f) == "nmParallel.h") next   # the one sanctioned mapping
      .l <- readLines(.f, warn = FALSE)
      .hit <- grep("\\b(getOrdId|nmOrdId|foceiOrdId)\\s*\\(", .l)
      .hit <- .hit[!grepl("^\\s*(//|\\*|/\\*)", .l[.hit])]
      # outerSolveFill() runs one body in two modes (all subjects, or a single
      # named one), so it keeps its own loop and calls the rule directly
      .hit <- .hit[!grepl("subject >= 0 ? subject : (nmOrdId(rx, i, nsub) - 1)",
                          .l[.hit], fixed = TRUE)]
      if (length(.hit)) {
        .bad <- c(.bad, sprintf("%s:%d: %s", basename(.f), .hit, trimws(.l[.hit])))
      }
    }
    expect_equal(.bad, character(0))
  })

  test_that("the mapping is refused when the loop does not cover every solve", {
    skip_if(!dir.exists(.src), "source tree not available (installed package)")
    .f <- file.path(.src, "nmParallel.h")
    skip_if(!file.exists(.f))
    .txt <- paste(readLines(.f, warn = FALSE), collapse = "\n")
    # rx->ordId is a permutation of the nsub*nsim solves, so a loop bounded by
    # anything else must fall back to the data order rather than read a subset
    expect_true(grepl("(n == getRxNsub(rxIn) * getRxNsim(rxIn)) ? getOrdId(rxIn, pos) : pos + 1",
                      .txt, fixed = TRUE))
  })

  test_that("the thread id is set once per thread, not once per subject", {
    skip_if(!dir.exists(.src), "source tree not available (installed package)")
    .f <- file.path(.src, "nmParallel.h")
    skip_if(!file.exists(.f))
    .l <- readLines(.f, warn = FALSE)
    # it cannot change within a thread inside one region, so the `parallel` and
    # the `for` are separated to hoist it out of the iteration
    .par <- grep("^#pragma omp parallel num_threads", .l)
    .for <- grep("^#pragma omp for", .l)
    .set <- grep("setRxThreadId(omp_get_thread_num())", .l, fixed = TRUE)
    expect_equal(length(.par), 1L)
    expect_equal(length(.for), 1L)
    expect_equal(length(.set), 1L)
    expect_true(.par < .set && .set < .for)
  })

  test_that("an escaping exception cannot reach the OpenMP boundary", {
    skip_if(!dir.exists(.src), "source tree not available (installed package)")
    .f <- file.path(.src, "nmParallel.h")
    skip_if(!file.exists(.f))
    .txt <- paste(readLines(.f, warn = FALSE), collapse = "\n")
    expect_true(grepl("catch (...) {", .txt, fixed = TRUE))
    # ... but only on the parallel path: serial, the caller can still see it
    expect_true(grepl("if (!par) {", .txt, fixed = TRUE))
  })
})
