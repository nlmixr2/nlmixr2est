nmTest({
  # rx->ordId maps a POSITION in a parallel loop to a SUBJECT ID; the loops here
  # then have to use that id for everything -- the solve AND the per-subject
  # arrays they write (inds_focei, ebes, Hb, okv, ...).  Using the position for
  # one and the id for the other mixes subjects, which is the bug class already
  # found on the rxode2 side, where `solveid` carried both meanings across ~100
  # ind_*() drivers (nlmixr2est#1020, guarded there by
  # test-ind-solve-subject-id.R).  This is the same convention's other half.
  #
  # It cannot be provoked from R while rx->ordId is the identity, and it is the
  # identity for any solve rxode2 sets up itself -- so the convention is
  # asserted against the sources, which is also where it regresses: the call
  # sites are ~20 near-copies of one line.

  .src <- file.path("..", "..", "src")

  test_that("every ordId walk goes through foceiOrdId(), not getOrdId()", {
    skip_if(!dir.exists(.src), "source tree not available (installed package)")
    .f <- file.path(.src, "inner.cpp")
    skip_if(!file.exists(.f))
    .l <- readLines(.f, warn = FALSE)
    # the single sanctioned raw call is inside foceiOrdId()'s own body
    .raw <- grep("\\bgetOrdId\\(", .l)
    .raw <- .raw[!grepl("getRxNsim(rxIn) == 1) ? getOrdId(rxIn, pos)", .l[.raw], fixed = TRUE)]
    .raw <- .raw[!grepl("^\\s*//", .l[.raw])]
    expect_equal(sprintf("%d: %s", .raw, trimws(.l[.raw])), character(0))
  })

  test_that("no ordId walk reuses the loop position after mapping it", {
    skip_if(!dir.exists(.src), "source tree not available (installed package)")
    .f <- file.path(.src, "inner.cpp")
    skip_if(!file.exists(.f))
    .l <- readLines(.f, warn = FALSE)
    .hit <- grep("foceiOrdId\\(\\s*rx", .l)
    .hit <- .hit[!grepl("static inline int foceiOrdId", .l[.hit])]
    expect_gt(length(.hit), 0L)
    .bad <- character(0)
    for (.h in .hit) {
      # `int <id> = ... foceiOrdId(rx, <pos>) ...`
      .m <- regmatches(.l[.h],
                       regexec("int\\s+(\\w+)\\s*=.*foceiOrdId\\(\\s*\\w+\\s*,\\s*(\\w+)\\s*\\)", .l[.h]))[[1]]
      if (length(.m) != 3L) {
        .bad <- c(.bad, sprintf("%d: unrecognized mapping form: %s", .h, trimws(.l[.h])))
        next
      }
      .pos <- .m[3]
      # body: to the end of the enclosing for loop, by brace depth from the line
      # holding the mapping (depth 0 there; < 0 closes the loop)
      .depth <- 0L
      .end <- length(.l)
      for (.j in seq(.h, min(.h + 150L, length(.l)))) {
        .depth <- .depth +
          nchar(gsub("[^{]", "", .l[.j])) - nchar(gsub("[^}]", "", .l[.j]))
        if (.j > .h && .depth < 0L) { .end <- .j; break }
      }
      if (.end <= .h + 1L) next
      .body <- .l[seq(.h + 1L, .end - 1L)]
      .lines <- seq(.h + 1L, .end - 1L)
      # strip comments and string literals before looking for the position
      .code <- sub("//.*$", "", .body)
      .code <- gsub('"[^"]*"', '""', .code)
      # a shadowing redeclaration re-uses the name legitimately; stop there
      .shadow <- grep(sprintf("\\bfor\\s*\\(\\s*(int|size_t|unsigned)\\s+%s\\b", .pos), .code)
      if (length(.shadow)) {
        .keep <- seq_len(.shadow[1] - 1L)
        .code <- .code[.keep]; .lines <- .lines[.keep]
      }
      .use <- grep(sprintf("\\b%s\\b", .pos), .code)
      if (length(.use)) {
        .bad <- c(.bad, sprintf("%d: position '%s' used after mapping: %s",
                                .lines[.use], .pos, trimws(.code[.use])))
      }
    }
    expect_equal(.bad, character(0))
  })

  test_that("every ordId walk is bounded by the subject count", {
    skip_if(!dir.exists(.src), "source tree not available (installed package)")
    .f <- file.path(.src, "inner.cpp")
    skip_if(!file.exists(.f))
    .l <- readLines(.f, warn = FALSE)
    .hit <- grep("foceiOrdId\\(\\s*rx", .l)
    .hit <- .hit[!grepl("static inline int foceiOrdId", .l[.hit])]
    .bad <- character(0)
    for (.h in .hit) {
      # nearest enclosing `for (int <v> = 0; <v> < <bound>;` above the mapping
      .bound <- NA_character_
      for (.j in seq(.h - 1L, max(1L, .h - 8L))) {
        .m <- regmatches(.l[.j],
                         regexec("for\\s*\\(\\s*int\\s+\\w+\\s*=\\s*0\\s*;\\s*\\w+\\s*<\\s*([^;]+);", .l[.j]))[[1]]
        if (length(.m) == 2L) { .bound <- trimws(.m[2]); break }
      }
      if (is.na(.bound)) {
        .bad <- c(.bad, sprintf("%d: no enclosing 0-based for loop found", .h))
        next
      }
      # the bound must be the subject count: `nsub`-named, or the single
      # documented exception (a per-subject matrix's row count), or the
      # single-subject short-circuit in outerSolveFill()
      .ok <- grepl("nsub", .bound, fixed = TRUE) ||
        grepl("subEta.n_rows", .bound, fixed = TRUE)
      if (!.ok) {
        .bad <- c(.bad, sprintf("%d: loop bound '%s' is not a subject count", .h, .bound))
      }
    }
    expect_equal(.bad, character(0))
  })

  test_that("foceiOrdId falls back to the identity when nsim != 1", {
    skip_if(!dir.exists(.src), "source tree not available (installed package)")
    .f <- file.path(.src, "inner.cpp")
    skip_if(!file.exists(.f))
    .txt <- paste(readLines(.f, warn = FALSE), collapse = "\n")
    # rx->ordId is a permutation of nsub*nsim solves; these loops walk nsub
    # subjects, so reading it is only meaningful at nsim == 1
    expect_true(grepl("return (getRxNsim(rxIn) == 1) ? getOrdId(rxIn, pos) : pos + 1;",
                      .txt, fixed = TRUE))
  })
})
