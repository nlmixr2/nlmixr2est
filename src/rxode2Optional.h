#ifndef NLMIXR2EST_RXODE2_OPTIONAL_H
#define NLMIXR2EST_RXODE2_OPTIONAL_H
// rxode2 entry points newer than the oldest rxode2 nlmixr2est supports (its
// DESCRIPTION).  Each is taken from rxode2's pointer table at load, by its
// slot and its name, only when the installed rxode2's table has it
// (iniRxode2Optional() in inner.cpp); otherwise it is NULL, and its callers
// keep the behaviour they have without it.  .rxode2Optional() says which ones
// the installed rxode2 has.  To add one: a typedef and an extern here, its
// definition and a row of rxode2OptionalSlots in inner.cpp.
#if defined(__cplusplus)
extern "C" {
#endif

  // setIndNewind() (slot 104): the newind the next calc_lhs of a subject reads.
  // Remove with test-dependency-workarounds.R's slot-104 block.
  typedef void (*nmSetIndNewind_t)(rx_solving_options_ind *ind, int newind);
  extern nmSetIndNewind_t nmSetIndNewind;

  // rxode2's output-record walk (slots 105-107) and keep columns / factor
  // templates (slots 108-112): the rows, dose columns and keep columns
  // rxSolve() lays out, for the pooled post-fit tables (tablePool.cpp).
  // Remove with test-dependency-workarounds.R's slot-105 block.
  typedef int (*nmDfWalkInit_t)(void *state, int stateSize, rx_solve *rx,
                                rx_solving_options_ind *ind, int doDose, int countOnly);
  extern nmDfWalkInit_t nmDfWalkInit;
  typedef int (*nmDfWalkStep_t)(void *state, int i, double *out, int nOut);
  extern nmDfWalkStep_t nmDfWalkStep;
  typedef void (*nmDfWalkDull_t)(const void *state, int *dull);
  extern nmDfWalkDull_t nmDfWalkDull;
  typedef int (*nmDfKeepSetup_t)(rx_solve *rx);
  extern nmDfKeepSetup_t nmDfKeepSetup;
  typedef int (*nmDfKeepFirst_t)(rx_solve *rx, int solveid);
  extern nmDfKeepFirst_t nmDfKeepFirst;
  typedef void (*nmDfKeepRow_t)(rx_solve *rx, rx_solving_options_ind *ind, int solveid,
                                int fid, int i, double t, double *out);
  extern nmDfKeepRow_t nmDfKeepRow;
  typedef SEXP (*nmDfKeepCols_t)(SEXP vals);
  extern nmDfKeepCols_t nmDfKeepCols;
  typedef SEXP (*nmDfLevels_t)(const char *item, rx_solve *rx, R_xlen_t nrow);
  extern nmDfLevels_t nmDfLevels;
  // rxDfSolveRefused() (slot 113): whether rxSolve() refuses the current solve
  // instead of returning it.  Remove with test-dependency-workarounds.R's
  // slot-113 block.
  typedef int (*nmDfSolveRefused_t)(rx_solve *rx);
  extern nmDfSolveRefused_t nmDfSolveRefused;

// rxode2.h's RXDF_* codes, for an rxode2 whose headers predate them (the
// slots above are then NULL and never called)
#ifndef RXDF_OUT_N
#define RXDF_BAD -2
#define RXDF_SKIP -1
#define RXDF_LHS 0
#define RXDF_ROW 1
#define RXDF_KNOT 2
#define RXDF_OUT_LHS_TIME 0
#define RXDF_OUT_TIME 1
#define RXDF_OUT_RESETNO 2
#define RXDF_OUT_EVID 3
#define RXDF_OUT_CMT 4
#define RXDF_OUT_SS 5
#define RXDF_OUT_AMT 6
#define RXDF_OUT_RATE 7
#define RXDF_OUT_DUR 8
#define RXDF_OUT_II 9
#define RXDF_OUT_SIMERR 10
#define RXDF_OUT_NEWIND 11
#define RXDF_OUT_N 12
#endif

#if defined(__cplusplus)
}
#endif
#endif
