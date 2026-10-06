#ifndef NLMIXR2EST_LBFGSB_TASK_H
#define NLMIXR2EST_LBFGSB_TASK_H
// L-BFGS-B exit code (lbfgsb3c's `fail`) -> message, matching lbfgsb3c's own
// table; lbfgsb3Cts ignores `msg`, so the callers report this instead.
static inline const char *lbfgsbTaskName(int itask) {
  static const char *const names[29] = {
    "NEW_X",
    "START",
    "STOP",
    "FG",
    "ABNORMAL_TERMINATION_IN_LNSRCH",
    "CONVERGENCE",
    "CONVERGENCE: NORM_OF_PROJECTED_GRADIENT_<=_PGTOL",
    "CONVERGENCE: REL_REDUCTION_OF_F_<=_FACTR*EPSMCH",
    "ERROR: FTOL .LT. ZERO",
    "ERROR: GTOL .LT. ZERO",
    "ERROR: INITIAL G .GE. ZERO",
    "ERROR: INVALID NBD",
    "ERROR: N .LE. 0",
    "ERROR: NO FEASIBLE SOLUTION",
    "ERROR: STP .GT. STPMAX",
    "ERROR: STP .LT. STPMIN",
    "ERROR: STPMAX .LT. STPMIN",
    "ERROR: STPMIN .LT. ZERO",
    "ERROR: XTOL .LT. ZERO",
    "FG_LNSRCH",
    "FG_START",
    "RESTART_FROM_LNSRCH",
    "WARNING: ROUNDING ERRORS PREVENT PROGRESS",
    "WARNING: STP .eq. STPMAX",
    "WARNING: STP .eq. STPMIN",
    "WARNING: XTOL TEST SATISFIED",
    "CONVERGENCE: Parameters differences below xtol",
    "Maximum number of iterations reached",
    "ERROR: INVALID LMM"};
  if (itask < 1 || itask > 29) return "UNKNOWN";
  return names[itask - 1];
}

// lbfgsb3c's convergence code for an exit code: 0 converged, 1 maxit,
// 51 warning, 52 error, NA otherwise.
static inline int lbfgsbConvergence(int itask) {
  switch (itask) {
  case 6: case 7: case 8: case 27: return 0;
  case 28: return 1;
  case 23: case 24: case 25: case 26: return 51;
  case 9: case 10: case 11: case 12: case 13: case 14: case 15: case 16:
  case 17: case 18: case 19: case 29: return 52;
  default: return NA_INTEGER;
  }
}
#endif
