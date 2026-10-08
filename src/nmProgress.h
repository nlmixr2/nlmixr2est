#ifndef __NMPROGRESS_H__
#define __NMPROGRESS_H__

#include <ctime>
#include <RcppArmadillo.h>
#include <rxode2ptr.h>
#include "utilc.h"

// par_progress() draws 100% only once per bar and clears that latch on a tick at 0,
// so a bar whose first tick is past 0 inherits the last bar's latch and stops short
// of 100% (#1171).  Start every bar here.
static inline int nmProgressStart(int totTick, clock_t t0) {
  return par_progress(0, totTick, 0, 1, t0, 0);
}

// Whether the bar redraws in place (\r) and so leaves the cursor on its own line;
// the line-printing style (non-interactive, RStudio) ends itself with a newline.
static inline bool nmProgressInPlace() {
  if (isRstudio()) return false;
  Rcpp::Environment rxNs = Rcpp::Environment::namespace_env("rxode2");
  Rcpp::Function getProgSupported = rxNs["getProgSupported"];
  return Rcpp::as<int>(getProgSupported()) == 1;
}

// Complete a bar started with nmProgressStart() and end its line.  Evaluates no R,
// so it is safe from a destructor.
static inline int nmProgressEnd(int totTick, int curTick, clock_t t0, bool inPlace) {
  curTick = par_progress(totTick, totTick, curTick, 1, t0, 0);
  if (inPlace) RSprintf("\n");
  return curTick;
}

#endif // __NMPROGRESS_H__
