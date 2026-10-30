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

#if defined(__cplusplus)
}
#endif
#endif
