#ifndef __FDHESS_H__
#define __FDHESS_H__
#if defined(__cplusplus)

// The objective a finite-difference Hessian differences: f(x) at the full
// parameter vector x, and restore(x0), which re-installs the base point after
// the last probe (nothing to do for an objective without state).
struct FdHessObj {
  virtual double f(double *x) = 0;
  virtual void restore(double *x0) {}
  virtual ~FdHessObj() {}
};

#endif
#endif
