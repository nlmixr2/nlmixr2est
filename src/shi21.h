#ifndef __SHI21_H__
#define __SHI21_H__
#if defined(__cplusplus)

using namespace arma;

typedef arma::vec (*shi21fn_type)(arma::vec &t, int id);

// Default bounds on a searched step.
#define shi21hMaxDefault 0.5
#define shi21hMinDefault 1e-4

double shi21Forward(shi21fn_type f, arma::vec &t, double &h,
                    arma::vec &f0, arma::vec &gr, int id, int idx,
                    double ef = 7e-7, double rl = 1.5, double ru = 6.0,
                    int maxiter=15, double hMax = shi21hMaxDefault,
                    double hMin = shi21hMinDefault);

double shi21Central(shi21fn_type f, arma::vec &t, double &h,
                    arma::vec &f0, arma::vec &gr, int id, int idx,
                    double ef = 7e-7, double rl = 1.5, double ru = 6.0,
                    double nu = 8.0,
                    int maxiter=15, double hMax = shi21hMaxDefault,
                    double hMin = shi21hMinDefault);

// Hessian by finite differences of the gradient `grad`, whose value at x is gr0:
// column k differences grad along coordinate k, then the matrix is symmetrized.
//
// hh[k] is coordinate k's step.  While it is <= 0, column k comes from a shi21Forward
// or shi21Central search (ef, maxiter, hMax and hMin[k] tune it; a NULL hMin means
// shi21hMinDefault), which stores the step it settles on in hh[k].  A negative hh[k]
// starts that search from -hh[k] instead of the default step.  With a step,
// column k is a forward or central difference; when one leg is non-finite it is the
// other one-sided difference, and when both are it stays 0.
//
// type is shi21HessForward or shi21HessCentral; any other value returns zeros without
// evaluating grad.  x is perturbed in place and put back exactly on every path.  Any
// other state grad writes is left at its last evaluation: restoring that is the
// caller's job.
#define shi21HessForward 1
#define shi21HessCentral 2
arma::mat shi21Hessian(shi21fn_type grad, arma::vec &x, arma::vec &gr0, int id,
                       int type, double *hh, double ef, int maxiter,
                       double hMax = shi21hMaxDefault, const double *hMin = NULL);


// 2/sqrt(3)
#define nm2divSqrt3 1.154700538379251684162 

#endif
#endif
