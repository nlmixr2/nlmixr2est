#ifndef __NEARPD_H__
#define __NEARPD_H__
#if defined(__cplusplus)

mat nmNearPD(mat x, bool corr = false, bool keepDiag = false,
             bool do2eigen = true, bool doSym = false, bool doDykstra = true, bool only_values = false,
             double eig_tol   = 1e-6, double conv_tol  = 1e-7, double posd_tol  = 1e-8,
             int maxit    = 100, bool trace = false // set to TRUE (or 1 ..) to trace iterations
             );
#endif
#endif
