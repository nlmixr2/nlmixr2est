#ifndef __NPDE_H__
#define __NPDE_H__

#if defined(__cplusplus)
extern "C" {
#endif

#define PHI(x) 0.5*(1.0+erf((x)/M_SQRT2))

#define min2( a , b )  ( (a) < (b) ? (a) : (b) )
#define max2( a , b )  ( (a) > (b) ? (a) : (b) )
#define innerOde(id) ind_solve(rx, id, rxInner.dydt_liblsoda, rxInner.dydt_lsoda_dum, rxInner.jdum_lsoda, rxInner.dydt, rxInner.update_inis, rxInner.global_jt)
#define predOde(id) ind_solve(rx, id, rxPred.dydt_liblsoda, rxPred.dydt_lsoda_dum, rxPred.jdum_lsoda, rxPred.dydt, rxPred.update_inis, rxPred.global_jt)
#define getCholOmegaInv() (as<arma::mat>(rxode2::rxSymInvCholEnvCalculate(_rxInv, "chol.omegaInv", R_NilValue)))
#define getOmega() (as<NumericMatrix>(rxode2::rxSymInvCholEnvCalculate(_rxInv, "omega", R_NilValue)))
#define getOmegaMat() (as<arma::mat>(rxode2::rxSymInvCholEnvCalculate(_rxInv, "omega", R_NilValue)))
#define getOmegaInv() (as<arma::mat>(rxode2::rxSymInvCholEnvCalculate(_rxInv, "omegaInv", R_NilValue)))
#define getOmegaDet() (as<double>(rxode2::rxSymInvCholEnvCalculate(_rxInv, "log.det.OMGAinv.5", R_NilValue)))
#define getOmegaN() as<int>(rxode2::rxSymInvCholEnvCalculate(_rxInv, "ntheta", R_NilValue))
#define getOmegaTheta() as<NumericVector>(rxode2::rxSymInvCholEnvCalculate(_rxInv, "theta", R_NilValue));
#define setOmegaTheta(x) rxode2::rxSymInvCholEnvCalculate(_rxInv, "theta", x)
#define tbs(x) _powerD(x,    ind->lambda, (int)(ind->yj), ind->logitLow, ind->logitHi)
#define tbsL(x) _powerL(x,   ind->lambda, (int)(ind->yj), ind->logitLow, ind->logitHi)
#define tbsDL(x) _powerDL(x, ind->lambda, (int)(ind->yj), ind->logitLow, ind->logitHi)
#define tbsD(x) _powerDD(x,  ind->lambda, (int)(ind->yj), ind->logitLow, ind->logitHi)
#define _safe_log(a) (((a) <= 0.0) ? log(DBL_EPSILON) : log(a))
// #define _safe_log(a) log(a)
#define _safe_zero(a) ((a) == 0 ? DBL_EPSILON : (a))
//#define _safe_zero(a) (a)
#define _safe_sqrt(a) ((a) <= 0 ? sqrt(DBL_EPSILON) : sqrt(a))
//#define _safe_sqrt(a) sqrt(a)
#define _as_dbleps(a) (fabs(a) < sqrt(DBL_EPSILON) ? ((a) < 0 ? -sqrt(DBL_EPSILON)  : sqrt(DBL_EPSILON)) : a)

#define expit(alpha, low, high) _powerDi(alpha, 1.0, 4, low, high)
#define probitInv(alpha, low, high) _powerDi(alpha, 1.0, 6, low, high)


  typedef void (*S2_fp) (int *, int *, double *, double *, double *, int *, float *,
                         double *, int *);
  typedef void (*n1qn1_fp)(S2_fp simul, int n[], double x[], double f[], double g[],
                           double var[], double eps[], int mode[], int niter[], int nsim[],
                           int imp[], double zm[], int izs[], float rzs[], double dzs[],
                           int id[]);

  extern n1qn1_fp n1qn1_;

  typedef double optimfn(int n, double *par, void *ex);

  typedef void optimgr(int n, double *par, double *gr, void *ex);

  void lbfgsbRX(int n, int lmm, double *x, double *lower,
                double *upper, int *nbd, double *Fmin, optimfn fn,
                optimgr gr, int *fail, void *ex, double factr,
                double pgtol, int *fncount, int *grcount,
                int maxit, char *msg, int trace, int nREPORT);


  typedef void (*ind_solve_t)(rx_solve *rx, unsigned int cid, t_dydt_liblsoda dydt_lls,
                              t_dydt_lsoda_dum dydt_lsoda, t_jdum_lsoda jdum,
                              t_dydt c_dydt, t_update_inis u_inis, int jt);
  extern ind_solve_t ind_solve;
  typedef int (*par_progress_t)(int c, int n, int d, int cores, clock_t t0, int stop);
  extern par_progress_t par_progress;
  typedef rx_solve* (*getRxSolve_t)();
  typedef int (*isRstudio_t)();
  extern isRstudio_t isRstudio;
  extern getRxSolve_t getRx;
  typedef const char *(*rxGetId_t)(int id);
  extern rxGetId_t rxGetId;
  typedef double (*getTime_t)(int idx, rx_solving_options_ind *ind);
  extern getTime_t getTimeF;
  typedef void (*sortIds_t)(rx_solve* rx, int ini);
  extern sortIds_t sortIdsF;

typedef int (*iniSubjectI_t)(int solveid, int inLhs, rx_solving_options_ind *ind, rx_solving_options *op, rx_solve *rx,
                             t_update_inis u_inis);

extern iniSubjectI_t iniSubjectI;


extern bool assignFn_;

extern void lin_cmt_stanC(double *obs_timeD, const int nobs, double *dose_timeD, const int ndose, double *doseD, double *TinfD,
                          double *paramsD, const int oral, const int infusion, const int ncmt, const int parameterization,
                          const int neta, double *fxD, double *dvdxD, double *fpD);


#if defined(__cplusplus)
}
extern List _rxInv;
extern Environment baseEnv;
extern Function doCall;
extern Function gillRfn_;


#endif

#endif
