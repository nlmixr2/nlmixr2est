// [[Rcpp::plugins(openmp)]]
#define ARMA_DONT_PRINT_ERRORS
#define STRICT_R_HEADER
#include "armahead.h"
#include "utilc.h"

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("nlmixr2est", String)
/* replace pkg as appropriate */
#else
#define _(String) (String)
#endif

#include "inner.h"

//[[Rcpp::export]]
RObject nlmixr2Hess_(RObject thetaT, RObject fT, RObject e,
                     RObject gillInfoT){
  par_progress = (par_progress_t) R_GetCCallable("rxode2", "par_progress");
  List par(1);
  NumericVector theta = as<NumericVector>(thetaT);
  Function f = as<Function>(fT);
  List gillInfo = as<List>(gillInfoT);
  arma::mat H(theta.size(), theta.size(), fill::zeros);
  double epsI, epsJ;
  NumericVector rEpsC = as<NumericVector>(gillInfo["rEpsC"]);
  NumericVector aEpsC = as<NumericVector>(gillInfo["aEpsC"]);
  NumericVector nF = as<NumericVector>(gillInfo["f"]);
  double lastOfv=nF[0];
  int n = theta.size();
  double f1,f2,f3,f4;
  double ti, tj;
  int i, j;
  int totTick= 4*n +2*n*(n-1);
  int cur = 0, curTick=0;
  clock_t t0=clock();
  for (i=n; i--;){
    epsI = (std::fabs(theta[i])*rEpsC[i] + aEpsC[i]);
    ti = theta[i];
    theta[i] = ti + 2*epsI;
    par[0]=theta;
    f1 = as<double>(doCall(_["what"] = f, _["args"]=par, _["envir"]=e));
    cur++;
    curTick = par_progress(cur, totTick, curTick, 1, t0, 0);
    theta[i] = ti + epsI;
    par[0]=theta;
    f2 = as<double>(doCall(_["what"] = f, _["args"]=par, _["envir"]=e));
    cur++;
    curTick = par_progress(cur, totTick, curTick, 1, t0, 0);
    theta[i] = ti - epsI;
    par[0]=theta;
    f3 = as<double>(doCall(_["what"] = f, _["args"]=par, _["envir"]=e));
    cur++;
    curTick = par_progress(cur, totTick, curTick, 1, t0, 0);
    theta[i] = ti - 2*epsI;
    par[0]=theta;
    f4 = as<double>(doCall(_["what"] = f, _["args"]=par, _["envir"]=e));
    cur++;
    curTick = par_progress(cur, totTick, curTick, 1, t0, 0);
    theta[i] = ti;
    H(i,i)=(-f1+16*f2-30*lastOfv+16*f3-f4)/(12*epsI*epsI);
    for (j = i; j--;){
      epsJ = (std::fabs(theta[j])*rEpsC[j] + aEpsC[j]);
      // eps = sqrt(epsI*epsJ);// 0.5*epsI+0.5*epsJ;
      // epsI = eps;
      // epsJ = eps;
      tj = theta[j];
      theta[i] = ti + epsI;
      theta[j] = tj + epsJ;
      par[0]=theta;
      f1 = as<double>(doCall(_["what"] = f, _["args"]=par, _["envir"]=e));
      cur++;
      curTick = par_progress(cur, totTick, curTick, 1, t0, 0);
      theta[i] = ti + epsI;
      theta[j] = tj - epsJ;
      par[0]=theta;
      f2 = as<double>(doCall(_["what"] = f, _["args"]=par, _["envir"]=e));
      cur++;
      curTick = par_progress(cur, totTick, curTick, 1, t0, 0);
      theta[i] = ti - epsI;
      theta[j] = tj + epsJ;
      par[0]=theta;
      f3 = as<double>(doCall(_["what"] = f, _["args"]=par, _["envir"]=e));
      cur++;
      curTick = par_progress(cur, totTick, curTick, 1, t0, 0);
      theta[i] = ti - epsI;
      theta[j] = tj - epsJ;
      par[0]=theta;
      f4 = as<double>(doCall(_["what"] = f, _["args"]=par, _["envir"]=e));
      cur++;
      curTick = par_progress(cur, totTick, curTick, 1, t0, 0);
      H(i,j)= (f1-f2-f3+f4)/(4*epsI*epsJ);
      H(j,i) = H(i,j);
      theta[i] = ti;
      theta[j] = tj;
    }
  }
  par_progress(totTick, totTick, cur, 1, t0, 0);
  if (isRstudio()){
    RSprintf("\n");
  } else {
    RSprintf("\r                                                                                \r");
  }
  return wrap(H);
}
