#define STRICT_R_HEADER
#include "shrink.h"

//[[Rcpp::export]]
List nlmixr2Parameters(NumericVector theta, DataFrame eta) {
  DataFrame eta2 = eta;
  eta2.erase(0);
  NumericVector pred(theta.size()+eta2.ncol());
  CharacterVector predn(pred.size());
  NumericMatrix ipred(eta.nrow(),pred.size());
  unsigned int j;
  for (j=theta.size(); j--;){
    pred[j] = theta[j];
    std::fill_n(ipred.begin()+eta.nrow()*j,eta.nrow(),theta[j]);
    predn[j] = "THETA[" + std::to_string(j+1) + "]";
  }
  unsigned int k = theta.size();
  unsigned int n, n1;
  double M1, M2, M3, M4, delta, delta_n, delta_n2, term1;
  // Create data frame for ETA/EPS mean, sd, variance, kurtosis and skewness
  List etadf(eta2.ncol()+1);
  NumericVector cur1(9);
  etadf[eta2.ncol()]=cur1;
  for (j=eta2.ncol();j--;){
    pred[k+j]=0;
    predn[k+j] = "ETA[" + std::to_string(j+1) + "]";
    NumericVector cur = as<NumericVector>(eta2[j]);
    std::copy(cur.begin(),cur.end(),ipred.begin()+(j+k)*eta.nrow());
    // This uses welford's method as described by
    // https://www.johndcook.com/blog/skewness_kurtosis/ adapted to
    // this problem.
    n =0; n1 = 0;
    M1 =  M2 = M3 = M4 =0;
    NumericVector stat(9);
    for (unsigned int i = eta.nrow(); i--;){
      n1=n; n++;
      delta = cur[i] - M1;
      delta_n = delta / n;
      delta_n2 = delta_n * delta_n;
      term1 = delta * delta_n * n1;
      M1 += delta_n;
      M4 += term1 * delta_n2 * (n*n - 3*n + 3) + 6 * delta_n2 * M2 - 4 * delta_n * M3;
      M3 += term1 * delta_n * (n - 2) - 3 * delta_n * M2;
      M2 += term1;
    }
    stat[0] = M1;
    stat[1] = M2/((double)n1);
    stat[2] = sqrt((double)stat[1]);
    stat[3] = sqrt((double)(n)) * M3/ pow(M2, 1.5);
    stat[4] = (double)(n)*M4 / (M2*M2) - 3.0;
    etadf[j] = stat;
  }
  pred.names() = predn;
  ipred.attr("dimnames") = List::create(R_NilValue,predn);
  return List::create(_["pred"]=pred, _["ipred"]=ipred,_["eta.lst"]=etadf);
}

void calcShrinkFinalize(arma::mat &omegaMat, int &nid, List& etaLst, arma::vec &iwres, arma::Col<int> &evid,
			CharacterVector &etaNames, int doIwres) {
  double tc = sqrt((double)nid);
  double om;
  int j;
  unsigned int neta = omegaMat.n_rows;
  for (j = neta;j--;){
    NumericVector cur = etaLst[j];
    // Finalize by adding shrinkage.
    om =omegaMat(j,j);
    cur[5]= (1.0-cur[1]/om)*100.0;
    cur[6]= (1.0-cur[2]/sqrt(om))*100.0;
    cur[7] = cur[0]*tc/(cur[2]);
    cur[8] = 2*Rf_pt(-fabs(cur[7]),(double)(nid-1),1,0);
  }
  unsigned int n =0, n1 = 0;
  double M1 =  0, M2 = 0, M3 = 0, M4 =0, term1, delta, delta_n, delta_n2;
  CharacterVector dimN2(etaNames.size()+1);
  std::copy(etaNames.begin(),etaNames.end(),dimN2.begin());
  dimN2[etaNames.size()] = "IWRES";
  if (doIwres) {
    NumericVector stat=etaLst[neta];
    for (unsigned int i = iwres.size(); i--;){
      if (evid[i] == 0 && !ISNA(iwres[i])){
        n1=n; n++;
        delta = iwres[i] - M1;
        delta_n = delta / n;
        delta_n2 = delta_n * delta_n;
        term1 = delta * delta_n * n1;
        M1 += delta_n;
        M4 += term1 * delta_n2 * (n*n - 3*n + 3) + 6 * delta_n2 * M2 - 4 * delta_n * M3;
        M3 += term1 * delta_n * (n - 2) - 3 * delta_n * M2;
        M2 += term1;
      }
    }
    stat[0] = M1;
    stat[1] = M2/((double)n1);
    stat[2] = sqrt(stat[1]);
    stat[3] = sqrt((double)(n)) * M3/ pow(M2, 1.5);
    stat[4] = (double)(n)*M4 / (M2*M2) - 3.0;
    stat[5] = (1-stat[1])*100;
    stat[6] = (1-stat[2])*100;
    stat[7] = M1*sqrt((double)n)/stat[2];
    stat[8] = 2*Rf_pt(stat[7],(double)n1,1,0);
  } else {
    NumericVector stat=etaLst[neta];
    stat[0] = NA_REAL;
    stat[1] = NA_REAL;
    stat[2] = NA_REAL;
    stat[3] = NA_REAL;
    stat[4] = NA_REAL;
    stat[5] = NA_REAL;
    stat[6] = NA_REAL;
    stat[7] = NA_REAL;
    stat[8] = NA_REAL;
  }
  etaLst.attr("names") = dimN2;
  etaLst.attr("row.names") = CharacterVector::create("mean","var","sd","skewness", "kurtosis","var shrinkage (%)", "sd shrinkage (%)","t statistic","p-value");
  etaLst.attr("class") = "data.frame";
}

extern "C" SEXP _nlmixr2est_calcShrinkOnly(SEXP omegaMatSEXP, SEXP etaLstSEXP, SEXP nidSEXP) {
BEGIN_RCPP
  // These are not needed because IWRES shrinkage isn't calculated
 arma::vec iwres;
 arma::Col<int> evid;
 arma::mat omegaMat = as<arma::mat>(omegaMatSEXP);
 int nid = INTEGER(nidSEXP)[0];
 CharacterVector etaNames = VECTOR_ELT(Rf_getAttrib(omegaMatSEXP, R_DimNamesSymbol), 1);
 List etaLst = as<List>(etaLstSEXP);
 calcShrinkFinalize(omegaMat, nid, etaLst, iwres, evid,  etaNames, 0);
 return wrap(etaLst);
END_RCPP
}
