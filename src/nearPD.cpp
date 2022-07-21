//https://stackoverflow.com/questions/51490499/results-for-calculating-nearest-positive-definite-matrix-are-different-in-r-func/51492402#51492402
// Modifications by M. Fidler
// Contributors Ding Li, Ralf Stubner
#define ARMA_DONT_PRINT_ERRORS
#define STRICT_R_HEADER
#include "armahead.h"
#include "nearPD.h"

using namespace arma;
using namespace Rcpp;

vec nmRepEach(const vec& x, const int each) {
  std::size_t n=x.n_elem;
  std::size_t n_out=n*each;
  vec res(n_out);
  auto begin = res.begin();
  for (std::size_t i = 0, ind = 0; i < n; ind += each, ++i) {
    auto start = begin + ind;
    auto end = start + each;
    std::fill(start, end, x[i]);
  }
  return res;
}

mat nmMatVecSameLen(mat mt1, vec v1){
  //do not check the input...
  int t=0;
  for(unsigned int i=0;i<mt1.n_cols;i++){
    for(unsigned int j=0;j<mt1.n_rows;j++){
      mt1(j,i)=mt1(j,i)*v1(t);
      t++;
    }
  }
  return(mt1);
}

vec nmPmaxC(double a, vec b){
  vec c(b.n_elem);
  for(unsigned int i=0;i<b.n_elem;i++){
    c(i)=std::max(a,b(i));
  }
  return c;
}

mat nmNearPD(mat x, 
             bool corr //= false
             , bool keepDiag// = false
             , bool do2eigen// = true  // if TRUE do a sfsmisc::posdefify() eigen step
             , bool doSym// = false // symmetrize after tcrossprod()
             , bool doDykstra// = true // do use Dykstra's correction
             , bool only_values// = false // if TRUE simply return lambda[j].
             , double eig_tol//   = 1e-6 // defines relative positiveness of eigenvalues compared to largest
             , double conv_tol//  = 1e-7 // convergence tolerance for algorithm
             , double posd_tol//  = 1e-8 // tolerance for enforcing positive definiteness
             , int maxit//    = 100 // maximum number of iterations allowed
             , bool trace// = false // set to TRUE (or 1 ..) to trace iterations
             ){

  int n = x.n_cols;
  vec diagX0;
  if(keepDiag) {
    diagX0 = x.diag();
  }
  mat D_S;
  if(doDykstra) {
    //D_S should be like x, but filled with '0' -- following also works for 'Matrix':
    D_S = x;
    D_S.zeros(); //set all element
  }

  mat X = x;
  int iter = 0 ;
  bool converged = false; 
  double conv = R_PosInf;


  mat Y;
  mat R;
  mat B;
  while (iter < maxit && !converged) {
    Y = X;
    if(doDykstra){
      R = Y - D_S;
    }

    vec d;
    mat Q;
    if(doDykstra){
      B=R;
    }else{
      B=Y;
    }

    eig_sym(d, Q, B);

    // create mask from relative positive eigenvalues
    uvec p= (d>eig_tol*d[1]);
    if(sum(p)==0){
      //stop("Matrix seems negative semi-definite")
      break;
    }

    // use p mask to only compute 'positive' part
    uvec p_indexes(sum(p));

    int p_i_i=0;
    for(unsigned int i=0;i<p.n_elem;i++){
      if(p(i)) {
        p_indexes(p_i_i)=i;
        p_i_i++;
      }
    }


    Q=Q.cols(p_indexes);

    X=nmMatVecSameLen(Q,nmRepEach(d.elem(p_indexes),Q.n_rows))*Q.t();

    // update Dykstra's correction D_S = \Delta S_k           
    if(doDykstra){
      D_S = X - R;
    }

    // project onto symmetric and possibly 'given diag' matrices:
    if(doSym){
      X = (X + X.t())/2;
    }

    if(corr){

      X.diag().ones(); //set diagnols as ones
    } 
    else if(keepDiag){
      X.diag() = diagX0;
    } 

    conv = norm(Y-X,"inf")/norm(Y,"inf");

    iter = iter + 1;
    if (trace){
      // cat(sprintf("iter %3d : #{p}=%d, ||Y-X|| / ||Y||= %11g\n",
      // iter, sum(p), conv))
      Rcpp::Rcout << "iter " << iter <<" : #{p}= "<< sum(p) << std::endl;
    }

    converged = (conv <= conv_tol);       

    // force symmetry is *NEVER* needed, we have symmetric X here!
    //X <- (X + t(X))/2
    if(do2eigen || only_values) {
      // begin from posdefify(sfsmisc)

      eig_sym(d, Q, X);

      double Eps = posd_tol * std::abs(d[n-1]);
      if (d(0) < Eps) {
        uvec d_comp = d < Eps;
        for(unsigned int i=0; i < sum(d_comp); i++){
          if(d_comp(i)){
            d(i)=Eps;
          }
        }

        // d[d < Eps] = Eps; //how to assign values likes this?
        if(!only_values) {

          vec o_diag = X.diag();
          X = Q * (d *Q.t());
          vec D = sqrt(nmPmaxC(Eps, o_diag)/X.diag());
          x=D * X * nmRepEach(D,  n);
        }
      }
      if(only_values) return(d);

      // unneeded(?!): X <- (X + t(X))/2
      if(corr) {
        X.diag().ones(); //set diag as ones
      }
      else if(keepDiag){
        X.diag()= diagX0;
      } 

    } //end from posdefify(sfsmisc)

  }

  // if(!converged){ //not converged
  //   Rcpp::Rcout << "did not converge! " <<std::endl;
  // }

  return X;
}

//[[Rxpp::export]]
RObject nmNearPD_(RObject x, 
                  bool corr = false
                  , bool keepDiag = false
                  , bool do2eigen = true  // if TRUE do a sfsmisc::posdefify() eigen step
                  , bool doSym = false // symmetrize after tcrossprod()
                  , bool doDykstra = true // do use Dykstra's correction
                  , bool only_values = false // if TRUE simply return lambda[j].
                  , double eig_tol   = 1e-6 // defines relative positiveness of eigenvalues compared to largest
                  , double conv_tol  = 1e-7 // convergence tolerance for algorithm
                  , double posd_tol  = 1e-8 // tolerance for enforcing positive definiteness
                  , int maxit    = 100 // maximum number of iterations allowed
                  , bool trace = false // set to TRUE (or 1 ..) to trace iterations
                  ){
  return wrap(nmNearPD(as<arma::mat>(x), corr, keepDiag, do2eigen, doSym,
                       doDykstra, only_values, eig_tol, conv_tol, posd_tol, maxit, trace));
}
