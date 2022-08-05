#ifndef __SHI21_H__
#define __SHI21_H__
#if defined(__cplusplus)

using namespace arma;

typedef arma::vec (*shi21fn_type)(arma::vec &t, int id);

double shi21Forward(shi21fn_type f, arma::vec &t, double &h,
                    arma::vec &f0, arma::vec &gr, int id, int idx,
                    double ef = 7e-7, double rl = 1.5, double ru = 6.0,
                    int maxiter=15);

double shi21Central(shi21fn_type f, arma::vec &t, double &h,
                    arma::vec &f0, arma::vec &gr, int id, int idx,
                    double ef = 7e-7, double rl = 1.5, double ru = 6.0,
                    double nu = 8.0,
                    int maxiter=15);

double shi21Stencil(shi21fn_type f, arma::vec &t, double &h,
                    arma::vec &f0, arma::vec &gr, int id, int idx,
                    double ef, double rl, double ru, double nu,
                    int maxiter);


// 2/sqrt(3)
#define nm2divSqrt3 1.154700538379251684162 

#endif
#endif
