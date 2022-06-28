#ifndef __LLIK_H__
#define __LLIK_H__

#if defined(__cplusplus)

typedef struct llikFxJ {
  Eigen::VectorXd fx;
  Eigen::Matrix<double, -1, -1> J;
} llikFxJ;


llikFxJ llik_binomial_c(Eigen::Map<Eigen::VectorXd> y,
                        Eigen::Map<Eigen::VectorXd> N,
                        Eigen::Map<Eigen::VectorXd> params);

llikFxJ llik_normal(Eigen::Map<Eigen::VectorXd> y,
                    Eigen::Map<Eigen::VectorXd> params);

llikFxJ llik_betabinomial(Eigen::Map<Eigen::VectorXd> y,
                          Eigen::Map<Eigen::VectorXd> N,
                          Eigen::Map<Eigen::VectorXd> params);

llikFxJ llik_student_t(Eigen::Map<Eigen::VectorXd> y,
                       Eigen::Map<Eigen::VectorXd> params);

llikFxJ llik_beta(Eigen::Map<Eigen::VectorXd> y,
                  Eigen::Map<Eigen::VectorXd> params);

llikFxJ llik_neg_binomial(Eigen::Map<Eigen::VectorXd> y,
                          Eigen::Map<Eigen::VectorXd> params);

#endif // c++
#endif // __LLIK_H__
