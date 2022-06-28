#define STRICT_R_HEADER
#include <stan/math/prim/mat/fun/Eigen.hpp> // must come before #include <RcppEigen.h>
#include "../inst/include/nlmixr2est_types.h"
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;

#include <vector>
#include <stan/math/rev/core.hpp>
#include <stan/math.hpp>
#include "llik.h"

//===============================================================
struct binomial_llik {
  const Eigen::VectorXd y_, N_;
  binomial_llik(const Eigen::VectorXd& y, const Eigen::VectorXd& N) : y_(y), N_(N) { }
  
  template <typename T>
  Eigen::Matrix<T, -1, 1> operator()(const Eigen::Matrix<T, -1, 1>& theta) const {
    Eigen::Matrix<T, -1, 1> lp(y_.size());
    for (int n = 0; n < y_.size(); ++n)
      lp[n] = binomial_log(y_[n], N_[n], theta[n]);
    return lp;
  }
};

llikFxJ llik_binomial_c(Eigen::Map<Eigen::VectorXd> y,
			  Eigen::Map<Eigen::VectorXd> N,
			  Eigen::Map<Eigen::VectorXd> params) {
    int i;
    for (i=0; i<params.size(); ++i) {
		if (params[i] > .99999) params[i] = .99999;
		if (params[i] < .00001) params[i] = .00001;
	}

	binomial_llik f(y, N);
	Eigen::VectorXd fx;
	Eigen::Matrix<double, -1, -1> J;
	stan::math::jacobian(f, params, fx, J);
  llikFxJ ret;
  ret.fx = fx;
  ret.J = J;
  return ret;
}

//===============================================================
struct poisson_llik {
  const Eigen::VectorXd y_;
  poisson_llik(const Eigen::VectorXd& y) : y_(y) { }
  
  template <typename T>
  Eigen::Matrix<T, -1, 1> operator()(const Eigen::Matrix<T, -1, 1>& theta) const {
    Eigen::Matrix<T, -1, 1> lp(y_.size());
    for (int n = 0; n < y_.size(); ++n)
      lp[n] = poisson_log(y_[n], theta[n]);
    return lp;
  }
};

llikFxJ llik_poisson(Eigen::Map<Eigen::VectorXd> y, Eigen::Map<Eigen::VectorXd> params) {
  poisson_llik f(y);
  Eigen::VectorXd fx;
  Eigen::Matrix<double, -1, -1> J;
  stan::math::jacobian(f, params, fx, J);

  llikFxJ ret;
  ret.fx = fx;
  ret.J = J;
  
  return ret;
}


//===============================================================
struct normal_llik {
  const Eigen::VectorXd y_;
  normal_llik(const Eigen::VectorXd& y) : y_(y) { }

  template <typename T>
  Eigen::Matrix<T, -1, 1> operator()(const Eigen::Matrix<T, -1, 1>& theta) const {
    T mu = theta[0];
    T sigma = theta[1];
		
    if (sigma <= 0) {
      // These are not thread safe
      Rcpp::Rcout << "Warning: sigma <= 0" <<std::endl;
      sigma = 1.0e-12;
    }
		
    Eigen::Matrix<T, -1, 1> lp(y_.size());
    for (int n = 0; n < y_.size(); ++n)
      lp[n] = normal_log(y_[n], mu, sigma);
    return lp;
  }
};

llikFxJ llik_normal(Eigen::Map<Eigen::VectorXd> y, Eigen::Map<Eigen::VectorXd> params) {
  normal_llik f(y);
  Eigen::VectorXd fx;
  Eigen::Matrix<double, -1, -1> J;
  stan::math::jacobian(f, params, fx, J);

  llikFxJ ret;
  ret.fx = fx;
  ret.J = J;
  return ret;
}


//===============================================================
struct betabinomial_llik {
	const Eigen::VectorXd y_, N_;
	betabinomial_llik(const Eigen::VectorXd& y, const Eigen::VectorXd& N) : y_(y), N_(N) { }

	template <typename T>
	Eigen::Matrix<T, -1, 1> operator()(const Eigen::Matrix<T, -1, 1>& theta) const {
		T alpha = theta[0];
		T beta  = theta[1];
		if (alpha <= 0) {
			Rcpp::Rcout << "Warning: alpha <= 0" <<std::endl;
			alpha = 1.0e-12;
		}
		if (beta <= 0) {
			Rcpp::Rcout << "Warning: beta <= 0" <<std::endl;
			beta = 1.0e-12;
		}

		Eigen::Matrix<T, -1, 1> lp(y_.size());
		for (int n = 0; n < y_.size(); ++n)
		lp[n] = beta_binomial_log(y_[n], N_[n], alpha, beta);
		return lp;
	}
};

llikFxJ llik_betabinomial(Eigen::Map<Eigen::VectorXd> y,
                          Eigen::Map<Eigen::VectorXd> N,
                          Eigen::Map<Eigen::VectorXd> params) {
  betabinomial_llik f(y, N);
  Eigen::VectorXd fx;
  Eigen::Matrix<double, -1, -1> J;
  stan::math::jacobian(f, params, fx, J);

  llikFxJ ret;
  ret.fx = fx;
  ret.J = J;
  return ret;
}

//===============================================================
struct student_t_llik {
	const Eigen::VectorXd y_;
	student_t_llik(const Eigen::VectorXd& y) : y_(y) { }

	template <typename T>
	Eigen::Matrix<T, -1, 1> operator()(const Eigen::Matrix<T, -1, 1>& theta) const {
		T nu = theta[0];
		T mu = theta[1];
		T sigma = theta[2];

		if (nu <= 0) {
			Rcpp::Rcout << "Warning: nu <= 0" <<std::endl;
			nu = 1.0e-12;	//FIXME
		}
		if (sigma <= 0) {
			Rcpp::Rcout << "Warning: sigma <= 0" <<std::endl;
			sigma = 1.0e-12;
		}

		Eigen::Matrix<T, -1, 1> lp(y_.size());
		for (int n = 0; n < y_.size(); ++n)
		lp[n] = student_t_log(y_[n], nu, mu, sigma);
		return lp;
	}
};

llikFxJ llik_student_t(Eigen::Map<Eigen::VectorXd> y,
                       Eigen::Map<Eigen::VectorXd> params) {
  student_t_llik f(y);
  Eigen::VectorXd fx;
  Eigen::Matrix<double, -1, -1> J;
  stan::math::jacobian(f, params, fx, J);

  llikFxJ ret;
  ret.fx = fx;
  ret.J = J;
  return ret;
}

//===============================================================
struct beta_llik {
	const Eigen::VectorXd y_;
	beta_llik(const Eigen::VectorXd& y) : y_(y) { }

	template <typename T>
	Eigen::Matrix<T, -1, 1> operator()(const Eigen::Matrix<T, -1, 1>& theta) const {
		T alpha = theta[0];
		T beta = theta[1];

		if (alpha <= 0) {
			Rcpp::Rcout << "Warning: alpha <= 0" <<std::endl;
			alpha = 1.0e-12;
		}
		if (beta <= 0) {
			Rcpp::Rcout << "Warning: beta <= 0" <<std::endl;
			beta = 1.0e-12;
		}

		Eigen::Matrix<T, -1, 1> lp(y_.size());
		for (int n = 0; n < y_.size(); ++n)
		lp[n] = beta_log(y_[n], alpha, beta);
		return lp;
	}
};

llikFxJ llik_beta(Eigen::Map<Eigen::VectorXd> y,
                  Eigen::Map<Eigen::VectorXd> params) {
  beta_llik f(y);
  Eigen::VectorXd fx;
  Eigen::Matrix<double, -1, -1> J;
  stan::math::jacobian(f, params, fx, J);

  llikFxJ ret;
  ret.fx = fx;
  ret.J = J;
  return ret;  
}


//===============================================================
struct neg_binomial_llik {
	const Eigen::VectorXd y_;
	neg_binomial_llik(const Eigen::VectorXd& y) : y_(y) { }

	template <typename T>
	Eigen::Matrix<T, -1, 1> operator()(const Eigen::Matrix<T, -1, 1>& theta) const {
		T alpha = theta[0];
		T beta = theta[1];

		if (alpha <= 0) {
			Rcpp::Rcout << "Warning: alpha <= 0" <<std::endl;
			alpha = 1.0e-12;
		}
		if (beta <= 0) {
			Rcpp::Rcout << "Warning: beta <= 0" <<std::endl;
			beta = 1.0e-12;
		}

		Eigen::Matrix<T, -1, 1> lp(y_.size());
		for (int n = 0; n < y_.size(); ++n)
		lp[n] = neg_binomial_log(y_[n], alpha, beta);
		return lp;
	}
};

llikFxJ llik_neg_binomial(Eigen::Map<Eigen::VectorXd> y,
                          Eigen::Map<Eigen::VectorXd> params) {
  neg_binomial_llik f(y);
  Eigen::VectorXd fx;
  Eigen::Matrix<double, -1, -1> J;
  stan::math::jacobian(f, params, fx, J);

  llikFxJ ret;
  ret.fx = fx;
  ret.J = J;
  return ret;
}
