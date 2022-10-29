#ifndef __CENSEST_H__
#define __CENSEST_H__

#define STRICT_R_HEADER
#include <stdio.h>
#include <stdarg.h>
#include <thread>
#include <chrono>
#include <R_ext/Rdynload.h>
#include <RcppArmadillo.h>


// Cumulative distribution function of a standardized normal distribution (mean
// of zero, standard deviation of 1); x = (value - mean)/standard deviation
#define PHI(x) 0.5*(1.0+erf((x)/M_SQRT2))

#define hasFiniteLimit(lim) (R_FINITE(lim) && !ISNA(lim))
#define isM2(cens, lim) (cens == 0.0 && hasFiniteLimit(lim))
#define isM3orM4(cens) (cens == 1.0 || cens == -1.0)

#define _safe_zero(a) ((a) == 0 ? DBL_EPSILON : (a))
#define _safe_sqrt(a) ((a) <= 0 ? sqrt(DBL_EPSILON) : sqrt(a))
#define _as_dbleps(a) (fabs(a) < sqrt(DBL_EPSILON) ? ((a) < 0 ? -sqrt(DBL_EPSILON)  : sqrt(DBL_EPSILON)) : a)


#define censFlagM2 1
#define censFlagM3 2
#define censFlagM4 4

extern int globalCensFlag;

static inline void updateCensFlag(int censMethod) {
  if ((globalCensFlag & censMethod) == 0) {
    globalCensFlag += censMethod;
  }
}

static inline void resetCensFlag() {
  globalCensFlag=0;
}

// This is the core censoring likelihood function for nlmixr2
//
// @param cens Censoring input value, can be -1, 0, or 1, based on the CENS column
// @param limDv represents the limit captured in the DV column
// @param lim represents the upper or lower limit based on the LIMIT column
// @param ll is the log-liklihood if this was not censored
// @param f is the prediction
// @param r is the variance estimate
// @param adjLik is the likelihood adjustment flag, if 1, remove the log(2*pi) to match NONMEM's focei
//
// @return the log-likelihood considering any adjustments from censoring for a normal distribution
//
// @noRd
static inline double doCensNormal1(double cens, double limDv, double lim, double ll,
                                   double f, double r, int adjLik) {
  // Subtract log(2*pi)/2 since this is the likelihood instead of the liklihood missing 0.5*log(2*pi)
  double adj = M_LN_SQRT_2PI*((double)adjLik);
  if (isM2(cens, lim)) {
    // M2 adds likelihood even when the observation is defined
    // Note the ll is the numerator of equation 1 in Beal 2001, Ways to fit a PK model with some data below the quantification limit
    // So the negative represents the fact this is the denominator
    // The M2 can also be applied for quantitation above the limit.
    updateCensFlag(censFlagM2);
    return ll - log(1.0 - PHI(((lim<f)*2.0 - 1.0)*(lim - f)/sqrt(r))) - adj;
  } else if (isM3orM4(cens)) {
    if (hasFiniteLimit(lim)) {
      // M4 method
      double sd = sqrt(r);
      double cum1 = PHI(cens*(limDv-f)/sd);
      double cum2 = PHI(cens*(lim - f)/sd);
      updateCensFlag(censFlagM4);
      return log(cum1-cum2) - log(1.0 - cum2) - adj;
    } else {
      // M3 method
      updateCensFlag(censFlagM3);
      return log(PHI(cens*(limDv-f)/sqrt(r))) - adj;
    }
  }
  return ll;
}
// Calculate the censoring derivative
//
// @param cens Censoring input value, can be -1, 0, or 1, based on the CENS column
// @param limDv represents the limit captured in the DV column
// @param lim represents the upper or lower limit based on the LIMIT column
// @param dll is the derivative of the log-liklihood if this was not censored
// @param f is the prediction
// @param r is the variance estimate
// @param df is the prediction derivative
// @param dr is the variance derivative
// @return the derivatives of  the log-likelihood considering any adjustments from censoring for a normal distribution
//
// @noRd
static inline double dCensNormal1(double cens, double limDv, double lim,
                                  double dll,
                                  double f, double r,
                                  double df, double dr) {
  // M2 has the contribution of the dll as well as the censoring
  if (isM2(cens, lim)) { /*
  #log(1.0 - PHI(((lim<f)*2.0 - 1.0)*(lim - f)/sqrt(r))) - adj;
# Case 1 lim < f
# > rxToSE("log(1.0 - phi((lim - f)/sqrt(r)))")
# Then changed to 
# D(S("log(1-0.5*(1+erf(((lim-f(x))/sqrt(r(x)))/sqrt(2))))"),"x")
(Mul)	-1.0*exp((-1/2)*(lim - f(x))^2/r(x))*((-1/2)*sqrt(2)*Derivative(f(x), x)/sqrt(r(x)) + (-1/4)*sqrt(2)*Derivative(r(x), x)*(lim - f(x))/r(x)^(3/2))/(sqrt(pi)*(1 - 0.5*(1 + erf((1/2)*sqrt(2)*(lim - f(x))/sqrt(r(x))))))

ll=-1.0*exp((-1/2)*(lim - f)^2/r)*((-1/2)*sqrt(2)*fpm/sqrt(r) + (-1/4)*sqrt(2)*dr*(lim - f)/r^(3/2))/(sqrt(pi)*(1 - 0.5*(1 + erf((1/2)*sqrt(2)*(lim - f)/sqrt(r)))))

tmp <- rxode2(rxOptExpr("ll=-1.0*exp((-1/2)*(lim - f)^2/r)*((-1/2)*sqrt(2)*fpm/sqrt(r) + (-1/4)*sqrt(2)*dr*(lim - f)/r^(3/2))/(sqrt(pi)*(1 - 0.5*(1 + erf((1/2)*sqrt(2)*(lim - f)/sqrt(r)))))"))
 */
    if (lim < f) {
      double rx_expr_0 =lim-f;
      double rx_expr_1 = M_SQRT2;
      double rx_expr_2 =_safe_sqrt(r);
 return dll -exp((-0.5)*((rx_expr_0)*(rx_expr_0))/_safe_zero(r))*((-0.5)*rx_expr_1*df/_safe_zero(rx_expr_2)+(-0.25)*rx_expr_1*dr*(rx_expr_0)/_safe_zero(R_pow(_as_dbleps(r),(1.5))))/_safe_zero((M_SQRT_PI*(1-0.5*(1+erf((0.5)*rx_expr_1*(rx_expr_0)/_safe_zero(rx_expr_2))))));
    } else {
      /*
Case 2:

D(S("log(1-0.5*(1+erf(((f(x)-lim)/sqrt(r(x)))/sqrt(2))))"),"x")
(Mul)	-1.0*exp((-1/2)*(-lim + f)^2/r)*((1/2)*sqrt(2)*fpm/sqrt(r) + (-1/4)*sqrt(2)*dr*(-lim + f)/r^(3/2))/(sqrt(pi)*(1 - 0.5*(1 + erf((1/2)*sqrt(2)*(-lim + f)/sqrt(r)))))
*/
      double rx_expr_0 =M_SQRT2;
      double rx_expr_1 =_safe_sqrt(r);
      double rx_expr_2 =(0.5)*rx_expr_0;
  return dll -exp((-0.5)*((-lim+f)*(-lim+f))/_safe_zero(r))*(rx_expr_2*df/_safe_zero(rx_expr_1)+(-0.25)*rx_expr_0*dr*(-lim+f)/_safe_zero(R_pow(_as_dbleps(r),(1.5))))/_safe_zero((M_SQRT_PI*(1-0.5*(1+erf(rx_expr_2*(-lim+f)/_safe_zero(rx_expr_1))))));
    }
  } else if (isM3orM4(cens)) {
    // M3 and M4 has no contribution based on the "normal" likelihood slope
    if (hasFiniteLimit(lim)) {                          /*
> rxToSE("log(phi(cens*(dv-f)/sqrt(g))-phi(cens*(lim-f)/sqrt(g)))-log(1-phi(cens*(lim-f)/sqrt(g))")

Then changed to

D(S("log(0.5*(1+erf((cens*(dv-f(x))/sqrt(g(x)))/sqrt(2)))-0.5*(1+erf((cens*(lim-f(x))/sqrt(g(x)))/sqrt(2))))-log(1-0.5*(1+erf((cens*(lim-f(x))/sqrt(g(x)))/sqrt(2))))"),"x")
(Add)	(-1.0*exp((-1/2)*cens^2*(lim - f(x))^2/g(x))*((-1/2)*sqrt(2)*Derivative(f(x), x)*cens/sqrt(g(x)) + (-1/4)*sqrt(2)*Derivative(g(x), x)*cens*(lim - f(x))/g(x)^(3/2))/sqrt(pi) + 1.0*exp((-1/2)*(dv - f(x))^2*cens^2/g(x))*((-1/2)*sqrt(2)*Derivative(f(x), x)*cens/sqrt(g(x)) + (-1/4)*sqrt(2)*Derivative(g(x), x)*(dv - f(x))*cens/g(x)^(3/2))/sqrt(pi))/(-0.5*(1 + erf((1/2)*sqrt(2)*cens*(lim - f(x))/sqrt(g(x)))) + 0.5*(1 + erf((1/2)*sqrt(2)*(dv - f(x))*cens/sqrt(g(x))))) + 1.0*exp((-1/2)*cens^2*(lim - f(x))^2/g(x))*((-1/2)*sqrt(2)*Derivative(f(x), x)*cens/sqrt(g(x)) + (-1/4)*sqrt(2)*Derivative(g(x), x)*cens*(lim - f(x))/g(x)^(3/2))/(sqrt(pi)*(1 - 0.5*(1 + erf((1/2)*sqrt(2)*cens*(lim - f(x))/sqrt(g(x))))))


# Note that:
#  - Derivative(g(x), x) = dr
#  - Derivative(f(x), x) = df
#  - f(x)                = f
#  - g(x)                = r
#  - cens^2 = 1 (removed)

	(-1.0*exp((-1/2)*(lim - f)^2/r)*((-1/2)*sqrt(2)*df*cens/sqrt(r) + (-1/4)*sqrt(2)*dr*cens*(lim - f)/r^(3/2))/sqrt(pi) + 1.0*exp((-1/2)*(dv - f)^2/r)*((-1/2)*sqrt(2)*df*cens/sqrt(r) + (-1/4)*sqrt(2)*dr*(dv - f)*cens/r^(3/2))/sqrt(pi))/(-0.5*(1 + erf((1/2)*sqrt(2)*cens*(lim - f)/sqrt(r))) + 0.5*(1 + erf((1/2)*sqrt(2)*(dv - f)*cens/sqrt(r)))) + 1.0*exp((-1/2)*(lim - f)^2/r)*((-1/2)*sqrt(2)*df*cens/sqrt(r) + (-1/4)*sqrt(2)*dr*cens*(lim - f)/r^(3/2))/(sqrt(pi)*(1 - 0.5*(1 + erf((1/2)*sqrt(2)*cens*(lim - f)/sqrt(r)))))

rxOptExpr(
rxode2({ll = (-1.0*exp((-1/2)*(lim - f)^2/r)*((-1/2)*sqrt(2)*df*cens/sqrt(r) + (-1/4)*sqrt(2)*dr*cens*(lim - f)/r^(3/2))/sqrt(pi) + 1.0*exp((-1/2)*(dv - f)^2/r)*((-1/2)*sqrt(2)*df*cens/sqrt(r) + (-1/4)*sqrt(2)*dr*(dv - f)*cens/r^(3/2))/sqrt(pi))/(-0.5*(1 + erf((1/2)*sqrt(2)*cens*(lim - f)/sqrt(r))) + 0.5*(1 + erf((1/2)*sqrt(2)*(dv - f)*cens/sqrt(r)))) + 1.0*exp((-1/2)*(lim - f)^2/r)*((-1/2)*sqrt(2)*df*cens/sqrt(r) + (-1/4)*sqrt(2)*dr*cens*(lim - f)/r^(3/2))/(sqrt(pi)*(1 - 0.5*(1 + erf((1/2)*sqrt(2)*cens*(lim - f)/sqrt(r)))))}))

rx_expr_0~dv-f
rx_expr_1~lim-f
rx_expr_2~sqrt(2)
rx_expr_3~sqrt(r)
rx_expr_4~r^(1.5)
rx_expr_5~sqrt(pi)
rx_expr_6~(0.5)*rx_expr_2
rx_expr_7~(-0.5)*rx_expr_2
rx_expr_8~(-0.25)*rx_expr_2
rx_expr_9~rx_expr_8*dr
rx_expr_10~rx_expr_7*df
rx_expr_11~rx_expr_6*cens
rx_expr_12~rx_expr_9*cens
rx_expr_13~rx_expr_10*cens
rx_expr_15~rx_expr_11*(rx_expr_1)
rx_expr_16~rx_expr_13/rx_expr_3
rx_expr_18~rx_expr_12*(rx_expr_1)
rx_expr_20~rx_expr_15/rx_expr_3
rx_expr_22~rx_expr_18/rx_expr_4
rx_expr_23~erf(rx_expr_20)
rx_expr_24~1+rx_expr_23
ll=(-exp((-0.5)*((rx_expr_1)*(rx_expr_1))/r)*(rx_expr_16+rx_expr_22)/rx_expr_5+exp((-0.5)*((rx_expr_0)*(rx_expr_0))/r)*(rx_expr_16+rx_expr_9*(rx_expr_0)*cens/rx_expr_4)/rx_expr_5)/(-0.5*(rx_expr_24)+0.5*(1+erf(rx_expr_6*(rx_expr_0)*cens/rx_expr_3)))+exp((-0.5)*((rx_expr_1)*(rx_expr_1))/r)*(rx_expr_16+rx_expr_22)/(rx_expr_5*(1-0.5*(rx_expr_24)))

rxOptExpr(x)
                    */
                    double rx_expr_0 =limDv-f;
                    double rx_expr_1 =lim-f;
                    double rx_expr_2 =M_SQRT2;
                    double rx_expr_3 =sqrt(r);
                    double rx_expr_4 =R_pow(_as_dbleps(r),(1.5));
                    double rx_expr_5 =M_SQRT_PI;
                    double rx_expr_6 =(0.5)*rx_expr_2;
                    double rx_expr_7 =(-0.5)*rx_expr_2;
                    double rx_expr_8 =(-0.25)*rx_expr_2;
                    double rx_expr_9 =rx_expr_8*dr;
                    double rx_expr_10 =rx_expr_7*df;
                    double rx_expr_11 =rx_expr_6*cens;
                    double rx_expr_12 =rx_expr_9*cens;
                    double rx_expr_13 =rx_expr_10*cens;
                    double rx_expr_15 =rx_expr_11*(rx_expr_1);
                    double rx_expr_16 =rx_expr_13/_safe_zero(rx_expr_3);
                    double rx_expr_18 =rx_expr_12*(rx_expr_1);
                    double rx_expr_20 =rx_expr_15/_safe_zero(rx_expr_3);
                    double rx_expr_22 =rx_expr_18/_safe_zero(rx_expr_4);
                    double rx_expr_23 =erf(rx_expr_20);
                    double rx_expr_24 =1+rx_expr_23;
    return  (-exp((-0.5)*((rx_expr_1)*(rx_expr_1))/_safe_zero(r))*(rx_expr_16+rx_expr_22)/_safe_zero(rx_expr_5)+exp((-0.5)*((rx_expr_0)*(rx_expr_0))/_safe_zero(r))*(rx_expr_16+rx_expr_9*(rx_expr_0)*cens/_safe_zero(rx_expr_4))/_safe_zero(rx_expr_5))/_safe_zero((-0.5*(rx_expr_24)+0.5*(1+erf(rx_expr_6*(rx_expr_0)*cens/_safe_zero(rx_expr_3)))))+exp((-0.5)*((rx_expr_1)*(rx_expr_1))/_safe_zero(r))*(rx_expr_16+rx_expr_22)/_safe_zero((rx_expr_5*(1-0.5*(rx_expr_24))));
    } else {
      // M3 method
      // logLik = log(phi((QL-f(x)/sqrt(g(x)))))
      // logLik = log(1/2*(1+erf((cens*(dv-f(x)))/sqrt(g(x))))/sqrt(2))
      // > D(S("log(1/2*(1+erf((cens*(dv-f(x)))/sqrt(g(x))/M_SQRT2)))"),"x")
      // (Mul) 2*exp(-(dv - f(x))^2*cens^2/(g(x)*M_SQRT2^2))*(-Derivative(f(x), x)*cens/(sqrt(g(x))*M_SQRT2) + (-1/2)*Derivative(g(x), x)*(limDv - f(x))*cens/(g(x)^(3/2)*M_SQRT2))/(sqrt(pi)*(1 + erf((limDv - f(x))*cens/(sqrt(g(x))*M_SQRT2))))
      // 2*exp(-(limDv - f)^2*cens^2/(r*M_SQRT2^2))*(-df*cens/(sqrt(r)*M_SQRT2) + (-0.5)*dr*(limDv - f)*cens/(r^(1.5)*M_SQRT2))/(sqrt(pi)*(1 + erf((limDv - f)*cens/(sqrt(r)*M_SQRT2))))
      double rx_expr_0=limDv-f;
      double rx_expr_1=_safe_sqrt(r);
      double rx_expr_2=rx_expr_1*M_SQRT2;
      return  2*exp(-(rx_expr_0*rx_expr_0)/_safe_zero(r*2))*(-df*cens/_safe_zero(rx_expr_2)-0.5*dr*rx_expr_0*cens/_safe_zero(R_pow(r,1.5)*M_SQRT2))/(M_SQRT_PI*(1+erf(rx_expr_0*cens/_safe_zero(rx_expr_2))));
    }
  }
  return dll;
}

#undef hasFiniteLimit
#undef isM2
#undef isM3orM4
#undef PHI
#undef _safe_zero
#undef _safe_sqrt
#undef _as_dbleps

#undef censFlagM2
#undef censFlagM3
#undef censFlagM4

SEXP censEstGetFactor();

#endif // __CENSEST_H__
