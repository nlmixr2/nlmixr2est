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
#pragma omp atomic update
  globalCensFlag |= censMethod;
}

static inline void resetCensFlag() {
  globalCensFlag=0;
}

static inline double logspace_sub(double logx, double logy) {
  if (logy > logx) {
    double tmp = logx;
    logx = logy;
    logy = tmp;
  }
  if (!R_FINITE(logy)) {
    return logx;
  }
  return logx + log1p(-exp(logy - logx));
}

static inline double finite_log_prob(double logp) {
  if (R_FINITE(logp)) {
    return logp;
  }
  return log(DBL_MIN);
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
    double z = ((lim<f)*2.0 - 1.0)*(lim - f)/sqrt(r);
    return ll - finite_log_prob(R::pnorm5(z, 0.0, 1.0, 0, 1)) - adj;
  } else if (isM3orM4(cens)) {
    if (hasFiniteLimit(lim)) {
      // M4 method
      double sd = sqrt(r);
      double logCum1 = R::pnorm5(cens*(limDv-f)/sd, 0.0, 1.0, 1, 1);
      double logCum2 = R::pnorm5(cens*(lim - f)/sd, 0.0, 1.0, 1, 1);
      double logTail = R::pnorm5(cens*(lim - f)/sd, 0.0, 1.0, 0, 1);
      updateCensFlag(censFlagM4);
      return finite_log_prob(logspace_sub(logCum1, logCum2)) -
        finite_log_prob(logTail) - adj;
    } else {
      // M3 method
      updateCensFlag(censFlagM3);
      return finite_log_prob(R::pnorm5(cens*(limDv-f)/sqrt(r), 0.0, 1.0, 1, 1)) - adj;
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

// Analytic partials of the censored normal objective rho = -logLik-per-obs (kernel
// orientation, transformed scale) w.r.t. the prediction f and variance r, for M2/M3/M4.
// out[0..8] = rho_f, rho_r, rho_ff, rho_fr, rho_rr, rho_fff, rho_ffr, rho_frr, rho_rrr
// (order<3 fills only 0..4).  These are the second/third derivatives the analytic FOCEI
// gradient/covariance (and the exact censored inner Laplace Hessian) need -- the pieces
// dCensNormal1 (first derivative) does not provide.  cens/limDv/lim/f/r follow doCensNormal1.
//
// DERIVED + VERIFIED by tools/genCensPartials.R (symengine D of rho, rxOptExpr CSE, then
// translated to C); every partial validated vs numDeriv of doCensNormal1 in f,r to ~1e-6
// for |z|<3.5 (deep tails share doCensNormal1's erf-form cancellation; _safe_zero guards
// keep them finite).  Regenerate with that script; do NOT hand-edit the body below.
static inline void censNormalPartials(double cens, double limDv, double lim,
                                      double f, double r, int order, double* out) {
  double dv = limDv;
  int hasFin = (R_FINITE(lim) && !ISNA(lim));
  if (cens == 0.0 && hasFin) {            // M2
    if (lim < f) {
      double rx_expr_0 = (dv - f);
      double rx_expr_2 = 1.414213562373095;
      double rx_expr_3 = sqrt(r);
      double rx_expr_6 = sqrt(M_PI);
      double rx_expr_7 = rx_expr_2;
      double rx_expr_10 = (0.500000000000000 * rx_expr_2);
      double rx_expr_15 = (rx_expr_10 * (((-f) + lim)));
      double rx_expr_16 = (((((-f) + lim)) * (((-f) + lim))));
      double rx_expr_17 = (-0.500000000000000 * rx_expr_16);
      double rx_expr_18 = (rx_expr_15 / _safe_zero(rx_expr_3));
      double rx_expr_20 = (rx_expr_17 / _safe_zero(r));
      double rx_expr_21 = erf(rx_expr_18);
      double rx_expr_22 = exp(rx_expr_20);
      double rx_expr_23 = (1.000000000000000 + rx_expr_21);
      double rx_expr_24 = (rx_expr_7 * rx_expr_22);
      out[0] = (((-(rx_expr_0)) / _safe_zero(r)) - (rx_expr_24 / _safe_zero((((rx_expr_3 * rx_expr_6) * (rx_expr_23))))));
      double rx_expr_1 = R_pow(r, -1.000000000000000);
      double rx_expr_4 = R_pow(r, 1.500000000000000);
      double rx_expr_8 = (((r) * (r)));
      double rx_expr_9 = (0.500000000000000 * rx_expr_2);
      double rx_expr_12 = (rx_expr_4 * rx_expr_6);
      double rx_expr_25 = (rx_expr_9 * rx_expr_22);
      double rx_expr_26 = (rx_expr_12 * (rx_expr_23));
      out[1] = ((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) / _safe_zero(rx_expr_8)) - ((rx_expr_25 * (((-f) + lim))) / _safe_zero((rx_expr_26)))) + (0.500000000000000 * rx_expr_1));
      double rx_expr_19 = exp(((-rx_expr_16) / _safe_zero(r)));
      out[2] = ((((-2.000000000000000 * rx_expr_19) / _safe_zero((((r * M_PI) * (((rx_expr_23) * (rx_expr_23))))))) - ((rx_expr_24 * (((-f) + lim))) / _safe_zero((rx_expr_26)))) + rx_expr_1);
      double rx_expr_5 = R_pow(r, 2.500000000000000);
      double rx_expr_13 = (rx_expr_5 * rx_expr_6);
      double rx_expr_27 = (rx_expr_13 * (rx_expr_23));
      out[3] = (((((rx_expr_0) / _safe_zero(rx_expr_8)) + (rx_expr_25 / _safe_zero((rx_expr_26)))) - ((rx_expr_19 * (((-f) + lim))) / _safe_zero((((rx_expr_8 * M_PI) * (((rx_expr_23) * (rx_expr_23)))))))) - ((rx_expr_25 * rx_expr_16) / _safe_zero((rx_expr_27))));
      double rx_expr_11 = ((((r) * (r)) * (r)));
      out[4] = ((((((((rx_expr_0) * (rx_expr_0))) / _safe_zero(rx_expr_11)) - (((0.500000000000000 * rx_expr_19) * rx_expr_16) / _safe_zero((((rx_expr_11 * M_PI) * (((rx_expr_23) * (rx_expr_23)))))))) - ((((0.250000000000000 * rx_expr_2) * rx_expr_22) * ((((((-f) + lim)) * (((-f) + lim))) * (((-f) + lim))))) / _safe_zero((((R_pow(r, 3.500000000000000) * rx_expr_6) * (rx_expr_23)))))) + ((((0.750000000000000 * rx_expr_2) * rx_expr_22) * (((-f) + lim))) / _safe_zero((rx_expr_27)))) - (0.500000000000000 * R_pow(r, -2.000000000000000)));
      if (order >= 3) {
        double rx_expr_1 = 1.414213562373095;
        double rx_expr_2 = R_pow(r, 1.500000000000000);
        double rx_expr_3 = sqrt(r);
        double rx_expr_4 = R_pow(r, 2.500000000000000);
        double rx_expr_7 = R_pow(M_PI, 1.500000000000000);
        double rx_expr_8 = sqrt(M_PI);
        double rx_expr_9 = rx_expr_1;
        double rx_expr_10 = (((r) * (r)));
        double rx_expr_13 = (rx_expr_10 * M_PI);
        double rx_expr_14 = (0.500000000000000 * rx_expr_1);
        double rx_expr_16 = (rx_expr_4 * rx_expr_8);
        double rx_expr_22 = (rx_expr_14 * (((-f) + lim)));
        double rx_expr_23 = (((((-f) + lim)) * (((-f) + lim))));
        double rx_expr_24 = (-1.500000000000000 * rx_expr_23);
        double rx_expr_25 = (rx_expr_22 / _safe_zero(rx_expr_3));
        double rx_expr_26 = (-0.500000000000000 * rx_expr_23);
        double rx_expr_27 = exp(((-rx_expr_23) / _safe_zero(r)));
        double rx_expr_28 = (rx_expr_24 / _safe_zero(r));
        double rx_expr_29 = (rx_expr_26 / _safe_zero(r));
        double rx_expr_31 = erf(rx_expr_25);
        double rx_expr_32 = exp(rx_expr_28);
        double rx_expr_33 = (1.000000000000000 + rx_expr_31);
        double rx_expr_34 = exp(rx_expr_29);
        double rx_expr_36 = (rx_expr_9 * rx_expr_34);
        double rx_expr_38 = (rx_expr_16 * (rx_expr_33));
        out[5] = ((((((-4.000000000000000 * rx_expr_1) * rx_expr_32) / _safe_zero((((rx_expr_2 * rx_expr_7) * ((((rx_expr_33) * (rx_expr_33)) * (rx_expr_33))))))) + (rx_expr_36 / _safe_zero((((rx_expr_2 * rx_expr_8) * (rx_expr_33)))))) - (((6.000000000000000 * rx_expr_27) * (((-f) + lim))) / _safe_zero(((rx_expr_13 * (((rx_expr_33) * (rx_expr_33)))))))) - ((rx_expr_36 * rx_expr_23) / _safe_zero((rx_expr_38))));
        double rx_expr_5 = R_pow(r, 3.500000000000000);
        double rx_expr_11 = (0.500000000000000 * rx_expr_1);
        double rx_expr_12 = (1.500000000000000 * rx_expr_1);
        double rx_expr_15 = ((((r) * (r)) * (r)));
        double rx_expr_17 = (rx_expr_15 * M_PI);
        double rx_expr_18 = (rx_expr_5 * rx_expr_8);
        double rx_expr_30 = ((((((-f) + lim)) * (((-f) + lim))) * (((-f) + lim))));
        double rx_expr_37 = (rx_expr_12 * rx_expr_34);
        double rx_expr_39 = (rx_expr_18 * (rx_expr_33));
        out[6] = (((((((2.000000000000000 * rx_expr_27) / _safe_zero(((rx_expr_13 * (((rx_expr_33) * (rx_expr_33))))))) - (((3.000000000000000 * rx_expr_27) * rx_expr_23) / _safe_zero(((rx_expr_17 * (((rx_expr_33) * (rx_expr_33)))))))) - (((rx_expr_11 * rx_expr_34) * rx_expr_30) / _safe_zero((rx_expr_39)))) - ((((2.000000000000000 * rx_expr_1) * rx_expr_32) * (((-f) + lim))) / _safe_zero((((rx_expr_4 * rx_expr_7) * ((((rx_expr_33) * (rx_expr_33)) * (rx_expr_33)))))))) + ((rx_expr_37 * (((-f) + lim))) / _safe_zero((rx_expr_38)))) - R_pow(r, -2.000000000000000));
        double rx_expr_0 = (dv - f);
        double rx_expr_6 = R_pow(r, 4.500000000000000);
        double rx_expr_19 = (rx_expr_6 * rx_expr_8);
        double rx_expr_20 = (((((r) * (r)) * (r)) * (r)));
        double rx_expr_21 = (rx_expr_20 * M_PI);
        double rx_expr_35 = (((((((-f) + lim)) * (((-f) + lim))) * (((-f) + lim))) * (((-f) + lim))));
        double rx_expr_40 = (rx_expr_19 * (rx_expr_33));
        out[7] = ((((((((-2.000000000000000 * (rx_expr_0)) / _safe_zero(rx_expr_15)) - (((0.750000000000000 * rx_expr_1) * rx_expr_34) / _safe_zero((rx_expr_38)))) - (((1.500000000000000 * rx_expr_27) * rx_expr_30) / _safe_zero(((rx_expr_21 * (((rx_expr_33) * (rx_expr_33)))))))) + (((2.500000000000000 * rx_expr_27) * (((-f) + lim))) / _safe_zero(((rx_expr_17 * (((rx_expr_33) * (rx_expr_33)))))))) - ((((0.250000000000000 * rx_expr_1) * rx_expr_34) * rx_expr_35) / _safe_zero((rx_expr_40)))) - (((rx_expr_9 * rx_expr_32) * rx_expr_23) / _safe_zero((((rx_expr_5 * rx_expr_7) * ((((rx_expr_33) * (rx_expr_33)) * (rx_expr_33)))))))) + ((rx_expr_37 * rx_expr_23) / _safe_zero((rx_expr_39))));
        out[8] = (((((((((-3.000000000000000 * (((rx_expr_0) * (rx_expr_0)))) / _safe_zero(rx_expr_20)) - (((0.750000000000000 * rx_expr_27) * rx_expr_35) / _safe_zero((((((((((r) * (r)) * (r)) * (r)) * (r))) * M_PI) * (((rx_expr_33) * (rx_expr_33)))))))) + (((2.250000000000000 * rx_expr_27) * rx_expr_23) / _safe_zero(((rx_expr_21 * (((rx_expr_33) * (rx_expr_33)))))))) - ((((0.125000000000000 * rx_expr_1) * rx_expr_34) * ((((((((-f) + lim)) * (((-f) + lim))) * (((-f) + lim))) * (((-f) + lim))) * (((-f) + lim))))) / _safe_zero((((R_pow(r, 5.500000000000000) * rx_expr_8) * (rx_expr_33)))))) - (((rx_expr_11 * rx_expr_32) * rx_expr_30) / _safe_zero((((rx_expr_6 * rx_expr_7) * ((((rx_expr_33) * (rx_expr_33)) * (rx_expr_33)))))))) + ((((1.250000000000000 * rx_expr_1) * rx_expr_34) * rx_expr_30) / _safe_zero((rx_expr_40)))) - ((((1.875000000000000 * rx_expr_1) * rx_expr_34) * (((-f) + lim))) / _safe_zero((rx_expr_39)))) + R_pow(r, -3.000000000000000));
      }
    } else {
      double rx_expr_0 = (dv - f);
      double rx_expr_2 = 1.414213562373095;
      double rx_expr_3 = sqrt(r);
      double rx_expr_6 = sqrt(M_PI);
      double rx_expr_7 = rx_expr_2;
      double rx_expr_10 = (0.500000000000000 * rx_expr_2);
      double rx_expr_15 = (rx_expr_10 * (((-f) + lim)));
      double rx_expr_16 = (((((-f) + lim)) * (((-f) + lim))));
      double rx_expr_17 = (-0.500000000000000 * rx_expr_16);
      double rx_expr_18 = (rx_expr_15 / _safe_zero(rx_expr_3));
      double rx_expr_20 = (rx_expr_17 / _safe_zero(r));
      double rx_expr_21 = erf(rx_expr_18);
      double rx_expr_22 = exp(rx_expr_20);
      double rx_expr_23 = (1.000000000000000 - rx_expr_21);
      double rx_expr_24 = (rx_expr_7 * rx_expr_22);
      out[0] = (((-(rx_expr_0)) / _safe_zero(r)) + (rx_expr_24 / _safe_zero((((rx_expr_3 * rx_expr_6) * (rx_expr_23))))));
      double rx_expr_1 = R_pow(r, -1.000000000000000);
      double rx_expr_4 = R_pow(r, 1.500000000000000);
      double rx_expr_8 = (((r) * (r)));
      double rx_expr_9 = (0.500000000000000 * rx_expr_2);
      double rx_expr_12 = (rx_expr_4 * rx_expr_6);
      double rx_expr_25 = (rx_expr_9 * rx_expr_22);
      double rx_expr_26 = (rx_expr_12 * (rx_expr_23));
      out[1] = ((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) / _safe_zero(rx_expr_8)) + ((rx_expr_25 * (((-f) + lim))) / _safe_zero((rx_expr_26)))) + (0.500000000000000 * rx_expr_1));
      double rx_expr_19 = exp(((-rx_expr_16) / _safe_zero(r)));
      out[2] = ((((-2.000000000000000 * rx_expr_19) / _safe_zero((((r * M_PI) * (((rx_expr_23) * (rx_expr_23))))))) + ((rx_expr_24 * (((-f) + lim))) / _safe_zero((rx_expr_26)))) + rx_expr_1);
      double rx_expr_5 = R_pow(r, 2.500000000000000);
      double rx_expr_13 = (rx_expr_5 * rx_expr_6);
      double rx_expr_27 = (rx_expr_13 * (rx_expr_23));
      out[3] = (((((rx_expr_0) / _safe_zero(rx_expr_8)) - (rx_expr_25 / _safe_zero((rx_expr_26)))) - ((rx_expr_19 * (((-f) + lim))) / _safe_zero((((rx_expr_8 * M_PI) * (((rx_expr_23) * (rx_expr_23)))))))) + ((rx_expr_25 * rx_expr_16) / _safe_zero((rx_expr_27))));
      double rx_expr_11 = ((((r) * (r)) * (r)));
      out[4] = ((((((((rx_expr_0) * (rx_expr_0))) / _safe_zero(rx_expr_11)) - (((0.500000000000000 * rx_expr_19) * rx_expr_16) / _safe_zero((((rx_expr_11 * M_PI) * (((rx_expr_23) * (rx_expr_23)))))))) + ((((0.250000000000000 * rx_expr_2) * rx_expr_22) * ((((((-f) + lim)) * (((-f) + lim))) * (((-f) + lim))))) / _safe_zero((((R_pow(r, 3.500000000000000) * rx_expr_6) * (rx_expr_23)))))) - ((((0.750000000000000 * rx_expr_2) * rx_expr_22) * (((-f) + lim))) / _safe_zero((rx_expr_27)))) - (0.500000000000000 * R_pow(r, -2.000000000000000)));
      if (order >= 3) {
        double rx_expr_1 = 1.414213562373095;
        double rx_expr_2 = R_pow(r, 1.500000000000000);
        double rx_expr_3 = sqrt(r);
        double rx_expr_4 = R_pow(r, 2.500000000000000);
        double rx_expr_7 = R_pow(M_PI, 1.500000000000000);
        double rx_expr_8 = sqrt(M_PI);
        double rx_expr_9 = rx_expr_1;
        double rx_expr_10 = (((r) * (r)));
        double rx_expr_13 = (rx_expr_10 * M_PI);
        double rx_expr_14 = (0.500000000000000 * rx_expr_1);
        double rx_expr_16 = (rx_expr_4 * rx_expr_8);
        double rx_expr_22 = (rx_expr_14 * (((-f) + lim)));
        double rx_expr_23 = (((((-f) + lim)) * (((-f) + lim))));
        double rx_expr_24 = (-1.500000000000000 * rx_expr_23);
        double rx_expr_25 = (rx_expr_22 / _safe_zero(rx_expr_3));
        double rx_expr_26 = (-0.500000000000000 * rx_expr_23);
        double rx_expr_27 = exp(((-rx_expr_23) / _safe_zero(r)));
        double rx_expr_28 = (rx_expr_24 / _safe_zero(r));
        double rx_expr_29 = (rx_expr_26 / _safe_zero(r));
        double rx_expr_31 = erf(rx_expr_25);
        double rx_expr_32 = exp(rx_expr_28);
        double rx_expr_33 = (1.000000000000000 - rx_expr_31);
        double rx_expr_34 = exp(rx_expr_29);
        double rx_expr_36 = (rx_expr_9 * rx_expr_34);
        double rx_expr_38 = (rx_expr_16 * (rx_expr_33));
        out[5] = ((((((4.000000000000000 * rx_expr_1) * rx_expr_32) / _safe_zero((((rx_expr_2 * rx_expr_7) * ((((rx_expr_33) * (rx_expr_33)) * (rx_expr_33))))))) - (rx_expr_36 / _safe_zero((((rx_expr_2 * rx_expr_8) * (rx_expr_33)))))) - (((6.000000000000000 * rx_expr_27) * (((-f) + lim))) / _safe_zero(((rx_expr_13 * (((rx_expr_33) * (rx_expr_33)))))))) + ((rx_expr_36 * rx_expr_23) / _safe_zero((rx_expr_38))));
        double rx_expr_5 = R_pow(r, 3.500000000000000);
        double rx_expr_11 = (0.500000000000000 * rx_expr_1);
        double rx_expr_12 = (1.500000000000000 * rx_expr_1);
        double rx_expr_15 = ((((r) * (r)) * (r)));
        double rx_expr_17 = (rx_expr_15 * M_PI);
        double rx_expr_18 = (rx_expr_5 * rx_expr_8);
        double rx_expr_30 = ((((((-f) + lim)) * (((-f) + lim))) * (((-f) + lim))));
        double rx_expr_37 = (rx_expr_12 * rx_expr_34);
        double rx_expr_39 = (rx_expr_18 * (rx_expr_33));
        out[6] = (((((((2.000000000000000 * rx_expr_27) / _safe_zero(((rx_expr_13 * (((rx_expr_33) * (rx_expr_33))))))) - (((3.000000000000000 * rx_expr_27) * rx_expr_23) / _safe_zero(((rx_expr_17 * (((rx_expr_33) * (rx_expr_33)))))))) + (((rx_expr_11 * rx_expr_34) * rx_expr_30) / _safe_zero((rx_expr_39)))) + ((((2.000000000000000 * rx_expr_1) * rx_expr_32) * (((-f) + lim))) / _safe_zero((((rx_expr_4 * rx_expr_7) * ((((rx_expr_33) * (rx_expr_33)) * (rx_expr_33)))))))) - ((rx_expr_37 * (((-f) + lim))) / _safe_zero((rx_expr_38)))) - R_pow(r, -2.000000000000000));
        double rx_expr_0 = (dv - f);
        double rx_expr_6 = R_pow(r, 4.500000000000000);
        double rx_expr_19 = (rx_expr_6 * rx_expr_8);
        double rx_expr_20 = (((((r) * (r)) * (r)) * (r)));
        double rx_expr_21 = (rx_expr_20 * M_PI);
        double rx_expr_35 = (((((((-f) + lim)) * (((-f) + lim))) * (((-f) + lim))) * (((-f) + lim))));
        double rx_expr_40 = (rx_expr_19 * (rx_expr_33));
        out[7] = ((((((((-2.000000000000000 * (rx_expr_0)) / _safe_zero(rx_expr_15)) + (((0.750000000000000 * rx_expr_1) * rx_expr_34) / _safe_zero((rx_expr_38)))) - (((1.500000000000000 * rx_expr_27) * rx_expr_30) / _safe_zero(((rx_expr_21 * (((rx_expr_33) * (rx_expr_33)))))))) + (((2.500000000000000 * rx_expr_27) * (((-f) + lim))) / _safe_zero(((rx_expr_17 * (((rx_expr_33) * (rx_expr_33)))))))) + ((((0.250000000000000 * rx_expr_1) * rx_expr_34) * rx_expr_35) / _safe_zero((rx_expr_40)))) + (((rx_expr_9 * rx_expr_32) * rx_expr_23) / _safe_zero((((rx_expr_5 * rx_expr_7) * ((((rx_expr_33) * (rx_expr_33)) * (rx_expr_33)))))))) - ((rx_expr_37 * rx_expr_23) / _safe_zero((rx_expr_39))));
        out[8] = (((((((((-3.000000000000000 * (((rx_expr_0) * (rx_expr_0)))) / _safe_zero(rx_expr_20)) - (((0.750000000000000 * rx_expr_27) * rx_expr_35) / _safe_zero((((((((((r) * (r)) * (r)) * (r)) * (r))) * M_PI) * (((rx_expr_33) * (rx_expr_33)))))))) + (((2.250000000000000 * rx_expr_27) * rx_expr_23) / _safe_zero(((rx_expr_21 * (((rx_expr_33) * (rx_expr_33)))))))) + ((((0.125000000000000 * rx_expr_1) * rx_expr_34) * ((((((((-f) + lim)) * (((-f) + lim))) * (((-f) + lim))) * (((-f) + lim))) * (((-f) + lim))))) / _safe_zero((((R_pow(r, 5.500000000000000) * rx_expr_8) * (rx_expr_33)))))) + (((rx_expr_11 * rx_expr_32) * rx_expr_30) / _safe_zero((((rx_expr_6 * rx_expr_7) * ((((rx_expr_33) * (rx_expr_33)) * (rx_expr_33)))))))) - ((((1.250000000000000 * rx_expr_1) * rx_expr_34) * rx_expr_30) / _safe_zero((rx_expr_40)))) + ((((1.875000000000000 * rx_expr_1) * rx_expr_34) * (((-f) + lim))) / _safe_zero((rx_expr_39)))) + R_pow(r, -3.000000000000000));
      }
    }
  } else if (hasFin) {                    // M4
    double rx_expr_0 = (dv - f);
    double rx_expr_1 = 1.414213562373095;
    double rx_expr_2 = sqrt(r);
    double rx_expr_6 = sqrt(M_PI);
    double rx_expr_7 = (0.500000000000000 * rx_expr_1);
    double rx_expr_9 = (0.500000000000000 * rx_expr_1);
    double rx_expr_12 = (((cens) * (cens)));
    double rx_expr_13 = (rx_expr_2 * rx_expr_6);
    double rx_expr_17 = (rx_expr_9 * cens);
    double rx_expr_19 = (rx_expr_9 * (rx_expr_0));
    double rx_expr_20 = (-0.500000000000000 * rx_expr_12);
    double rx_expr_22 = (((((-f) + lim)) * (((-f) + lim))));
    double rx_expr_23 = (rx_expr_19 * cens);
    double rx_expr_25 = (rx_expr_17 * (((-f) + lim)));
    double rx_expr_26 = (rx_expr_23 / _safe_zero(rx_expr_2));
    double rx_expr_28 = (rx_expr_25 / _safe_zero(rx_expr_2));
    double rx_expr_29 = erf(rx_expr_26);
    double rx_expr_30 = erf(rx_expr_28);
    double rx_expr_31 = (1.000000000000000 + rx_expr_29);
    double rx_expr_33 = (1.000000000000000 + rx_expr_30);
    double rx_expr_35 = (rx_expr_20 * rx_expr_22);
    double rx_expr_36 = (0.500000000000000 * (rx_expr_31));
    double rx_expr_38 = (rx_expr_35 / _safe_zero(r));
    double rx_expr_39 = (0.500000000000000 * (rx_expr_33));
    double rx_expr_41 = (1.000000000000000 - rx_expr_39);
    double rx_expr_42 = exp(rx_expr_38);
    double rx_expr_45 = (rx_expr_7 * rx_expr_42);
    double rx_expr_53 = (rx_expr_45 * cens);
    double rx_expr_61 = (rx_expr_53 / _safe_zero((rx_expr_13)));
    out[0] = (-(((((rx_expr_61 - (((rx_expr_7 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_12) / _safe_zero(r)))) * cens) / _safe_zero((rx_expr_13))))) / _safe_zero((((-0.500000000000000 * (rx_expr_33)) + rx_expr_36)))) - (rx_expr_53 / _safe_zero(((rx_expr_13 * (rx_expr_41))))))));
    double rx_expr_3 = R_pow(r, 1.500000000000000);
    double rx_expr_8 = (0.250000000000000 * rx_expr_1);
    double rx_expr_14 = (rx_expr_3 * rx_expr_6);
    double rx_expr_47 = (rx_expr_8 * rx_expr_42);
    double rx_expr_51 = (rx_expr_14 * (rx_expr_41));
    double rx_expr_54 = (rx_expr_47 * cens);
    double rx_expr_57 = (rx_expr_54 * (((-f) + lim)));
    double rx_expr_66 = (rx_expr_57 / _safe_zero((rx_expr_14)));
    out[1] = (-(((((rx_expr_66 - ((((rx_expr_8 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_12) / _safe_zero(r)))) * (rx_expr_0)) * cens) / _safe_zero((rx_expr_14))))) / _safe_zero((((-0.500000000000000 * (rx_expr_33)) + rx_expr_36)))) - (rx_expr_57 / _safe_zero((rx_expr_51))))));
    double rx_expr_21 = ((((cens) * (cens)) * (cens)));
    double rx_expr_37 = exp((((-rx_expr_12) * rx_expr_22) / _safe_zero(r)));
    double rx_expr_60 = (rx_expr_45 * rx_expr_21);
    double rx_expr_65 = (rx_expr_60 * (((-f) + lim)));
    out[2] = (-((((((-((((rx_expr_61 - (((rx_expr_7 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_12) / _safe_zero(r)))) * cens) / _safe_zero((rx_expr_13))))) * ((rx_expr_61 - (((rx_expr_7 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_12) / _safe_zero(r)))) * cens) / _safe_zero((rx_expr_13)))))))) / _safe_zero((((((-0.500000000000000 * (rx_expr_33)) + rx_expr_36)) * (((-0.500000000000000 * (rx_expr_33)) + rx_expr_36)))))) + ((((rx_expr_65 / _safe_zero((rx_expr_14))) - ((((rx_expr_7 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_12) / _safe_zero(r)))) * (rx_expr_0)) * rx_expr_21) / _safe_zero((rx_expr_14))))) / _safe_zero((((-0.500000000000000 * (rx_expr_33)) + rx_expr_36))))) + (((0.500000000000000 * rx_expr_37) * rx_expr_12) / _safe_zero((((r * M_PI) * (((rx_expr_41) * (rx_expr_41)))))))) - (rx_expr_65 / _safe_zero((rx_expr_51))))));
    double rx_expr_4 = R_pow(r, 2.500000000000000);
    double rx_expr_15 = (rx_expr_4 * rx_expr_6);
    double rx_expr_52 = (rx_expr_15 * (rx_expr_41));
    double rx_expr_62 = (rx_expr_47 * rx_expr_21);
    out[3] = (-((((((((((((((-0.250000000000000 * rx_expr_1) * rx_expr_42) * cens) / _safe_zero((rx_expr_14))) + (((rx_expr_8 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_12) / _safe_zero(r)))) * cens) / _safe_zero((rx_expr_14)))) + ((rx_expr_62 * rx_expr_22) / _safe_zero((rx_expr_15)))) - ((((rx_expr_8 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_12) / _safe_zero(r)))) * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_21) / _safe_zero((rx_expr_15))))) / _safe_zero((((-0.500000000000000 * (rx_expr_33)) + rx_expr_36)))) - ((((rx_expr_66 - ((((rx_expr_8 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_12) / _safe_zero(r)))) * (rx_expr_0)) * cens) / _safe_zero((rx_expr_14))))) * ((rx_expr_61 - (((rx_expr_7 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_12) / _safe_zero(r)))) * cens) / _safe_zero((rx_expr_13)))))) / _safe_zero((((((-0.500000000000000 * (rx_expr_33)) + rx_expr_36)) * (((-0.500000000000000 * (rx_expr_33)) + rx_expr_36))))))) + (rx_expr_54 / _safe_zero((rx_expr_51)))) + ((((0.250000000000000 * rx_expr_37) * rx_expr_12) * (((-f) + lim))) / _safe_zero(((((((r) * (r))) * M_PI) * (((rx_expr_41) * (rx_expr_41)))))))) - ((rx_expr_62 * rx_expr_22) / _safe_zero((rx_expr_52))))));
    double rx_expr_5 = R_pow(r, 3.500000000000000);
    double rx_expr_10 = (0.125000000000000 * rx_expr_1);
    double rx_expr_11 = (0.375000000000000 * rx_expr_1);
    double rx_expr_16 = (rx_expr_5 * rx_expr_6);
    double rx_expr_27 = ((((((-f) + lim)) * (((-f) + lim))) * (((-f) + lim))));
    double rx_expr_48 = (rx_expr_10 * rx_expr_42);
    double rx_expr_49 = (rx_expr_11 * rx_expr_42);
    double rx_expr_55 = (rx_expr_49 * cens);
    double rx_expr_58 = (rx_expr_55 * (((-f) + lim)));
    double rx_expr_63 = (rx_expr_48 * rx_expr_21);
    double rx_expr_69 = (rx_expr_63 * rx_expr_27);
    out[4] = (-(((((((-((((rx_expr_66 - ((((rx_expr_8 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_12) / _safe_zero(r)))) * (rx_expr_0)) * cens) / _safe_zero((rx_expr_14))))) * ((rx_expr_66 - ((((rx_expr_8 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_12) / _safe_zero(r)))) * (rx_expr_0)) * cens) / _safe_zero((rx_expr_14)))))))) / _safe_zero((((((-0.500000000000000 * (rx_expr_33)) + rx_expr_36)) * (((-0.500000000000000 * (rx_expr_33)) + rx_expr_36)))))) + ((((((rx_expr_69 / _safe_zero((rx_expr_16))) - ((((rx_expr_10 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_12) / _safe_zero(r)))) * ((((rx_expr_0) * (rx_expr_0)) * (rx_expr_0)))) * rx_expr_21) / _safe_zero((rx_expr_16)))) - (rx_expr_58 / _safe_zero((rx_expr_15)))) + ((((rx_expr_11 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_12) / _safe_zero(r)))) * (rx_expr_0)) * cens) / _safe_zero((rx_expr_15))))) / _safe_zero((((-0.500000000000000 * (rx_expr_33)) + rx_expr_36))))) + ((((0.125000000000000 * rx_expr_37) * rx_expr_12) * rx_expr_22) / _safe_zero((((((((r) * (r)) * (r))) * M_PI) * (((rx_expr_41) * (rx_expr_41)))))))) - (rx_expr_69 / _safe_zero(((rx_expr_16 * (rx_expr_41)))))) + (rx_expr_58 / _safe_zero((rx_expr_52))))));
    if (order >= 3) {
      double rx_expr_0 = (dv - f);
      double rx_expr_1 = 1.414213562373095;
      double rx_expr_2 = sqrt(r);
      double rx_expr_3 = R_pow(r, 1.500000000000000);
      double rx_expr_4 = R_pow(r, 2.500000000000000);
      double rx_expr_7 = sqrt(M_PI);
      double rx_expr_8 = R_pow(M_PI, 1.500000000000000);
      double rx_expr_10 = (((r) * (r)));
      double rx_expr_11 = (0.500000000000000 * rx_expr_1);
      double rx_expr_12 = (rx_expr_10 * M_PI);
      double rx_expr_15 = (0.500000000000000 * rx_expr_1);
      double rx_expr_22 = (((cens) * (cens)));
      double rx_expr_23 = (rx_expr_2 * rx_expr_7);
      double rx_expr_24 = (rx_expr_3 * rx_expr_7);
      double rx_expr_25 = (rx_expr_4 * rx_expr_7);
      double rx_expr_31 = (rx_expr_15 * cens);
      double rx_expr_33 = (rx_expr_15 * (rx_expr_0));
      double rx_expr_35 = (-0.500000000000000 * rx_expr_22);
      double rx_expr_36 = ((((cens) * (cens)) * (cens)));
      double rx_expr_37 = (-1.500000000000000 * rx_expr_22);
      double rx_expr_38 = (((((-f) + lim)) * (((-f) + lim))));
      double rx_expr_39 = (rx_expr_33 * cens);
      double rx_expr_41 = (rx_expr_31 * (((-f) + lim)));
      double rx_expr_43 = (((((cens) * (cens)) * (cens)) * (cens)));
      double rx_expr_44 = (rx_expr_39 / _safe_zero(rx_expr_2));
      double rx_expr_46 = (rx_expr_41 / _safe_zero(rx_expr_2));
      double rx_expr_47 = ((((((cens) * (cens)) * (cens)) * (cens)) * (cens)));
      double rx_expr_48 = erf(rx_expr_44);
      double rx_expr_49 = erf(rx_expr_46);
      double rx_expr_50 = (1.000000000000000 + rx_expr_48);
      double rx_expr_52 = (1.000000000000000 + rx_expr_49);
      double rx_expr_55 = (rx_expr_35 * rx_expr_38);
      double rx_expr_56 = (0.500000000000000 * (rx_expr_50));
      double rx_expr_57 = (rx_expr_37 * rx_expr_38);
      double rx_expr_58 = exp((((-rx_expr_22) * rx_expr_38) / _safe_zero(r)));
      double rx_expr_59 = (rx_expr_55 / _safe_zero(r));
      double rx_expr_60 = (rx_expr_57 / _safe_zero(r));
      double rx_expr_61 = (0.500000000000000 * (rx_expr_52));
      double rx_expr_63 = (1.000000000000000 - rx_expr_61);
      double rx_expr_64 = exp(rx_expr_59);
      double rx_expr_65 = exp(rx_expr_60);
      double rx_expr_72 = (rx_expr_11 * rx_expr_64);
      double rx_expr_83 = (rx_expr_25 * (rx_expr_63));
      double rx_expr_87 = (rx_expr_72 * cens);
      double rx_expr_100 = (rx_expr_72 * rx_expr_36);
      double rx_expr_101 = (rx_expr_87 / _safe_zero((rx_expr_23)));
      double rx_expr_110 = (rx_expr_100 * (((-f) + lim)));
      double rx_expr_114 = (rx_expr_72 * rx_expr_47);
      double rx_expr_124 = (rx_expr_110 / _safe_zero((rx_expr_24)));
      out[5] = (-(((((((((2.000000000000000 * (((((rx_expr_101 - (((rx_expr_11 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_22) / _safe_zero(r)))) * cens) / _safe_zero((rx_expr_23))))) * ((rx_expr_101 - (((rx_expr_11 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_22) / _safe_zero(r)))) * cens) / _safe_zero((rx_expr_23)))))) * ((rx_expr_101 - (((rx_expr_11 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_22) / _safe_zero(r)))) * cens) / _safe_zero((rx_expr_23)))))))) / _safe_zero(((((((-0.500000000000000 * (rx_expr_52)) + rx_expr_56)) * (((-0.500000000000000 * (rx_expr_52)) + rx_expr_56))) * (((-0.500000000000000 * (rx_expr_52)) + rx_expr_56)))))) + (((((((((-0.500000000000000 * rx_expr_1) * rx_expr_64) * rx_expr_36) / _safe_zero((rx_expr_24))) + (((rx_expr_11 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_22) / _safe_zero(r)))) * rx_expr_36) / _safe_zero((rx_expr_24)))) + ((rx_expr_114 * rx_expr_38) / _safe_zero((rx_expr_25)))) - ((((rx_expr_11 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_22) / _safe_zero(r)))) * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_47) / _safe_zero((rx_expr_25))))) / _safe_zero((((-0.500000000000000 * (rx_expr_52)) + rx_expr_56))))) - (((3.000000000000000 * ((rx_expr_124 - ((((rx_expr_11 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_22) / _safe_zero(r)))) * (rx_expr_0)) * rx_expr_36) / _safe_zero((rx_expr_24)))))) * ((rx_expr_101 - (((rx_expr_11 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_22) / _safe_zero(r)))) * cens) / _safe_zero((rx_expr_23)))))) / _safe_zero((((((-0.500000000000000 * (rx_expr_52)) + rx_expr_56)) * (((-0.500000000000000 * (rx_expr_52)) + rx_expr_56))))))) - (((rx_expr_11 * rx_expr_65) * rx_expr_36) / _safe_zero((((rx_expr_3 * rx_expr_8) * ((((rx_expr_63) * (rx_expr_63)) * (rx_expr_63)))))))) + (rx_expr_100 / _safe_zero(((rx_expr_24 * (rx_expr_63)))))) + ((((1.500000000000000 * rx_expr_58) * rx_expr_43) * (((-f) + lim))) / _safe_zero(((rx_expr_12 * (((rx_expr_63) * (rx_expr_63)))))))) - ((rx_expr_114 * rx_expr_38) / _safe_zero((rx_expr_83))))));
      double rx_expr_5 = R_pow(r, 3.500000000000000);
      double rx_expr_13 = (0.250000000000000 * rx_expr_1);
      double rx_expr_14 = (0.750000000000000 * rx_expr_1);
      double rx_expr_16 = ((((r) * (r)) * (r)));
      double rx_expr_26 = (rx_expr_5 * rx_expr_7);
      double rx_expr_27 = (rx_expr_16 * M_PI);
      double rx_expr_45 = ((((((-f) + lim)) * (((-f) + lim))) * (((-f) + lim))));
      double rx_expr_74 = (rx_expr_13 * rx_expr_64);
      double rx_expr_75 = (rx_expr_14 * rx_expr_64);
      double rx_expr_85 = (rx_expr_26 * (rx_expr_63));
      double rx_expr_89 = (rx_expr_74 * cens);
      double rx_expr_94 = (rx_expr_89 * (((-f) + lim)));
      double rx_expr_102 = (rx_expr_75 * rx_expr_36);
      double rx_expr_103 = (rx_expr_74 * rx_expr_36);
      double rx_expr_111 = (rx_expr_102 * (((-f) + lim)));
      double rx_expr_112 = (rx_expr_94 / _safe_zero((rx_expr_24)));
      double rx_expr_115 = (rx_expr_74 * rx_expr_47);
      double rx_expr_134 = (rx_expr_115 * rx_expr_45);
      out[6] = (-(((((((((((((((rx_expr_134 / _safe_zero((rx_expr_26))) - ((((rx_expr_13 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_22) / _safe_zero(r)))) * ((((rx_expr_0) * (rx_expr_0)) * (rx_expr_0)))) * rx_expr_47) / _safe_zero((rx_expr_26)))) - (rx_expr_111 / _safe_zero((rx_expr_25)))) + ((((rx_expr_14 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_22) / _safe_zero(r)))) * (rx_expr_0)) * rx_expr_36) / _safe_zero((rx_expr_25))))) / _safe_zero((((-0.500000000000000 * (rx_expr_52)) + rx_expr_56)))) + (((2.000000000000000 * ((rx_expr_112 - ((((rx_expr_13 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_22) / _safe_zero(r)))) * (rx_expr_0)) * cens) / _safe_zero((rx_expr_24)))))) * ((((rx_expr_101 - (((rx_expr_11 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_22) / _safe_zero(r)))) * cens) / _safe_zero((rx_expr_23))))) * ((rx_expr_101 - (((rx_expr_11 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_22) / _safe_zero(r)))) * cens) / _safe_zero((rx_expr_23)))))))) / _safe_zero(((((((-0.500000000000000 * (rx_expr_52)) + rx_expr_56)) * (((-0.500000000000000 * (rx_expr_52)) + rx_expr_56))) * (((-0.500000000000000 * (rx_expr_52)) + rx_expr_56))))))) - ((((rx_expr_112 - ((((rx_expr_13 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_22) / _safe_zero(r)))) * (rx_expr_0)) * cens) / _safe_zero((rx_expr_24))))) * ((rx_expr_124 - ((((rx_expr_11 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_22) / _safe_zero(r)))) * (rx_expr_0)) * rx_expr_36) / _safe_zero((rx_expr_24)))))) / _safe_zero((((((-0.500000000000000 * (rx_expr_52)) + rx_expr_56)) * (((-0.500000000000000 * (rx_expr_52)) + rx_expr_56))))))) - (((2.000000000000000 * ((((((((-0.250000000000000 * rx_expr_1) * rx_expr_64) * cens) / _safe_zero((rx_expr_24))) + (((rx_expr_13 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_22) / _safe_zero(r)))) * cens) / _safe_zero((rx_expr_24)))) + ((rx_expr_103 * rx_expr_38) / _safe_zero((rx_expr_25)))) - ((((rx_expr_13 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_22) / _safe_zero(r)))) * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_36) / _safe_zero((rx_expr_25)))))) * ((rx_expr_101 - (((rx_expr_11 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_22) / _safe_zero(r)))) * cens) / _safe_zero((rx_expr_23)))))) / _safe_zero((((((-0.500000000000000 * (rx_expr_52)) + rx_expr_56)) * (((-0.500000000000000 * (rx_expr_52)) + rx_expr_56))))))) - (((0.500000000000000 * rx_expr_58) * rx_expr_22) / _safe_zero(((rx_expr_12 * (((rx_expr_63) * (rx_expr_63)))))))) + ((((0.750000000000000 * rx_expr_58) * rx_expr_43) * rx_expr_38) / _safe_zero(((rx_expr_27 * (((rx_expr_63) * (rx_expr_63)))))))) - (rx_expr_134 / _safe_zero((rx_expr_85)))) - ((((rx_expr_13 * rx_expr_65) * rx_expr_36) * (((-f) + lim))) / _safe_zero((((rx_expr_4 * rx_expr_8) * ((((rx_expr_63) * (rx_expr_63)) * (rx_expr_63)))))))) + (rx_expr_111 / _safe_zero((rx_expr_83))))));
      double rx_expr_6 = R_pow(r, 4.500000000000000);
      double rx_expr_17 = (0.375000000000000 * rx_expr_1);
      double rx_expr_18 = (0.125000000000000 * rx_expr_1);
      double rx_expr_28 = (rx_expr_6 * rx_expr_7);
      double rx_expr_29 = (((((r) * (r)) * (r)) * (r)));
      double rx_expr_34 = (rx_expr_29 * M_PI);
      double rx_expr_54 = (((((((-f) + lim)) * (((-f) + lim))) * (((-f) + lim))) * (((-f) + lim))));
      double rx_expr_77 = (rx_expr_17 * rx_expr_64);
      double rx_expr_78 = (rx_expr_18 * rx_expr_64);
      double rx_expr_86 = (rx_expr_28 * (rx_expr_63));
      double rx_expr_90 = (rx_expr_77 * cens);
      double rx_expr_95 = (rx_expr_90 * (((-f) + lim)));
      double rx_expr_104 = (rx_expr_78 * rx_expr_36);
      double rx_expr_113 = (rx_expr_95 / _safe_zero((rx_expr_25)));
      double rx_expr_116 = (rx_expr_78 * rx_expr_47);
      double rx_expr_127 = (rx_expr_104 * rx_expr_45);
      double rx_expr_135 = (rx_expr_127 / _safe_zero((rx_expr_26)));
      double rx_expr_136 = (rx_expr_116 * rx_expr_54);
      out[7] = (-((((((((((((((((((rx_expr_90 / _safe_zero((rx_expr_25))) - (((rx_expr_17 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_22) / _safe_zero(r)))) * cens) / _safe_zero((rx_expr_25)))) + (rx_expr_136 / _safe_zero((rx_expr_28)))) - ((((rx_expr_18 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_22) / _safe_zero(r)))) * (((((rx_expr_0) * (rx_expr_0)) * (rx_expr_0)) * (rx_expr_0)))) * rx_expr_47) / _safe_zero((rx_expr_28)))) - ((rx_expr_102 * rx_expr_38) / _safe_zero((rx_expr_26)))) + ((((rx_expr_14 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_22) / _safe_zero(r)))) * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_36) / _safe_zero((rx_expr_26))))) / _safe_zero((((-0.500000000000000 * (rx_expr_52)) + rx_expr_56)))) + (((2.000000000000000 * ((((rx_expr_112 - ((((rx_expr_13 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_22) / _safe_zero(r)))) * (rx_expr_0)) * cens) / _safe_zero((rx_expr_24))))) * ((rx_expr_112 - ((((rx_expr_13 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_22) / _safe_zero(r)))) * (rx_expr_0)) * cens) / _safe_zero((rx_expr_24)))))))) * ((rx_expr_101 - (((rx_expr_11 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_22) / _safe_zero(r)))) * cens) / _safe_zero((rx_expr_23)))))) / _safe_zero(((((((-0.500000000000000 * (rx_expr_52)) + rx_expr_56)) * (((-0.500000000000000 * (rx_expr_52)) + rx_expr_56))) * (((-0.500000000000000 * (rx_expr_52)) + rx_expr_56))))))) - (((2.000000000000000 * ((((((((-0.250000000000000 * rx_expr_1) * rx_expr_64) * cens) / _safe_zero((rx_expr_24))) + (((rx_expr_13 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_22) / _safe_zero(r)))) * cens) / _safe_zero((rx_expr_24)))) + ((rx_expr_103 * rx_expr_38) / _safe_zero((rx_expr_25)))) - ((((rx_expr_13 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_22) / _safe_zero(r)))) * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_36) / _safe_zero((rx_expr_25)))))) * ((rx_expr_112 - ((((rx_expr_13 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_22) / _safe_zero(r)))) * (rx_expr_0)) * cens) / _safe_zero((rx_expr_24)))))) / _safe_zero((((((-0.500000000000000 * (rx_expr_52)) + rx_expr_56)) * (((-0.500000000000000 * (rx_expr_52)) + rx_expr_56))))))) - ((((((rx_expr_135 - ((((rx_expr_18 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_22) / _safe_zero(r)))) * ((((rx_expr_0) * (rx_expr_0)) * (rx_expr_0)))) * rx_expr_36) / _safe_zero((rx_expr_26)))) - rx_expr_113) + ((((rx_expr_17 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_22) / _safe_zero(r)))) * (rx_expr_0)) * cens) / _safe_zero((rx_expr_25))))) * ((rx_expr_101 - (((rx_expr_11 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_22) / _safe_zero(r)))) * cens) / _safe_zero((rx_expr_23)))))) / _safe_zero((((((-0.500000000000000 * (rx_expr_52)) + rx_expr_56)) * (((-0.500000000000000 * (rx_expr_52)) + rx_expr_56))))))) - (rx_expr_90 / _safe_zero((rx_expr_83)))) + ((((0.375000000000000 * rx_expr_58) * rx_expr_43) * rx_expr_45) / _safe_zero(((rx_expr_34 * (((rx_expr_63) * (rx_expr_63)))))))) - ((((0.625000000000000 * rx_expr_58) * rx_expr_22) * (((-f) + lim))) / _safe_zero(((rx_expr_27 * (((rx_expr_63) * (rx_expr_63)))))))) - (rx_expr_136 / _safe_zero((rx_expr_86)))) - ((((rx_expr_18 * rx_expr_65) * rx_expr_36) * rx_expr_38) / _safe_zero((((rx_expr_5 * rx_expr_8) * ((((rx_expr_63) * (rx_expr_63)) * (rx_expr_63)))))))) + ((rx_expr_102 * rx_expr_38) / _safe_zero((rx_expr_85))))));
      double rx_expr_9 = R_pow(r, 5.500000000000000);
      double rx_expr_19 = (0.625000000000000 * rx_expr_1);
      double rx_expr_20 = (0.062500000000000 * rx_expr_1);
      double rx_expr_21 = (0.937500000000000 * rx_expr_1);
      double rx_expr_30 = (rx_expr_9 * rx_expr_7);
      double rx_expr_66 = ((((((((-f) + lim)) * (((-f) + lim))) * (((-f) + lim))) * (((-f) + lim))) * (((-f) + lim))));
      double rx_expr_79 = (rx_expr_19 * rx_expr_64);
      double rx_expr_81 = (rx_expr_20 * rx_expr_64);
      double rx_expr_82 = (rx_expr_21 * rx_expr_64);
      double rx_expr_91 = (rx_expr_82 * cens);
      double rx_expr_97 = (rx_expr_91 * (((-f) + lim)));
      double rx_expr_105 = (rx_expr_79 * rx_expr_36);
      double rx_expr_118 = (rx_expr_81 * rx_expr_47);
      double rx_expr_128 = (rx_expr_105 * rx_expr_45);
      double rx_expr_137 = (rx_expr_118 * rx_expr_66);
      out[8] = (-(((((((((((2.000000000000000 * (((((rx_expr_112 - ((((rx_expr_13 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_22) / _safe_zero(r)))) * (rx_expr_0)) * cens) / _safe_zero((rx_expr_24))))) * ((rx_expr_112 - ((((rx_expr_13 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_22) / _safe_zero(r)))) * (rx_expr_0)) * cens) / _safe_zero((rx_expr_24)))))) * ((rx_expr_112 - ((((rx_expr_13 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_22) / _safe_zero(r)))) * (rx_expr_0)) * cens) / _safe_zero((rx_expr_24)))))))) / _safe_zero(((((((-0.500000000000000 * (rx_expr_52)) + rx_expr_56)) * (((-0.500000000000000 * (rx_expr_52)) + rx_expr_56))) * (((-0.500000000000000 * (rx_expr_52)) + rx_expr_56)))))) + ((((((((rx_expr_137 / _safe_zero((rx_expr_30))) - ((((rx_expr_20 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_22) / _safe_zero(r)))) * ((((((rx_expr_0) * (rx_expr_0)) * (rx_expr_0)) * (rx_expr_0)) * (rx_expr_0)))) * rx_expr_47) / _safe_zero((rx_expr_30)))) - (rx_expr_128 / _safe_zero((rx_expr_28)))) + ((((rx_expr_19 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_22) / _safe_zero(r)))) * ((((rx_expr_0) * (rx_expr_0)) * (rx_expr_0)))) * rx_expr_36) / _safe_zero((rx_expr_28)))) + (rx_expr_97 / _safe_zero((rx_expr_26)))) - ((((rx_expr_21 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_22) / _safe_zero(r)))) * (rx_expr_0)) * cens) / _safe_zero((rx_expr_26))))) / _safe_zero((((-0.500000000000000 * (rx_expr_52)) + rx_expr_56))))) - (((3.000000000000000 * ((((rx_expr_135 - ((((rx_expr_18 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_22) / _safe_zero(r)))) * ((((rx_expr_0) * (rx_expr_0)) * (rx_expr_0)))) * rx_expr_36) / _safe_zero((rx_expr_26)))) - rx_expr_113) + ((((rx_expr_17 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_22) / _safe_zero(r)))) * (rx_expr_0)) * cens) / _safe_zero((rx_expr_25)))))) * ((rx_expr_112 - ((((rx_expr_13 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_22) / _safe_zero(r)))) * (rx_expr_0)) * cens) / _safe_zero((rx_expr_24)))))) / _safe_zero((((((-0.500000000000000 * (rx_expr_52)) + rx_expr_56)) * (((-0.500000000000000 * (rx_expr_52)) + rx_expr_56))))))) + ((((0.187500000000000 * rx_expr_58) * rx_expr_43) * rx_expr_54) / _safe_zero((((((((((r) * (r)) * (r)) * (r)) * (r))) * M_PI) * (((rx_expr_63) * (rx_expr_63)))))))) - ((((0.562500000000000 * rx_expr_58) * rx_expr_22) * rx_expr_38) / _safe_zero(((rx_expr_34 * (((rx_expr_63) * (rx_expr_63)))))))) - (rx_expr_137 / _safe_zero(((rx_expr_30 * (rx_expr_63)))))) - ((((rx_expr_20 * rx_expr_65) * rx_expr_36) * rx_expr_45) / _safe_zero((((rx_expr_6 * rx_expr_8) * ((((rx_expr_63) * (rx_expr_63)) * (rx_expr_63)))))))) + (rx_expr_128 / _safe_zero((rx_expr_86)))) - (rx_expr_97 / _safe_zero((rx_expr_85))))));
    }
  } else {                                // M3
    double rx_expr_0 = (dv - f);
    double rx_expr_1 = 1.414213562373095;
    double rx_expr_2 = sqrt(r);
    double rx_expr_5 = sqrt(M_PI);
    double rx_expr_6 = rx_expr_1;
    double rx_expr_8 = (0.500000000000000 * rx_expr_1);
    double rx_expr_9 = (((cens) * (cens)));
    double rx_expr_13 = (rx_expr_8 * (rx_expr_0));
    double rx_expr_15 = (rx_expr_13 * cens);
    double rx_expr_17 = (rx_expr_15 / _safe_zero(rx_expr_2));
    double rx_expr_18 = erf(rx_expr_17);
    double rx_expr_19 = (1.000000000000000 + rx_expr_18);
    out[0] = (((rx_expr_6 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_9) / _safe_zero(r)))) * cens) / _safe_zero((((rx_expr_2 * rx_expr_5) * (rx_expr_19)))));
    double rx_expr_3 = R_pow(r, 1.500000000000000);
    double rx_expr_7 = (0.500000000000000 * rx_expr_1);
    double rx_expr_10 = (rx_expr_3 * rx_expr_5);
    double rx_expr_25 = (rx_expr_10 * (rx_expr_19));
    out[1] = ((((rx_expr_7 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_9) / _safe_zero(r)))) * (rx_expr_0)) * cens) / _safe_zero((rx_expr_25)));
    double rx_expr_14 = ((((cens) * (cens)) * (cens)));
    out[2] = ((((2.000000000000000 * exp((((-(((rx_expr_0) * (rx_expr_0)))) * rx_expr_9) / _safe_zero(r)))) * rx_expr_9) / _safe_zero((((r * M_PI) * (((rx_expr_19) * (rx_expr_19))))))) + ((((rx_expr_6 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_9) / _safe_zero(r)))) * (rx_expr_0)) * rx_expr_14) / _safe_zero((rx_expr_25))));
    double rx_expr_4 = R_pow(r, 2.500000000000000);
    double rx_expr_11 = (rx_expr_4 * rx_expr_5);
    double rx_expr_26 = (rx_expr_11 * (rx_expr_19));
    out[3] = ((((((-0.500000000000000 * rx_expr_1) * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_9) / _safe_zero(r)))) * cens) / _safe_zero((rx_expr_25))) + (((exp((((-(((rx_expr_0) * (rx_expr_0)))) * rx_expr_9) / _safe_zero(r))) * (rx_expr_0)) * rx_expr_9) / _safe_zero(((((((r) * (r))) * M_PI) * (((rx_expr_19) * (rx_expr_19)))))))) + ((((rx_expr_7 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_9) / _safe_zero(r)))) * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_14) / _safe_zero((rx_expr_26))));
    out[4] = ((((((0.500000000000000 * exp((((-(((rx_expr_0) * (rx_expr_0)))) * rx_expr_9) / _safe_zero(r)))) * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_9) / _safe_zero((((((((r) * (r)) * (r))) * M_PI) * (((rx_expr_19) * (rx_expr_19))))))) + (((((0.250000000000000 * rx_expr_1) * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_9) / _safe_zero(r)))) * ((((rx_expr_0) * (rx_expr_0)) * (rx_expr_0)))) * rx_expr_14) / _safe_zero((((R_pow(r, 3.500000000000000) * rx_expr_5) * (rx_expr_19)))))) - (((((0.750000000000000 * rx_expr_1) * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_9) / _safe_zero(r)))) * (rx_expr_0)) * cens) / _safe_zero((rx_expr_26))));
    if (order >= 3) {
      double rx_expr_0 = (dv - f);
      double rx_expr_1 = 1.414213562373095;
      double rx_expr_2 = R_pow(r, 1.500000000000000);
      double rx_expr_3 = sqrt(r);
      double rx_expr_4 = R_pow(r, 2.500000000000000);
      double rx_expr_7 = R_pow(M_PI, 1.500000000000000);
      double rx_expr_8 = sqrt(M_PI);
      double rx_expr_9 = rx_expr_1;
      double rx_expr_10 = (((r) * (r)));
      double rx_expr_13 = (rx_expr_10 * M_PI);
      double rx_expr_14 = (0.500000000000000 * rx_expr_1);
      double rx_expr_16 = (((cens) * (cens)));
      double rx_expr_17 = (rx_expr_4 * rx_expr_8);
      double rx_expr_23 = (rx_expr_14 * (rx_expr_0));
      double rx_expr_25 = ((((cens) * (cens)) * (cens)));
      double rx_expr_26 = (rx_expr_23 * cens);
      double rx_expr_30 = (((((cens) * (cens)) * (cens)) * (cens)));
      double rx_expr_31 = (rx_expr_26 / _safe_zero(rx_expr_3));
      double rx_expr_32 = ((((((cens) * (cens)) * (cens)) * (cens)) * (cens)));
      double rx_expr_34 = erf(rx_expr_31);
      double rx_expr_35 = (1.000000000000000 + rx_expr_34);
      double rx_expr_44 = (rx_expr_17 * (rx_expr_35));
      out[5] = (((((((4.000000000000000 * rx_expr_1) * exp((((-1.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_16) / _safe_zero(r)))) * rx_expr_25) / _safe_zero((((rx_expr_2 * rx_expr_7) * ((((rx_expr_35) * (rx_expr_35)) * (rx_expr_35))))))) - (((rx_expr_9 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_16) / _safe_zero(r)))) * rx_expr_25) / _safe_zero((((rx_expr_2 * rx_expr_8) * (rx_expr_35)))))) + ((((6.000000000000000 * exp((((-(((rx_expr_0) * (rx_expr_0)))) * rx_expr_16) / _safe_zero(r)))) * (rx_expr_0)) * rx_expr_30) / _safe_zero(((rx_expr_13 * (((rx_expr_35) * (rx_expr_35)))))))) + ((((rx_expr_9 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_16) / _safe_zero(r)))) * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_32) / _safe_zero((rx_expr_44))));
      double rx_expr_5 = R_pow(r, 3.500000000000000);
      double rx_expr_11 = (0.500000000000000 * rx_expr_1);
      double rx_expr_12 = (1.500000000000000 * rx_expr_1);
      double rx_expr_15 = ((((r) * (r)) * (r)));
      double rx_expr_18 = (rx_expr_15 * M_PI);
      double rx_expr_19 = (rx_expr_5 * rx_expr_8);
      double rx_expr_45 = (rx_expr_19 * (rx_expr_35));
      out[6] = (((((((-2.000000000000000 * exp((((-(((rx_expr_0) * (rx_expr_0)))) * rx_expr_16) / _safe_zero(r)))) * rx_expr_16) / _safe_zero(((rx_expr_13 * (((rx_expr_35) * (rx_expr_35))))))) + ((((3.000000000000000 * exp((((-(((rx_expr_0) * (rx_expr_0)))) * rx_expr_16) / _safe_zero(r)))) * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_30) / _safe_zero(((rx_expr_18 * (((rx_expr_35) * (rx_expr_35)))))))) + ((((rx_expr_11 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_16) / _safe_zero(r)))) * ((((rx_expr_0) * (rx_expr_0)) * (rx_expr_0)))) * rx_expr_32) / _safe_zero((rx_expr_45)))) + (((((2.000000000000000 * rx_expr_1) * exp((((-1.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_16) / _safe_zero(r)))) * (rx_expr_0)) * rx_expr_25) / _safe_zero((((rx_expr_4 * rx_expr_7) * ((((rx_expr_35) * (rx_expr_35)) * (rx_expr_35)))))))) - ((((rx_expr_12 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_16) / _safe_zero(r)))) * (rx_expr_0)) * rx_expr_25) / _safe_zero((rx_expr_44))));
      double rx_expr_6 = R_pow(r, 4.500000000000000);
      double rx_expr_20 = (rx_expr_6 * rx_expr_8);
      double rx_expr_21 = (((((r) * (r)) * (r)) * (r)));
      double rx_expr_24 = (rx_expr_21 * M_PI);
      double rx_expr_46 = (rx_expr_20 * (rx_expr_35));
      out[7] = (((((((((0.750000000000000 * rx_expr_1) * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_16) / _safe_zero(r)))) * cens) / _safe_zero((rx_expr_44))) + ((((1.500000000000000 * exp((((-(((rx_expr_0) * (rx_expr_0)))) * rx_expr_16) / _safe_zero(r)))) * ((((rx_expr_0) * (rx_expr_0)) * (rx_expr_0)))) * rx_expr_30) / _safe_zero(((rx_expr_24 * (((rx_expr_35) * (rx_expr_35)))))))) - ((((2.500000000000000 * exp((((-(((rx_expr_0) * (rx_expr_0)))) * rx_expr_16) / _safe_zero(r)))) * (rx_expr_0)) * rx_expr_16) / _safe_zero(((rx_expr_18 * (((rx_expr_35) * (rx_expr_35)))))))) + (((((0.250000000000000 * rx_expr_1) * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_16) / _safe_zero(r)))) * (((((rx_expr_0) * (rx_expr_0)) * (rx_expr_0)) * (rx_expr_0)))) * rx_expr_32) / _safe_zero((rx_expr_46)))) + ((((rx_expr_9 * exp((((-1.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_16) / _safe_zero(r)))) * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_25) / _safe_zero((((rx_expr_5 * rx_expr_7) * ((((rx_expr_35) * (rx_expr_35)) * (rx_expr_35)))))))) - ((((rx_expr_12 * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_16) / _safe_zero(r)))) * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_25) / _safe_zero((rx_expr_45))));
      out[8] = (((((((((0.750000000000000 * exp((((-(((rx_expr_0) * (rx_expr_0)))) * rx_expr_16) / _safe_zero(r)))) * (((((rx_expr_0) * (rx_expr_0)) * (rx_expr_0)) * (rx_expr_0)))) * rx_expr_30) / _safe_zero((((((((((r) * (r)) * (r)) * (r)) * (r))) * M_PI) * (((rx_expr_35) * (rx_expr_35))))))) - ((((2.250000000000000 * exp((((-(((rx_expr_0) * (rx_expr_0)))) * rx_expr_16) / _safe_zero(r)))) * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_16) / _safe_zero(((rx_expr_24 * (((rx_expr_35) * (rx_expr_35)))))))) + (((((0.125000000000000 * rx_expr_1) * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_16) / _safe_zero(r)))) * ((((((rx_expr_0) * (rx_expr_0)) * (rx_expr_0)) * (rx_expr_0)) * (rx_expr_0)))) * rx_expr_32) / _safe_zero((((R_pow(r, 5.500000000000000) * rx_expr_8) * (rx_expr_35)))))) + ((((rx_expr_11 * exp((((-1.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_16) / _safe_zero(r)))) * ((((rx_expr_0) * (rx_expr_0)) * (rx_expr_0)))) * rx_expr_25) / _safe_zero((((rx_expr_6 * rx_expr_7) * ((((rx_expr_35) * (rx_expr_35)) * (rx_expr_35)))))))) - (((((1.250000000000000 * rx_expr_1) * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_16) / _safe_zero(r)))) * ((((rx_expr_0) * (rx_expr_0)) * (rx_expr_0)))) * rx_expr_25) / _safe_zero((rx_expr_46)))) + (((((1.875000000000000 * rx_expr_1) * exp((((-0.500000000000000 * (((rx_expr_0) * (rx_expr_0)))) * rx_expr_16) / _safe_zero(r)))) * (rx_expr_0)) * cens) / _safe_zero((rx_expr_45))));
    }
  }
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
