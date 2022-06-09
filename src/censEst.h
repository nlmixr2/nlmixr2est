#ifndef __CENSEST_H__
#define __CENSEST_H__

// Cumulative distribution function of a standardized normal distribution (mean
// of zero, standard deviation of 1); x = (value - mean)/standard deviation
#define PHI(x) 0.5*(1.0+erf((x)/M_SQRT2))

static inline double doCensNormal1(double cens, double limDv, double lim, double ll,
                                   double f, double r, int adjLik){
  // Subtract log(2*pi)/2 since this is the likelihood instead of the liklihood missing 0.5*log(2*pi)
  double adj = 0.918938533204672669541*((double)adjLik);
  if (cens == 0.0) {
    // M2 adds likelihood even when the observation is defined
    if (R_FINITE(lim) && !ISNA(lim)) {
      return ll - log(1.0 - PHI(((lim<f)*2.0 - 1.0)*(lim - f)/sqrt(r))) - adj;
    }
  } else if (cens == 1.0 || cens == -1.0) {
    if (R_FINITE(lim) && !ISNA(lim)) {
      // M4 method
      double sd = sqrt(r);
      double cum1 = PHI(cens*(limDv-fj)/sd);
      double cum2 = PHI(cens*(lim - f)/sd);
      return log(cum1-cum2) - log(1.0 - cum2) - adj;
    } else {
      // M3 method
      return log(PHI(cens*(limDv-fc)/sqrt(r))) - adj;
    }
  }
  return NA_REAL;
}



#endif // __CENSEST_H__
