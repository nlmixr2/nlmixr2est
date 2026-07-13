// vaeEncoder.h -- arma-native VAE LSTM encoder core, shared by the
// vaeEncoderFwdBwd R export (vaeEncoder.cpp) and the C++ training loop
// (vaeTrainCpp_ in inner.cpp).
#ifndef __VAEENCODER_H__
#define __VAEENCODER_H__

#include <RcppArmadillo.h>

// Single-layer unidirectional LSTM (PyTorch i,f,g,o gate convention) + FC head
// -> (mu, logSigma, L) with reparameterization z = mu + L eps.  When
// `backward` is true, runs the exact analytic backward (LSTM BPTT + FC +
// reparam) accumulating the parameter gradients into gWih..gFcB (which must be
// pre-sized and are zeroed here); when false the gradient outputs are left
// untouched.  Outputs mu/logSigma/Lout/zOut must be pre-sized
// ([N,zDim], [N,zDim], [zDim,zDim,N], [N,zDim]).
void vaeEncoderFwdBwdCore(const arma::cube& dataIn,     // [N, Tmax, xDim]
                          const arma::ivec& lengths,    // [N]
                          const arma::mat& covIn,       // [N, nCov]
                          const arma::mat& eps,         // [N, zDim]
                          const arma::mat& Wih,         // [4h, xDim]
                          const arma::mat& Whh,         // [4h, h]
                          const arma::vec& bih,         // [4h]
                          const arma::vec& bhh,         // [4h]
                          const arma::mat& fcW,         // [outDim, h + nCov]
                          const arma::vec& fcB,         // [outDim]
                          int zDim,
                          const arma::mat& gZ,              // [N, zDim]
                          const arma::mat& gLogSigmaDirect, // [N, zDim]
                          bool backward,
                          arma::mat& mu, arma::mat& logSigma,
                          arma::cube& Lout, arma::mat& zOut,
                          arma::mat& gWih, arma::mat& gWhh,
                          arma::vec& gbih, arma::vec& gbhh,
                          arma::mat& gFcW, arma::vec& gFcB);

#endif // __VAEENCODER_H__
