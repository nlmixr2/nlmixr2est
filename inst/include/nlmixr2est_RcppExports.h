// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#ifndef RCPP_nlmixr2est_RCPPEXPORTS_H_GEN_
#define RCPP_nlmixr2est_RCPPEXPORTS_H_GEN_

#include "nlmixr2est_types.h"
#include <RcppEigen.h>
#include <Rcpp.h>

namespace nlmixr2est {

    using namespace Rcpp;

    namespace {
        void validateSignature(const char* sig) {
            Rcpp::Function require = Rcpp::Environment::base_env()["require"];
            require("nlmixr2est", Rcpp::Named("quietly") = true);
            typedef int(*Ptr_validate)(const char*);
            static Ptr_validate p_validate = (Ptr_validate)
                R_GetCCallable("nlmixr2est", "_nlmixr2est_RcppExport_validate");
            if (!p_validate(sig)) {
                throw Rcpp::function_not_exported(
                    "C++ function with signature '" + std::string(sig) + "' not found in nlmixr2est");
            }
        }
    }

    inline List nlmixrExpandFdParNlme_(CharacterVector state, CharacterVector vars) {
        typedef SEXP(*Ptr_nlmixrExpandFdParNlme_)(SEXP,SEXP);
        static Ptr_nlmixrExpandFdParNlme_ p_nlmixrExpandFdParNlme_ = NULL;
        if (p_nlmixrExpandFdParNlme_ == NULL) {
            validateSignature("List(*nlmixrExpandFdParNlme_)(CharacterVector,CharacterVector)");
            p_nlmixrExpandFdParNlme_ = (Ptr_nlmixrExpandFdParNlme_)R_GetCCallable("nlmixr2est", "_nlmixr2est_nlmixrExpandFdParNlme_");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_nlmixrExpandFdParNlme_(Shield<SEXP>(Rcpp::wrap(state)), Shield<SEXP>(Rcpp::wrap(vars)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<List >(rcpp_result_gen);
    }

}

#endif // RCPP_nlmixr2est_RCPPEXPORTS_H_GEN_
