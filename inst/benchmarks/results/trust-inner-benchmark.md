# innerOpt="trust" vs innerOpt="n1qn1" benchmark

Total corpus entries: 112 (both converged: 46, skipped/failed: 66)

Median speedup (n1qn1 time / trust time): 1.07x
Median |objf diff|: 0.004034
Median max |param diff|: 0.003452

| source | model | n1qn1 ok | n1qn1 time | n1qn1 objf | trust ok | trust time | trust objf | dObjf | speedup | maxParDiff |
|---|---|---|---|---|---|---|---|---|---|---|
| discussion-924 | d924-theo-1cmt | TRUE | 4.74 | 116.8041 | TRUE | 0.735 | 116.8053 | 0.001127051 | 6.44898 | 0.002429226 |
| discussion-924 | d924-theo-block-omega | TRUE | 2.452 | 106.6603 | TRUE | 0.935 | 107.472 | 0.811765 | 2.62246 | 0.003717699 |
| discussion-924 | d924-oral-2cmt | TRUE | 2.677 | 113.1409 | TRUE | 1.115 | 113.1413 | 0.0004665952 | 2.400897 | 0.01357849 |
| discussion-924 | d924-oral-1cmt-mm | TRUE | 3.547 | 117.3518 | TRUE | 2.004 | 117.9007 | 0.5489103 | 1.76996 | 0.4411873 |
| discussion-924 | d924-pheno-sparse | TRUE | 3.901 | 725.8398 | TRUE | 1.823 | 725.6129 | -0.2268526 | 2.139879 | 0.1318147 |
| vignette:addingCovariances | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| vignette:babelmixr2-external-engines | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| vignette:imp-impmap-qrpem | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| vignette:linearized-quadrature-ladder | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| vignette:mixture-models | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| vignette:nlm-family-optimizers | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| vignette:nlmixr2-hooks | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| vignette:nlmixr2-team-and-advisory-committee | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| vignette:nonparametric-npag-npb | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| vignette:saem | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| vignette:simulate-titrated-dosing | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| vignette:vaeNeonatal | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| vignette:variational-inference | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| vignette:broom | vignette:broom#1 | TRUE | 2.724 | 1255.203 | TRUE | 4.02 | 688.742 | -566.4606 | 0.6776119 | 4.517704 |
| vignette:censoring | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| vignette:citations | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| vignette:delays | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| vignette:mavoglurant | vignette:mavoglurant#1 | TRUE | 11.008 | 1957.387 | TRUE | 8.401 | 1957.368 | -0.01877788 | 1.31032 | 0.08942586 |
| vignette:mavoglurant | vignette:mavoglurant#2 | TRUE | 10.666 | 1963.991 | TRUE | 7.381 | 1965.786 | 1.794811 | 1.445062 | 0.07856663 |
| vignette:modelPiping | vignette:modelPiping#1 | TRUE | 1.157 | 116.8127 | TRUE | 1.125 | 116.8175 | 0.004836594 | 1.028444 | 0.007581369 |
| vignette:modelPiping | vignette:modelPiping#2 | TRUE | 1.309 | 116.8432 | TRUE | 1.119 | 116.8432 | 2.283664e-05 | 1.169794 | 0.0004403274 |
| vignette:modelPiping | vignette:modelPiping#3 | TRUE | 1.68 | 116.8107 | TRUE | 1.071 | 116.8099 | -0.0008603793 | 1.568627 | 0.0004673686 |
| vignette:modelPiping | vignette:modelPiping#4 | TRUE | 1.022 | 116.8226 | TRUE | 1.034 | 116.8349 | 0.01235811 | 0.9883946 | 0.008581612 |
| vignette:modelPiping | vignette:modelPiping#5 | TRUE | 0.967 | 176.5783 | TRUE | 0.827 | 176.5989 | 0.02060642 | 1.169287 | 0.003799629 |
| vignette:modelPiping | vignette:modelPiping#6 | TRUE | 1.247 | 116.8738 | TRUE | 1.07 | 116.8484 | -0.02534256 | 1.165421 | 0.007751424 |
| vignette:modelPiping | vignette:modelPiping#7 | TRUE | 1.254 | 116.1757 | TRUE | 1.125 | 116.1789 | 0.003210156 | 1.114667 | 0.00167244 |
| vignette:modelPiping | vignette:modelPiping#8 | TRUE | 1.539 | 104.5666 | TRUE | 1.39 | 104.3809 | -0.1856818 | 1.107194 | 0.013516 |
| vignette:multiple-endpoints | vignette:multiple-endpoints#1 | TRUE | 43.351 | 2430.013 | TRUE | 34.297 | 1331.012 | -1099.002 | 1.263988 | 6.874621 |
| vignette:nimo | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| vignette:precompute-articles | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| vignette:residualErrors | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| vignette:running_nlmixr | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| vignette:wbc | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| vignette:xgxr-nlmixr-ggpmx | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-1 | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-bobyqa-retry | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-cens-t-fit | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-cens-t | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-cens | test:test-focei-cens#1 | TRUE | 0.364 | 80.28091 | TRUE | 0.39 | 80.07174 | -0.2091661 | 0.9333333 | 0.0009734435 |
| test:test-focei-char | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-cores-fresh-rx | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-eta-only | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-eta-reset-path-dependence | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-factr-control | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-factr-fit | test:test-focei-factr-fit#1 | TRUE | 10.662 | 19592.63 | TRUE | 8.579 | 19592.56 | -0.07107747 | 1.242802 | 0.008905675 |
| test:test-focei-factr-fit | test:test-focei-factr-fit#2 | TRUE | 32.754 | 19592.53 | TRUE | 24.125 | 19592.54 | 0.003230679 | 1.357679 | 0.006559009 |
| test:test-focei-family-control | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-fast-grad | test:test-focei-fast-grad#1 | TRUE | 0.752 | 154.5595 | TRUE | 0.703 | 154.5597 | 0.0001911382 | 1.069701 | 0 |
| test:test-focei-fast-grad | test:test-focei-fast-grad#2 | TRUE | 1.187 | 76.99156 | TRUE | 1.171 | 76.99307 | 0.001507704 | 1.013664 | 0 |
| test:test-focei-fast-grad | test:test-focei-fast-grad#3 | TRUE | 1.33 | 53695.56 | TRUE | 1.305 | 53695.56 | 0.0001147799 | 1.019157 | 0 |
| test:test-focei-fast-grad | test:test-focei-fast-grad#4 | TRUE | 2.517 | 13407.83 | TRUE | 1.974 | 13407.54 | -0.2884417 | 1.275076 | 0 |
| test:test-focei-fast-grad | test:test-focei-fast-grad#5 | TRUE | 1.012 | 161.8901 | TRUE | 0.896 | 161.8898 | -0.0002326117 | 1.129464 | 0 |
| test:test-focei-fast-grad | test:test-focei-fast-grad#6 | TRUE | 0.412 | 187.288 | TRUE | 0.389 | 188.3979 | 1.109921 | 1.059126 | 0.05485484 |
| test:test-focei-fast-methods-fit | test:test-focei-fast-methods-fit#1 | TRUE | 0.82 | 116.8076 | TRUE | 0.974 | 116.8079 | 0.000357685 | 0.8418891 | 0.0005122744 |
| test:test-focei-fast-methods-fit | test:test-focei-fast-methods-fit#2 | TRUE | 0.825 | 116.8076 | TRUE | 0.944 | 116.8079 | 0.000357685 | 0.8739407 | 0.0005122744 |
| test:test-focei-fast-methods-fit | test:test-focei-fast-methods-fit#3 | TRUE | 1.002 | 116.8198 | TRUE | 0.88 | 116.8223 | 0.002458825 | 1.138636 | 0.003186513 |
| test:test-focei-fast-methods | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-fd-fallbacks | test:test-focei-fd-fallbacks#1 | TRUE | 0.721 | 119.2472 | TRUE | 0.706 | 124.5619 | 5.314643 | 1.021246 | 0.06344021 |
| test:test-focei-fd-fallbacks | test:test-focei-fd-fallbacks#2 | TRUE | 0.678 | 118.2311 | TRUE | 0.742 | 122.268 | 4.036918 | 0.9137466 | 0.06206696 |
| test:test-focei-foce-plus | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-hessian-etastep | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-hessian-method | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-inner | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-lincmt-adr-default | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-lincmt-alag-sens | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-lincmt-carry-eligible-gates | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-lincmt-carry-eligible | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-lincmt-carry-emit-fd | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-lincmt-carry-emit | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-lincmt-carry-fit-fallback | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-lincmt-carry-fit | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-lincmt-carry-jump-fit | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-lincmt-carry-jump | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-lincmt-carry-ll-fit | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-lincmt-carry-ll | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-lincmt-carry-theta-fit | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-lincmt-carry-theta | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-lincmt-carry-trans | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-lincmt-phi-engage | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-ll-fast-grad-fit | test:test-focei-ll-fast-grad-fit#1 | TRUE | 3.778 | 53697.53 | TRUE | 3.658 | 53697.53 | 1.018454e-05 | 1.032805 | 0 |
| test:test-focei-ll-fast-grad-fit | test:test-focei-ll-fast-grad-fit#2 | TRUE | 1.303 | 53695.56 | TRUE | 1.436 | 53695.56 | 0.0001147799 | 0.9073816 | 0 |
| test:test-focei-ll-fast-grad-fit | test:test-focei-ll-fast-grad-fit#3 | TRUE | 3.957 | 53697.53 | TRUE | 3.793 | 53697.53 | 1.018454e-05 | 1.043238 | 0 |
| test:test-focei-ll-fast-grad | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-llik | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-mceta-empty-cube | test:test-focei-mceta-empty-cube#1 | TRUE | 1.136 | 219.2103 | TRUE | 1.14 | 219.2103 | 0 | 0.9964912 | 0 |
| test:test-focei-outer-fd-fallback | test:test-focei-outer-fd-fallback#1 | TRUE | 0.68 | 123.5806 | TRUE | 0.703 | 123.5924 | 0.01181243 | 0.9672831 | 0.00425447 |
| test:test-focei-outer-fd-fallback | test:test-focei-outer-fd-fallback#2 | TRUE | 0.935 | 123.5965 | TRUE | 0.901 | 123.5977 | 0.001199802 | 1.037736 | 0.004467526 |
| test:test-focei-outer-fd-fallback | test:test-focei-outer-fd-fallback#3 | TRUE | 0.948 | 133.3747 | TRUE | 0.973 | 133.4471 | 0.07238877 | 0.9743063 | 0.06597879 |
| test:test-focei-outer-fd-fallback | test:test-focei-outer-fd-fallback#4 | TRUE | 0.793 | 117.5506 | TRUE | 0.777 | 117.5382 | -0.0123998 | 1.020592 | 0.0007914417 |
| test:test-focei-parallel | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-pooled-solve-args | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-preprocess | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-prior | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-ptrs | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-reproducible | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-shi21-bounds | test:test-focei-shi21-bounds#1 | TRUE | 0.841 | 117.2652 | TRUE | 0.628 | 135.0752 | 17.81003 | 1.339172 | 0.205212 |
| test:test-focei-shi21-bounds | test:test-focei-shi21-bounds#2 | TRUE | 0.944 | 134.0206 | TRUE | 0.65 | 137.5915 | 3.570931 | 1.452308 | 0.1751681 |
| test:test-focei-theta-reset-bounds | test:test-focei-theta-reset-bounds#1 | FALSE | NA | NA | FALSE | NA | NA | NA | NA | NA |
| test:test-focei-theta-reset | test:test-focei-theta-reset#1 | FALSE | NA | NA | FALSE | NA | NA | NA | NA | NA |
| test:test-focei-theta-reset | test:test-focei-theta-reset#2 | FALSE | NA | NA | FALSE | NA | NA | NA | NA | NA |
| test:test-focei-vid-alloc-1010 | test:test-focei-vid-alloc-1010#1 | TRUE | 0.339 | 92688.51 | TRUE | 0.318 | 92688.51 | 3.064508e-05 | 1.066038 | 0 |
| test:test-focei-vid-alloc-1010 | test:test-focei-vid-alloc-1010#2 | TRUE | 0.227 | -75.55361 | TRUE | 0.228 | -75.55361 | 0 | 0.995614 | 0 |
| test:test-focei-vid-alloc-1010 | test:test-focei-vid-alloc-1010#3 | TRUE | 0.233 | 17816.12 | TRUE | 0.196 | 17816.12 | 0 | 1.188776 | 0 |
| test:test-focei-warm | test:test-focei-warm#1 | TRUE | 0.423 | 133.5751 | TRUE | 0.448 | 133.5765 | 0.001388782 | 0.9441964 | 0 |
| test:test-focei-warm | test:test-focei-warm#2 | TRUE | 0.404 | 133.5751 | TRUE | 0.401 | 133.5765 | 0.001374886 | 1.007481 | 0 |
| test:test-focei-zero-init-scale | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-zero-theta | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
