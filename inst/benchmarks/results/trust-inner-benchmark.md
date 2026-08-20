# innerOpt="trust" vs innerOpt="n1qn1" benchmark

Total corpus entries: 85 (both converged: 48, skipped/failed: 37)

Median speedup (n1qn1 time / trust time): 2.59x
Median |objf diff|: 0.01116
Median max |param diff|: 0.004909

| source | model | n1qn1 ok | n1qn1 time | n1qn1 objf | trust ok | trust time | trust objf | dObjf | speedup | maxParDiff |
|---|---|---|---|---|---|---|---|---|---|---|
| discussion-924 | d924-theo-1cmt | TRUE | 4.783 | 116.8041 | TRUE | 0.626 | 116.834 | 0.02993734 | 7.640575 | 0.008659051 |
| discussion-924 | d924-theo-block-omega | TRUE | 3.228 | 106.6109 | TRUE | 0.587 | 108.4691 | 1.858167 | 5.499148 | 0.01541409 |
| discussion-924 | d924-oral-2cmt | TRUE | 3.585 | 113.1409 | TRUE | 0.918 | 113.1594 | 0.01851589 | 3.905229 | 0.02328601 |
| discussion-924 | d924-oral-1cmt-mm | TRUE | 4.66 | 117.3518 | TRUE | 1.19 | 118.5838 | 1.232056 | 3.915966 | 0.9559803 |
| discussion-924 | d924-pheno-sparse | TRUE | 4.234 | 725.8398 | TRUE | 1.375 | 726.2155 | 0.3757322 | 3.079273 | 0.03724781 |
| vignette:addingCovariances | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| vignette:babelmixr2-external-engines | vignette:babelmixr2-external-engines#1 | TRUE | 3.938 | 116.8127 | TRUE | 0.77 | 116.8234 | 0.01070937 | 5.114286 | 0.008142089 |
| vignette:imp-impmap-qrpem | vignette:imp-impmap-qrpem#1 | TRUE | 1.387 | 116.8127 | TRUE | 0.739 | 116.8234 | 0.01070937 | 1.876861 | 0.008142089 |
| vignette:linearized-quadrature-ladder | vignette:linearized-quadrature-ladder#1 | TRUE | 1.113 | 116.8127 | TRUE | 0.95 | 116.8234 | 0.01070937 | 1.171579 | 0.008142089 |
| vignette:mixture-models | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| vignette:nlm-family-optimizers | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| vignette:nlmixr2-hooks | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| vignette:nlmixr2-team-and-advisory-committee | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| vignette:nonparametric-npag-npb | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| vignette:saem | vignette:saem#1 | TRUE | 0.977 | 116.8127 | TRUE | 0.741 | 116.8234 | 0.01070937 | 1.318489 | 0.008142089 |
| vignette:simulate-titrated-dosing | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| vignette:vaeNeonatal | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| vignette:variational-inference | vignette:variational-inference#1 | TRUE | 3.815 | 130.295 | TRUE | 0.836 | 130.2899 | -0.005149068 | 4.563397 | 0.003028254 |
| vignette:broom | vignette:broom#1 | TRUE | 3.367 | 1255.203 | TRUE | 7.01 | 1175.857 | -79.34537 | 0.4803138 | 3.685358 |
| vignette:censoring | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| vignette:citations | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| vignette:delays | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| vignette:mavoglurant | vignette:mavoglurant#1 | TRUE | 12.513 | 1957.387 | TRUE | 6.101 | 1957.349 | -0.03806424 | 2.050975 | 0.01618362 |
| vignette:mavoglurant | vignette:mavoglurant#2 | TRUE | 12.441 | 1963.339 | TRUE | 4.801 | 1963.433 | 0.09341492 | 2.591335 | 0.02763026 |
| vignette:modelPiping | vignette:modelPiping#1 | TRUE | 4.6 | 116.8127 | TRUE | 1.259 | 116.8234 | 0.01070937 | 3.653693 | 0.008142089 |
| vignette:modelPiping | vignette:modelPiping#2 | TRUE | 6.1 | 116.8438 | TRUE | 1.032 | 116.8434 | -0.0003425377 | 5.910853 | 2.220446e-15 |
| vignette:modelPiping | vignette:modelPiping#3 | TRUE | 6.273 | 116.8097 | TRUE | 1.041 | 116.8099 | 0.0001659137 | 6.025937 | 0.001641903 |
| vignette:modelPiping | vignette:modelPiping#4 | TRUE | 2.003 | 116.8226 | TRUE | 1.088 | 116.8331 | 0.0105012 | 1.840993 | 0.003297857 |
| vignette:modelPiping | vignette:modelPiping#5 | TRUE | 5.96 | 176.5975 | TRUE | 1.099 | 176.5965 | -0.001066129 | 5.423112 | 0.001817615 |
| vignette:modelPiping | vignette:modelPiping#6 | TRUE | 6.889 | 116.832 | TRUE | 1.498 | 116.8448 | 0.01287024 | 4.598798 | 0.02367671 |
| vignette:modelPiping | vignette:modelPiping#7 | TRUE | 6.797 | 116.1704 | TRUE | 1.213 | 116.1733 | 0.002898836 | 5.603462 | 0.001094617 |
| vignette:modelPiping | vignette:modelPiping#8 | TRUE | 7.859 | 104.4766 | TRUE | 1.493 | 104.4378 | -0.03880686 | 5.263898 | 0.008408933 |
| vignette:multiple-endpoints | vignette:multiple-endpoints#1 | TRUE | 55.43 | 2267.142 | TRUE | 23.494 | 1332.598 | -934.5437 | 2.359326 | 5.922689 |
| vignette:nimo | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| vignette:precompute-articles | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| vignette:residualErrors | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| vignette:running_nlmixr | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| vignette:wbc | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| vignette:xgxr-nlmixr-ggpmx | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-1 | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-cens | test:test-focei-cens#1 | TRUE | 0.422 | 80.28091 | TRUE | 0.371 | 80.07175 | -0.2091607 | 1.137466 | 0.001001241 |
| test:test-focei-char | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-cores-fresh-rx | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-eta-only | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-eta-reset-path-dependence | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-factr-control | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-factr-fit | test:test-focei-factr-fit#1 | TRUE | 15.199 | 19592.63 | TRUE | 6.79 | 19592.88 | 0.2523923 | 2.238439 | 0.04314825 |
| test:test-focei-factr-fit | test:test-focei-factr-fit#2 | TRUE | 50.461 | 19592.53 | TRUE | 19.508 | 19592.44 | -0.09166179 | 2.586682 | 0.005749821 |
| test:test-focei-family-control | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-fast-grad | test:test-focei-fast-grad#1 | TRUE | 6.378 | 154.5595 | TRUE | 0.433 | 154.5611 | 0.001636159 | 14.72979 | 0 |
| test:test-focei-fast-grad | test:test-focei-fast-grad#2 | TRUE | 13.051 | 76.99156 | TRUE | 0.445 | 76.99308 | 0.00151846 | 29.32809 | 0 |
| test:test-focei-fast-grad | test:test-focei-fast-grad#3 | TRUE | 19.722 | 53695.56 | TRUE | 0.595 | 53695.56 | 0.0002067363 | 33.14622 | 0 |
| test:test-focei-fast-grad | test:test-focei-fast-grad#4 | TRUE | 10.047 | 13407.83 | TRUE | 1.535 | 13407.54 | -0.2886197 | 6.545277 | 0 |
| test:test-focei-fast-grad | test:test-focei-fast-grad#5 | TRUE | 10.907 | 161.8901 | TRUE | 0.512 | 161.8888 | -0.001248319 | 21.30273 | 0 |
| test:test-focei-fast-grad | test:test-focei-fast-grad#6 | TRUE | 3.058 | 187.288 | TRUE | 0.442 | 248.1239 | 60.83598 | 6.918552 | 0.3792874 |
| test:test-focei-fast-methods-fit | test:test-focei-fast-methods-fit#1 | TRUE | 0.677 | 116.8076 | TRUE | 0.618 | 116.8148 | 0.007201311 | 1.095469 | 0.004067583 |
| test:test-focei-fast-methods-fit | test:test-focei-fast-methods-fit#2 | TRUE | 0.701 | 116.8076 | TRUE | 0.635 | 116.8148 | 0.007201311 | 1.103937 | 0.004067583 |
| test:test-focei-fast-methods-fit | test:test-focei-fast-methods-fit#3 | TRUE | 8.167 | 116.8198 | TRUE | 0.805 | 116.8219 | 0.002109196 | 10.14534 | 0.001054478 |
| test:test-focei-fast-methods | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-fd-fallbacks | test:test-focei-fd-fallbacks#1 | TRUE | 0.704 | 119.2472 | TRUE | 1.155 | 120.3894 | 1.142182 | 0.6095238 | 0.04590297 |
| test:test-focei-fd-fallbacks | test:test-focei-fd-fallbacks#2 | TRUE | 0.627 | 118.2311 | TRUE | 1.014 | 118.3819 | 0.1507552 | 0.6183432 | 0.07884764 |
| test:test-focei-foce-plus | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-inner | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-lincmt-alag-sens | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-ll-fast-grad-fit | test:test-focei-ll-fast-grad-fit#1 | TRUE | 48.802 | 53697.53 | TRUE | 0.556 | 53697.53 | 1.018506e-05 | 87.77338 | 0 |
| test:test-focei-ll-fast-grad-fit | test:test-focei-ll-fast-grad-fit#2 | TRUE | 0.541 | 53695.56 | TRUE | 0.416 | 53695.56 | 0.0002067363 | 1.300481 | 0 |
| test:test-focei-ll-fast-grad-fit | test:test-focei-ll-fast-grad-fit#3 | TRUE | 0.631 | 53697.53 | TRUE | 0.624 | 53697.53 | 1.018506e-05 | 1.011218 | 0 |
| test:test-focei-ll-fast-grad | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-llik | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-mceta-empty-cube | test:test-focei-mceta-empty-cube#1 | TRUE | 1.3 | 219.2103 | TRUE | 1.3 | 219.2103 | 0 | 1 | 0 |
| test:test-focei-outer-fd-fallback | test:test-focei-outer-fd-fallback#1 | TRUE | 0.465 | 123.5806 | TRUE | 0.456 | 123.575 | -0.005595646 | 1.019737 | 0.0007502844 |
| test:test-focei-outer-fd-fallback | test:test-focei-outer-fd-fallback#2 | TRUE | 0.715 | 123.5965 | TRUE | 0.6 | 123.5819 | -0.01462968 | 1.191667 | 0.0008523617 |
| test:test-focei-outer-fd-fallback | test:test-focei-outer-fd-fallback#3 | TRUE | 1.042 | 133.3747 | TRUE | 1.021 | 133.4864 | 0.1117324 | 1.020568 | 0.06710815 |
| test:test-focei-outer-fd-fallback | test:test-focei-outer-fd-fallback#4 | TRUE | 7.959 | 117.5506 | TRUE | 0.592 | 117.55 | -0.0006331308 | 13.44426 | 0.0006663919 |
| test:test-focei-parallel | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-pooled-solve-args | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-preprocess | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-ptrs | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-shi21-bounds | test:test-focei-shi21-bounds#1 | TRUE | 1.152 | 117.2652 | TRUE | 0.542 | 137.5873 | 20.32216 | 2.125461 | 0.1855818 |
| test:test-focei-shi21-bounds | test:test-focei-shi21-bounds#2 | TRUE | 1.021 | 134.0206 | TRUE | 0.582 | 137.5925 | 3.571897 | 1.754296 | 0.1751682 |
| test:test-focei-theta-reset-bounds | test:test-focei-theta-reset-bounds#1 | TRUE | 3.461 | 2400 | FALSE | NA | NA | NA | NA | NA |
| test:test-focei-warm | test:test-focei-warm#1 | TRUE | 0.489 | 133.5751 | TRUE | 0.408 | 133.5867 | 0.01162049 | 1.198529 | 0 |
| test:test-focei-warm | test:test-focei-warm#2 | TRUE | 0.466 | 133.5751 | TRUE | 0.382 | 133.5867 | 0.01160659 | 1.219895 | 0 |
| test:test-focei-zero-init-scale | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-zero-theta | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
