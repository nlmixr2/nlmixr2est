# innerOpt="trust" vs innerOpt="n1qn1" benchmark

Total corpus entries: 88 (both converged: 43, skipped/failed: 45)

Median speedup (n1qn1 time / trust time): 1.14x
Median |objf diff|: 0.01162
Median max |param diff|: 0.004068

| source | model | n1qn1 ok | n1qn1 time | n1qn1 objf | trust ok | trust time | trust objf | dObjf | speedup | maxParDiff |
|---|---|---|---|---|---|---|---|---|---|---|
| discussion-924 | d924-theo-1cmt | TRUE | 6.413 | 116.8041 | TRUE | 0.721 | 116.8335 | 0.02942797 | 8.894591 | 0.008659051 |
| discussion-924 | d924-theo-block-omega | TRUE | 5.613 | 106.6109 | TRUE | 0.972 | 108.4687 | 1.857736 | 5.774691 | 0.01501075 |
| discussion-924 | d924-oral-2cmt | TRUE | 4.864 | 113.1409 | TRUE | 1.629 | 113.5058 | 0.36492 | 2.985881 | 0.4433785 |
| discussion-924 | d924-oral-1cmt-mm | TRUE | 7.191 | 117.3518 | TRUE | 2.597 | 117.2997 | -0.05212468 | 2.768964 | 0.08702107 |
| discussion-924 | d924-pheno-sparse | TRUE | 7.216 | 725.8398 | TRUE | 2.69 | 726.5347 | 0.6949065 | 2.682528 | 0.08495316 |
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
| vignette:broom | vignette:broom#1 | TRUE | 4.506 | 1255.203 | TRUE | 5.875 | 1010.257 | -244.946 | 0.7669787 | 2.373538 |
| vignette:censoring | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| vignette:citations | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| vignette:delays | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| vignette:mavoglurant | vignette:mavoglurant#1 | TRUE | 20.811 | 1957.387 | TRUE | 10.482 | 1957.35 | -0.03744922 | 1.985404 | 0.01618362 |
| vignette:mavoglurant | vignette:mavoglurant#2 | TRUE | 21.4 | 1963.361 | TRUE | 8.896 | 1964.279 | 0.9183189 | 2.405576 | 0.09426949 |
| vignette:modelPiping | vignette:modelPiping#1 | TRUE | 1.738 | 116.8127 | TRUE | 1.266 | 116.8234 | 0.01070936 | 1.372828 | 0.008142084 |
| vignette:modelPiping | vignette:modelPiping#2 | TRUE | 1.712 | 116.8457 | TRUE | 1.201 | 116.8482 | 0.002520586 | 1.425479 | 0.001805649 |
| vignette:modelPiping | vignette:modelPiping#3 | TRUE | 1.677 | 116.8111 | TRUE | 1.167 | 116.8103 | -0.000816895 | 1.437018 | 0.0004828706 |
| vignette:modelPiping | vignette:modelPiping#4 | TRUE | 1.444 | 116.8226 | TRUE | 1.283 | 116.834 | 0.01145704 | 1.125487 | 0.003213825 |
| vignette:modelPiping | vignette:modelPiping#5 | TRUE | 1.388 | 176.5784 | TRUE | 1.132 | 176.5793 | 0.0009077998 | 1.226148 | 0.001020729 |
| vignette:modelPiping | vignette:modelPiping#6 | TRUE | 1.828 | 116.8758 | TRUE | 1.417 | 116.8701 | -0.005719779 | 1.290049 | 0.02119987 |
| vignette:modelPiping | vignette:modelPiping#7 | TRUE | 1.6 | 116.1705 | TRUE | 1.635 | 116.178 | 0.007417781 | 0.9785933 | 0.00117362 |
| vignette:modelPiping | vignette:modelPiping#8 | TRUE | 2.409 | 104.4474 | TRUE | 1.45 | 104.6133 | 0.1659244 | 1.661379 | 0.02899287 |
| vignette:multiple-endpoints | vignette:multiple-endpoints#1 | TRUE | 89.248 | 2267.062 | TRUE | 37.595 | 1346.193 | -920.8697 | 2.373933 | 6.078567 |
| vignette:nimo | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| vignette:precompute-articles | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| vignette:residualErrors | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| vignette:running_nlmixr | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| vignette:wbc | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| vignette:xgxr-nlmixr-ggpmx | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-1 | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-bobyqa-retry | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-cens | test:test-focei-cens#1 | TRUE | 0.372 | 80.28091 | TRUE | 0.381 | 80.07175 | -0.2091607 | 0.976378 | 0.001001241 |
| test:test-focei-char | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-cores-fresh-rx | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-eta-only | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-eta-reset-path-dependence | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-factr-control | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-factr-fit | test:test-focei-factr-fit#1 | TRUE | 23.071 | 19592.63 | TRUE | 13.806 | 19592.75 | 0.1166556 | 1.671085 | 0.01340993 |
| test:test-focei-factr-fit | test:test-focei-factr-fit#2 | TRUE | 37.329 | 19592.53 | TRUE | 21.208 | 19592.65 | 0.1163041 | 1.760138 | 0.004281016 |
| test:test-focei-family-control | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-fast-grad | test:test-focei-fast-grad#1 | TRUE | 1.183 | 154.5595 | TRUE | 0.309 | 154.5611 | 0.001636159 | 3.828479 | 0 |
| test:test-focei-fast-grad | test:test-focei-fast-grad#2 | TRUE | 0.657 | 76.99156 | TRUE | 0.575 | 76.99308 | 0.00151846 | 1.142609 | 0 |
| test:test-focei-fast-grad | test:test-focei-fast-grad#3 | TRUE | 0.501 | 53695.56 | TRUE | 0.507 | 53695.56 | 0.0002067363 | 0.9881657 | 0 |
| test:test-focei-fast-grad | test:test-focei-fast-grad#4 | TRUE | 4.359 | 13407.83 | TRUE | 2.507 | 13407.54 | -0.2886197 | 1.738732 | 0 |
| test:test-focei-fast-grad | test:test-focei-fast-grad#5 | TRUE | 0.714 | 161.8901 | TRUE | 0.788 | 161.8888 | -0.001248319 | 0.9060914 | 0 |
| test:test-focei-fast-grad | test:test-focei-fast-grad#6 | TRUE | 0.494 | 187.288 | TRUE | 0.496 | 248.1239 | 60.83598 | 0.9959677 | 0.3792874 |
| test:test-focei-fast-methods-fit | test:test-focei-fast-methods-fit#1 | TRUE | 0.872 | 116.8076 | TRUE | 0.879 | 116.8148 | 0.007201299 | 0.9920364 | 0.004067578 |
| test:test-focei-fast-methods-fit | test:test-focei-fast-methods-fit#2 | TRUE | 0.882 | 116.8076 | TRUE | 0.882 | 116.8148 | 0.007201299 | 1 | 0.004067578 |
| test:test-focei-fast-methods-fit | test:test-focei-fast-methods-fit#3 | TRUE | 1.285 | 116.8198 | TRUE | 1.241 | 116.8219 | 0.002109177 | 1.035455 | 0.001054489 |
| test:test-focei-fast-methods | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-fd-fallbacks | test:test-focei-fd-fallbacks#1 | TRUE | 0.941 | 119.2472 | TRUE | 1.48 | 123.6053 | 4.358045 | 0.6358108 | 0.02722319 |
| test:test-focei-fd-fallbacks | test:test-focei-fd-fallbacks#2 | TRUE | 1.142 | 118.2311 | TRUE | 1.584 | 122.8792 | 4.648035 | 0.7209596 | 0.01893422 |
| test:test-focei-foce-plus | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-hessian-method | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-inner | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-lincmt-alag-sens | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-ll-fast-grad-fit | test:test-focei-ll-fast-grad-fit#1 | TRUE | 0.729 | 53697.53 | TRUE | 0.762 | 53697.53 | 1.018506e-05 | 0.9566929 | 0 |
| test:test-focei-ll-fast-grad-fit | test:test-focei-ll-fast-grad-fit#2 | TRUE | 0.509 | 53695.56 | TRUE | 0.584 | 53695.56 | 0.0002067363 | 0.8715753 | 0 |
| test:test-focei-ll-fast-grad-fit | test:test-focei-ll-fast-grad-fit#3 | TRUE | 0.676 | 53697.53 | TRUE | 0.945 | 53697.53 | 1.018506e-05 | 0.7153439 | 0 |
| test:test-focei-ll-fast-grad | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-llik | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-mceta-empty-cube | test:test-focei-mceta-empty-cube#1 | TRUE | 1.772 | 219.2103 | TRUE | 2.053 | 219.2103 | 0 | 0.8631271 | 0 |
| test:test-focei-outer-fd-fallback | test:test-focei-outer-fd-fallback#1 | TRUE | 0.748 | 123.5806 | TRUE | 0.787 | 123.575 | -0.005595604 | 0.9504447 | 0.0007502843 |
| test:test-focei-outer-fd-fallback | test:test-focei-outer-fd-fallback#2 | TRUE | 0.96 | 123.5965 | TRUE | 0.969 | 123.5819 | -0.01462961 | 0.9907121 | 0.0008523612 |
| test:test-focei-outer-fd-fallback | test:test-focei-outer-fd-fallback#3 | TRUE | 1.086 | 133.3747 | TRUE | 2.036 | 133.5397 | 0.1650066 | 0.5333988 | 0.06702624 |
| test:test-focei-outer-fd-fallback | test:test-focei-outer-fd-fallback#4 | TRUE | 0.795 | 117.5506 | TRUE | 0.941 | 117.55 | -0.0006330954 | 0.8448459 | 0.000666393 |
| test:test-focei-parallel | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-pooled-solve-args | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-preprocess | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-prior | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-ptrs | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-shi21-bounds | test:test-focei-shi21-bounds#1 | TRUE | 1.38 | 117.2652 | TRUE | 0.691 | 137.5873 | 20.32216 | 1.997106 | 0.1855818 |
| test:test-focei-shi21-bounds | test:test-focei-shi21-bounds#2 | TRUE | 1.725 | 134.0206 | TRUE | 0.94 | 137.5925 | 3.571901 | 1.835106 | 0.1751682 |
| test:test-focei-theta-reset-bounds | test:test-focei-theta-reset-bounds#1 | FALSE | NA | NA | FALSE | NA | NA | NA | NA | NA |
| test:test-focei-warm | test:test-focei-warm#1 | TRUE | 0.832 | 133.5751 | TRUE | 0.512 | 133.5867 | 0.01162049 | 1.625 | 0 |
| test:test-focei-warm | test:test-focei-warm#2 | TRUE | 0.441 | 133.5751 | TRUE | 0.429 | 133.5867 | 0.01160659 | 1.027972 | 0 |
| test:test-focei-zero-init-scale | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-zero-theta | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
