# innerOpt="trust" vs innerOpt="n1qn1" benchmark

Total corpus entries: 85 (both converged: 43, skipped/failed: 42)

Median speedup (n1qn1 time / trust time): 2.64x
Median |objf diff|: 0.01162
Median max |param diff|: 0.004068

| source | model | n1qn1 ok | n1qn1 time | n1qn1 objf | trust ok | trust time | trust objf | dObjf | speedup | maxParDiff |
|---|---|---|---|---|---|---|---|---|---|---|
| discussion-924 | d924-theo-1cmt | TRUE | 4.352 | 116.8041 | TRUE | 0.531 | 116.8335 | 0.02942797 | 8.195857 | 0.008659051 |
| discussion-924 | d924-theo-block-omega | TRUE | 3.238 | 106.6109 | TRUE | 0.565 | 108.4687 | 1.857736 | 5.730973 | 0.01501075 |
| discussion-924 | d924-oral-2cmt | TRUE | 2.874 | 113.1409 | TRUE | 0.894 | 113.5058 | 0.36492 | 3.214765 | 0.4433785 |
| discussion-924 | d924-oral-1cmt-mm | TRUE | 4.404 | 117.3518 | TRUE | 1.611 | 117.2997 | -0.05212468 | 2.733706 | 0.08702107 |
| discussion-924 | d924-pheno-sparse | TRUE | 4.031 | 725.8398 | TRUE | 1.42 | 726.5347 | 0.6949065 | 2.838732 | 0.08495316 |
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
| vignette:broom | vignette:broom#1 | TRUE | 2.965 | 1255.203 | TRUE | 2.694 | 688.9636 | -566.239 | 1.100594 | 4.52051 |
| vignette:censoring | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| vignette:citations | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| vignette:delays | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| vignette:mavoglurant | vignette:mavoglurant#1 | TRUE | 12.261 | 1957.342 | TRUE | 5.901 | 1957.4 | 0.05750313 | 2.077783 | 0.128241 |
| vignette:mavoglurant | vignette:mavoglurant#2 | TRUE | 10.56 | 1965.677 | TRUE | 4.964 | 1962.971 | -2.705454 | 2.127317 | 0.2373275 |
| vignette:modelPiping | vignette:modelPiping#1 | TRUE | 4.652 | 116.8127 | TRUE | 0.772 | 116.8234 | 0.01070936 | 6.025907 | 0.008142084 |
| vignette:modelPiping | vignette:modelPiping#2 | TRUE | 4.486 | 116.8438 | TRUE | 0.709 | 116.8434 | -0.0003425376 | 6.327221 | 2.220446e-15 |
| vignette:modelPiping | vignette:modelPiping#3 | TRUE | 4.502 | 116.8097 | TRUE | 0.63 | 116.8099 | 0.0001659066 | 7.146032 | 0.00164191 |
| vignette:modelPiping | vignette:modelPiping#4 | TRUE | 1.249 | 116.8226 | TRUE | 0.781 | 116.834 | 0.01145704 | 1.599232 | 0.003213825 |
| vignette:modelPiping | vignette:modelPiping#5 | TRUE | 4.483 | 176.5975 | TRUE | 0.695 | 176.5965 | -0.001028291 | 6.45036 | 0.001063984 |
| vignette:modelPiping | vignette:modelPiping#6 | TRUE | 4.848 | 116.832 | TRUE | 0.986 | 116.8333 | 0.001318554 | 4.916836 | 0.01918816 |
| vignette:modelPiping | vignette:modelPiping#7 | TRUE | 4.895 | 116.1704 | TRUE | 0.825 | 116.1732 | 0.002889219 | 5.933333 | 0.001094615 |
| vignette:modelPiping | vignette:modelPiping#8 | TRUE | 5.623 | 104.4766 | TRUE | 1.087 | 104.4378 | -0.03882549 | 5.172953 | 0.008404756 |
| vignette:multiple-endpoints | vignette:multiple-endpoints#1 | TRUE | 43.956 | 2267.523 | TRUE | 20.55 | 1333.794 | -933.7295 | 2.138978 | 5.153331 |
| vignette:nimo | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| vignette:precompute-articles | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| vignette:residualErrors | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| vignette:running_nlmixr | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| vignette:wbc | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| vignette:xgxr-nlmixr-ggpmx | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-1 | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-cens | test:test-focei-cens#1 | TRUE | 0.288 | 80.28091 | TRUE | 0.267 | 80.07175 | -0.2091607 | 1.078652 | 0.001001241 |
| test:test-focei-char | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-cores-fresh-rx | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-eta-only | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-eta-reset-path-dependence | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-factr-control | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-factr-fit | test:test-focei-factr-fit#1 | TRUE | 14.489 | 19592.63 | TRUE | 6.106 | 19592.75 | 0.1166556 | 2.372912 | 0.01340993 |
| test:test-focei-factr-fit | test:test-focei-factr-fit#2 | TRUE | 45.155 | 19592.53 | TRUE | 17.088 | 19592.65 | 0.1163041 | 2.642498 | 0.004281016 |
| test:test-focei-family-control | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-fast-grad | test:test-focei-fast-grad#1 | TRUE | 5.533 | 154.5595 | TRUE | 0.398 | 154.5611 | 0.001636159 | 13.90201 | 0 |
| test:test-focei-fast-grad | test:test-focei-fast-grad#2 | TRUE | 11.374 | 76.99156 | TRUE | 0.346 | 76.99308 | 0.00151846 | 32.87283 | 0 |
| test:test-focei-fast-grad | test:test-focei-fast-grad#3 | TRUE | 17.36 | 53695.56 | TRUE | 0.302 | 53695.56 | 0.0002067363 | 57.48344 | 0 |
| test:test-focei-fast-grad | test:test-focei-fast-grad#4 | TRUE | 9.761 | 13407.83 | TRUE | 1.255 | 13407.54 | -0.2886197 | 7.777689 | 0 |
| test:test-focei-fast-grad | test:test-focei-fast-grad#5 | TRUE | 10 | 161.8901 | TRUE | 0.478 | 161.8888 | -0.001248319 | 20.9205 | 0 |
| test:test-focei-fast-grad | test:test-focei-fast-grad#6 | TRUE | 2.824 | 187.288 | TRUE | 0.344 | 248.1239 | 60.83598 | 8.209302 | 0.3792874 |
| test:test-focei-fast-methods-fit | test:test-focei-fast-methods-fit#1 | TRUE | 0.638 | 116.8076 | TRUE | 0.601 | 116.8148 | 0.007201299 | 1.061564 | 0.004067578 |
| test:test-focei-fast-methods-fit | test:test-focei-fast-methods-fit#2 | TRUE | 0.599 | 116.8076 | TRUE | 0.507 | 116.8148 | 0.007201299 | 1.18146 | 0.004067578 |
| test:test-focei-fast-methods-fit | test:test-focei-fast-methods-fit#3 | TRUE | 8.26 | 116.8198 | TRUE | 0.685 | 116.8219 | 0.002109177 | 12.05839 | 0.001054489 |
| test:test-focei-fast-methods | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-fd-fallbacks | test:test-focei-fd-fallbacks#1 | TRUE | 0.477 | 119.2472 | TRUE | 0.741 | 123.6053 | 4.358045 | 0.6437247 | 0.02722319 |
| test:test-focei-fd-fallbacks | test:test-focei-fd-fallbacks#2 | TRUE | 0.523 | 118.2311 | TRUE | 0.777 | 122.8792 | 4.648035 | 0.6731017 | 0.01893422 |
| test:test-focei-foce-plus | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-inner | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-lincmt-alag-sens | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-ll-fast-grad-fit | test:test-focei-ll-fast-grad-fit#1 | TRUE | 44.333 | 53697.53 | TRUE | 0.463 | 53697.53 | 1.018506e-05 | 95.75162 | 0 |
| test:test-focei-ll-fast-grad-fit | test:test-focei-ll-fast-grad-fit#2 | TRUE | 0.344 | 53695.56 | TRUE | 0.299 | 53695.56 | 0.0002067363 | 1.150502 | 0 |
| test:test-focei-ll-fast-grad-fit | test:test-focei-ll-fast-grad-fit#3 | TRUE | 0.442 | 53697.53 | TRUE | 0.428 | 53697.53 | 1.018506e-05 | 1.03271 | 0 |
| test:test-focei-ll-fast-grad | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-llik | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-mceta-empty-cube | test:test-focei-mceta-empty-cube#1 | TRUE | 1.137 | 219.2103 | TRUE | 1.102 | 219.2103 | 0 | 1.03176 | 0 |
| test:test-focei-outer-fd-fallback | test:test-focei-outer-fd-fallback#1 | TRUE | 0.326 | 123.5806 | TRUE | 0.35 | 123.575 | -0.005595604 | 0.9314286 | 0.0007502843 |
| test:test-focei-outer-fd-fallback | test:test-focei-outer-fd-fallback#2 | TRUE | 0.523 | 123.5965 | TRUE | 0.516 | 123.5819 | -0.01462961 | 1.013566 | 0.0008523612 |
| test:test-focei-outer-fd-fallback | test:test-focei-outer-fd-fallback#3 | TRUE | 0.55 | 133.3747 | TRUE | 1.114 | 133.5397 | 0.1650066 | 0.4937163 | 0.06702624 |
| test:test-focei-outer-fd-fallback | test:test-focei-outer-fd-fallback#4 | TRUE | 7.138 | 117.5506 | TRUE | 0.448 | 117.55 | -0.0006330954 | 15.93304 | 0.000666393 |
| test:test-focei-parallel | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-pooled-solve-args | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-preprocess | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-ptrs | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-shi21-bounds | test:test-focei-shi21-bounds#1 | TRUE | 0.783 | 117.2652 | TRUE | 0.388 | 137.5873 | 20.32216 | 2.018041 | 0.1855818 |
| test:test-focei-shi21-bounds | test:test-focei-shi21-bounds#2 | TRUE | 0.82 | 134.0206 | TRUE | 0.481 | 137.5925 | 3.571901 | 1.704782 | 0.1751682 |
| test:test-focei-theta-reset-bounds | test:test-focei-theta-reset-bounds#1 | TRUE | 3.039 | 2400 | FALSE | NA | NA | NA | NA | NA |
| test:test-focei-warm | test:test-focei-warm#1 | TRUE | 0.321 | 133.5751 | TRUE | 0.334 | 133.5867 | 0.01162049 | 0.9610778 | 0 |
| test:test-focei-warm | test:test-focei-warm#2 | TRUE | 0.272 | 133.5751 | TRUE | 0.274 | 133.5867 | 0.01160659 | 0.9927007 | 0 |
| test:test-focei-zero-init-scale | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-focei-zero-theta | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA |
