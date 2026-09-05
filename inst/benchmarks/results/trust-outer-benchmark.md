# est="trust" hessianMethod= vs est="bobyqa" benchmark

Total corpus entries: 23

## hessianMethod="fd" (both converged: 21/23)

- Median speedup (bobyqa time / trust time): 0.66x
- Median |objf diff|: 1.551
- Median max |param diff|: 0.05681
- Entries with |objf diff| > 1 (likely a different local optimum, not just numeric noise): 11

## hessianMethod="bfgs" (both converged: 21/23)

- Median speedup (bobyqa time / trust time): 0.83x
- Median |objf diff|: 1.298
- Median max |param diff|: 0.1979
- Entries with |objf diff| > 1 (likely a different local optimum, not just numeric noise): 11

## hessianMethod="sr1" (both converged: 21/23)

- Median speedup (bobyqa time / trust time): 1.30x
- Median |objf diff|: 1.545
- Median max |param diff|: 0.09383
- Entries with |objf diff| > 1 (likely a different local optimum, not just numeric noise): 11

## hessianMethod="bofill" (both converged: 21/23)

- Median speedup (bobyqa time / trust time): 1.27x
- Median |objf diff|: 1.551
- Median max |param diff|: 0.09385
- Entries with |objf diff| > 1 (likely a different local optimum, not just numeric noise): 11

| source | model | bobyqa ok | bobyqa time | bobyqa objf | trust[fd] ok | trust[fd] time | trust[fd] objf | trust[fd] dObjf | trust[fd] speedup | trust[fd] maxParDiff | trust[bfgs] ok | trust[bfgs] time | trust[bfgs] objf | trust[bfgs] dObjf | trust[bfgs] speedup | trust[bfgs] maxParDiff | trust[sr1] ok | trust[sr1] time | trust[sr1] objf | trust[sr1] dObjf | trust[sr1] speedup | trust[sr1] maxParDiff | trust[bofill] ok | trust[bofill] time | trust[bofill] objf | trust[bofill] dObjf | trust[bofill] speedup | trust[bofill] maxParDiff |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| test:test-cens-dist-warn | test:test-cens-dist-warn#1 | TRUE | 0.884 | 258.1839 | TRUE | 0.591 | 258.1843 | 0.0004609432 | 1.49577 | 0.0279133 | TRUE | 0.55 | 258.2757 | 0.09182249 | 1.607273 | 0.6986012 | TRUE | 0.458 | 258.1843 | 0.000451142 | 1.930131 | 1.261659 | TRUE | 0.487 | 258.1843 | 0.0004556921 | 1.815195 | 0.4232278 |
| test:test-cens-dist-warn | test:test-cens-dist-warn#2 | TRUE | 0.989 | 22.28706 | TRUE | 1.509 | 22.25837 | -0.02868386 | 0.6554009 | 0.04975605 | TRUE | 0.631 | 22.23087 | -0.05619005 | 1.567353 | 0.04023307 | TRUE | 0.637 | 22.30889 | 0.02183063 | 1.55259 | 0.06458114 | TRUE | 0.632 | 22.22968 | -0.05737305 | 1.564873 | 0.04074211 |
| test:test-lincmt-ode-fit | test:test-lincmt-ode-fit#1 | TRUE | 1.668 | 536.3186 | TRUE | 17.84 | 497.2813 | -39.03738 | 0.09349776 | 1.036766 | TRUE | 7.294 | 497.2813 | -39.03738 | 0.2286811 | 1.036816 | TRUE | 6.224 | 497.2813 | -39.03737 | 0.2679949 | 1.036807 | TRUE | 14.321 | 501.1979 | -35.12073 | 0.1164723 | 0.8421035 |
| test:test-lincmt-ode-fit | test:test-lincmt-ode-fit#2 | TRUE | 2.26 | 536.3186 | TRUE | 18.674 | 497.2813 | -39.03738 | 0.1210239 | 1.036766 | TRUE | 7.154 | 497.2813 | -39.03738 | 0.3159072 | 1.036816 | TRUE | 6.159 | 497.2813 | -39.03737 | 0.3669427 | 1.036807 | TRUE | 14.118 | 501.1979 | -35.12073 | 0.1600793 | 0.8421035 |
| test:test-matexp | test:test-matexp#1 | TRUE | 1.385 | -140.4195 | TRUE | 1.437 | -141.9705 | -1.551089 | 0.9638135 | 0.6339321 | TRUE | 2.396 | -141.7172 | -1.297794 | 0.5780467 | 0.3750075 | TRUE | 1.063 | -141.9644 | -1.544926 | 1.302916 | 0.5929592 | TRUE | 1.088 | -141.9705 | -1.551089 | 1.272978 | 0.6338801 |
| test:test-matexp | test:test-matexp#2 | TRUE | 1.265 | -140.4195 | TRUE | 1.324 | -141.9705 | -1.551089 | 0.9554381 | 0.6339321 | TRUE | 2.087 | -141.7172 | -1.297794 | 0.6061332 | 0.3750075 | TRUE | 0.896 | -141.9644 | -1.544926 | 1.41183 | 0.5929592 | TRUE | 0.887 | -141.9705 | -1.551089 | 1.426156 | 0.6338801 |
| test:test-nlm-cens-t | test:test-nlm-cens-t#1 | TRUE | 0.618 | 42.80645 | TRUE | 1.647 | 48.2984 | 5.491957 | 0.3752277 | 0.1668793 | TRUE | 0.742 | 48.2984 | 5.491957 | 0.8328841 | 0.1668793 | TRUE | 0.78 | 48.2984 | 5.491957 | 0.7923077 | 0.1668793 | TRUE | 0.746 | 48.2984 | 5.491957 | 0.8284182 | 0.1668793 |
| test:test-nlm-cens-t | test:test-nlm-cens-t#2 | TRUE | 0.605 | 31.12565 | TRUE | 1.139 | 31.44066 | 0.3150148 | 0.5311677 | 0.09371813 | TRUE | 0.873 | 31.44179 | 0.3161355 | 0.6930126 | 0.09379709 | TRUE | 0.877 | 31.44213 | 0.3164787 | 0.6898518 | 0.0938289 | TRUE | 0.79 | 31.44278 | 0.3171301 | 0.7658228 | 0.09384942 |
| test:test-nlm-cens-t | test:test-nlm-cens-t#3 | TRUE | 0.536 | 31.19994 | TRUE | 1.412 | 38.1996 | 6.999666 | 0.3796034 | 0.2666448 | TRUE | 0.856 | 38.1996 | 6.999666 | 0.6261682 | 0.2666448 | TRUE | 0.806 | 38.1996 | 6.999666 | 0.6650124 | 0.2666448 | TRUE | 0.77 | 38.1996 | 6.999666 | 0.6961039 | 0.2666448 |
| test:test-nlm-cens-t | test:test-nlm-cens-t#4 | TRUE | 1.206 | 43.36694 | TRUE | 1.138 | 36.92497 | -6.441965 | 1.059754 | 0.05680591 | TRUE | 0.992 | 36.92497 | -6.441965 | 1.215726 | 0.05680591 | TRUE | 0.734 | 36.92497 | -6.441965 | 1.643052 | 0.05680591 | TRUE | 0.663 | 36.92497 | -6.441965 | 1.819005 | 0.05680591 |
| test:test-nlm-cens-t | test:test-nlm-cens-t#5 | TRUE | 1.085 | 43.36694 | TRUE | 1.175 | 36.92497 | -6.441965 | 0.9234043 | 0.05680591 | TRUE | 0.823 | 36.92497 | -6.441965 | 1.318348 | 0.05680591 | TRUE | 0.759 | 36.92497 | -6.441965 | 1.429513 | 0.05680591 | TRUE | 0.7 | 36.92497 | -6.441965 | 1.55 | 0.05680591 |
| test:test-nlm-cens-t | test:test-nlm-cens-t#6 | TRUE | 0.432 | 43.36694 | TRUE | 1.094 | 36.92497 | -6.441965 | 0.3948812 | 0.05680591 | TRUE | 0.707 | 36.92497 | -6.441965 | 0.6110325 | 0.05680591 | TRUE | 0.666 | 36.92497 | -6.441965 | 0.6486486 | 0.05680591 | TRUE | 0.699 | 36.92497 | -6.441965 | 0.6180258 | 0.05680591 |
| test:test-nlm-lik-contrib | test:test-nlm-lik-contrib#1 | TRUE | 1.263 | 216.1453 | TRUE | 0.982 | 216.1514 | 0.006166986 | 1.286151 | 0.00165489 | TRUE | 0.683 | 216.1515 | 0.006168388 | 1.849195 | 0.001741454 | TRUE | 0.686 | 216.1514 | 0.006167211 | 1.841108 | 0.001652988 | TRUE | 0.721 | 216.1516 | 0.006273663 | 1.751734 | 0.002404527 |
| test:test-nlm | test:test-nlm#1 | TRUE | 1.089 | -728.2062 | TRUE | 0.524 | -728.2812 | -0.07498424 | 2.078244 | 1.929441 | TRUE | 0.428 | -728.2812 | -0.07498375 | 2.544393 | 1.929296 | TRUE | 0.431 | -728.2812 | -0.07498425 | 2.526682 | 1.929975 | TRUE | 0.434 | -728.2811 | -0.0749447 | 2.509217 | 1.888426 |
| test:test-nlm | test:test-nlm#2 | TRUE | 2.606 | 496.285 | TRUE | 1.181 | 496.2426 | -0.04242127 | 2.206605 | 0.01143937 | TRUE | 0.861 | 496.2426 | -0.04239962 | 3.026713 | 0.01144973 | TRUE | 0.892 | 496.2426 | -0.04242141 | 2.921525 | 0.01146594 | TRUE | 0.844 | 496.2427 | -0.04234763 | 3.087678 | 0.01107343 |
| test:test-nlm | test:test-nlm#3 | TRUE | 0.573 | -161.754 | TRUE | 1.093 | -161.8738 | -0.1198215 | 0.5242452 | 0.05341999 | TRUE | 1.003 | -161.6431 | 0.1109331 | 0.5712861 | 0.1979262 | TRUE | 0.802 | -161.8738 | -0.1198214 | 0.7144638 | 0.05344711 | TRUE | 0.851 | -161.8738 | -0.1198202 | 0.6733255 | 0.05361187 |
| test:test-nlm | test:test-nlm#4 | TRUE | 0.587 | -161.754 | TRUE | 1.224 | -161.8738 | -0.1198215 | 0.4795752 | 0.05341999 | TRUE | 1.099 | -161.6431 | 0.1109331 | 0.5341219 | 0.1979262 | TRUE | 0.679 | -161.8738 | -0.1198214 | 0.8645066 | 0.05344711 | TRUE | 0.714 | -161.8738 | -0.1198202 | 0.8221289 | 0.05361187 |
| test:test-nlm | test:test-nlm#5 | TRUE | 1.629 | 98.4812 | TRUE | 1.635 | 98.48119 | -1.043081e-05 | 0.9963303 | 0.0003183686 | TRUE | 0.333 | 98.48119 | -9.363996e-06 | 4.891892 | 0.0003679814 | TRUE | 0.356 | 98.48119 | -9.800433e-06 | 4.575843 | 0.0003530292 | TRUE | 0.39 | 98.48119 | -1.043081e-05 | 4.176923 | 0.0003185158 |
| test:test-nlmsetup-fresh-rx | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-optexpression-ebe | test:test-optexpression-ebe#1 | TRUE | 0.556 | 342.5069 | TRUE | 1.084 | 216.1514 | -126.3555 | 0.5129151 | 0.5636927 | TRUE | 0.733 | 216.1515 | -126.3554 | 0.7585266 | 0.5637793 | TRUE | 0.788 | 216.1514 | -126.3555 | 0.7055838 | 0.5636908 | TRUE | 0.728 | 216.1516 | -126.3553 | 0.7637363 | 0.5644424 |
| test:test-optexpression-ebe | test:test-optexpression-ebe#2 | TRUE | 2.436 | 342.5069 | TRUE | 0.973 | 216.1514 | -126.3555 | 2.503597 | 0.5636927 | TRUE | 0.755 | 216.1515 | -126.3554 | 3.22649 | 0.5637793 | TRUE | 0.794 | 216.1514 | -126.3555 | 3.06801 | 0.5636908 | TRUE | 0.761 | 216.1516 | -126.3553 | 3.201051 | 0.5644424 |
| test:test-optim | test:test-optim#1 | TRUE | 0.481 | -649.7053 | TRUE | 1.718 | -649.7053 | -1.158339e-05 | 0.2799767 | 0.004068166 | TRUE | 0.543 | -649.7053 | -1.10036e-05 | 0.8858195 | 0.00501942 | TRUE | 0.54 | -649.7053 | -1.158242e-05 | 0.8907407 | 0.004102418 | TRUE | 0.495 | -649.7052 | 0.0001287166 | 0.9717172 | 0.01793185 |
| test:test-splitbolus-interp | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA |
