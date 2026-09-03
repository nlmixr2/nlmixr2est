# est="trust" hessianMethod= vs est="bobyqa" benchmark

Total corpus entries: 23

## hessianMethod="fd" (both converged: 21/23)

- Median speedup (bobyqa time / trust time): 0.85x
- Median |objf diff|: 1.551
- Median max |param diff|: 0.07977
- Entries with |objf diff| > 1 (likely a different local optimum, not just numeric noise): 11

## hessianMethod="bfgs" (both converged: 21/23)

- Median speedup (bobyqa time / trust time): 1.04x
- Median |objf diff|: 0.4343
- Median max |param diff|: 0.08164
- Entries with |objf diff| > 1 (likely a different local optimum, not just numeric noise): 9

## hessianMethod="sr1" (both converged: 21/23)

- Median speedup (bobyqa time / trust time): 1.50x
- Median |objf diff|: 1.532
- Median max |param diff|: 0.0706
- Entries with |objf diff| > 1 (likely a different local optimum, not just numeric noise): 11

## hessianMethod="bofill" (both converged: 21/23)

- Median speedup (bobyqa time / trust time): 1.54x
- Median |objf diff|: 1.551
- Median max |param diff|: 0.09381
- Entries with |objf diff| > 1 (likely a different local optimum, not just numeric noise): 11

| source | model | bobyqa ok | bobyqa time | bobyqa objf | trust[fd] ok | trust[fd] time | trust[fd] objf | trust[fd] dObjf | trust[fd] speedup | trust[fd] maxParDiff | trust[bfgs] ok | trust[bfgs] time | trust[bfgs] objf | trust[bfgs] dObjf | trust[bfgs] speedup | trust[bfgs] maxParDiff | trust[sr1] ok | trust[sr1] time | trust[sr1] objf | trust[sr1] dObjf | trust[sr1] speedup | trust[sr1] maxParDiff | trust[bofill] ok | trust[bofill] time | trust[bofill] objf | trust[bofill] dObjf | trust[bofill] speedup | trust[bofill] maxParDiff |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| test:test-cens-dist-warn | test:test-cens-dist-warn#1 | TRUE | 0.869 | 258.1839 | TRUE | 0.504 | 258.1847 | 0.0008755444 | 1.724206 | 0.9544515 | TRUE | 0.436 | 258.7913 | 0.6074574 | 1.993119 | 0.7369447 | TRUE | 0.415 | 258.1846 | 0.0007097233 | 2.093976 | 1.606828 | TRUE | 0.436 | 258.1845 | 0.0006446797 | 1.993119 | 0.4204929 |
| test:test-cens-dist-warn | test:test-cens-dist-warn#2 | TRUE | 0.921 | 22.28706 | TRUE | 1.078 | 22.2633 | -0.02375948 | 0.8543599 | 0.05169591 | TRUE | 0.603 | 22.23087 | -0.05619005 | 1.527363 | 0.04023307 | TRUE | 0.597 | 22.30889 | 0.02183063 | 1.542714 | 0.06458114 | TRUE | 0.58 | 22.22968 | -0.05737305 | 1.587931 | 0.04074211 |
| test:test-lincmt-ode-fit | test:test-lincmt-ode-fit#1 | TRUE | 1.619 | 536.3186 | TRUE | 14.869 | 497.2813 | -39.03738 | 0.1088843 | 1.036766 | TRUE | 5.366 | 497.2813 | -39.03734 | 0.3017145 | 1.037014 | TRUE | 5.413 | 497.2813 | -39.03737 | 0.2990948 | 1.036816 | TRUE | 6.114 | 519.6119 | -16.70678 | 0.2648021 | 2.542535 |
| test:test-lincmt-ode-fit | test:test-lincmt-ode-fit#2 | TRUE | 1.594 | 536.3186 | TRUE | 14.189 | 497.2813 | -39.03738 | 0.1123405 | 1.036766 | TRUE | 5.543 | 497.2813 | -39.03734 | 0.2875699 | 1.037014 | TRUE | 5.51 | 497.2813 | -39.03737 | 0.2892922 | 1.036816 | TRUE | 6.288 | 519.6119 | -16.70678 | 0.2534987 | 2.542535 |
| test:test-matexp | test:test-matexp#1 | TRUE | 0.975 | -140.4195 | TRUE | 1.01 | -141.9705 | -1.551088 | 0.9653465 | 0.6334603 | TRUE | 1.052 | -140.8181 | -0.398622 | 0.9268061 | 0.08164248 | TRUE | 0.625 | -141.9517 | -1.532214 | 1.56 | 0.5626262 | TRUE | 0.613 | -141.9705 | -1.551081 | 1.590538 | 0.6334126 |
| test:test-matexp | test:test-matexp#2 | TRUE | 1.034 | -140.4195 | TRUE | 0.918 | -141.9705 | -1.551088 | 1.126362 | 0.6334603 | TRUE | 0.999 | -140.8181 | -0.398622 | 1.035035 | 0.08164248 | TRUE | 0.591 | -141.9517 | -1.532214 | 1.749577 | 0.5626262 | TRUE | 0.6 | -141.9705 | -1.551081 | 1.723333 | 0.6334126 |
| test:test-nlm-cens-t | test:test-nlm-cens-t#1 | TRUE | 0.433 | 42.80645 | TRUE | 1.083 | 48.2984 | 5.491957 | 0.3998153 | 0.1668793 | TRUE | 0.679 | 48.2984 | 5.491957 | 0.6377025 | 0.1668793 | TRUE | 0.655 | 48.2984 | 5.491957 | 0.6610687 | 0.1668793 | TRUE | 0.658 | 48.2984 | 5.491957 | 0.6580547 | 0.1668793 |
| test:test-nlm-cens-t | test:test-nlm-cens-t#2 | TRUE | 0.468 | 31.12565 | TRUE | 0.773 | 31.44066 | 0.3150148 | 0.6054334 | 0.09371813 | TRUE | 0.625 | 31.44861 | 0.3229593 | 0.7488 | 0.09391893 | TRUE | 0.676 | 31.44631 | 0.320662 | 0.6923077 | 0.09379141 | TRUE | 0.648 | 31.44742 | 0.3217702 | 0.7222222 | 0.09381115 |
| test:test-nlm-cens-t | test:test-nlm-cens-t#3 | TRUE | 0.454 | 31.19994 | TRUE | 1.161 | 38.1996 | 6.999666 | 0.3910422 | 0.2666448 | TRUE | 0.665 | 38.1996 | 6.999666 | 0.6827068 | 0.2666448 | TRUE | 0.658 | 38.1996 | 6.999666 | 0.6899696 | 0.2666448 | TRUE | 0.67 | 38.1996 | 6.999666 | 0.6776119 | 0.2666448 |
| test:test-nlm-cens-t | test:test-nlm-cens-t#4 | TRUE | 0.938 | 43.36694 | TRUE | 0.865 | 36.92497 | -6.441965 | 1.084393 | 0.05680591 | TRUE | 0.839 | 36.92497 | -6.441965 | 1.117998 | 0.05680591 | TRUE | 0.582 | 36.92497 | -6.441965 | 1.611684 | 0.05680591 | TRUE | 0.611 | 36.92497 | -6.441965 | 1.535188 | 0.05680591 |
| test:test-nlm-cens-t | test:test-nlm-cens-t#5 | TRUE | 0.949 | 43.36694 | TRUE | 0.936 | 36.92497 | -6.441965 | 1.013889 | 0.05680591 | TRUE | 0.66 | 36.92497 | -6.441965 | 1.437879 | 0.05680591 | TRUE | 0.632 | 36.92497 | -6.441965 | 1.501582 | 0.05680591 | TRUE | 0.618 | 36.92497 | -6.441965 | 1.535599 | 0.05680591 |
| test:test-nlm-cens-t | test:test-nlm-cens-t#6 | TRUE | 0.365 | 43.36694 | TRUE | 0.887 | 36.92497 | -6.441965 | 0.4114994 | 0.05680591 | TRUE | 0.603 | 36.92497 | -6.441965 | 0.6053068 | 0.05680591 | TRUE | 0.612 | 36.92497 | -6.441965 | 0.5964052 | 0.05680591 | TRUE | 0.596 | 36.92497 | -6.441965 | 0.6124161 | 0.05680591 |
| test:test-nlm-lik-contrib | test:test-nlm-lik-contrib#1 | TRUE | 1.129 | 216.1453 | TRUE | 0.773 | 216.1514 | 0.006166986 | 1.460543 | 0.00165489 | TRUE | 0.57 | 216.1517 | 0.006433385 | 1.980702 | 0.001932282 | TRUE | 0.559 | 216.1516 | 0.006330244 | 2.019678 | 0.001807752 | TRUE | 0.573 | 216.1518 | 0.006484259 | 1.970332 | 0.00194956 |
| test:test-nlm | test:test-nlm#1 | TRUE | 0.902 | -688.9902 | TRUE | 0.448 | -688.9902 | -4.851295e-05 | 2.013393 | 0.002507636 | TRUE | 0.404 | -688.9901 | 0.0001128133 | 2.232673 | 0.003118019 | TRUE | 0.37 | -688.9902 | 3.129912e-05 | 2.437838 | 0.004602835 | TRUE | 0.375 | -688.9899 | 0.000325373 | 2.405333 | 0.005000613 |
| test:test-nlm | test:test-nlm#2 | TRUE | 2.078 | 496.285 | TRUE | 0.923 | 496.2426 | -0.04242127 | 2.251354 | 0.01143937 | TRUE | 0.736 | 496.2429 | -0.04209633 | 2.82337 | 0.01247407 | TRUE | 0.728 | 496.2427 | -0.0423753 | 2.854396 | 0.01146277 | TRUE | 0.714 | 496.2449 | -0.04011282 | 2.910364 | 0.01399659 |
| test:test-nlm | test:test-nlm#3 | TRUE | 0.473 | -161.754 | TRUE | 0.885 | -161.8738 | -0.1198215 | 0.5344633 | 0.05341999 | TRUE | 0.666 | -161.3197 | 0.4342762 | 0.7102102 | 0.2782779 | TRUE | 0.671 | -161.8738 | -0.1198178 | 0.704918 | 0.05299995 | TRUE | 0.681 | -161.8738 | -0.1198048 | 0.6945668 | 0.05460171 |
| test:test-nlm | test:test-nlm#4 | TRUE | 0.444 | -161.754 | TRUE | 0.9 | -161.8738 | -0.1198215 | 0.4933333 | 0.05341999 | TRUE | 0.588 | -161.3197 | 0.4342762 | 0.755102 | 0.2782779 | TRUE | 0.571 | -161.8738 | -0.1198178 | 0.7775832 | 0.05299995 | TRUE | 0.619 | -161.8738 | -0.1198048 | 0.7172859 | 0.05460171 |
| test:test-nlm | test:test-nlm#5 | TRUE | 1.347 | 98.4812 | TRUE | 1.417 | 98.48119 | -1.042989e-05 | 0.9505999 | 0.00031254 | TRUE | 0.321 | 98.48119 | -9.363996e-06 | 4.196262 | 0.0003679814 | TRUE | 0.3 | 98.48119 | -9.800433e-06 | 4.49 | 0.0003530292 | TRUE | 0.317 | 98.48299 | 0.001785315 | 4.249211 | 0.007202662 |
| test:test-nlmsetup-fresh-rx | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-optexpression-ebe | test:test-optexpression-ebe#1 | TRUE | 0.484 | 342.5069 | TRUE | 0.811 | 216.1514 | -126.3555 | 0.5967941 | 0.5636927 | TRUE | 0.624 | 216.1517 | -126.3552 | 0.775641 | 0.5638879 | TRUE | 0.62 | 216.1516 | -126.3553 | 0.7806452 | 0.5638456 | TRUE | 0.641 | 216.1518 | -126.3551 | 0.7550702 | 0.5639268 |
| test:test-optexpression-ebe | test:test-optexpression-ebe#2 | TRUE | 2.108 | 342.5069 | TRUE | 0.823 | 216.1514 | -126.3555 | 2.561361 | 0.5636927 | TRUE | 0.646 | 216.1517 | -126.3552 | 3.263158 | 0.5638879 | TRUE | 0.654 | 216.1516 | -126.3553 | 3.223242 | 0.5638456 | TRUE | 0.621 | 216.1518 | -126.3551 | 3.394525 | 0.5639268 |
| test:test-optim | test:test-optim#1 | TRUE | 0.43 | -696.5537 | TRUE | 1.579 | -696.5539 | -0.0002146411 | 0.2723243 | 0.07976865 | TRUE | 0.402 | -696.5535 | 0.0002186757 | 1.069652 | 0.02795079 | TRUE | 0.407 | -696.5539 | -0.0002122423 | 1.056511 | 0.07060465 | TRUE | 0.409 | -696.5475 | 0.006205034 | 1.051345 | 0.1594278 |
| test:test-splitbolus-interp | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA |
