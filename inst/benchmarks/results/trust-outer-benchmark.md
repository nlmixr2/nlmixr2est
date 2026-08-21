# est="trust" hessianMethod= vs est="bobyqa" benchmark

Total corpus entries: 23

## hessianMethod="fd" (both converged: 21/23)

- Median speedup (bobyqa time / trust time): 0.60x
- Median |objf diff|: 1.551
- Median max |param diff|: 0.07084
- Entries with |objf diff| > 1 (likely a different local optimum, not just numeric noise): 11

## hessianMethod="bfgs" (both converged: 21/23)

- Median speedup (bobyqa time / trust time): 0.70x
- Median |objf diff|: 0.4343
- Median max |param diff|: 0.08697
- Entries with |objf diff| > 1 (likely a different local optimum, not just numeric noise): 9

## hessianMethod="sr1" (both converged: 21/23)

- Median speedup (bobyqa time / trust time): 0.80x
- Median |objf diff|: 1.532
- Median max |param diff|: 0.07055
- Entries with |objf diff| > 1 (likely a different local optimum, not just numeric noise): 11

## hessianMethod="bofill" (both converged: 21/23)

- Median speedup (bobyqa time / trust time): 0.77x
- Median |objf diff|: 1.551
- Median max |param diff|: 0.07621
- Entries with |objf diff| > 1 (likely a different local optimum, not just numeric noise): 11

| source | model | bobyqa ok | bobyqa time | bobyqa objf | trust[fd] ok | trust[fd] time | trust[fd] objf | trust[fd] dObjf | trust[fd] speedup | trust[fd] maxParDiff | trust[bfgs] ok | trust[bfgs] time | trust[bfgs] objf | trust[bfgs] dObjf | trust[bfgs] speedup | trust[bfgs] maxParDiff | trust[sr1] ok | trust[sr1] time | trust[sr1] objf | trust[sr1] dObjf | trust[sr1] speedup | trust[sr1] maxParDiff | trust[bofill] ok | trust[bofill] time | trust[bofill] objf | trust[bofill] dObjf | trust[bofill] speedup | trust[bofill] maxParDiff |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| test:test-cens-dist-warn | test:test-cens-dist-warn#1 | TRUE | 1.066 | 258.1839 | TRUE | 1.299 | 258.1847 | 0.0008755444 | 0.8206313 | 0.9544515 | TRUE | 1.221 | 258.7913 | 0.6074574 | 0.8730549 | 0.7369447 | TRUE | 1.23 | 258.1846 | 0.0007097233 | 0.8666667 | 1.606828 | TRUE | 1.233 | 258.1845 | 0.0006446797 | 0.864558 | 0.4204929 |
| test:test-cens-dist-warn | test:test-cens-dist-warn#2 | TRUE | 1.356 | 22.28706 | TRUE | 2.089 | 22.2633 | -0.02375948 | 0.6491144 | 0.05169591 | TRUE | 1.639 | 22.23087 | -0.05619005 | 0.8273337 | 0.04023307 | TRUE | 1.687 | 22.30889 | 0.02183063 | 0.8037937 | 0.06458114 | TRUE | 1.753 | 22.22968 | -0.05737305 | 0.7735311 | 0.04074211 |
| test:test-lincmt-ode-fit | test:test-lincmt-ode-fit#1 | TRUE | 3.131 | 536.3186 | TRUE | 26.939 | 179350.7 | 178814.4 | 0.1162255 | 1.845656 | TRUE | 16.841 | 179350.7 | 178814.4 | 0.1859153 | 1.845656 | TRUE | 16.724 | 179350.7 | 178814.4 | 0.187216 | 1.845656 | TRUE | 17.717 | 179350.7 | 178814.4 | 0.1767229 | 1.845656 |
| test:test-lincmt-ode-fit | test:test-lincmt-ode-fit#2 | TRUE | 3.149 | 536.3186 | TRUE | 23.544 | 497.2813 | -39.03738 | 0.1337496 | 1.036766 | TRUE | 13.985 | 497.2813 | -39.03734 | 0.2251698 | 1.037014 | TRUE | 14.424 | 497.2813 | -39.03737 | 0.2183167 | 1.036816 | TRUE | 15.575 | 519.6119 | -16.70678 | 0.202183 | 2.542535 |
| test:test-matexp | test:test-matexp#1 | TRUE | 1.497 | -140.4195 | TRUE | 2.483 | -141.9705 | -1.551088 | 0.6028997 | 0.6334603 | TRUE | 2.498 | -140.8181 | -0.398622 | 0.5992794 | 0.08164248 | TRUE | 2.157 | -141.9517 | -1.532214 | 0.6940195 | 0.5626262 | TRUE | 2.185 | -141.9705 | -1.551081 | 0.6851259 | 0.6334126 |
| test:test-matexp | test:test-matexp#2 | TRUE | 1.644 | -140.4195 | TRUE | 2.206 | -141.9705 | -1.551088 | 0.7452403 | 0.6334603 | TRUE | 2.35 | -140.8181 | -0.398622 | 0.6995745 | 0.08164248 | TRUE | 1.877 | -141.9517 | -1.532214 | 0.8758657 | 0.5626262 | TRUE | 1.838 | -141.9705 | -1.551081 | 0.8944505 | 0.6334126 |
| test:test-nlm-cens-t | test:test-nlm-cens-t#1 | TRUE | 0.963 | 42.80645 | TRUE | 2.079 | 48.2984 | 5.491957 | 0.4632035 | 0.1668793 | TRUE | 1.761 | 48.2984 | 5.491957 | 0.5468484 | 0.1668793 | TRUE | 1.779 | 48.2984 | 5.491957 | 0.5413153 | 0.1668793 | TRUE | 1.781 | 48.2984 | 5.491957 | 0.5407075 | 0.1668793 |
| test:test-nlm-cens-t | test:test-nlm-cens-t#2 | TRUE | 1.016 | 31.12565 | TRUE | 2.033 | 31.44066 | 0.3150148 | 0.4997541 | 0.09371813 | TRUE | 1.826 | 31.44861 | 0.3229593 | 0.5564074 | 0.09391893 | TRUE | 1.826 | 31.44631 | 0.320662 | 0.5564074 | 0.09379141 | TRUE | 1.89 | 31.44742 | 0.3217702 | 0.5375661 | 0.09381115 |
| test:test-nlm-cens-t | test:test-nlm-cens-t#3 | TRUE | 1.088 | 31.19994 | TRUE | 2.423 | 38.1996 | 6.999666 | 0.4490301 | 0.2666448 | TRUE | 1.99 | 38.1996 | 6.999666 | 0.5467337 | 0.2666448 | TRUE | 2.042 | 38.1996 | 6.999666 | 0.532811 | 0.2666448 | TRUE | 2.075 | 38.1996 | 6.999666 | 0.5243373 | 0.2666448 |
| test:test-nlm-cens-t | test:test-nlm-cens-t#4 | TRUE | 1.729 | 43.36694 | TRUE | 2.427 | 36.92497 | -6.441965 | 0.7124021 | 0.05680591 | TRUE | 2.097 | 36.92497 | -6.441965 | 0.8245112 | 0.05680591 | TRUE | 2.119 | 36.92497 | -6.441965 | 0.8159509 | 0.05680591 | TRUE | 2.129 | 36.92497 | -6.441965 | 0.8121184 | 0.05680591 |
| test:test-nlm-cens-t | test:test-nlm-cens-t#5 | TRUE | 1.857 | 43.36694 | TRUE | 3.079 | 36.92497 | -6.441965 | 0.6031179 | 0.05680591 | TRUE | 2.424 | 36.92497 | -6.441965 | 0.7660891 | 0.05680591 | TRUE | 1.903 | 36.92497 | -6.441965 | 0.9758276 | 0.05680591 | TRUE | 1.908 | 36.92497 | -6.441965 | 0.9732704 | 0.05680591 |
| test:test-nlm-cens-t | test:test-nlm-cens-t#6 | TRUE | 1.055 | 43.36694 | TRUE | 2.205 | 36.92497 | -6.441965 | 0.478458 | 0.05680591 | TRUE | 1.931 | 36.92497 | -6.441965 | 0.546349 | 0.05680591 | TRUE | 1.96 | 36.92497 | -6.441965 | 0.5382653 | 0.05680591 | TRUE | 2.126 | 36.92497 | -6.441965 | 0.4962371 | 0.05680591 |
| test:test-nlm-lik-contrib | test:test-nlm-lik-contrib#1 | TRUE | 2.04 | 216.1453 | TRUE | 2.483 | 216.1514 | 0.006166986 | 0.8215868 | 0.00165489 | TRUE | 2.437 | 216.1517 | 0.006433385 | 0.8370948 | 0.001932282 | TRUE | 2.469 | 216.1516 | 0.006330244 | 0.8262454 | 0.001807752 | TRUE | 2.446 | 216.1518 | 0.006484259 | 0.8340147 | 0.00194956 |
| test:test-nlm | test:test-nlm#1 | TRUE | 1.297 | -672.9741 | TRUE | 1.009 | -672.9756 | -0.001550284 | 1.285431 | 0.07084283 | TRUE | 0.972 | -672.9755 | -0.001358496 | 1.334362 | 0.08697369 | TRUE | 0.981 | -672.9756 | -0.001549998 | 1.32212 | 0.07055283 | TRUE | 0.958 | -672.9753 | -0.001238407 | 1.353862 | 0.07621495 |
| test:test-nlm | test:test-nlm#2 | TRUE | 2.448 | 496.285 | TRUE | 1.748 | 496.2426 | -0.04242127 | 1.400458 | 0.01143937 | TRUE | 1.506 | 496.2429 | -0.04209633 | 1.625498 | 0.01247407 | TRUE | 1.572 | 496.2427 | -0.0423753 | 1.557252 | 0.01146277 | TRUE | 1.6 | 496.2449 | -0.04011282 | 1.53 | 0.01399659 |
| test:test-nlm | test:test-nlm#3 | TRUE | 1.05 | -161.754 | TRUE | 2.496 | -161.8738 | -0.1198215 | 0.4206731 | 0.05341999 | TRUE | 2.179 | -161.3197 | 0.4342762 | 0.4818724 | 0.2782779 | TRUE | 2.233 | -161.8738 | -0.1198178 | 0.4702194 | 0.05299995 | TRUE | 2.243 | -161.8738 | -0.1198048 | 0.468123 | 0.05460171 |
| test:test-nlm | test:test-nlm#4 | TRUE | 1.047 | -161.754 | TRUE | 2.049 | -161.8738 | -0.1198215 | 0.510981 | 0.05341999 | TRUE | 1.787 | -161.3197 | 0.4342762 | 0.5858982 | 0.2782779 | TRUE | 1.815 | -161.8738 | -0.1198178 | 0.5768595 | 0.05299995 | TRUE | 1.805 | -161.8738 | -0.1198048 | 0.5800554 | 0.05460171 |
| test:test-nlm | test:test-nlm#5 | TRUE | 1.686 | 98.4812 | TRUE | 1.749 | 98.48119 | -1.042989e-05 | 0.9639794 | 0.00031254 | TRUE | 0.753 | 98.48119 | -9.363996e-06 | 2.239044 | 0.0003679814 | TRUE | 0.799 | 98.48119 | -9.800433e-06 | 2.110138 | 0.0003530292 | TRUE | 0.771 | 98.48299 | 0.001785315 | 2.18677 | 0.007202662 |
| test:test-nlmsetup-fresh-rx | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-optexpression-ebe | test:test-optexpression-ebe#1 | TRUE | 1.144 | 342.5069 | TRUE | 2.194 | 216.1514 | -126.3555 | 0.5214221 | 0.5636927 | TRUE | 2.254 | 216.1517 | -126.3552 | 0.5075421 | 0.5638879 | TRUE | 1.673 | 216.1516 | -126.3553 | 0.6838016 | 0.5638456 | TRUE | 1.709 | 216.1518 | -126.3551 | 0.6693973 | 0.5639268 |
| test:test-optexpression-ebe | test:test-optexpression-ebe#2 | TRUE | 2.482 | 342.5069 | TRUE | 2.039 | 216.1514 | -126.3555 | 1.217263 | 0.5636927 | TRUE | 1.753 | 216.1517 | -126.3552 | 1.415859 | 0.5638879 | TRUE | 1.791 | 216.1516 | -126.3553 | 1.385818 | 0.5638456 | TRUE | 1.777 | 216.1518 | -126.3551 | 1.396736 | 0.5639268 |
| test:test-optim | test:test-optim#1 | TRUE | 0.726 | -695.5321 | TRUE | 1.902 | -695.5322 | -2.388019e-05 | 0.3817035 | 0.004732219 | TRUE | 0.858 | -695.5314 | 0.0007521613 | 0.8461538 | 0.007899038 | TRUE | 0.845 | -695.5322 | -2.359206e-05 | 0.8591716 | 0.005327788 | TRUE | 0.827 | -695.5319 | 0.0002170529 | 0.8778718 | 0.009432982 |
| test:test-splitbolus-interp | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA |
