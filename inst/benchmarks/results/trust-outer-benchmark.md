# est="trust" hessianMethod= vs est="bobyqa" benchmark

Total corpus entries: 23

## hessianMethod="fd" (both converged: 21/23)

- Median speedup (bobyqa time / trust time): 0.61x
- Median |objf diff|: 1.551
- Median max |param diff|: 0.09372
- Entries with |objf diff| > 1 (likely a different local optimum, not just numeric noise): 11

## hessianMethod="bfgs" (both converged: 21/23)

- Median speedup (bobyqa time / trust time): 0.73x
- Median |objf diff|: 0.4343
- Median max |param diff|: 0.09392
- Entries with |objf diff| > 1 (likely a different local optimum, not just numeric noise): 9

## hessianMethod="sr1" (both converged: 21/23)

- Median speedup (bobyqa time / trust time): 0.77x
- Median |objf diff|: 1.532
- Median max |param diff|: 0.09379
- Entries with |objf diff| > 1 (likely a different local optimum, not just numeric noise): 11

## hessianMethod="bofill" (both converged: 21/23)

- Median speedup (bobyqa time / trust time): 0.76x
- Median |objf diff|: 1.551
- Median max |param diff|: 0.09381
- Entries with |objf diff| > 1 (likely a different local optimum, not just numeric noise): 11

| source | model | bobyqa ok | bobyqa time | bobyqa objf | trust[fd] ok | trust[fd] time | trust[fd] objf | trust[fd] dObjf | trust[fd] speedup | trust[fd] maxParDiff | trust[bfgs] ok | trust[bfgs] time | trust[bfgs] objf | trust[bfgs] dObjf | trust[bfgs] speedup | trust[bfgs] maxParDiff | trust[sr1] ok | trust[sr1] time | trust[sr1] objf | trust[sr1] dObjf | trust[sr1] speedup | trust[sr1] maxParDiff | trust[bofill] ok | trust[bofill] time | trust[bofill] objf | trust[bofill] dObjf | trust[bofill] speedup | trust[bofill] maxParDiff |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| test:test-cens-dist-warn | test:test-cens-dist-warn#1 | TRUE | 1.222 | 258.1839 | TRUE | 1.349 | 258.1847 | 0.0008755444 | 0.9058562 | 0.9544515 | TRUE | 1.237 | 258.7913 | 0.6074574 | 0.9878739 | 0.7369447 | TRUE | 1.237 | 258.1846 | 0.0007097233 | 0.9878739 | 1.606828 | TRUE | 1.327 | 258.1845 | 0.0006446797 | 0.9208742 | 0.4204929 |
| test:test-cens-dist-warn | test:test-cens-dist-warn#2 | TRUE | 1.454 | 22.28706 | TRUE | 2.399 | 22.2633 | -0.02375948 | 0.6060859 | 0.05169591 | TRUE | 1.789 | 22.23087 | -0.05619005 | 0.8127446 | 0.04023307 | TRUE | 1.821 | 22.30889 | 0.02183063 | 0.7984624 | 0.06458114 | TRUE | 1.792 | 22.22968 | -0.05737305 | 0.8113839 | 0.04074211 |
| test:test-lincmt-ode-fit | test:test-lincmt-ode-fit#1 | TRUE | 3.195 | 536.3186 | TRUE | 29.331 | 179350.7 | 178814.4 | 0.1089291 | 1.845656 | TRUE | 20.149 | 179350.7 | 178814.4 | 0.1585687 | 1.845656 | TRUE | 21.615 | 179350.7 | 178814.4 | 0.147814 | 1.845656 | TRUE | 24.181 | 179350.7 | 178814.4 | 0.1321285 | 1.845656 |
| test:test-lincmt-ode-fit | test:test-lincmt-ode-fit#2 | TRUE | 4.314 | 536.3186 | TRUE | 31.373 | 497.2813 | -39.03738 | 0.1375068 | 1.036766 | TRUE | 16.202 | 497.2813 | -39.03734 | 0.2662634 | 1.037014 | TRUE | 16.428 | 497.2813 | -39.03737 | 0.2626004 | 1.036816 | TRUE | 16.481 | 519.6119 | -16.70678 | 0.261756 | 2.542535 |
| test:test-matexp | test:test-matexp#1 | TRUE | 1.54 | -140.4195 | TRUE | 2.468 | -141.9705 | -1.551088 | 0.623987 | 0.6334603 | TRUE | 2.587 | -140.8181 | -0.398622 | 0.5952841 | 0.08164248 | TRUE | 2.112 | -141.9517 | -1.532214 | 0.7291667 | 0.5626262 | TRUE | 2.12 | -141.9705 | -1.551081 | 0.7264151 | 0.6334126 |
| test:test-matexp | test:test-matexp#2 | TRUE | 1.593 | -140.4195 | TRUE | 2.079 | -141.9705 | -1.551088 | 0.7662338 | 0.6334603 | TRUE | 2.172 | -140.8181 | -0.398622 | 0.7334254 | 0.08164248 | TRUE | 1.729 | -141.9517 | -1.532214 | 0.9213418 | 0.5626262 | TRUE | 1.699 | -141.9705 | -1.551081 | 0.9376104 | 0.6334126 |
| test:test-nlm-cens-t | test:test-nlm-cens-t#1 | TRUE | 1.085 | 42.80645 | TRUE | 2.459 | 48.2984 | 5.491957 | 0.4412363 | 0.1668793 | TRUE | 2.076 | 48.2984 | 5.491957 | 0.5226397 | 0.1668793 | TRUE | 2.107 | 48.2984 | 5.491957 | 0.5149502 | 0.1668793 | TRUE | 2.136 | 48.2984 | 5.491957 | 0.5079588 | 0.1668793 |
| test:test-nlm-cens-t | test:test-nlm-cens-t#2 | TRUE | 1.206 | 31.12565 | TRUE | 2.286 | 31.44066 | 0.3150148 | 0.5275591 | 0.09371813 | TRUE | 1.825 | 31.44861 | 0.3229593 | 0.6608219 | 0.09391893 | TRUE | 1.738 | 31.44631 | 0.320662 | 0.693901 | 0.09379141 | TRUE | 1.744 | 31.44742 | 0.3217702 | 0.6915138 | 0.09381115 |
| test:test-nlm-cens-t | test:test-nlm-cens-t#3 | TRUE | 1.03 | 31.19994 | TRUE | 2.339 | 38.1996 | 6.999666 | 0.4403591 | 0.2666448 | TRUE | 1.891 | 38.1996 | 6.999666 | 0.5446854 | 0.2666448 | TRUE | 1.865 | 38.1996 | 6.999666 | 0.5522788 | 0.2666448 | TRUE | 1.94 | 38.1996 | 6.999666 | 0.5309278 | 0.2666448 |
| test:test-nlm-cens-t | test:test-nlm-cens-t#4 | TRUE | 1.491 | 43.36694 | TRUE | 2.151 | 36.92497 | -6.441965 | 0.693166 | 0.05680591 | TRUE | 1.973 | 36.92497 | -6.441965 | 0.755702 | 0.05680591 | TRUE | 1.935 | 36.92497 | -6.441965 | 0.7705426 | 0.05680591 | TRUE | 1.987 | 36.92497 | -6.441965 | 0.7503775 | 0.05680591 |
| test:test-nlm-cens-t | test:test-nlm-cens-t#5 | TRUE | 1.625 | 43.36694 | TRUE | 2.402 | 36.92497 | -6.441965 | 0.6765196 | 0.05680591 | TRUE | 2.101 | 36.92497 | -6.441965 | 0.7734412 | 0.05680591 | TRUE | 2.105 | 36.92497 | -6.441965 | 0.7719715 | 0.05680591 | TRUE | 2.104 | 36.92497 | -6.441965 | 0.7723384 | 0.05680591 |
| test:test-nlm-cens-t | test:test-nlm-cens-t#6 | TRUE | 1.294 | 43.36694 | TRUE | 2.573 | 36.92497 | -6.441965 | 0.5029149 | 0.05680591 | TRUE | 2.271 | 36.92497 | -6.441965 | 0.569793 | 0.05680591 | TRUE | 2.317 | 36.92497 | -6.441965 | 0.5584808 | 0.05680591 | TRUE | 1.74 | 36.92497 | -6.441965 | 0.7436782 | 0.05680591 |
| test:test-nlm-lik-contrib | test:test-nlm-lik-contrib#1 | TRUE | 1.659 | 216.1453 | TRUE | 1.951 | 216.1514 | 0.006166986 | 0.8503332 | 0.00165489 | TRUE | 1.757 | 216.1517 | 0.006433385 | 0.9442231 | 0.001932282 | TRUE | 1.751 | 216.1516 | 0.006330244 | 0.9474586 | 0.001807752 | TRUE | 1.755 | 216.1518 | 0.006484259 | 0.9452991 | 0.00194956 |
| test:test-nlm | test:test-nlm#1 | TRUE | 1.236 | -723.5896 | TRUE | 0.947 | -723.5932 | -0.003529778 | 1.305174 | 0.1836563 | TRUE | 0.901 | -723.5931 | -0.003509126 | 1.371809 | 0.1993866 | TRUE | 0.882 | -723.593 | -0.003331791 | 1.401361 | 0.1398983 | TRUE | 0.926 | -723.5847 | 0.004893633 | 1.334773 | 0.4302022 |
| test:test-nlm | test:test-nlm#2 | TRUE | 2.63 | 496.285 | TRUE | 1.979 | 496.2426 | -0.04242127 | 1.328954 | 0.01143937 | TRUE | 1.774 | 496.2429 | -0.04209633 | 1.482525 | 0.01247407 | TRUE | 1.791 | 496.2427 | -0.0423753 | 1.468453 | 0.01146277 | TRUE | 1.708 | 496.2449 | -0.04011282 | 1.539813 | 0.01399659 |
| test:test-nlm | test:test-nlm#3 | TRUE | 1.207 | -161.754 | TRUE | 2.69 | -161.8738 | -0.1198215 | 0.4486989 | 0.05341999 | TRUE | 2.397 | -161.3197 | 0.4342762 | 0.5035461 | 0.2782779 | TRUE | 2.41 | -161.8738 | -0.1198178 | 0.5008299 | 0.05299995 | TRUE | 2.488 | -161.8738 | -0.1198048 | 0.4851286 | 0.05460171 |
| test:test-nlm | test:test-nlm#4 | TRUE | 1.197 | -161.754 | TRUE | 2.271 | -161.8738 | -0.1198215 | 0.5270806 | 0.05341999 | TRUE | 1.951 | -161.3197 | 0.4342762 | 0.6135315 | 0.2782779 | TRUE | 1.607 | -161.8738 | -0.1198178 | 0.7448662 | 0.05299995 | TRUE | 1.576 | -161.8738 | -0.1198048 | 0.7595178 | 0.05460171 |
| test:test-nlm | test:test-nlm#5 | TRUE | 1.526 | 98.4812 | TRUE | 1.638 | 98.48119 | -1.042989e-05 | 0.9316239 | 0.00031254 | TRUE | 0.641 | 98.48119 | -9.363996e-06 | 2.380655 | 0.0003679814 | TRUE | 0.642 | 98.48119 | -9.800433e-06 | 2.376947 | 0.0003530292 | TRUE | 0.717 | 98.48299 | 0.001785315 | 2.128312 | 0.007202662 |
| test:test-nlmsetup-fresh-rx | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-optexpression-ebe | test:test-optexpression-ebe#1 | TRUE | 1.014 | 342.5069 | TRUE | 2.01 | 216.1514 | -126.3555 | 0.5044776 | 0.5636927 | TRUE | 1.763 | 216.1517 | -126.3552 | 0.575156 | 0.5638879 | TRUE | 1.817 | 216.1516 | -126.3553 | 0.5580627 | 0.5638456 | TRUE | 1.839 | 216.1518 | -126.3551 | 0.5513866 | 0.5639268 |
| test:test-optexpression-ebe | test:test-optexpression-ebe#2 | TRUE | 2.525 | 342.5069 | TRUE | 2.195 | 216.1514 | -126.3555 | 1.150342 | 0.5636927 | TRUE | 1.836 | 216.1517 | -126.3552 | 1.375272 | 0.5638879 | TRUE | 1.953 | 216.1516 | -126.3553 | 1.292883 | 0.5638456 | TRUE | 1.959 | 216.1518 | -126.3551 | 1.288923 | 0.5639268 |
| test:test-optim | test:test-optim#1 | TRUE | 0.735 | -688.5949 | TRUE | 1.938 | -688.595 | -0.0001177126 | 0.379257 | 0.01070052 | TRUE | 0.888 | -688.5946 | 0.0002557761 | 0.8277027 | 0.0208239 | TRUE | 0.941 | -688.595 | -7.287138e-05 | 0.781084 | 0.01091715 | TRUE | 0.876 | -688.5948 | 7.805241e-05 | 0.8390411 | 0.007214459 |
| test:test-splitbolus-interp | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA |
