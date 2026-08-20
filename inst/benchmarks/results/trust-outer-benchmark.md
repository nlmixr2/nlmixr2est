# est="trust" hessianMethod= vs est="bobyqa" benchmark

Total corpus entries: 23

## hessianMethod="fd" (both converged: 21/23)

- Median speedup (bobyqa time / trust time): 0.45x
- Median |objf diff|: 1.551
- Median max |param diff|: 0.2666
- Entries with |objf diff| > 1 (likely a different local optimum, not just numeric noise): 11

## hessianMethod="bfgs" (both converged: 21/23)

- Median speedup (bobyqa time / trust time): 0.75x
- Median |objf diff|: 5.492
- Median max |param diff|: 0.2666
- Entries with |objf diff| > 1 (likely a different local optimum, not just numeric noise): 11

## hessianMethod="sr1" (both converged: 21/23)

- Median speedup (bobyqa time / trust time): 0.76x
- Median |objf diff|: 1.525
- Median max |param diff|: 0.2666
- Entries with |objf diff| > 1 (likely a different local optimum, not just numeric noise): 11

## hessianMethod="bofill" (both converged: 21/23)

- Median speedup (bobyqa time / trust time): 0.73x
- Median |objf diff|: 1.55
- Median max |param diff|: 0.2666
- Entries with |objf diff| > 1 (likely a different local optimum, not just numeric noise): 11

| source | model | bobyqa ok | bobyqa time | bobyqa objf | trust[fd] ok | trust[fd] time | trust[fd] objf | trust[fd] dObjf | trust[fd] speedup | trust[fd] maxParDiff | trust[bfgs] ok | trust[bfgs] time | trust[bfgs] objf | trust[bfgs] dObjf | trust[bfgs] speedup | trust[bfgs] maxParDiff | trust[sr1] ok | trust[sr1] time | trust[sr1] objf | trust[sr1] dObjf | trust[sr1] speedup | trust[sr1] maxParDiff | trust[bofill] ok | trust[bofill] time | trust[bofill] objf | trust[bofill] dObjf | trust[bofill] speedup | trust[bofill] maxParDiff |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| test:test-cens-dist-warn | test:test-cens-dist-warn#1 | TRUE | 1.157 | 258.1839 | TRUE | 1.808 | 258.1843 | 0.0004511889 | 0.6399336 | 0.5102856 | TRUE | 1.245 | 258.1988 | 0.01497269 | 0.9293173 | 0.5386631 | TRUE | 1.224 | 258.1843 | 0.0004513506 | 0.9452614 | 0.5280679 | TRUE | 1.258 | 258.185 | 0.001188884 | 0.9197138 | 0.5113524 |
| test:test-cens-dist-warn | test:test-cens-dist-warn#2 | TRUE | 1.359 | 22.28706 | TRUE | 2.242 | 38.33546 | 16.04841 | 0.6061552 | 0.4896932 | TRUE | 1.708 | 38.33546 | 16.04841 | 0.7956674 | 0.4896932 | TRUE | 1.746 | 38.33546 | 16.04841 | 0.7783505 | 0.4896932 | TRUE | 1.719 | 38.33546 | 16.04841 | 0.7905759 | 0.4896932 |
| test:test-lincmt-ode-fit | test:test-lincmt-ode-fit#1 | TRUE | 3.215 | 536.3186 | TRUE | 21.604 | -352.8724 | -889.191 | 0.148815 | 883.8711 | TRUE | 18.055 | -352.8724 | -889.191 | 0.178067 | 883.8711 | TRUE | 18.123 | -352.8724 | -889.191 | 0.1773989 | 883.8711 | TRUE | 18.999 | -352.8724 | -889.191 | 0.1692194 | 883.8711 |
| test:test-lincmt-ode-fit | test:test-lincmt-ode-fit#2 | TRUE | 3.245 | 536.3186 | TRUE | 509.249 | 2104.291 | 1567.972 | 0.006372128 | 1.423501 | TRUE | 20.328 | 2104.291 | 1567.972 | 0.159632 | 1.423501 | TRUE | 20.752 | 2104.29 | 1567.971 | 0.1563705 | 1.423501 | TRUE | 21.131 | 2104.291 | 1567.972 | 0.1535659 | 1.423501 |
| test:test-matexp | test:test-matexp#1 | TRUE | 1.551 | -140.4195 | TRUE | 2.517 | -141.9705 | -1.551089 | 0.6162098 | 0.6339014 | TRUE | 2.072 | -62.38832 | 78.03113 | 0.7485521 | 0.3596372 | TRUE | 2.176 | -141.9605 | -1.541028 | 0.7127757 | 0.5993832 | TRUE | 2.22 | -141.9697 | -1.550222 | 0.6986486 | 0.6310022 |
| test:test-matexp | test:test-matexp#2 | TRUE | 1.517 | -140.4195 | TRUE | 2.182 | -141.9705 | -1.551089 | 0.6952337 | 0.6339014 | TRUE | 1.685 | -62.38832 | 78.03113 | 0.9002967 | 0.3596372 | TRUE | 1.777 | -141.9448 | -1.525383 | 0.853686 | 0.5539603 | TRUE | 1.818 | -141.9697 | -1.550222 | 0.8344334 | 0.6310022 |
| test:test-nlm-cens-t | test:test-nlm-cens-t#1 | TRUE | 0.962 | 42.80645 | TRUE | 2.168 | 48.2984 | 5.491957 | 0.4437269 | 0.1668793 | TRUE | 1.765 | 48.2984 | 5.491957 | 0.5450425 | 0.1668793 | TRUE | 1.755 | 48.2984 | 5.491957 | 0.5481481 | 0.1668793 | TRUE | 1.771 | 48.2984 | 5.491957 | 0.5431959 | 0.1668793 |
| test:test-nlm-cens-t | test:test-nlm-cens-t#2 | TRUE | 1.011 | 31.12565 | TRUE | 2.206 | 38.13989 | 7.014235 | 0.4582956 | 0.2738045 | TRUE | 1.844 | 38.13989 | 7.014235 | 0.5482646 | 0.2738045 | TRUE | 1.867 | 38.13989 | 7.014235 | 0.5415104 | 0.2738045 | TRUE | 1.922 | 38.13989 | 7.014235 | 0.5260146 | 0.2738045 |
| test:test-nlm-cens-t | test:test-nlm-cens-t#3 | TRUE | 1.083 | 31.19994 | TRUE | 2.528 | 38.1996 | 6.999666 | 0.4284019 | 0.2666448 | TRUE | 2.023 | 38.1996 | 6.999666 | 0.5353435 | 0.2666448 | TRUE | 1.982 | 38.1996 | 6.999666 | 0.5464178 | 0.2666448 | TRUE | 1.985 | 38.1996 | 6.999666 | 0.5455919 | 0.2666448 |
| test:test-nlm-cens-t | test:test-nlm-cens-t#4 | TRUE | 1.685 | 43.36694 | TRUE | 2.408 | 43.36157 | -0.005366348 | 0.6997508 | 0.1629693 | TRUE | 2.097 | 43.36157 | -0.005366348 | 0.8035289 | 0.1629693 | TRUE | 2.077 | 43.36157 | -0.005366348 | 0.8112662 | 0.1629693 | TRUE | 2.064 | 43.36157 | -0.005366348 | 0.816376 | 0.1629693 |
| test:test-nlm-cens-t | test:test-nlm-cens-t#5 | TRUE | 1.506 | 43.36694 | TRUE | 2.111 | 43.36157 | -0.005366348 | 0.713406 | 0.1629693 | TRUE | 1.827 | 43.36157 | -0.005366348 | 0.8243021 | 0.1629693 | TRUE | 1.793 | 43.36157 | -0.005366348 | 0.8399331 | 0.1629693 | TRUE | 1.786 | 43.36157 | -0.005366348 | 0.8432251 | 0.1629693 |
| test:test-nlm-cens-t | test:test-nlm-cens-t#6 | TRUE | 0.987 | 43.36694 | TRUE | 2.198 | 43.36157 | -0.005366348 | 0.4490446 | 0.1629693 | TRUE | 1.826 | 43.36157 | -0.005366348 | 0.5405257 | 0.1629693 | TRUE | 1.877 | 43.36157 | -0.005366348 | 0.5258391 | 0.1629693 | TRUE | 1.88 | 43.36157 | -0.005366348 | 0.525 | 0.1629693 |
| test:test-nlm-lik-contrib | test:test-nlm-lik-contrib#1 | TRUE | 1.753 | 216.1453 | TRUE | 5.511 | 216.1514 | 0.006167006 | 0.3180911 | 0.001646493 | TRUE | 2.017 | 216.1514 | 0.006167452 | 0.8691125 | 0.001596452 | TRUE | 2.017 | 216.1514 | 0.00616745 | 0.8691125 | 0.001596552 | TRUE | 2.014 | 216.1514 | 0.006167454 | 0.8704071 | 0.001596326 |
| test:test-nlm | test:test-nlm#1 | TRUE | 1.35 | -711.1335 | TRUE | 2.992 | -711.1336 | -6.307744e-05 | 0.4512032 | 0.0068795 | TRUE | 1.038 | -711.1336 | -6.630461e-05 | 1.300578 | 0.01422672 | TRUE | 1.032 | -711.1335 | -2.175261e-05 | 1.30814 | 0.01936846 | TRUE | 1.048 | -711.1329 | 0.0006462618 | 1.288168 | 0.02186664 |
| test:test-nlm | test:test-nlm#2 | TRUE | 2.473 | 496.285 | TRUE | 2.222 | 1148.513 | 652.2276 | 1.112961 | 0.8410139 | TRUE | 1.521 | 1148.513 | 652.2276 | 1.625904 | 0.8410139 | TRUE | 1.557 | 1148.513 | 652.2275 | 1.588311 | 0.8410139 | TRUE | 1.578 | 1148.513 | 652.2276 | 1.567174 | 0.8410139 |
| test:test-nlm | test:test-nlm#3 | TRUE | 1.045 | -161.754 | TRUE | 2.576 | -161.8738 | -0.1198215 | 0.4056677 | 0.05342663 | TRUE | 2.162 | -161.3917 | 0.3622692 | 0.4833488 | 0.2135402 | TRUE | 2.209 | -161.8738 | -0.1198206 | 0.4730647 | 0.05376127 | TRUE | 2.248 | -161.8736 | -0.1195617 | 0.4648577 | 0.05215243 |
| test:test-nlm | test:test-nlm#4 | TRUE | 1.06 | -161.754 | TRUE | 2.136 | -161.8738 | -0.1198215 | 0.4962547 | 0.05342663 | TRUE | 1.755 | -161.3917 | 0.3622692 | 0.6039886 | 0.2135402 | TRUE | 1.875 | -161.8738 | -0.1198206 | 0.5653333 | 0.05376127 | TRUE | 1.734 | -161.8736 | -0.1195617 | 0.6113033 | 0.05215243 |
| test:test-nlm | test:test-nlm#5 | TRUE | 1.635 | 98.4812 | TRUE | 1.737 | 98.48119 | -1.042678e-05 | 0.9412781 | 0.0003082133 | TRUE | 0.718 | 98.4812 | 3.394639e-06 | 2.277159 | 0.0003426962 | TRUE | 0.783 | 98.48119 | -1.006477e-05 | 2.088123 | 0.0002709622 | TRUE | 0.725 | 98.48119 | -9.287834e-06 | 2.255172 | 0.0002488461 |
| test:test-nlmsetup-fresh-rx | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| test:test-optexpression-ebe | test:test-optexpression-ebe#1 | TRUE | 1.111 | 342.5069 | TRUE | 5.861 | 216.1514 | -126.3555 | 0.1895581 | 0.5636843 | TRUE | 2.372 | 216.1514 | -126.3554 | 0.4683811 | 0.5636343 | TRUE | 1.808 | 216.1514 | -126.3554 | 0.6144912 | 0.5636344 | TRUE | 1.82 | 216.1514 | -126.3554 | 0.6104396 | 0.5636342 |
| test:test-optexpression-ebe | test:test-optexpression-ebe#2 | TRUE | 2.444 | 342.5069 | TRUE | 6.077 | 216.1514 | -126.3555 | 0.4021721 | 0.5636843 | TRUE | 1.836 | 216.1514 | -126.3554 | 1.331155 | 0.5636343 | TRUE | 1.836 | 216.1514 | -126.3554 | 1.331155 | 0.5636344 | TRUE | 1.89 | 216.1514 | -126.3554 | 1.293122 | 0.5636342 |
| test:test-optim | test:test-optim#1 | TRUE | 0.659 | -706.9558 | TRUE | 3.873 | -706.9558 | -5.109187e-06 | 0.1701523 | 0.0006024971 | TRUE | 0.874 | -706.9558 | 3.62903e-05 | 0.7540046 | 0.003419314 | TRUE | 0.865 | -706.9557 | 7.487037e-05 | 0.7618497 | 0.00399977 | TRUE | 0.904 | -706.9556 | 0.0002447424 | 0.7289823 | 0.0100113 |
| test:test-splitbolus-interp | (whole file) | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA |
