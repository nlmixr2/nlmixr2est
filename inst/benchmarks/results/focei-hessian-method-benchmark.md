# foceiControl(hessianMethod=) benchmark (non-normal-endpoint models)

Total models: 2

## hessianMethod="fd" (converged: 2/2)

- Median time: 3.248s
- Median |objf diff| vs fd: 0
- Total nHessianQN calls: 0

## hessianMethod="bfgs" (converged: 2/2)

- Median time: 0.598s
- Median |objf diff| vs fd: 0.0261
- Total nHessianQN calls: 22140

## hessianMethod="sr1" (converged: 2/2)

- Median time: 0.945s
- Median |objf diff| vs fd: 0.06753
- Total nHessianQN calls: 24061

## hessianMethod="bofill" (converged: 2/2)

- Median time: 0.502s
- Median |objf diff| vs fd: 0.07856
- Total nHessianQN calls: 22476

| model | fd ok | fd time | fd objf | fd dObjf | bfgs ok | bfgs time | bfgs objf | bfgs dObjf | sr1 ok | sr1 time | sr1 objf | sr1 dObjf | bofill ok | bofill time | bofill objf | bofill dObjf |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| poisson-e0 | TRUE | 2.276 | 318.3235 | 0 | TRUE | 0.293 | 318.3635 | 0.04000911 | TRUE | 0.163 | 318.3635 | 0.04000911 | TRUE | 0.157 | 318.3635 | 0.04000911 |
| ll-gaussian-equivalent | TRUE | 4.22 | 130.9893 | 0 | TRUE | 0.902 | 130.9771 | -0.01218806 | TRUE | 1.727 | 131.0844 | 0.09504399 | TRUE | 0.848 | 131.1064 | 0.1171078 |
