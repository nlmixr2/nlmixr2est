# foceiControl(hessianMethod=) benchmark (non-normal-endpoint models)

Total models: 2

## hessianMethod="fd" (converged: 2/2)

- Median time: 4.002s
- Median |objf diff| vs fd: 0
- Total nHessianQN calls: 0

## hessianMethod="bfgs" (converged: 2/2)

- Median time: 0.732s
- Median |objf diff| vs fd: 0.1597
- Total nHessianQN calls: 28005

## hessianMethod="sr1" (converged: 2/2)

- Median time: 0.815s
- Median |objf diff| vs fd: 0.06688
- Total nHessianQN calls: 30281

## hessianMethod="bofill" (converged: 2/2)

- Median time: 0.648s
- Median |objf diff| vs fd: 0.06069
- Total nHessianQN calls: 25869

| model | fd ok | fd time | fd objf | fd dObjf | bfgs ok | bfgs time | bfgs objf | bfgs dObjf | sr1 ok | sr1 time | sr1 objf | sr1 dObjf | bofill ok | bofill time | bofill objf | bofill dObjf |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| poisson-e0 | TRUE | 4.252 | 318.3239 | 0 | TRUE | 0.41 | 318.2181 | -0.1057629 | TRUE | 0.222 | 318.3221 | -0.001770251 | TRUE | 0.348 | 318.3221 | -0.001770251 |
| ll-gaussian-equivalent | TRUE | 3.752 | 130.9302 | 0 | TRUE | 1.055 | 130.7165 | -0.213719 | TRUE | 1.407 | 131.0622 | 0.1319926 | TRUE | 0.949 | 131.0498 | 0.119606 |
