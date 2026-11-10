# SAEM readback in the workers: tested, no wall-time benefit

Experiment of 2026-11-10, kept as a record. Nothing from it is on main.

## What was changed

After each bulk `par_solve()`, `user_function()` in `src/saem.cpp` reads every
subject's predictions back in one serial loop: `saem_lhs()` (the model's
`calc_lhs`) at each record, and one `g` row per observation. With several
threads the solve is parallel but this readback runs on one.

The experiment moved that loop into the per-subject workers
(`nmForEachSubject()`), as `saemReadRowsPooled()` already does for the pooled
solve. Each subject's rows start at a fixed offset, the number of observations
of the subjects before it, so `g` does not depend on the thread count or the
schedule. `do_mcmc()`'s copies were not part of the experiment.

## A race, and its fix

The first version read `ind->lhs` before `iniSubjectE(..., inLhs = 1, ...)`.
That call points `ind->lhs` at the calling thread's slice, so two threads
could write to the same slice. At 4 threads the estimates then changed from
run to run. Reading `getIndLhs(ind)` after `iniSubjectE()` fixed it.

## Results

With the fix, the fit is bit-identical to main at 1, 2 and 4 threads. The
check compared theta, omega and the etas through their `%a` representations.
`test-saem-readback-threads.R` on the park branch pins this on `theo_sd`.

The readback was 1.5% of the single-thread instructions in the profile that
proposed the change, so the most it could save is about 1% of the
single-thread time.

Wall time on a shared 32-core Linux machine, 4 threads, 5 runs per build,
interleaved, with `bench.R` (a 2-compartment oral model with proportional
error, 3 states, 200 simulated subjects with 3 doses and 12 samples each;
`saem` with `nBurn = 60`, `nEm = 60`, `covMethod = ""`,
`calcTables = FALSE`):

| build | median | range |
|---|---|---|
| main | 116.5 s | 106.8-127.8 s |
| readback in the workers | 115.8 s | 106.0-117.8 s |

The difference is within the noise. `covMethod = ""` skips the covariance
step (default `linFim`), which took about 685 s of the 200-subject fit on one
thread, so it is out of the comparison.

## Conclusion

The change gives no wall-time benefit, so it was not merged. It is parked on
the nlmixr2est branch `park/saem-readback-workers` at
160d485bd557d8a378d7ca2f4718e7fbf9dd0f8b, together with its thread test.

## Running the benchmark

```sh
# once per build, with the build to time first on the library path
Rscript design/saem-readback/bench.R main 4 1 200 60 60
```

This runs one fit (`nrep` = 1). The table above came from 5 such runs per
build, alternating the builds one fit at a time; compare the medians. The
fingerprint printed with each run must be the same for every build and every
thread count.
