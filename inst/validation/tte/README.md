# Time-to-event validation scenarios

Reference scenarios for `tte()` estimation (nlmixr2/nlmixr2est#1461).
Each scenario is simulated with rxode2's `tte()` from known parameters
and is fitted by nlmixr2est, NONMEM and Monolix from the same data set.

| Scenario      | Model                                                              |
|---------------|--------------------------------------------------------------------|
| `exponential` | constant hazard, single event, right-censored at 20 (no random effect) |
| `weibull`     | Weibull hazard, random effect on the scale, single event           |
| `rtte`        | repeated events, constant hazard with a subject frailty            |
| `gapTime`     | repeated events, Weibull hazard in the time since the last event (`tteTime()`) |
| `pkFixed`     | hazard inhibited by a one-compartment PK model, individual PK known (`CLI`, `VI`) |
| `pkJoint`     | the same hazard with the PK estimated jointly from a PK endpoint   |
| `tvCov`       | repeated events, a covariate `X` switching from 0 to 1 at `TSW`; nested model without it |
| `interval`    | interval-censored single event (inspections every 4)               |
| `lateEntry`   | Weibull single event, observation starting between 0 and 10        |
| `dropout`     | single event with dropout (`dropout = TRUE`) as a competing process |

## Files

- `scenarios.R` -- the scenarios: the nlmixr2 model (true values in
  `ini()`), the starting values every tool uses, and the simulation
  design.  `tteValidationData(name, dir)` reads a data set for nlmixr2.
- `simulate.R` -- regenerates `data/` (`Rscript simulate.R`).
- `data/<scenario>.csv` -- the data sets.

## Data sets

Every data set has NONMEM-style columns, so all tools read the same file:

| Column | Meaning |
|--------|---------|
| `ID`, `TIME` | subject and time |
| `EVID` | 0 observation, 1 dose, 2 covariate change |
| `AMT`  | dose amount |
| `CMT`  | NONMEM compartment: the dosing compartment, or the cumulative-hazard (or central) compartment of an observation |
| `DV`   | 0 or 1 for an event endpoint; the concentration for a PK endpoint |
| `MDV`  | 1 for a dose or covariate record |
| `DVID` | endpoint number, in the model's endpoint order (0 for a dose) |

The event records follow rxode2/Monolix conventions: each subject's
first record of an endpoint (`DV = 0`) starts its observation, `DV = 1`
marks an event (for `interval`, the end of an inspection interval that
holds an event) and a final `DV = 0` record right-censors it.

Conventions that differ between the tools are carried as extra columns:

- `tvCov$XNM`, `gapTime$TLAST` -- NONMEM uses the data items of the
  record that *ends* an integration interval, so these hold the value in
  force *before* each record (the covariate on the previous record; the
  last event time before the record).  rxode2 and Monolix carry `X`
  forward from the record where it changes.
- `dropout` -- a dropout event ends every endpoint's observation, but
  rxode2's simulation leaves no record of the other endpoint there.  The
  data set adds a right-censoring `DV = 0` record for `ev` at each
  dropout, without which `ev`'s hazard would be over-estimated.
