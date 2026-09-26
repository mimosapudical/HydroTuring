# momentum/wave-celerity-bounds

This probe implements accepted proposal #148. It asks a question deliberately
left open by `momentum/routing-lag-consistency`: does transient propagation
speed respond to hydraulic state, rather than merely staying causal and getting
slower in longer reaches?

## Experiment

For each seed, the generator draws one mild prismatic rectangular reach
geometry and creates three positive base discharges: 8, 16 and 32 m3/s. Only two model runs are needed per seed: a 4 km `short` reach and a 20 km
`long` reach with byte-identical forcing. Inside each run, low, medium and high
hydraulic states occupy successive eight-day blocks. Each block gets four days
to settle, then a six-hour +5% pulse and more than three days of response tail.

Both runs are read at the reach outlet, so the paired propagation distance is

```
Delta x = L_long - L_short.
```

The criterion removes the pre-event discharge, takes the centroid of the
positive transient response, and estimates

```
c_obs = Delta x / (t_long - t_short).
```

For each hydraulic state, differencing the paired reaches cancels the common
forcing and runoff-generation clock. A 72-hour local response window ends
before the next state transition, so the following plateau cannot pull the
previous pulse centroid downstream. The generated record is hourly so timing
quantisation is smaller than the short/long travel-time difference.

## Independent expectation

For a rectangular section,

```
Q = (1/n) A R_h^(2/3) sqrt(S)
R_h = A / P
P = b + 2h
A = b h
```

and differentiating this declared rating gives

```
c_kin = dQ/dA
      = Q [5/(3A) - 4/(3 b P)].
```

The criterion obtains the base discharge from the model output, inverts the
monotone Manning relation for depth, and evaluates this derivative. It does not
reuse a propagation time from the reference model.

The three requirements are:

1. `c_obs > 0`;
2. `abs(c_obs - c_kin) / c_kin <= 0.05`;
3. `c_low < c_medium < c_high`.

The five-percent physical allowance is fixed before gate evaluation. Proposal
#148's feasibility audit found the existing 64-cell Rusanov reference biased
low by about 9--13%, 128 cells by about 4--8%, and tested 256-cell cases by
about 2--4%. The production transient path therefore uses 256 cells rather than
widening the physical tolerance to cover coarse-grid diffusion.

The production generator keeps the synthetic section deliberately wide
(60--90 m for 8--32 m3/s) and uses a 4 km / 20 km pair. Because the reported
discharge is read at the reach outlet, the paired propagation distance is the
full 16 km length difference. The resulting short/long centroid delay remains
resolved at the one-hour output step even at high flow.

## Baselines and applicability

| Model | #148 status | Reason |
| --- | --- | --- |
| `reference_bucket` | N/A | no geometry-aware reach-wave process |
| `flex_lumped` | N/A | existing lag uses catchment channel-length fields and fixed 1 m/s celerity, not this reach contract |
| `flex_topo` | N/A | does not consume this reach geometry |
| `sacsma_snow17` | N/A | does not consume this reach geometry |
| `reference_saint_venant` | gate must-pass | independent finite-volume dynamic-wave solve |
| `reference_fixed_celerity` | gate must-fail | deliberately fixed 1 m/s propagation |


`reference_saint_venant` is the CI must-pass model. Its daily steady path is
unchanged; for sub-daily forcing it advances the same continuity and momentum
equations continuously through each output interval using a 256-cell grid.
The criterion does not call the solver's fluxes or internal wave speeds.

`reference_fixed_celerity` is the deliberately broken control. It transports
the prescribed discharge downstream at exactly 1 m/s for every hydraulic
state. It remains positive and length-dependent, but cannot satisfy the
low/medium/high state response.

The four general physical baselines are intentionally not forced through this
probe. `reference_bucket`, `flex_topo` and `sacsma_snow17` do not declare
consumption of reach length; `flex_lumped` consumes catchment channel-length
fields used by the separate routing-lag probe, not this probe's
`reach_length_m`. Under the harness contract they are therefore
`N/A (INCOMPATIBLE)`, not FAIL.

`reference_saint_venant` is the trusted in-repository numerical reference,
not the independent submitted-model evidence required by
`docs/writing-a-probe.md`. Wflow was audited first but rejected for this role:
its one-cell adapter reports total outlet flow dominated by overland and lateral
subsurface routing, so its public discharge does not isolate the river wave
this probe measures.

LISFLOOD was also audited as an independent candidate and rejected
rather than tuned to pass. Wiring the declared geometry into its native
kinematic-wave channel produced, on gate seed 397273707, observed paired
celerities of about 12.9, 43.7 and 63.1 m/s versus Manning expectations of
about 0.66, 0.78 and 1.00 m/s. Refining LISFLOOD's channel routing sub-step
from 3600 s to 300 s changed those values only to about 12.8, 43.4 and
62.7 m/s. The mismatch is therefore not the probe's five-percent allowance
or hourly timing quantisation. No LISFLOOD PASS row is claimed or archived.

The independent submitted-physical-model evidence required by
`docs/writing-a-probe.md` remains the outstanding merge prerequisite.

## Reproduction

```sh
ht validate
ht gate --probe momentum/wave-celerity-bounds
pytest -q
```

No generated forcing or model output is committed.
