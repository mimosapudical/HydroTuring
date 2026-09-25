# momentum/wave-celerity-bounds

This probe implements accepted proposal #148. It asks a question deliberately
left open by `momentum/routing-lag-consistency`: does transient propagation
speed respond to hydraulic state, rather than merely staying causal and getting
slower in longer reaches?

## Experiment

For each seed, the generator draws one mild prismatic rectangular reach
geometry and creates three positive base discharges: 8, 16 and 32 m3/s. Each
state is run twice with byte-identical forcing, once at 4 km and once at 12 km.
After twenty days of spinup and one scored day at the base state, a six-hour
5% inflow pulse is applied. The case also supplies a small 80 mm soil capacity
and 1.5 mm canopy capacity so a full hydrologic model can settle before the
routing transient is timed, rather than leaving a slowly filling catchment
store inside the response centroid.

The gauge is at the reach centre, so the extra propagation distance is

```
Delta x = (L_long - L_short) / 2.
```

The criterion removes the pre-event discharge, takes the centroid of the
positive transient response, and estimates

```
c_obs = Delta x / (t_long - t_short).
```

Differencing the paired reaches cancels the common forcing and runoff-generation
clock. The generated record is hourly so timing quantisation is much smaller
than the travel-time difference.

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

## Baselines and applicability

| Model | #148 status | Reason |
| --- | --- | --- |
| `reference_bucket` | N/A | no geometry-aware reach-wave process |
| `flex_lumped` | N/A | existing lag uses catchment channel-length fields and fixed 1 m/s celerity, not this reach contract |
| `flex_topo` | N/A | does not consume this reach geometry |
| `sacsma_snow17` | N/A | does not consume this reach geometry |
| `reference_saint_venant` | gate must-pass | independent finite-volume dynamic-wave solve |
| `reference_fixed_celerity` | gate must-fail | deliberately fixed 1 m/s propagation |
| `wflow_sbm` | gate must-pass, pending Docker evidence | submitted Wflow.jl model; its river kinematic wave consumes the explicit reach geometry |


`reference_saint_venant` and `wflow_sbm` are the two must-pass models. The Saint-Venant reference's daily steady path is
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

`reference_saint_venant` is the trusted in-repository numerical reference.
`wflow_sbm` is the independent submitted physical baseline required by
`docs/writing-a-probe.md`. Its adapter genuinely maps explicit
`reach_length_m`, `width_m`, `slope` and `manning_n` into Wflow's river
static maps and lets Wflow's own kinematic-wave routing respond to them.

Wflow's public `dis` output is still the model's total outlet flow converted
from `mrro`; it is not a hidden river-only diagnostic. Therefore the Docker
gate is load-bearing: it must demonstrate that, under this controlled case, the
paired short/long timing remains identifiable in the reported outlet response.
No PASS row is written in advance. After the real Docker gate passes, the
`wflow_sbm` row for this probe must be appended to `models/result.csv`.

## Reproduction

```sh
ht validate
ht gate --probe momentum/wave-celerity-bounds
pytest -q
```

No generated forcing or model output is committed.
