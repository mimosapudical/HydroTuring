# momentum/wave-celerity-bounds

This probe implements accepted proposal #148. It asks a question deliberately
left open by `momentum/routing-lag-consistency`: does transient propagation
speed respond to hydraulic state, rather than merely staying causal and getting
slower in longer reaches?

## Experiment

For each seed, the generator draws one mild prismatic rectangular reach
geometry and creates three positive base discharges: 8, 16 and 32 m3/s. Each
state is run twice with byte-identical forcing, once at 4 km and once at 12 km.
After two days of spinup and one scored day at the base state, a six-hour 5%
inflow pulse is applied.

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

`reference_saint_venant` is the must-pass model. Its daily steady path is
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

`reference_saint_venant` is the CI must-pass reference, but it is not by itself
the independent submitted-model evidence required by the exception in
`docs/writing-a-probe.md`. Before merge, a submitted physical model that
actually implements reach-wave propagation must also pass and have that row
archived in `models/result.csv`. `wflow_sbm` is the first candidate to audit,
because its manifest explicitly describes kinematic-wave river routing; its
adapter must not be declared compatible until it genuinely consumes this
probe's reach geometry and the archived run passes.

## Reproduction

```sh
ht validate
ht gate --probe momentum/wave-celerity-bounds
pytest -q
```

No generated forcing or model output is committed.
